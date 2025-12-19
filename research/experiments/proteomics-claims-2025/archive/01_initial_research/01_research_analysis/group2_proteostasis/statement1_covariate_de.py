"""
Finding Group 2 - Statement 1: Covariate-Controlled Differential Expression
Claim: A covariate-controlled differential expression analysis (age, PMI, and PatientID)
       across 5,853 proteins (BH-FDR) identified 2,115 proteins (36.14%) significantly
       altered between tau-positive and tau-negative neurons.

This is a comprehensive differential expression analysis with confounding variable control.
"""

import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

# ============================================
# THEORETICAL BACKGROUND
# ============================================

"""
Covariate-Controlled Analysis:
- Removes confounding effects from age, PMI (post-mortem interval), PatientID
- Uses linear models: Expression ~ TauStatus + Age + PMI + PatientID
- More accurate than simple t-tests when confounders present

Statistical Framework:
- Linear model: Y = β₀ + β₁·TauStatus + β₂·Age + β₃·PMI + β₄·PatientID + ε
- β₁ represents the tau effect after controlling for covariates
- FDR correction handles multiple testing across 5,853 proteins

References:
- Smyth (2004) Linear models and empirical Bayes methods
- Benjamini & Hochberg (1995) Controlling the false discovery rate
- ISLP Chapter 3: Linear Regression
"""

# ============================================
# STEP 1: Data Preparation and Validation
# ============================================

def prepare_data_for_de(adata):
    """
    Prepare data and validate covariates for differential expression
    """

    print("="*60)
    print("DATA PREPARATION")
    print("="*60)

    # Check dimensions
    print(f"Dataset dimensions: {adata.shape}")
    print(f"Number of cells: {adata.n_obs}")
    print(f"Number of proteins: {adata.n_vars}")

    # Check for required columns
    required_cols = ['tau_status', 'age', 'PMI', 'PatientID']
    available_cols = []
    missing_cols = []

    for col in required_cols:
        if col in adata.obs.columns:
            available_cols.append(col)
        else:
            # Try case-insensitive search
            matches = [c for c in adata.obs.columns if col.lower() in c.lower()]
            if matches:
                print(f"Found alternative for {col}: {matches[0]}")
                adata.obs[col] = adata.obs[matches[0]]
                available_cols.append(col)
            else:
                missing_cols.append(col)

    print(f"\nAvailable covariates: {available_cols}")
    print(f"Missing covariates: {missing_cols}")

    # Analyze tau status distribution
    if 'tau_status' in adata.obs.columns:
        tau_counts = adata.obs['tau_status'].value_counts()
        print(f"\nTau status distribution:")
        print(tau_counts)

    # Check covariate distributions
    for col in available_cols:
        if col == 'tau_status':
            continue
        print(f"\n{col} statistics:")
        if adata.obs[col].dtype in ['float64', 'int64']:
            print(f"  Mean: {adata.obs[col].mean():.2f}")
            print(f"  Std: {adata.obs[col].std():.2f}")
            print(f"  Range: [{adata.obs[col].min():.2f}, {adata.obs[col].max():.2f}]")
        else:
            print(f"  Unique values: {adata.obs[col].nunique()}")

    return available_cols, missing_cols

# ============================================
# STEP 2: Simple Differential Expression (Baseline)
# ============================================

def simple_differential_expression(adata, n_proteins=None):
    """
    Perform simple t-test based DE for comparison with covariate-controlled
    """

    print("\n" + "="*60)
    print("SIMPLE DIFFERENTIAL EXPRESSION (Baseline)")
    print("="*60)

    tau_pos = adata[adata.obs['tau_status'] == 'positive']
    tau_neg = adata[adata.obs['tau_status'] == 'negative']

    if n_proteins is None:
        n_proteins = adata.n_vars

    results = []

    print(f"Analyzing {n_proteins} proteins...")
    for i, gene in enumerate(adata.var_names[:n_proteins]):
        if i % 1000 == 0:
            print(f"  Progress: {i}/{n_proteins}")

        try:
            pos_expr = tau_pos[:, gene].X.flatten()
            neg_expr = tau_neg[:, gene].X.flatten()

            # Calculate statistics
            t_stat, p_val = stats.ttest_ind(pos_expr, neg_expr)
            log2fc = np.mean(pos_expr) - np.mean(neg_expr)

            results.append({
                'protein': gene,
                'log2FC_simple': log2fc,
                'pvalue_simple': p_val,
                't_statistic': t_stat
            })
        except:
            results.append({
                'protein': gene,
                'log2FC_simple': 0,
                'pvalue_simple': 1,
                't_statistic': 0
            })

    results_df = pd.DataFrame(results)

    # Apply FDR correction
    rejected, pvals_corrected, _, _ = multipletests(
        results_df['pvalue_simple'],
        method='fdr_bh',
        alpha=0.05
    )
    results_df['FDR_simple'] = pvals_corrected
    results_df['significant_simple'] = rejected

    n_sig = sum(results_df['significant_simple'])
    pct_sig = 100 * n_sig / len(results_df)

    print(f"\nSimple DE Results:")
    print(f"  Significant proteins (FDR < 0.05): {n_sig}/{len(results_df)} ({pct_sig:.2f}%)")

    return results_df

# ============================================
# STEP 3: Covariate-Controlled DE (Main Analysis)
# ============================================

def covariate_controlled_de(adata, available_covariates, n_proteins=None):
    """
    Perform differential expression with covariate control using linear models

    This is the key analysis that accounts for confounding variables
    """

    print("\n" + "="*60)
    print("COVARIATE-CONTROLLED DIFFERENTIAL EXPRESSION")
    print("="*60)

    if n_proteins is None:
        n_proteins = adata.n_vars

    results = []

    # Prepare the data frame for modeling
    obs_df = adata.obs.copy()

    # Ensure tau_status is categorical
    obs_df['tau_status'] = pd.Categorical(obs_df['tau_status'])

    # Build the formula based on available covariates
    formula_parts = ['expression ~ tau_status']

    if 'age' in available_covariates:
        formula_parts.append('age')
    if 'PMI' in available_covariates:
        formula_parts.append('PMI')
    if 'PatientID' in available_covariates:
        formula_parts.append('C(PatientID)')  # Categorical

    formula = ' + '.join(formula_parts)
    print(f"Model formula: {formula}")

    print(f"\nAnalyzing {n_proteins} proteins with covariates...")

    # Progress bar for long computation
    from tqdm import tqdm

    for gene in tqdm(adata.var_names[:n_proteins], desc="Processing proteins"):
        try:
            # Extract expression for this gene
            expr = adata[:, gene].X.flatten()
            obs_df['expression'] = expr

            # Fit linear model
            model = ols(formula, data=obs_df)
            results_fit = model.fit()

            # Extract tau effect
            # Look for the tau_status coefficient
            tau_coef = None
            tau_pval = None

            for param_name in results_fit.params.index:
                if 'tau_status' in param_name and 'positive' in param_name:
                    tau_coef = results_fit.params[param_name]
                    tau_pval = results_fit.pvalues[param_name]
                    break

            if tau_coef is None:
                # Try alternative parameterization
                for param_name in results_fit.params.index:
                    if 'tau_status' in param_name:
                        tau_coef = results_fit.params[param_name]
                        tau_pval = results_fit.pvalues[param_name]
                        break

            # Also calculate simple fold change for comparison
            tau_pos = adata[adata.obs['tau_status'] == 'positive'][:, gene].X.mean()
            tau_neg = adata[adata.obs['tau_status'] == 'negative'][:, gene].X.mean()
            simple_fc = tau_pos - tau_neg

            results.append({
                'protein': gene,
                'log2FC_adjusted': tau_coef if tau_coef is not None else 0,
                'pvalue_adjusted': tau_pval if tau_pval is not None else 1,
                'log2FC_simple': simple_fc,
                'r_squared': results_fit.rsquared,
                'aic': results_fit.aic
            })

        except Exception as e:
            # Handle any errors in model fitting
            results.append({
                'protein': gene,
                'log2FC_adjusted': 0,
                'pvalue_adjusted': 1,
                'log2FC_simple': 0,
                'r_squared': 0,
                'aic': np.inf
            })

    results_df = pd.DataFrame(results)

    # Apply FDR correction
    rejected, pvals_corrected, _, _ = multipletests(
        results_df['pvalue_adjusted'],
        method='fdr_bh',
        alpha=0.05
    )
    results_df['FDR_adjusted'] = pvals_corrected
    results_df['significant_adjusted'] = rejected

    n_sig = sum(results_df['significant_adjusted'])
    pct_sig = 100 * n_sig / len(results_df)

    print(f"\nCovariate-Controlled DE Results:")
    print(f"  Significant proteins (FDR < 0.05): {n_sig}/{len(results_df)} ({pct_sig:.2f}%)")
    print(f"  Expected: 2,115/5,853 (36.14%)")

    # Show top hits
    top_up = results_df.nlargest(5, 'log2FC_adjusted')
    top_down = results_df.nsmallest(5, 'log2FC_adjusted')

    print(f"\nTop 5 upregulated proteins:")
    print(top_up[['protein', 'log2FC_adjusted', 'FDR_adjusted']])

    print(f"\nTop 5 downregulated proteins:")
    print(top_down[['protein', 'log2FC_adjusted', 'FDR_adjusted']])

    return results_df

# ============================================
# STEP 4: Compare Simple vs Adjusted Results
# ============================================

def compare_de_methods(simple_df, adjusted_df):
    """
    Compare results from simple and covariate-adjusted analyses
    """

    print("\n" + "="*60)
    print("COMPARISON: SIMPLE vs COVARIATE-ADJUSTED")
    print("="*60)

    # Merge results
    merged = pd.merge(simple_df, adjusted_df, on='protein', suffixes=('_simple', '_adj'))

    # Compare number of significant proteins
    n_sig_simple = sum(simple_df['significant_simple'])
    n_sig_adjusted = sum(adjusted_df['significant_adjusted'])

    print(f"Significant proteins (FDR < 0.05):")
    print(f"  Simple t-test: {n_sig_simple}")
    print(f"  Covariate-adjusted: {n_sig_adjusted}")
    print(f"  Difference: {n_sig_adjusted - n_sig_simple}")

    # Correlation between fold changes
    fc_corr = np.corrcoef(merged['log2FC_simple_simple'], merged['log2FC_adjusted'])[0, 1]
    print(f"\nCorrelation between fold changes: {fc_corr:.3f}")

    # Proteins significant in one but not the other
    sig_simple_only = merged[merged['significant_simple'] & ~merged['significant_adjusted']]
    sig_adjusted_only = merged[~merged['significant_simple'] & merged['significant_adjusted']]

    print(f"\nDiscordant results:")
    print(f"  Significant in simple only: {len(sig_simple_only)}")
    print(f"  Significant in adjusted only: {len(sig_adjusted_only)}")

    # Effect of adjustment on fold changes
    fc_diff = merged['log2FC_adjusted'] - merged['log2FC_simple_simple']
    print(f"\nEffect of covariate adjustment on fold changes:")
    print(f"  Mean difference: {fc_diff.mean():.3f}")
    print(f"  Std difference: {fc_diff.std():.3f}")

    return merged

# ============================================
# STEP 5: Visualization
# ============================================

def visualize_de_results(simple_df, adjusted_df, merged_df):
    """
    Create comprehensive visualizations of DE results
    """

    fig = plt.figure(figsize=(16, 12))

    # 1. Volcano plot - Simple DE
    ax1 = plt.subplot(3, 3, 1)
    x = simple_df['log2FC_simple']
    y = -np.log10(simple_df['pvalue_simple'])
    colors = ['red' if sig else 'gray' for sig in simple_df['significant_simple']]

    ax1.scatter(x, y, c=colors, alpha=0.5, s=10)
    ax1.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.3)
    ax1.axvline(x=0, color='black', linestyle='-', alpha=0.3)
    ax1.set_xlabel('Log2 Fold Change')
    ax1.set_ylabel('-Log10(p-value)')
    ax1.set_title(f'Simple DE (n={sum(simple_df["significant_simple"])} sig)')

    # 2. Volcano plot - Adjusted DE
    ax2 = plt.subplot(3, 3, 2)
    x = adjusted_df['log2FC_adjusted']
    y = -np.log10(adjusted_df['pvalue_adjusted'])
    colors = ['red' if sig else 'gray' for sig in adjusted_df['significant_adjusted']]

    ax2.scatter(x, y, c=colors, alpha=0.5, s=10)
    ax2.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.3)
    ax2.axvline(x=0, color='black', linestyle='-', alpha=0.3)
    ax2.set_xlabel('Log2 Fold Change')
    ax2.set_ylabel('-Log10(p-value)')
    ax2.set_title(f'Covariate-Adjusted DE (n={sum(adjusted_df["significant_adjusted"])} sig)')

    # 3. Fold change comparison
    ax3 = plt.subplot(3, 3, 3)
    ax3.scatter(merged_df['log2FC_simple_simple'], merged_df['log2FC_adjusted'],
               alpha=0.5, s=10)
    ax3.plot([-5, 5], [-5, 5], 'r--', alpha=0.5)  # Identity line
    ax3.set_xlabel('Simple Log2FC')
    ax3.set_ylabel('Adjusted Log2FC')
    ax3.set_title('Fold Change Comparison')
    ax3.set_xlim([-5, 5])
    ax3.set_ylim([-5, 5])

    # 4. P-value comparison
    ax4 = plt.subplot(3, 3, 4)
    ax4.scatter(-np.log10(merged_df['pvalue_simple']),
               -np.log10(merged_df['pvalue_adjusted']),
               alpha=0.5, s=10)
    ax4.plot([0, 10], [0, 10], 'r--', alpha=0.5)
    ax4.set_xlabel('Simple -Log10(p)')
    ax4.set_ylabel('Adjusted -Log10(p)')
    ax4.set_title('P-value Comparison')

    # 5. MA plot - Adjusted results
    ax5 = plt.subplot(3, 3, 5)
    A = (adjusted_df['log2FC_adjusted'] + adjusted_df['log2FC_simple']) / 2  # Average
    M = adjusted_df['log2FC_adjusted']  # Difference
    colors = ['red' if sig else 'gray' for sig in adjusted_df['significant_adjusted']]

    ax5.scatter(A, M, c=colors, alpha=0.5, s=10)
    ax5.axhline(y=0, color='black', linestyle='-', alpha=0.3)
    ax5.set_xlabel('Average Expression')
    ax5.set_ylabel('Log2 Fold Change')
    ax5.set_title('MA Plot (Adjusted)')

    # 6. FDR distribution
    ax6 = plt.subplot(3, 3, 6)
    ax6.hist(adjusted_df['FDR_adjusted'], bins=50, edgecolor='black', alpha=0.7)
    ax6.axvline(x=0.05, color='red', linestyle='--', label='FDR = 0.05')
    ax6.set_xlabel('FDR')
    ax6.set_ylabel('Count')
    ax6.set_title('FDR Distribution')
    ax6.legend()

    # 7. Effect size distribution
    ax7 = plt.subplot(3, 3, 7)
    sig_fc = adjusted_df[adjusted_df['significant_adjusted']]['log2FC_adjusted']
    nonsig_fc = adjusted_df[~adjusted_df['significant_adjusted']]['log2FC_adjusted']

    ax7.hist([nonsig_fc, sig_fc], bins=30, label=['Non-sig', 'Significant'],
            color=['gray', 'red'], alpha=0.7, edgecolor='black')
    ax7.set_xlabel('Log2 Fold Change')
    ax7.set_ylabel('Count')
    ax7.set_title('Fold Change Distribution')
    ax7.legend()

    # 8. R-squared distribution (model fit)
    ax8 = plt.subplot(3, 3, 8)
    if 'r_squared' in adjusted_df.columns:
        ax8.hist(adjusted_df['r_squared'], bins=50, edgecolor='black', alpha=0.7)
        ax8.set_xlabel('R²')
        ax8.set_ylabel('Count')
        ax8.set_title('Model Fit (R²) Distribution')

    # 9. Summary statistics
    ax9 = plt.subplot(3, 3, 9)
    ax9.axis('off')

    summary_text = f"""
    SUMMARY STATISTICS

    Simple DE:
    - Significant: {sum(simple_df['significant_simple'])}
    - Percentage: {100*sum(simple_df['significant_simple'])/len(simple_df):.2f}%

    Covariate-Adjusted DE:
    - Significant: {sum(adjusted_df['significant_adjusted'])}
    - Percentage: {100*sum(adjusted_df['significant_adjusted'])/len(adjusted_df):.2f}%

    Expected:
    - Significant: 2,115
    - Percentage: 36.14%

    Difference from Expected:
    - Count: {sum(adjusted_df['significant_adjusted']) - 2115}
    - Percentage: {100*sum(adjusted_df['significant_adjusted'])/len(adjusted_df) - 36.14:.2f}%
    """

    ax9.text(0.1, 0.5, summary_text, fontsize=10, verticalalignment='center')

    plt.tight_layout()
    plt.show()

# ============================================
# STEP 6: Evaluation
# ============================================

def evaluate_statement(adjusted_df):
    """
    Evaluate whether the statement is supported
    """

    print("\n" + "="*60)
    print("STATEMENT EVALUATION")
    print("="*60)

    # Expected values from statement
    expected_total = 5853
    expected_significant = 2115
    expected_percentage = 36.14

    # Observed values
    observed_total = len(adjusted_df)
    observed_significant = sum(adjusted_df['significant_adjusted'])
    observed_percentage = 100 * observed_significant / observed_total

    print(f"Expected vs Observed:")
    print(f"  Total proteins: {expected_total} vs {observed_total}")
    print(f"  Significant: {expected_significant} vs {observed_significant}")
    print(f"  Percentage: {expected_percentage:.2f}% vs {observed_percentage:.2f}%")

    # Determine if the statement is supported
    # Allow for some tolerance in the numbers
    percentage_diff = abs(observed_percentage - expected_percentage)
    count_diff = abs(observed_significant - expected_significant)

    if percentage_diff < 5 and count_diff < 200:
        evaluation = "SUPPORTED"
        explanation = f"Found {observed_significant}/{observed_total} ({observed_percentage:.2f}%) significant proteins, closely matching expected 36.14%"
    elif percentage_diff < 10:
        evaluation = "PARTIALLY SUPPORTED"
        explanation = f"Found {observed_percentage:.2f}% significant proteins, somewhat different from expected 36.14%"
    else:
        evaluation = "REFUTED"
        explanation = f"Found {observed_percentage:.2f}% significant proteins, substantially different from expected 36.14%"

    print(f"\nEVALUATION: {evaluation}")
    print(f"EXPLANATION: {explanation}")

    return evaluation, explanation

# ============================================
# MAIN PIPELINE
# ============================================

def main():
    """
    Complete pipeline for covariate-controlled differential expression
    """

    # Load data
    print("Loading data...")
    adata = sc.read_h5ad('data/pool_processed_v2.h5ad')

    # Step 1: Prepare data
    available_covariates, missing_covariates = prepare_data_for_de(adata)

    # For testing, use subset of proteins (remove this for full analysis)
    n_proteins = min(5853, adata.n_vars)  # Use all available proteins up to 5853

    # Step 2: Simple DE (baseline)
    simple_df = simple_differential_expression(adata, n_proteins)

    # Step 3: Covariate-controlled DE
    adjusted_df = covariate_controlled_de(adata, available_covariates, n_proteins)

    # Step 4: Compare methods
    merged_df = compare_de_methods(simple_df, adjusted_df)

    # Step 5: Visualize
    visualize_de_results(simple_df, adjusted_df, merged_df)

    # Step 6: Evaluate
    evaluation, explanation = evaluate_statement(adjusted_df)

    # Save results
    adjusted_df.to_csv('covariate_adjusted_de_results.csv', index=False)
    print(f"\nResults saved to covariate_adjusted_de_results.csv")

    return evaluation, explanation

# ============================================
# KEY CONCEPTS AND NOTES
# ============================================

"""
CRITICAL CONCEPTS:

1. Why Covariate Control Matters:
   - Age affects protein expression
   - PMI (post-mortem interval) causes degradation
   - PatientID accounts for individual differences
   - Without control, these confound the tau effect

2. Linear Model Interpretation:
   - Coefficient β₁ = pure tau effect
   - Removes variance explained by covariates
   - More accurate than simple comparison

3. FDR vs P-value:
   - 5,853 tests = massive multiple testing problem
   - FDR controls expected false discovery rate
   - BH-FDR is standard for genomics/proteomics

4. ISLP Connections:
   - Chapter 3: Linear Regression (covariate adjustment)
   - Chapter 13: Multiple Testing
   - Chapter 6: Linear Model Selection (which covariates to include)

5. Expected Result (36.14%):
   - High percentage indicates major proteomic changes
   - Consistent with neurodegeneration
   - Validates the biological relevance

TROUBLESHOOTING:
- Missing covariates: May need to impute or exclude
- Model convergence: Check for multicollinearity
- Memory issues: Process proteins in batches
- Outliers: Consider robust regression methods
"""

if __name__ == "__main__":
    evaluation, explanation = main()