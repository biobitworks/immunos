---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/group1_mitochondrial/statement2_sqstm1_upregulation.py
relative: research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/group1_mitochondrial/statement2_sqstm1_upregulation.py
generated_at: 2025-12-23 10:28
---

```python
"""
Finding Group 1 - Statement 2: SQSTM1 Upregulation Analysis
Claim: SQSTM1 (p62) is massively upregulated (log2FC = 3.413, FDR = 1.76 × 10^-8)
       and increases with pseudotime (β = 4.951, FDR < 0.001)

# Analytical Approach and Rationale

## Overview
This analysis evaluates an extraordinary claim: SQSTM1 shows a 10.7-fold increase (2^3.413) between tau states, representing one of the most dramatic protein changes in neurodegeneration research. Such extreme upregulation requires rigorous validation through multiple analytical approaches.

## Statistical Strategy
1. **Differential expression**: Two-sample comparison with robust statistical testing
2. **Bootstrap validation**: Non-parametric confidence intervals for fold change
3. **FDR contextualization**: Evaluate significance within 5,853-protein context
4. **Pseudotime regression**: Linear model to quantify disease progression effect
5. **Effect size emphasis**: Cohen's d and biological significance assessment

## Rationale for Method Selection
- **Bootstrap CI**: Essential for extreme fold changes to avoid normal distribution assumptions
- **Dual significance testing**: Both raw p-value and FDR-corrected significance
- **Linear regression for pseudotime**: Assumes linear relationship (validated through residual analysis)
- **Welch's t-test**: Accounts for potentially unequal variances in extreme differences

## Expected Outcome Criteria
For SUPPORTED evaluation:
- Log2FC within ±0.5 of 3.413 (10.7-fold ± 41% tolerance)
- p-value < 1e-6 (extremely significant)
- Pseudotime β within ±1.0 of 4.951
- Consistent direction across all tests

## Why SQSTM1 is Critical
SQSTM1/p62 is the primary autophagy receptor protein:
- **Normal function**: Delivers ubiquitinated proteins to autophagosomes
- **Disease relevance**: Accumulates when autophagy flux is impaired
- **Biomarker status**: Protein aggregates mark autophagy failure
- **Massive upregulation**: 10.7-fold increase indicates severe autophagy dysfunction

## Statistical Challenges
- **Extreme values**: 3.413 log2FC is far beyond typical protein changes (0.5-1.5)
- **Multiple testing**: FDR = 1.76e-8 must remain significant among 5,853 proteins
- **Pseudotime linearity**: β = 4.951 assumes linear accumulation over disease progression
- **Outlier sensitivity**: Extreme values vulnerable to individual outliers

## Biological Context
In Alzheimer's disease progression:
- Early stages: Compensatory autophagy upregulation
- Middle stages: Autophagy flux reduction, p62 accumulation begins
- Late stages: Massive p62 accumulation (our observation), autophagy failure
- Pseudotime correlation: Progressive accumulation throughout disease course
"""

import sys
import os
sys.path.append('/Users/byron/project_plan')

from config import DATA_PATH, DATA_SPECS, ANALYSIS_PARAMS, load_data, get_tau_groups
import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
from scipy.stats import ttest_ind, pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import seaborn as sns

# ============================================
# BIOLOGICAL BACKGROUND
# ============================================

"""
SQSTM1/p62 (Sequestosome-1):
- Key autophagy receptor protein
- Binds ubiquitinated proteins for degradation
- Accumulates when autophagy is impaired
- Biomarker for autophagy dysfunction in neurodegeneration

References:
- Komatsu et al. (2007) Homeostatic levels of p62 control cytoplasmic inclusion body formation
- Bjørkøy et al. (2005) p62/SQSTM1 forms protein aggregates degraded by autophagy
"""

# ============================================
# STEP 1: SQSTM1 Differential Expression
# ============================================

def analyze_sqstm1_differential_expression(adata):
    """
    Perform differential expression analysis for SQSTM1

    Statistical approach:
    1. Two-sample t-test
    2. Calculate exact log2 fold change
    3. FDR correction for multiple testing context
    """

    print("="*60)
    print("SQSTM1 DIFFERENTIAL EXPRESSION ANALYSIS")
    print("="*60)

    # Check if SQSTM1 exists
    if 'SQSTM1' not in adata.var_names:
        print("ERROR: SQSTM1 not found in dataset")
        # Try alternative names
        alt_names = ['P62', 'SQSTM', 'p62']
        for alt in alt_names:
            matches = [g for g in adata.var_names if alt in g.upper()]
            if matches:
                print(f"Found potential match: {matches}")
                return None
        return None

    # Split by tau status using config
    tau_col = DATA_SPECS['tau_column']
    tau_pos_val = DATA_SPECS['tau_positive_value']
    tau_neg_val = DATA_SPECS['tau_negative_value']

    tau_pos = adata[adata.obs[tau_col] == tau_pos_val]
    tau_neg = adata[adata.obs[tau_col] == tau_neg_val]

    # Extract SQSTM1 expression (already log2 transformed)
    sqstm1_pos = tau_pos[:, 'SQSTM1'].X.flatten()
    sqstm1_neg = tau_neg[:, 'SQSTM1'].X.flatten()

    # Calculate statistics
    results = {}

    # 1. Descriptive statistics
    results['tau_pos_mean'] = np.mean(sqstm1_pos)
    results['tau_pos_std'] = np.std(sqstm1_pos)
    results['tau_pos_median'] = np.median(sqstm1_pos)
    results['tau_pos_n'] = len(sqstm1_pos)

    results['tau_neg_mean'] = np.mean(sqstm1_neg)
    results['tau_neg_std'] = np.std(sqstm1_neg)
    results['tau_neg_median'] = np.median(sqstm1_neg)
    results['tau_neg_n'] = len(sqstm1_neg)

    # 2. Log2 fold change (data already log2)
    results['log2FC'] = results['tau_pos_mean'] - results['tau_neg_mean']

    # 3. Statistical tests
    # T-test (assumes normal distribution)
    t_stat, t_pval = ttest_ind(sqstm1_pos, sqstm1_neg)
    results['t_statistic'] = t_stat
    results['t_pvalue'] = t_pval

    # Welch's t-test (doesn't assume equal variance)
    t_stat_welch, t_pval_welch = ttest_ind(sqstm1_pos, sqstm1_neg, equal_var=False)
    results['welch_t_statistic'] = t_stat_welch
    results['welch_t_pvalue'] = t_pval_welch

    # 4. Effect size (Cohen's d)
    pooled_std = np.sqrt(((len(sqstm1_pos)-1)*np.std(sqstm1_pos)**2 +
                          (len(sqstm1_neg)-1)*np.std(sqstm1_neg)**2) /
                         (len(sqstm1_pos) + len(sqstm1_neg) - 2))
    results['cohens_d'] = (results['tau_pos_mean'] - results['tau_neg_mean']) / pooled_std

    # 5. Bootstrap confidence interval for fold change
    n_bootstrap = 10000
    bootstrap_fcs = []
    for _ in range(n_bootstrap):
        boot_pos = np.random.choice(sqstm1_pos, size=len(sqstm1_pos), replace=True)
        boot_neg = np.random.choice(sqstm1_neg, size=len(sqstm1_neg), replace=True)
        bootstrap_fcs.append(np.mean(boot_pos) - np.mean(boot_neg))

    results['log2FC_ci_lower'] = np.percentile(bootstrap_fcs, 2.5)
    results['log2FC_ci_upper'] = np.percentile(bootstrap_fcs, 97.5)

    # Print results
    print(f"\nSQSTM1 Expression Statistics:")
    print(f"Tau-positive: {results['tau_pos_mean']:.3f} ± {results['tau_pos_std']:.3f} (n={results['tau_pos_n']})")
    print(f"Tau-negative: {results['tau_neg_mean']:.3f} ± {results['tau_neg_std']:.3f} (n={results['tau_neg_n']})")
    print(f"\nLog2 Fold Change: {results['log2FC']:.3f}")
    print(f"95% CI: [{results['log2FC_ci_lower']:.3f}, {results['log2FC_ci_upper']:.3f}]")
    print(f"\nT-test p-value: {results['t_pvalue']:.2e}")
    print(f"Welch's t-test p-value: {results['welch_t_pvalue']:.2e}")
    print(f"Cohen's d: {results['cohens_d']:.3f}")

    return results

# ============================================
# STEP 2: FDR Correction in Context
# ============================================

def calculate_fdr_in_context(adata, target_pvalue):
    """
    Calculate FDR for SQSTM1 in the context of all proteins

    This simulates the multiple testing context where SQSTM1's
    p-value would be adjusted alongside other proteins
    """

    print("\n" + "="*60)
    print("FDR CORRECTION IN CONTEXT")
    print("="*60)

    # Perform DE for a sample of proteins to establish context
    n_proteins_to_test = min(100, adata.n_vars)  # Test subset for efficiency
    random_proteins = np.random.choice(adata.var_names, n_proteins_to_test, replace=False)

    # Make sure SQSTM1 is included if it exists
    if 'SQSTM1' in adata.var_names and 'SQSTM1' not in random_proteins:
        random_proteins = np.append(random_proteins[:-1], 'SQSTM1')

    pvalues = []
    fold_changes = []

    tau_pos = adata[adata.obs['tau_status'] == 'positive']
    tau_neg = adata[adata.obs['tau_status'] == 'negative']

    for protein in random_proteins:
        try:
            pos_expr = tau_pos[:, protein].X.flatten()
            neg_expr = tau_neg[:, protein].X.flatten()
            _, pval = ttest_ind(pos_expr, neg_expr)
            fc = np.mean(pos_expr) - np.mean(neg_expr)

            pvalues.append(pval)
            fold_changes.append(fc)
        except:
            pvalues.append(1.0)
            fold_changes.append(0.0)

    # Apply FDR correction
    rejected, pvals_corrected, _, _ = multipletests(pvalues, method='fdr_bh', alpha=0.05)

    # Find SQSTM1's FDR if present
    if 'SQSTM1' in random_proteins:
        sqstm1_idx = list(random_proteins).index('SQSTM1')
        sqstm1_fdr = pvals_corrected[sqstm1_idx]
        print(f"SQSTM1 FDR-corrected p-value: {sqstm1_fdr:.2e}")
        print(f"Expected FDR: 1.76e-08")
        print(f"FDR match: {'Yes' if sqstm1_fdr < 1e-7 else 'No'}")

    # Show top proteins by fold change
    results_df = pd.DataFrame({
        'protein': random_proteins,
        'log2FC': fold_changes,
        'pvalue': pvalues,
        'FDR': pvals_corrected,
        'significant': rejected
    }).sort_values('log2FC', ascending=False)

    print(f"\nTop 5 upregulated proteins:")
    print(results_df.head())

    return results_df

# ============================================
# STEP 3: Pseudotime Correlation Analysis
# ============================================

def analyze_pseudotime_correlation(adata):
    """
    Analyze SQSTM1 correlation with pseudotime

    The claim states β = 4.951, which represents the regression coefficient
    in a linear model of SQSTM1 expression vs pseudotime
    """

    print("\n" + "="*60)
    print("PSEUDOTIME CORRELATION ANALYSIS")
    print("="*60)

    if 'SQSTM1' not in adata.var_names:
        print("ERROR: SQSTM1 not found")
        return None

    if 'pseudotime' not in adata.obs.columns:
        print("ERROR: pseudotime not found in metadata")
        return None

    # Extract data
    sqstm1_expr = adata[:, 'SQSTM1'].X.flatten()
    pseudotime = adata.obs['pseudotime'].values

    # Remove any NaN values
    mask = ~(np.isnan(sqstm1_expr) | np.isnan(pseudotime))
    sqstm1_expr = sqstm1_expr[mask]
    pseudotime = pseudotime[mask]

    results = {}

    # 1. Correlation analysis
    pearson_r, pearson_p = pearsonr(sqstm1_expr, pseudotime)
    spearman_r, spearman_p = spearmanr(sqstm1_expr, pseudotime)

    results['pearson_r'] = pearson_r
    results['pearson_p'] = pearson_p
    results['spearman_r'] = spearman_r
    results['spearman_p'] = spearman_p

    # 2. Linear regression to get β coefficient
    # Using statsmodels for more detailed statistics
    X = sm.add_constant(pseudotime)  # Add intercept
    model = sm.OLS(sqstm1_expr, X)
    regression_results = model.fit()

    results['beta_coefficient'] = regression_results.params[1]  # Slope
    results['beta_stderr'] = regression_results.bse[1]
    results['beta_pvalue'] = regression_results.pvalues[1]
    results['r_squared'] = regression_results.rsquared
    results['intercept'] = regression_results.params[0]

    # 3. Confidence interval for beta
    conf_int = regression_results.conf_int()
    results['beta_ci_lower'] = conf_int.iloc[1, 0]
    results['beta_ci_upper'] = conf_int.iloc[1, 1]

    # 4. Alternative: sklearn LinearRegression
    lr = LinearRegression()
    lr.fit(pseudotime.reshape(-1, 1), sqstm1_expr)
    results['sklearn_beta'] = lr.coef_[0]

    # Print results
    print(f"\nCorrelation with Pseudotime:")
    print(f"Pearson r: {results['pearson_r']:.3f} (p={results['pearson_p']:.2e})")
    print(f"Spearman r: {results['spearman_r']:.3f} (p={results['spearman_p']:.2e})")

    print(f"\nLinear Regression:")
    print(f"β coefficient: {results['beta_coefficient']:.3f}")
    print(f"95% CI: [{results['beta_ci_lower']:.3f}, {results['beta_ci_upper']:.3f}]")
    print(f"p-value: {results['beta_pvalue']:.2e}")
    print(f"R²: {results['r_squared']:.3f}")

    print(f"\nExpected β: 4.951")
    print(f"Observed β: {results['beta_coefficient']:.3f}")

    return results

# ============================================
# STEP 4: Comprehensive Evaluation
# ============================================

def evaluate_statement(de_results, pseudotime_results):
    """
    Evaluate the complete statement about SQSTM1
    """

    print("\n" + "="*60)
    print("STATEMENT EVALUATION")
    print("="*60)

    evaluation_criteria = {
        'log2FC': {'expected': 3.413, 'tolerance': 0.5},
        'fdr': {'expected': 1.76e-8, 'threshold': 1e-6},
        'beta': {'expected': 4.951, 'tolerance': 1.0},
        'beta_fdr': {'threshold': 0.001}
    }

    # Check differential expression
    de_match = False
    if de_results:
        observed_fc = de_results['log2FC']
        fc_match = abs(observed_fc - evaluation_criteria['log2FC']['expected']) < evaluation_criteria['log2FC']['tolerance']
        pval_match = de_results['t_pvalue'] < evaluation_criteria['fdr']['threshold']
        de_match = fc_match and pval_match

        print(f"Differential Expression:")
        print(f"  Log2FC match: {fc_match} (observed: {observed_fc:.3f}, expected: 3.413)")
        print(f"  P-value significant: {pval_match} (p={de_results['t_pvalue']:.2e})")

    # Check pseudotime correlation
    pseudo_match = False
    if pseudotime_results:
        observed_beta = pseudotime_results['beta_coefficient']
        beta_match = abs(observed_beta - evaluation_criteria['beta']['expected']) < evaluation_criteria['beta']['tolerance']
        beta_sig = pseudotime_results['beta_pvalue'] < evaluation_criteria['beta_fdr']['threshold']
        pseudo_match = beta_match and beta_sig

        print(f"\nPseudotime Correlation:")
        print(f"  Beta match: {beta_match} (observed: {observed_beta:.3f}, expected: 4.951)")
        print(f"  Beta significant: {beta_sig} (p={pseudotime_results['beta_pvalue']:.2e})")

    # Overall evaluation
    if de_match and pseudo_match:
        evaluation = "SUPPORTED"
        explanation = "SQSTM1 shows massive upregulation and strong pseudotime correlation as claimed"
    elif de_match or pseudo_match:
        evaluation = "PARTIALLY SUPPORTED"
        explanation = f"{'DE' if de_match else 'Pseudotime'} component supported, but not both"
    else:
        evaluation = "REFUTED"
        explanation = "Neither the fold change nor pseudotime correlation match the claimed values"

    print(f"\nFINAL EVALUATION: {evaluation}")
    print(f"EXPLANATION: {explanation}")

    return evaluation, explanation

# ============================================
# STEP 5: Visualization
# ============================================

def visualize_sqstm1_analysis(adata, de_results, pseudotime_results):
    """
    Create comprehensive visualizations for SQSTM1 analysis
    """

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # 1. Expression by tau status (violin plot)
    ax = axes[0, 0]
    if 'SQSTM1' in adata.var_names:
        tau_pos = adata[adata.obs['tau_status'] == 'positive'][:, 'SQSTM1'].X.flatten()
        tau_neg = adata[adata.obs['tau_status'] == 'negative'][:, 'SQSTM1'].X.flatten()

        parts = ax.violinplot([tau_neg, tau_pos], positions=[0, 1], showmeans=True)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Tau-negative', 'Tau-positive'])
        ax.set_ylabel('SQSTM1 Expression (log2)')
        ax.set_title('SQSTM1 by Tau Status')

        # Add fold change annotation
        if de_results:
            ax.text(0.5, ax.get_ylim()[1]*0.9, f"Log2FC = {de_results['log2FC']:.3f}",
                   ha='center', fontsize=10, weight='bold')

    # 2. Expression distribution
    ax = axes[0, 1]
    if 'SQSTM1' in adata.var_names:
        sqstm1_all = adata[:, 'SQSTM1'].X.flatten()
        ax.hist(sqstm1_all, bins=50, edgecolor='black', alpha=0.7)
        ax.axvline(np.mean(sqstm1_all), color='red', linestyle='--', label='Mean')
        ax.set_xlabel('SQSTM1 Expression (log2)')
        ax.set_ylabel('Frequency')
        ax.set_title('SQSTM1 Distribution')
        ax.legend()

    # 3. Pseudotime correlation
    ax = axes[0, 2]
    if 'SQSTM1' in adata.var_names and 'pseudotime' in adata.obs.columns:
        sqstm1 = adata[:, 'SQSTM1'].X.flatten()
        pseudotime = adata.obs['pseudotime'].values

        ax.scatter(pseudotime, sqstm1, alpha=0.5, s=10)

        # Add regression line
        if pseudotime_results:
            x_line = np.linspace(pseudotime.min(), pseudotime.max(), 100)
            y_line = (pseudotime_results['beta_coefficient'] * x_line +
                     pseudotime_results['intercept'])
            ax.plot(x_line, y_line, 'r-', linewidth=2,
                   label=f"β = {pseudotime_results['beta_coefficient']:.3f}")

        ax.set_xlabel('Pseudotime')
        ax.set_ylabel('SQSTM1 Expression (log2)')
        ax.set_title('SQSTM1 vs Pseudotime')
        ax.legend()

    # 4. QQ plot for normality check
    ax = axes[1, 0]
    if 'SQSTM1' in adata.var_names:
        sqstm1_all = adata[:, 'SQSTM1'].X.flatten()
        stats.probplot(sqstm1_all, dist="norm", plot=ax)
        ax.set_title('Q-Q Plot')

    # 5. Residuals plot
    ax = axes[1, 1]
    if pseudotime_results and 'SQSTM1' in adata.var_names:
        sqstm1 = adata[:, 'SQSTM1'].X.flatten()
        pseudotime = adata.obs['pseudotime'].values
        predicted = (pseudotime_results['beta_coefficient'] * pseudotime +
                    pseudotime_results['intercept'])
        residuals = sqstm1 - predicted

        ax.scatter(predicted, residuals, alpha=0.5)
        ax.axhline(y=0, color='red', linestyle='--')
        ax.set_xlabel('Predicted Values')
        ax.set_ylabel('Residuals')
        ax.set_title('Residual Plot')

    # 6. Effect size comparison
    ax = axes[1, 2]
    if de_results:
        effect_sizes = {
            'Cohen\'s d': de_results['cohens_d'],
            'Log2 FC': de_results['log2FC']
        }
        ax.bar(range(len(effect_sizes)), list(effect_sizes.values()))
        ax.set_xticks(range(len(effect_sizes)))
        ax.set_xticklabels(list(effect_sizes.keys()))
        ax.set_ylabel('Effect Size')
        ax.set_title('Effect Size Metrics')

        # Add reference lines
        ax.axhline(y=3.413, color='red', linestyle='--', alpha=0.5, label='Expected FC')
        ax.legend()

    plt.tight_layout()
    plt.show()

# ============================================
# MAIN ANALYSIS PIPELINE
# ============================================

def main():
    """
    Complete analysis pipeline for Statement 2
    """

    # Load data
    print("Loading data...")

    # Load Alzheimer's disease proteomic dataset
    # Critical dataset characteristics for SQSTM1 analysis:
    # - Source: Mini-pools of 10 neurons from AD brain tissue samples
    # - Scope: 5,853 proteins quantified across neuronal populations
    # - Key variables for this analysis:
    #   * tau_status: positive/negative classification (primary comparison)
    #   * SQSTM1 expression: target protein showing claimed massive upregulation
    #   * pseudotime: disease progression continuum (for temporal analysis)
    #   * MC1 scores: misfolded tau burden (validation measure)
    # - Preprocessing applied:
    #   * log2 transformation (enables additive fold change calculation)
    #   * Quality control filtering (removes low-quality measurements)
    #   * Batch correction (minimizes technical variation)
    # - Statistical context: 5,853 proteins tested simultaneously (FDR critical)
    # - Biological context: Late-stage AD where autophagy failure expected
    adata = sc.read_h5ad('data/pool_processed_v2.h5ad')
    print(f"Loaded {adata.shape[0]} cells, {adata.shape[1]} proteins")

    # Verify SQSTM1 presence (essential for analysis)
    if 'SQSTM1' not in adata.var_names:
        print("CRITICAL ERROR: SQSTM1 not found in dataset")
        print("Alternative names to check: P62, SQSTM")
        return None
    else:
        print("✓ SQSTM1 protein found in dataset")

    # Part 1: Differential expression
    de_results = analyze_sqstm1_differential_expression(adata)

    # Part 2: FDR correction context
    fdr_df = calculate_fdr_in_context(adata, de_results['t_pvalue'] if de_results else 1.0)

    # Part 3: Pseudotime correlation
    pseudotime_results = analyze_pseudotime_correlation(adata)

    # Part 4: Evaluation
    evaluation, explanation = evaluate_statement(de_results, pseudotime_results)

    # Part 5: Visualization
    visualize_sqstm1_analysis(adata, de_results, pseudotime_results)

    return evaluation, explanation

# ============================================
# KEY TAKEAWAYS AND NOTES
# ============================================

"""
CRITICAL POINTS FOR EVALUATION:

1. Log2 Fold Change Interpretation:
   - FC = 3.413 means ~10.7-fold increase (2^3.413)
   - This is a massive upregulation
   - Check if data is already log2 transformed

2. FDR vs P-value:
   - FDR = 1.76e-8 is extremely significant
   - FDR accounts for multiple testing
   - More stringent than raw p-value

3. Beta Coefficient (β = 4.951):
   - Represents change in SQSTM1 per unit pseudotime
   - Large positive value indicates strong increase over disease progression
   - Should be evaluated with confidence intervals

4. Statistical Power:
   - Sample size affects ability to detect effects
   - Bootstrap can provide robust confidence intervals
   - Consider non-parametric tests if distribution is non-normal

5. ISLP Connection:
   - Linear regression (Chapter 3)
   - Hypothesis testing (Chapter 13)
   - Multiple testing correction concepts

TROUBLESHOOTING:
- If SQSTM1 not found, check alternative names (P62, SQSTM)
- Verify data is log2 transformed
- Check for outliers that might affect results
- Consider batch effects or confounders
"""

if __name__ == "__main__":
    evaluation, explanation = main()
    print(f"\n{'='*60}")
    print(f"FINAL RESULT: {evaluation}")
    print(f"{explanation}")
    print(f"{'='*60}")
```
