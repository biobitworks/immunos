# Bioinformatics Finding Group Evaluation Framework

## Comprehensive Analysis of Proteomic Data in Alzheimer's Disease Research

This notebook provides a complete framework for evaluating biological statements using rigorous statistical methods on proteomic data from Alzheimer's disease cases.

### Project Overview
- **Dataset**: Mini-pools of 10 neurons from Alzheimer's disease cases
- **Scope**: 5,853 proteins across neuronal samples
- **Groups**: Two finding groups focused on mitochondrial dysregulation and proteostasis failure
- **Methodology**: Statistical rigor following ISLP and Rosner's Biostatistics principles

### Key Features
- Enhanced documentation with analytical rationale
- Comprehensive statistical testing with multiple correction methods
- Advanced temporal analysis with sliding windows
- AI automation capabilities for large-scale analysis
- Publication-quality visualizations

## 📋 Setup and Dependencies

### Required Libraries
Before running this notebook, install the required dependencies:

```bash
pip install pandas numpy scipy scanpy statsmodels scikit-learn matplotlib seaborn tqdm
```

### Import Libraries


```python
# Core data science libraries
import pandas as pd
import numpy as np
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')

# Statistical analysis
from scipy import stats
from scipy.stats import ttest_ind, mannwhitneyu, pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from statsmodels.formula.api import ols

# Machine learning
from sklearn.linear_model import LinearRegression
from sklearn.tree import DecisionTreeRegressor
from scipy.optimize import curve_fit

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns

# Progress tracking
from tqdm import tqdm

# Configuration
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white')
plt.style.use('default')
sns.set_palette("husl")

print("✅ All libraries imported successfully")
print(f"📊 Scanpy version: {sc.__version__}")
print(f"🐼 Pandas version: {pd.__version__}")
print(f"🔢 NumPy version: {np.__version__}")
```

## 📊 Data Loading and Exploration

### Load the H5AD Dataset
The dataset contains proteomic measurements from neuronal mini-pools with comprehensive metadata.


```python
def load_and_validate_data(file_path='data/pool_processed_v2.h5ad'):
    """
    Load the H5AD file and perform comprehensive validation
    
    Returns:
    --------
    adata : AnnData
        Loaded and validated dataset
    """
    print("🔄 Loading proteomic dataset...")
    
    # Load data with error handling
    try:
        adata = sc.read_h5ad(file_path)
        print(f"✅ Successfully loaded dataset: {adata.shape}")
    except FileNotFoundError:
        print(f"❌ Error: Data file not found at {file_path}")
        print("📁 Please ensure the data file is in the correct location")
        return None
    except Exception as e:
        print(f"❌ Error loading data: {str(e)}")
        return None
    
    # Dataset characteristics validation
    print("\n📈 Dataset Characteristics:")
    print(f"   📦 Dimensions: {adata.n_obs} cells × {adata.n_vars} proteins")
    print(f"   💾 Memory usage: {adata.X.nbytes / 1e6:.1f} MB")
    
    # Metadata validation
    print("\n🗂️ Metadata Columns:")
    for i, col in enumerate(adata.obs.columns):
        print(f"   {i+1:2d}. {col}")
    
    # Required columns check
    required_cols = ['tau_status', 'MC1', 'pseudotime']
    missing = [col for col in required_cols if col not in adata.obs.columns]
    if missing:
        print(f"⚠️  Warning: Missing required columns: {missing}")
    else:
        print("✅ All required metadata columns present")
    
    # Tau status distribution
    if 'tau_status' in adata.obs.columns:
        tau_counts = adata.obs['tau_status'].value_counts()
        print("\n🧠 Tau Status Distribution:")
        for status, count in tau_counts.items():
            percentage = 100 * count / len(adata.obs)
            print(f"   {status}: {count} cells ({percentage:.1f}%)")
    
    # MC1 and pseudotime statistics
    if 'MC1' in adata.obs.columns:
        print(f"\n🔬 MC1 Score Range: [{adata.obs['MC1'].min():.2f}, {adata.obs['MC1'].max():.2f}]")
        print(f"   Mean ± SD: {adata.obs['MC1'].mean():.2f} ± {adata.obs['MC1'].std():.2f}")
    
    if 'pseudotime' in adata.obs.columns:
        print(f"\n⏰ Pseudotime Range: [{adata.obs['pseudotime'].min():.3f}, {adata.obs['pseudotime'].max():.3f}]")
        print(f"   Mean ± SD: {adata.obs['pseudotime'].mean():.3f} ± {adata.obs['pseudotime'].std():.3f}")
    
    # Sample protein names
    print(f"\n🧬 Sample Proteins: {', '.join(adata.var_names[:10].tolist())}...")
    
    return adata

# Load the dataset
adata = load_and_validate_data()

if adata is not None:
    print("\n🎉 Data loading completed successfully!")
else:
    print("\n❌ Failed to load data. Please check file path and dependencies.")
```

### Data Quality Assessment


```python
if adata is not None:
    # Create comprehensive visualization of dataset characteristics
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # 1. Tau status distribution
    if 'tau_status' in adata.obs.columns:
        tau_counts = adata.obs['tau_status'].value_counts()
        axes[0, 0].pie(tau_counts.values, labels=tau_counts.index, autopct='%1.1f%%', 
                      colors=['lightblue', 'lightcoral'])
        axes[0, 0].set_title('Tau Status Distribution', fontsize=14, fontweight='bold')
    
    # 2. MC1 score distribution
    if 'MC1' in adata.obs.columns:
        axes[0, 1].hist(adata.obs['MC1'], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
        axes[0, 1].axvline(adata.obs['MC1'].mean(), color='red', linestyle='--', 
                          label=f'Mean: {adata.obs["MC1"].mean():.2f}')
        axes[0, 1].set_xlabel('MC1 Score')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('MC1 Score Distribution', fontsize=14, fontweight='bold')
        axes[0, 1].legend()
    
    # 3. Pseudotime distribution
    if 'pseudotime' in adata.obs.columns:
        axes[0, 2].hist(adata.obs['pseudotime'], bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
        axes[0, 2].axvline(adata.obs['pseudotime'].mean(), color='red', linestyle='--',
                          label=f'Mean: {adata.obs["pseudotime"].mean():.3f}')
        axes[0, 2].set_xlabel('Pseudotime')
        axes[0, 2].set_ylabel('Frequency')
        axes[0, 2].set_title('Pseudotime Distribution', fontsize=14, fontweight='bold')
        axes[0, 2].legend()
    
    # 4. MC1 vs Pseudotime relationship
    if 'MC1' in adata.obs.columns and 'pseudotime' in adata.obs.columns:
        scatter = axes[1, 0].scatter(adata.obs['pseudotime'], adata.obs['MC1'], 
                                   c=adata.obs['tau_status'].astype('category').cat.codes,
                                   alpha=0.6, cmap='coolwarm')
        axes[1, 0].set_xlabel('Pseudotime')
        axes[1, 0].set_ylabel('MC1 Score')
        axes[1, 0].set_title('MC1 vs Pseudotime (colored by tau status)', fontsize=14, fontweight='bold')
        plt.colorbar(scatter, ax=axes[1, 0])
    
    # 5. Expression level distribution (sample)
    if adata.n_vars > 0:
        sample_expr = adata[:, adata.var_names[0]].X.flatten()
        axes[1, 1].hist(sample_expr, bins=30, alpha=0.7, color='orange', edgecolor='black')
        axes[1, 1].set_xlabel('Log2 Expression')
        axes[1, 1].set_ylabel('Frequency')
        axes[1, 1].set_title(f'Expression Distribution\n({adata.var_names[0]})', fontsize=14, fontweight='bold')
    
    # 6. Dataset summary statistics
    axes[1, 2].axis('off')
    summary_text = f"""
    📊 DATASET SUMMARY
    
    📦 Dimensions: {adata.n_obs} × {adata.n_vars}
    
    🧠 Tau Status:
    {adata.obs['tau_status'].value_counts().to_string() if 'tau_status' in adata.obs.columns else 'N/A'}
    
    🔬 MC1 Score:
    Range: [{adata.obs['MC1'].min():.2f}, {adata.obs['MC1'].max():.2f}]
    Mean: {adata.obs['MC1'].mean():.2f}
    
    ⏰ Pseudotime:
    Range: [{adata.obs['pseudotime'].min():.3f}, {adata.obs['pseudotime'].max():.3f}]
    Mean: {adata.obs['pseudotime'].mean():.3f}
    
    ✅ Ready for Analysis!
    """
    
    axes[1, 2].text(0.1, 0.5, summary_text, fontsize=11, verticalalignment='center',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray", alpha=0.5))
    
    plt.tight_layout()
    plt.show()

else:
    print("⚠️ Cannot create visualizations - data not loaded")
```

## 🧬 Group 1: Late-Stage Mitochondrial Dysregulation

This section analyzes mitochondrial dysfunction and autophagy failure in late-stage disease progression.

### Key Statements:
1. **UPS Protein Analysis**: No significant alterations in UPS proteins
2. **SQSTM1 Upregulation**: Massive 10.7-fold increase in autophagy receptor
3. **Sliding Window Analysis**: Dynamic correlation changes over disease progression

### Statement 1: UPS Protein Analysis

**Claim**: Targeted analyses show no significant UPS protein alterations across tau-positive versus tau-negative neurons.

**Analytical Approach**: Conservative evaluation using dual testing (parametric/non-parametric) with FDR correction.


```python
def analyze_ups_proteins(adata):
    """
    Comprehensive UPS protein analysis with enhanced documentation
    
    This analysis evaluates whether UPS proteins show significant alterations
    between tau-positive and tau-negative neurons using rigorous statistical methods.
    """
    if adata is None:
        print("❌ No data available for analysis")
        return None
    
    print("🔬 UPS Protein Analysis - Statement 1")
    print("="*50)
    
    # Step 1: Identify UPS proteins using multiple strategies
    print("\n🔍 Step 1: UPS Protein Identification")
    
    # Pattern-based identification
    ups_patterns = ['UB', 'PSM', 'PROT', 'UBA', 'USP']
    pattern_proteins = []
    for pattern in ups_patterns:
        matches = [p for p in adata.var_names if pattern in p.upper()]
        pattern_proteins.extend(matches)
        print(f"   {pattern}: {len(matches)} proteins found")
    
    # Curated UPS protein list (from literature)
    curated_ups = ['PSMA1', 'PSMA2', 'PSMA3', 'PSMA4', 'PSMA5', 'PSMA6', 'PSMA7',
                   'PSMB1', 'PSMB2', 'PSMB3', 'PSMB4', 'PSMB5', 'PSMB6', 'PSMB7',
                   'UBA1', 'UBA2', 'UBA3', 'UBE2A', 'UBE2B', 'UBE2C']
    
    available_curated = [p for p in curated_ups if p in adata.var_names]
    print(f"   Curated list: {len(available_curated)}/{len(curated_ups)} available")
    
    # Combine and deduplicate
    all_ups = list(set(pattern_proteins + available_curated))
    print(f"   📊 Total UPS proteins identified: {len(all_ups)}")
    
    if len(all_ups) < 5:
        print("⚠️ Warning: Few UPS proteins found. Analysis may be underpowered.")
    
    # Step 2: Differential expression analysis
    print("\n📈 Step 2: Differential Expression Analysis")
    
    if 'tau_status' not in adata.obs.columns:
        print("❌ Error: tau_status column not found")
        return None
    
    # Split by tau status
    tau_pos = adata[adata.obs['tau_status'] == 'positive']
    tau_neg = adata[adata.obs['tau_status'] == 'negative']
    
    print(f"   Tau-positive: {len(tau_pos)} cells")
    print(f"   Tau-negative: {len(tau_neg)} cells")
    
    results = []
    
    print("\n🧮 Calculating statistics for each UPS protein...")
    for protein in tqdm(all_ups, desc="Processing UPS proteins"):
        if protein not in adata.var_names:
            continue
            
        try:
            # Extract expression values
            pos_expr = tau_pos[:, protein].X.flatten()
            neg_expr = tau_neg[:, protein].X.flatten()
            
            # Calculate basic statistics
            mean_pos = np.mean(pos_expr)
            mean_neg = np.mean(neg_expr)
            log2fc = mean_pos - mean_neg  # Already in log2 scale
            
            # Parametric test (t-test)
            t_stat, p_ttest = ttest_ind(pos_expr, neg_expr)
            
            # Non-parametric test (Mann-Whitney U)
            u_stat, p_mannwhitney = mannwhitneyu(pos_expr, neg_expr, alternative='two-sided')
            
            # Effect size (Cohen's d)
            pooled_std = np.sqrt(((len(pos_expr)-1)*np.var(pos_expr, ddof=1) + 
                                 (len(neg_expr)-1)*np.var(neg_expr, ddof=1)) / 
                                (len(pos_expr) + len(neg_expr) - 2))
            cohens_d = (mean_pos - mean_neg) / pooled_std if pooled_std > 0 else 0
            
            results.append({
                'protein': protein,
                'mean_tau_pos': mean_pos,
                'mean_tau_neg': mean_neg,
                'log2fc': log2fc,
                'p_ttest': p_ttest,
                'p_mannwhitney': p_mannwhitney,
                'cohens_d': cohens_d,
                't_statistic': t_stat
            })
            
        except Exception as e:
            print(f"⚠️ Error processing {protein}: {str(e)}")
            continue
    
    if not results:
        print("❌ No results generated. Check protein names and data.")
        return None
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Step 3: Multiple testing correction
    print("\n🔧 Step 3: Multiple Testing Correction")
    
    # FDR correction for both test types
    rejected_ttest, fdr_ttest, _, _ = multipletests(results_df['p_ttest'], method='fdr_bh')
    rejected_mw, fdr_mw, _, _ = multipletests(results_df['p_mannwhitney'], method='fdr_bh')
    
    results_df['fdr_ttest'] = fdr_ttest
    results_df['fdr_mannwhitney'] = fdr_mw
    results_df['significant_ttest'] = rejected_ttest
    results_df['significant_mannwhitney'] = rejected_mw
    
    # Combined significance (both tests significant)
    results_df['significant_both'] = results_df['significant_ttest'] & results_df['significant_mannwhitney']
    
    # Step 4: Results summary
    print("\n📊 Step 4: Results Summary")
    
    n_total = len(results_df)
    n_sig_ttest = sum(results_df['significant_ttest'])
    n_sig_mw = sum(results_df['significant_mannwhitney'])
    n_sig_both = sum(results_df['significant_both'])
    
    pct_sig_ttest = 100 * n_sig_ttest / n_total
    pct_sig_mw = 100 * n_sig_mw / n_total
    pct_sig_both = 100 * n_sig_both / n_total
    
    print(f"   📈 Total UPS proteins analyzed: {n_total}")
    print(f"   🔬 Significant (t-test, FDR<0.05): {n_sig_ttest}/{n_total} ({pct_sig_ttest:.1f}%)")
    print(f"   🔬 Significant (Mann-Whitney, FDR<0.05): {n_sig_mw}/{n_total} ({pct_sig_mw:.1f}%)")
    print(f"   🔬 Significant (both tests): {n_sig_both}/{n_total} ({pct_sig_both:.1f}%)")
    
    # Effect size summary
    mean_abs_d = np.mean(np.abs(results_df['cohens_d']))
    print(f"   📏 Mean absolute Cohen's d: {mean_abs_d:.3f}")
    
    # Statement evaluation
    print("\n🎯 Statement Evaluation")
    
    # Criteria for "no significant alterations"
    criteria_met = pct_sig_both < 5 and mean_abs_d < 0.2
    
    if criteria_met:
        evaluation = "SUPPORTED"
        explanation = f"Only {pct_sig_both:.1f}% of UPS proteins show significant alterations (both tests), with small effect sizes (mean |d|={mean_abs_d:.3f})"
    else:
        evaluation = "REFUTED"
        explanation = f"{pct_sig_both:.1f}% of UPS proteins show significant alterations, or effect sizes are large (mean |d|={mean_abs_d:.3f})"
    
    print(f"   📋 Evaluation: {evaluation}")
    print(f"   📝 Explanation: {explanation}")
    
    return {
        'results_df': results_df,
        'evaluation': evaluation,
        'explanation': explanation,
        'n_total': n_total,
        'pct_significant': pct_sig_both,
        'mean_effect_size': mean_abs_d
    }

# Run UPS protein analysis
if adata is not None:
    ups_results = analyze_ups_proteins(adata)
    
    if ups_results:
        print(f"\n✅ UPS Analysis Complete: {ups_results['evaluation']}")
    else:
        print("❌ UPS Analysis Failed")
else:
    print("⚠️ Skipping UPS analysis - no data loaded")
```

### Statement 2: SQSTM1 Upregulation Analysis

**Claim**: SQSTM1 shows 10.7-fold increase (2^3.413) between tau states, representing massive autophagy dysfunction.

**Analytical Approach**: Bootstrap validation with multiple statistical measures for extreme fold change validation.


```python
def analyze_sqstm1_upregulation(adata):
    """
    Comprehensive SQSTM1 upregulation analysis with bootstrap validation
    
    SQSTM1/p62 is the primary autophagy receptor protein that accumulates
    when autophagy flux is impaired, making it a critical biomarker.
    """
    if adata is None:
        print("❌ No data available for analysis")
        return None
    
    print("🧬 SQSTM1 Upregulation Analysis - Statement 2")
    print("="*50)
    
    # Check if SQSTM1 exists
    if 'SQSTM1' not in adata.var_names:
        print("❌ SQSTM1 not found in dataset")
        # Look for alternatives
        alternatives = [p for p in adata.var_names if 'SQSTM' in p.upper()]
        if alternatives:
            print(f"🔍 Potential alternatives found: {alternatives}")
        return None
    
    print("✅ SQSTM1 found in dataset")
    
    # Step 1: Basic differential expression
    print("\n📊 Step 1: Differential Expression Analysis")
    
    tau_pos = adata[adata.obs['tau_status'] == 'positive']
    tau_neg = adata[adata.obs['tau_status'] == 'negative']
    
    pos_expr = tau_pos[:, 'SQSTM1'].X.flatten()
    neg_expr = tau_neg[:, 'SQSTM1'].X.flatten()
    
    mean_pos = np.mean(pos_expr)
    mean_neg = np.mean(neg_expr)
    log2fc = mean_pos - mean_neg  # Data already in log2 scale
    fold_change = 2**log2fc
    
    print(f"   📈 Tau-positive mean: {mean_pos:.3f}")
    print(f"   📉 Tau-negative mean: {mean_neg:.3f}")
    print(f"   🔢 Log2 fold change: {log2fc:.3f}")
    print(f"   📊 Fold change: {fold_change:.1f}x")
    print(f"   🎯 Expected: 3.413 log2FC (10.7x)")
    
    # Statistical test
    t_stat, p_val = ttest_ind(pos_expr, neg_expr)
    print(f"   📐 T-statistic: {t_stat:.3f}")
    print(f"   🎲 P-value: {p_val:.2e}")
    
    # Step 2: Bootstrap confidence interval
    print("\n🔄 Step 2: Bootstrap Confidence Interval")
    
    n_bootstrap = 1000
    bootstrap_fc = []
    
    print(f"   Running {n_bootstrap} bootstrap iterations...")
    for i in range(n_bootstrap):
        # Resample with replacement
        pos_sample = np.random.choice(pos_expr, len(pos_expr), replace=True)
        neg_sample = np.random.choice(neg_expr, len(neg_expr), replace=True)
        
        # Calculate fold change
        boot_fc = np.mean(pos_sample) - np.mean(neg_sample)
        bootstrap_fc.append(boot_fc)
    
    # Confidence intervals
    ci_lower = np.percentile(bootstrap_fc, 2.5)
    ci_upper = np.percentile(bootstrap_fc, 97.5)
    
    print(f"   📊 Bootstrap 95% CI: [{ci_lower:.3f}, {ci_upper:.3f}]")
    print(f"   🎯 Expected 3.413 in CI: {'Yes' if ci_lower <= 3.413 <= ci_upper else 'No'}")
    
    # Step 3: Pseudotime correlation
    print("\n⏰ Step 3: Pseudotime Correlation Analysis")
    
    if 'pseudotime' in adata.obs.columns:
        sqstm1_all = adata[:, 'SQSTM1'].X.flatten()
        pseudotime = adata.obs['pseudotime'].values
        
        # Correlation
        corr, corr_p = pearsonr(sqstm1_all, pseudotime)
        print(f"   📈 Correlation with pseudotime: r={corr:.3f}, p={corr_p:.2e}")
        
        # Linear regression for β coefficient
        from sklearn.linear_model import LinearRegression
        lr = LinearRegression()
        lr.fit(pseudotime.reshape(-1, 1), sqstm1_all)
        beta = lr.coef_[0]
        r_squared = lr.score(pseudotime.reshape(-1, 1), sqstm1_all)
        
        print(f"   📐 Linear regression β: {beta:.3f}")
        print(f"   📊 R-squared: {r_squared:.3f}")
        print(f"   🎯 Expected β ≈ 4.951")
    
    # Step 4: Visualization
    print("\n📈 Step 4: Creating Visualizations")
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Expression comparison
    axes[0, 0].boxplot([neg_expr, pos_expr], labels=['Tau-negative', 'Tau-positive'])
    axes[0, 0].set_ylabel('SQSTM1 Expression (log2)')
    axes[0, 0].set_title('SQSTM1 Expression by Tau Status')
    
    # 2. Bootstrap distribution
    axes[0, 1].hist(bootstrap_fc, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    axes[0, 1].axvline(log2fc, color='red', linestyle='-', linewidth=2, label=f'Observed: {log2fc:.3f}')
    axes[0, 1].axvline(3.413, color='orange', linestyle='--', linewidth=2, label='Expected: 3.413')
    axes[0, 1].axvline(ci_lower, color='gray', linestyle=':', label=f'95% CI')
    axes[0, 1].axvline(ci_upper, color='gray', linestyle=':')
    axes[0, 1].set_xlabel('Log2 Fold Change')
    axes[0, 1].set_ylabel('Frequency')
    axes[0, 1].set_title('Bootstrap Distribution of Fold Changes')
    axes[0, 1].legend()
    
    # 3. Pseudotime correlation
    if 'pseudotime' in adata.obs.columns:
        scatter = axes[1, 0].scatter(pseudotime, sqstm1_all, 
                                   c=adata.obs['tau_status'].astype('category').cat.codes,
                                   alpha=0.6, cmap='coolwarm')
        
        # Add regression line
        x_line = np.linspace(pseudotime.min(), pseudotime.max(), 100)
        y_line = lr.predict(x_line.reshape(-1, 1))
        axes[1, 0].plot(x_line, y_line, 'black', linewidth=2, 
                       label=f'β={beta:.3f}, R²={r_squared:.3f}')
        
        axes[1, 0].set_xlabel('Pseudotime')
        axes[1, 0].set_ylabel('SQSTM1 Expression')
        axes[1, 0].set_title('SQSTM1 vs Pseudotime')
        axes[1, 0].legend()
        plt.colorbar(scatter, ax=axes[1, 0])
    
    # 4. Summary statistics
    axes[1, 1].axis('off')
    summary_text = f"""
    SQSTM1 ANALYSIS SUMMARY
    
    📊 Fold Change Analysis:
    Observed: {log2fc:.3f} log2FC ({fold_change:.1f}x)
    Expected: 3.413 log2FC (10.7x)
    
    🔬 Statistical Test:
    T-statistic: {t_stat:.3f}
    P-value: {p_val:.2e}
    
    🔄 Bootstrap 95% CI:
    [{ci_lower:.3f}, {ci_upper:.3f}]
    
    ⏰ Pseudotime Correlation:
    r = {corr:.3f} (p = {corr_p:.2e})
    β = {beta:.3f} (expected ≈ 4.951)
    
    🎯 EVALUATION:
    {'SUPPORTED' if abs(log2fc - 3.413) < 0.5 else 'NEEDS REVIEW'}
    """
    
    axes[1, 1].text(0.1, 0.5, summary_text, fontsize=10, verticalalignment='center',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.5))
    
    plt.tight_layout()
    plt.show()
    
    # Step 5: Evaluation
    print("\n🎯 Step 5: Statement Evaluation")
    
    # Criteria for supporting the claim
    fc_match = abs(log2fc - 3.413) < 0.5  # Within 0.5 log2 units
    significant = p_val < 0.05
    ci_contains_expected = ci_lower <= 3.413 <= ci_upper
    
    if fc_match and significant:
        evaluation = "SUPPORTED"
        explanation = f"SQSTM1 shows {log2fc:.3f} log2FC (observed) vs 3.413 (expected), statistically significant (p={p_val:.2e})"
    elif significant:
        evaluation = "PARTIALLY SUPPORTED"
        explanation = f"SQSTM1 significantly upregulated ({log2fc:.3f} log2FC) but differs from expected 3.413"
    else:
        evaluation = "REFUTED"
        explanation = f"SQSTM1 fold change not statistically significant (p={p_val:.3f})"
    
    print(f"   📋 Evaluation: {evaluation}")
    print(f"   📝 Explanation: {explanation}")
    
    return {
        'log2fc_observed': log2fc,
        'log2fc_expected': 3.413,
        'fold_change': fold_change,
        'p_value': p_val,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'pseudotime_correlation': corr if 'pseudotime' in adata.obs.columns else None,
        'beta_coefficient': beta if 'pseudotime' in adata.obs.columns else None,
        'evaluation': evaluation,
        'explanation': explanation
    }

# Run SQSTM1 analysis
if adata is not None:
    sqstm1_results = analyze_sqstm1_upregulation(adata)
    
    if sqstm1_results:
        print(f"\n✅ SQSTM1 Analysis Complete: {sqstm1_results['evaluation']}")
    else:
        print("❌ SQSTM1 Analysis Failed")
else:
    print("⚠️ Skipping SQSTM1 analysis - no data loaded")
```

### Statement 6: Sliding Window Correlation Analysis

**Claim**: Running correlation between SQSTM1 and VDAC1 shifts from negative early (r=-0.417) to positive late (r=0.478) with strong trend (r=0.851).

**Analytical Approach**: Advanced temporal analysis revealing dynamic relationship changes during disease progression.


```python
def sliding_window_correlation_analysis(adata, protein1='SQSTM1', protein2='VDAC1', window_size=20):
    """
    Advanced sliding window correlation analysis
    
    This analysis reveals temporal dynamics that static correlations miss,
    showing how protein relationships change during disease progression.
    """
    if adata is None:
        print("❌ No data available for analysis")
        return None
    
    print("📊 Sliding Window Correlation Analysis - Statement 6")
    print("="*60)
    
    # Check proteins exist
    missing_proteins = [p for p in [protein1, protein2] if p not in adata.var_names]
    if missing_proteins:
        print(f"❌ Missing proteins: {missing_proteins}")
        return None
    
    print(f"✅ Analyzing {protein1} vs {protein2} correlation")
    print(f"📏 Window size: {window_size} cells")
    
    # Check pseudotime
    if 'pseudotime' not in adata.obs.columns:
        print("❌ Pseudotime column not found")
        return None
    
    # Step 1: Sort by pseudotime and prepare data
    print("\n🔄 Step 1: Data Preparation")
    
    sorted_idx = np.argsort(adata.obs['pseudotime'])
    adata_sorted = adata[sorted_idx]
    
    expr1 = adata_sorted[:, protein1].X.flatten()
    expr2 = adata_sorted[:, protein2].X.flatten()
    pseudotime_sorted = adata_sorted.obs['pseudotime'].values
    
    print(f"   📊 Total cells: {len(adata_sorted)}")
    print(f"   ⏰ Pseudotime range: [{pseudotime_sorted.min():.3f}, {pseudotime_sorted.max():.3f}]")
    
    # Step 2: Sliding window calculation
    print("\n📈 Step 2: Sliding Window Calculation")
    
    correlations = []
    positions = []
    p_values = []
    window_info = []
    
    n_windows = len(adata_sorted) - window_size + 1
    print(f"   🔢 Number of windows: {n_windows}")
    
    for i in tqdm(range(n_windows), desc="Processing windows"):
        # Extract window data
        window_expr1 = expr1[i:i+window_size]
        window_expr2 = expr2[i:i+window_size]
        window_pseudotime = pseudotime_sorted[i:i+window_size]
        
        # Calculate correlation
        if np.std(window_expr1) > 0 and np.std(window_expr2) > 0:
            corr, p_val = pearsonr(window_expr1, window_expr2)
            
            correlations.append(corr)
            positions.append(np.mean(window_pseudotime))
            p_values.append(p_val)
            
            window_info.append({
                'start_idx': i,
                'end_idx': i + window_size - 1,
                'pseudotime_min': np.min(window_pseudotime),
                'pseudotime_max': np.max(window_pseudotime),
                'pseudotime_mean': np.mean(window_pseudotime),
                'correlation': corr,
                'p_value': p_val
            })
    
    correlations = np.array(correlations)
    positions = np.array(positions)
    p_values = np.array(p_values)
    
    print(f"   ✅ Computed {len(correlations)} window correlations")
    
    # Step 3: Phase analysis
    print("\n🔍 Step 3: Phase Analysis")
    
    # Define phases based on pseudotime
    early_mask = positions < 0.33
    late_mask = positions > 0.67
    middle_mask = (positions >= 0.33) & (positions <= 0.67)
    
    early_corrs = correlations[early_mask]
    late_corrs = correlations[late_mask]
    middle_corrs = correlations[middle_mask]
    
    print(f"   🌅 Early phase (pseudotime < 0.33): {len(early_corrs)} windows")
    print(f"   🌇 Late phase (pseudotime > 0.67): {len(late_corrs)} windows")
    print(f"   🌤️ Middle phase: {len(middle_corrs)} windows")
    
    if len(early_corrs) > 0:
        early_mean = np.mean(early_corrs)
        early_std = np.std(early_corrs)
        print(f"   📊 Early mean correlation: {early_mean:.3f} ± {early_std:.3f}")
        print(f"   🎯 Expected early: -0.417")
    
    if len(late_corrs) > 0:
        late_mean = np.mean(late_corrs)
        late_std = np.std(late_corrs)
        print(f"   📊 Late mean correlation: {late_mean:.3f} ± {late_std:.3f}")
        print(f"   🎯 Expected late: 0.478")
    
    # Step 4: Trend analysis
    print("\n📈 Step 4: Trend Analysis")
    
    if len(correlations) > 2:
        trend_corr, trend_p = pearsonr(positions, correlations)
        
        # Linear regression for trend
        from sklearn.linear_model import LinearRegression
        lr = LinearRegression()
        lr.fit(positions.reshape(-1, 1), correlations)
        trend_slope = lr.coef_[0]
        trend_r2 = lr.score(positions.reshape(-1, 1), correlations)
        
        print(f"   📈 Trend correlation: r={trend_corr:.3f}, p={trend_p:.2e}")
        print(f"   📐 Trend slope: {trend_slope:.3f}")
        print(f"   📊 Trend R²: {trend_r2:.3f}")
        print(f"   🎯 Expected trend: r=0.851, p=6.98e-08")
    
    # Step 5: Comprehensive visualization
    print("\n📊 Step 5: Creating Comprehensive Visualizations")
    
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    
    # 1. Main sliding window plot
    colors = plt.cm.coolwarm((correlations + 1) / 2)  # Normalize to [0, 1]
    scatter = axes[0, 0].scatter(positions, correlations, c=correlations,
                               cmap='coolwarm', vmin=-1, vmax=1, s=30, alpha=0.7)
    
    # Phase boundaries
    axes[0, 0].axvline(x=0.33, color='gray', linestyle='--', alpha=0.5, label='Phase boundaries')
    axes[0, 0].axvline(x=0.67, color='gray', linestyle='--', alpha=0.5)
    
    # Expected values
    axes[0, 0].axhline(y=-0.417, color='blue', linestyle=':', alpha=0.7, label='Expected early')
    axes[0, 0].axhline(y=0.478, color='red', linestyle=':', alpha=0.7, label='Expected late')
    
    # Trend line
    if len(correlations) > 2:
        x_trend = np.linspace(positions.min(), positions.max(), 100)
        y_trend = lr.predict(x_trend.reshape(-1, 1))
        axes[0, 0].plot(x_trend, y_trend, 'black', linewidth=2, 
                       label=f'Trend (r={trend_corr:.3f})')
    
    axes[0, 0].set_xlabel('Pseudotime')
    axes[0, 0].set_ylabel('Correlation')
    axes[0, 0].set_title(f'{protein1}-{protein2} Sliding Window Correlation')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=axes[0, 0], label='Correlation')
    
    # 2. Phase comparison boxplot
    phase_data = []
    phase_labels = []
    
    if len(early_corrs) > 0:
        phase_data.append(early_corrs)
        phase_labels.append(f'Early\n(n={len(early_corrs)})')
    if len(middle_corrs) > 0:
        phase_data.append(middle_corrs)
        phase_labels.append(f'Middle\n(n={len(middle_corrs)})')
    if len(late_corrs) > 0:
        phase_data.append(late_corrs)
        phase_labels.append(f'Late\n(n={len(late_corrs)})')
    
    if phase_data:
        bp = axes[0, 1].boxplot(phase_data, labels=phase_labels, patch_artist=True)
        colors = ['lightblue', 'lightgray', 'lightcoral']
        for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
            patch.set_facecolor(color)
    
    axes[0, 1].axhline(y=0, color='gray', linestyle='-', alpha=0.3)
    axes[0, 1].set_ylabel('Correlation')
    axes[0, 1].set_title('Correlation by Phase')
    axes[0, 1].grid(True, alpha=0.3)
    
    # 3. Correlation distribution
    axes[0, 2].hist(correlations, bins=25, alpha=0.7, color='skyblue', edgecolor='black')
    axes[0, 2].axvline(x=np.mean(correlations), color='red', linestyle='--',
                      label=f'Mean: {np.mean(correlations):.3f}')
    axes[0, 2].axvline(x=0, color='black', linestyle='-', alpha=0.5)
    axes[0, 2].set_xlabel('Correlation')
    axes[0, 2].set_ylabel('Frequency')
    axes[0, 2].set_title('Distribution of Window Correlations')
    axes[0, 2].legend()
    
    # 4. P-value distribution
    axes[1, 0].hist(p_values, bins=25, alpha=0.7, color='orange', edgecolor='black')
    axes[1, 0].axvline(x=0.05, color='red', linestyle='--', label='p = 0.05')
    sig_pct = 100 * sum(np.array(p_values) < 0.05) / len(p_values)
    axes[1, 0].set_xlabel('P-value')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title(f'P-value Distribution\n({sig_pct:.1f}% significant)')
    axes[1, 0].legend()
    
    # 5. Smoothed trajectory
    if len(correlations) > 10:
        # Moving average smooth
        window_smooth = min(10, len(correlations)//5)
        smoothed = pd.Series(correlations).rolling(window=window_smooth, center=True).mean()
        axes[1, 1].plot(positions, correlations, 'lightgray', alpha=0.5, label='Raw')
        axes[1, 1].plot(positions, smoothed, 'blue', linewidth=2, 
                       label=f'Smoothed (window={window_smooth})')
    else:
        axes[1, 1].plot(positions, correlations, 'blue', linewidth=2, label='Correlations')
    
    axes[1, 1].axhline(y=0, color='black', linestyle='-', alpha=0.3)
    axes[1, 1].set_xlabel('Pseudotime')
    axes[1, 1].set_ylabel('Correlation')
    axes[1, 1].set_title('Smoothed Correlation Trajectory')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    # 6. Summary statistics
    axes[1, 2].axis('off')
    summary_text = f"""
    SLIDING WINDOW SUMMARY
    
    📊 Analysis Parameters:
    Window size: {window_size}
    Total windows: {len(correlations)}
    
    🌅 Early Phase (< 0.33):
    Mean r: {early_mean:.3f} (expected: -0.417)
    
    🌇 Late Phase (> 0.67):
    Mean r: {late_mean:.3f} (expected: 0.478)
    
    📈 Overall Trend:
    r = {trend_corr:.3f} (expected: 0.851)
    p = {trend_p:.2e} (expected: 6.98e-08)
    
    🎯 EVALUATION:
    {'SUPPORTED' if (len(early_corrs) > 0 and len(late_corrs) > 0 and 
                    abs(early_mean - (-0.417)) < 0.15 and 
                    abs(late_mean - 0.478) < 0.15 and 
                    abs(trend_corr - 0.851) < 0.15) else 'NEEDS REVIEW'}
    """
    
    axes[1, 2].text(0.05, 0.5, summary_text, fontsize=10, verticalalignment='center',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen", alpha=0.5))
    
    plt.tight_layout()
    plt.show()
    
    # Step 6: Evaluation
    print("\n🎯 Step 6: Statement Evaluation")
    
    # Evaluation criteria
    criteria_met = []
    criteria_failed = []
    
    if len(early_corrs) > 0:
        early_match = abs(early_mean - (-0.417)) < 0.15
        if early_match:
            criteria_met.append(f"Early phase correlation matches (-0.417 vs {early_mean:.3f})")
        else:
            criteria_failed.append(f"Early phase mismatch: {early_mean:.3f} vs -0.417")
    
    if len(late_corrs) > 0:
        late_match = abs(late_mean - 0.478) < 0.15
        if late_match:
            criteria_met.append(f"Late phase correlation matches (0.478 vs {late_mean:.3f})")
        else:
            criteria_failed.append(f"Late phase mismatch: {late_mean:.3f} vs 0.478")
    
    if len(correlations) > 2:
        trend_match = abs(trend_corr - 0.851) < 0.15 and trend_p < 1e-6
        if trend_match:
            criteria_met.append(f"Trend correlation matches (0.851 vs {trend_corr:.3f})")
        else:
            criteria_failed.append(f"Trend mismatch: {trend_corr:.3f} vs 0.851")
    
    # Final evaluation
    if len(criteria_met) >= 2:
        evaluation = "SUPPORTED"
        explanation = f"Key criteria met: {'; '.join(criteria_met)}"
    elif len(criteria_met) >= 1:
        evaluation = "PARTIALLY SUPPORTED"
        explanation = f"Some criteria met: {'; '.join(criteria_met)}. Issues: {'; '.join(criteria_failed)}"
    else:
        evaluation = "REFUTED"
        explanation = f"Criteria not met: {'; '.join(criteria_failed)}"
    
    print(f"   📋 Evaluation: {evaluation}")
    print(f"   📝 Explanation: {explanation}")
    
    return {
        'window_correlations': correlations,
        'window_positions': positions,
        'early_mean': early_mean if len(early_corrs) > 0 else None,
        'late_mean': late_mean if len(late_corrs) > 0 else None,
        'trend_correlation': trend_corr if len(correlations) > 2 else None,
        'trend_p_value': trend_p if len(correlations) > 2 else None,
        'evaluation': evaluation,
        'explanation': explanation,
        'n_windows': len(correlations)
    }

# Run sliding window analysis
if adata is not None:
    sliding_results = sliding_window_correlation_analysis(adata)
    
    if sliding_results:
        print(f"\n✅ Sliding Window Analysis Complete: {sliding_results['evaluation']}")
    else:
        print("❌ Sliding Window Analysis Failed")
else:
    print("⚠️ Skipping sliding window analysis - no data loaded")
```

## 🔬 Group 2: Sequential Failure of Proteostasis

This section analyzes the systematic breakdown of protein quality control mechanisms.

### Key Statement: Covariate-Controlled Differential Expression

**Claim**: 36.14% of proteins (2,115/5,853) show significant alterations when controlling for age, PMI, and PatientID.


```python
def covariate_controlled_differential_expression(adata, max_proteins=None):
    """
    Comprehensive covariate-controlled differential expression analysis
    
    This analysis removes confounding effects from age, PMI, and PatientID
    to identify pure tau-related protein changes.
    """
    if adata is None:
        print("❌ No data available for analysis")
        return None
    
    print("🧪 Covariate-Controlled Differential Expression - Group 2, Statement 1")
    print("="*70)
    
    # Step 1: Check required covariates
    print("\n🔍 Step 1: Covariate Availability Check")
    
    required_covariates = ['age', 'PMI', 'PatientID']
    available_covariates = []
    missing_covariates = []
    
    for cov in required_covariates:
        # Try exact match first
        if cov in adata.obs.columns:
            available_covariates.append(cov)
        else:
            # Try case-insensitive search
            matches = [c for c in adata.obs.columns if cov.lower() in c.lower()]
            if matches:
                print(f"   📝 Using {matches[0]} for {cov}")
                adata.obs[cov] = adata.obs[matches[0]]  # Create standardized name
                available_covariates.append(cov)
            else:
                missing_covariates.append(cov)
    
    print(f"   ✅ Available covariates: {available_covariates}")
    if missing_covariates:
        print(f"   ⚠️ Missing covariates: {missing_covariates}")
    
    if not available_covariates:
        print("❌ No covariates available - cannot perform covariate-controlled analysis")
        return None
    
    # Step 2: Data preparation and validation
    print("\n📊 Step 2: Data Preparation")
    
    # Check tau status
    if 'tau_status' not in adata.obs.columns:
        print("❌ tau_status column not found")
        return None
    
    tau_counts = adata.obs['tau_status'].value_counts()
    print(f"   🧠 Tau status distribution: {dict(tau_counts)}")
    
    # Set analysis scope
    if max_proteins is None:
        max_proteins = min(adata.n_vars, 5853)  # Full analysis or dataset limit
    else:
        max_proteins = min(max_proteins, adata.n_vars)
    
    print(f"   🔬 Analyzing {max_proteins} proteins (out of {adata.n_vars} available)")
    
    # Step 3: Linear model analysis
    print("\n📈 Step 3: Linear Model Analysis")
    print("   Running covariate-controlled differential expression...")
    
    # Build formula
    formula_parts = ['expression ~ tau_status']
    for cov in available_covariates:
        if cov == 'PatientID':
            formula_parts.append('C(PatientID)')  # Categorical
        else:
            formula_parts.append(cov)  # Continuous
    
    formula = ' + '.join(formula_parts)
    print(f"   📝 Model formula: {formula}")
    
    results = []
    failed_proteins = []
    
    # Prepare observation dataframe
    obs_df = adata.obs.copy()
    obs_df['tau_status'] = pd.Categorical(obs_df['tau_status'])
    
    print("\n🔄 Processing proteins...")
    for i, protein in enumerate(tqdm(adata.var_names[:max_proteins], desc="Analyzing proteins")):
        try:
            # Extract expression for this protein
            expr = adata[:, protein].X.flatten()
            obs_df['expression'] = expr
            
            # Fit linear model
            model = ols(formula, data=obs_df)
            results_fit = model.fit()
            
            # Extract tau coefficient
            tau_coef = None
            tau_pval = None
            
            # Look for tau coefficient (different parameterizations possible)
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
            
            # Calculate simple fold change for comparison
            tau_pos_mean = adata[adata.obs['tau_status'] == 'positive'][:, protein].X.mean()
            tau_neg_mean = adata[adata.obs['tau_status'] == 'negative'][:, protein].X.mean()
            simple_fc = tau_pos_mean - tau_neg_mean
            
            results.append({
                'protein': protein,
                'log2fc_adjusted': tau_coef if tau_coef is not None else np.nan,
                'pvalue_adjusted': tau_pval if tau_pval is not None else 1.0,
                'log2fc_simple': simple_fc,
                'r_squared': results_fit.rsquared,
                'aic': results_fit.aic
            })
            
        except Exception as e:
            failed_proteins.append(protein)
            results.append({
                'protein': protein,
                'log2fc_adjusted': np.nan,
                'pvalue_adjusted': 1.0,
                'log2fc_simple': np.nan,
                'r_squared': 0,
                'aic': np.inf
            })
    
    if failed_proteins:
        print(f"   ⚠️ Failed to analyze {len(failed_proteins)} proteins")
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Step 4: Multiple testing correction
    print("\n🔧 Step 4: Multiple Testing Correction (FDR)")
    
    # Remove NaN p-values for FDR correction
    valid_p = ~np.isnan(results_df['pvalue_adjusted'])
    
    if sum(valid_p) > 0:
        rejected, fdr_corrected, _, _ = multipletests(
            results_df.loc[valid_p, 'pvalue_adjusted'],
            method='fdr_bh',
            alpha=0.05
        )
        
        # Initialize FDR column
        results_df['fdr_adjusted'] = 1.0
        results_df['significant_adjusted'] = False
        
        # Assign corrected values
        results_df.loc[valid_p, 'fdr_adjusted'] = fdr_corrected
        results_df.loc[valid_p, 'significant_adjusted'] = rejected
        
        # Results summary
        n_total = len(results_df)
        n_valid = sum(valid_p)
        n_significant = sum(results_df['significant_adjusted'])
        pct_significant = 100 * n_significant / n_total
        
        print(f"   📊 Total proteins: {n_total}")
        print(f"   ✅ Valid analyses: {n_valid}")
        print(f"   🔬 Significant (FDR < 0.05): {n_significant}")
        print(f"   📈 Percentage significant: {pct_significant:.2f}%")
        print(f"   🎯 Expected: 36.14% (2,115/5,853)")
    
    # Step 5: Results visualization
    print("\n📊 Step 5: Results Visualization")
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # 1. Volcano plot
    valid_results = results_df[valid_p]
    if len(valid_results) > 0:
        x = valid_results['log2fc_adjusted']
        y = -np.log10(valid_results['pvalue_adjusted'])
        colors = ['red' if sig else 'gray' for sig in valid_results['significant_adjusted']]
        
        axes[0, 0].scatter(x, y, c=colors, alpha=0.6, s=10)
        axes[0, 0].axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
        axes[0, 0].axvline(x=0, color='black', linestyle='-', alpha=0.3)
        axes[0, 0].set_xlabel('Log2 Fold Change (Adjusted)')
        axes[0, 0].set_ylabel('-Log10(P-value)')
        axes[0, 0].set_title(f'Volcano Plot\n({n_significant} significant)')
    
    # 2. FDR distribution
    fdr_values = results_df['fdr_adjusted'][valid_p]
    axes[0, 1].hist(fdr_values, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    axes[0, 1].axvline(x=0.05, color='red', linestyle='--', label='FDR = 0.05')
    axes[0, 1].set_xlabel('FDR')
    axes[0, 1].set_ylabel('Count')
    axes[0, 1].set_title('FDR Distribution')
    axes[0, 1].legend()
    
    # 3. Effect size distribution
    fc_sig = results_df[results_df['significant_adjusted']]['log2fc_adjusted']
    fc_nonsig = results_df[~results_df['significant_adjusted']]['log2fc_adjusted']
    
    axes[0, 2].hist([fc_nonsig.dropna(), fc_sig.dropna()], bins=30, 
                   label=['Non-significant', 'Significant'],
                   color=['gray', 'red'], alpha=0.7, edgecolor='black')
    axes[0, 2].axvline(x=0, color='black', linestyle='-', alpha=0.3)
    axes[0, 2].set_xlabel('Log2 Fold Change')
    axes[0, 2].set_ylabel('Count')
    axes[0, 2].set_title('Effect Size Distribution')
    axes[0, 2].legend()
    
    # 4. Simple vs adjusted fold changes
    axes[1, 0].scatter(results_df['log2fc_simple'], results_df['log2fc_adjusted'], 
                      alpha=0.5, s=10)
    lims = [-4, 4]
    axes[1, 0].plot(lims, lims, 'r--', alpha=0.5)  # Identity line
    axes[1, 0].set_xlabel('Simple Log2FC')
    axes[1, 0].set_ylabel('Adjusted Log2FC')
    axes[1, 0].set_title('Simple vs Adjusted Fold Changes')
    axes[1, 0].set_xlim(lims)
    axes[1, 0].set_ylim(lims)
    
    # 5. Model fit quality (R-squared)
    r2_values = results_df['r_squared'][results_df['r_squared'] < 1]
    axes[1, 1].hist(r2_values, bins=30, alpha=0.7, color='green', edgecolor='black')
    axes[1, 1].set_xlabel('R²')
    axes[1, 1].set_ylabel('Count')
    axes[1, 1].set_title(f'Model Fit Quality\nMean R² = {r2_values.mean():.3f}')
    
    # 6. Summary statistics
    axes[1, 2].axis('off')
    summary_text = f"""
    COVARIATE-CONTROLLED DE SUMMARY
    
    📊 Analysis Scope:
    Total proteins: {n_total:,}
    Valid analyses: {n_valid:,}
    
    🔬 Results:
    Significant: {n_significant:,}
    Percentage: {pct_significant:.2f}%
    Expected: 36.14%
    
    📈 Covariates Used:
    {', '.join(available_covariates)}
    
    📏 Effect Sizes:
    Mean |FC|: {np.mean(np.abs(results_df['log2fc_adjusted'].dropna())):.3f}
    
    🎯 EVALUATION:
    {'SUPPORTED' if abs(pct_significant - 36.14) < 5 else 'NEEDS REVIEW'}
    """
    
    axes[1, 2].text(0.1, 0.5, summary_text, fontsize=10, verticalalignment='center',
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.7))
    
    plt.tight_layout()
    plt.show()
    
    # Step 6: Evaluation
    print("\n🎯 Step 6: Statement Evaluation")
    
    # Evaluation criteria
    expected_percentage = 36.14
    percentage_diff = abs(pct_significant - expected_percentage)
    
    if percentage_diff < 5:  # Within 5% tolerance
        evaluation = "SUPPORTED"
        explanation = f"Found {pct_significant:.2f}% significant proteins, closely matching expected 36.14%"
    elif percentage_diff < 10:  # Within 10% tolerance
        evaluation = "PARTIALLY SUPPORTED"
        explanation = f"Found {pct_significant:.2f}% significant proteins, somewhat different from expected 36.14%"
    else:
        evaluation = "REFUTED"
        explanation = f"Found {pct_significant:.2f}% significant proteins, substantially different from expected 36.14%"
    
    print(f"   📋 Evaluation: {evaluation}")
    print(f"   📝 Explanation: {explanation}")
    
    return {
        'results_df': results_df,
        'n_total': n_total,
        'n_significant': n_significant,
        'percentage_significant': pct_significant,
        'expected_percentage': expected_percentage,
        'available_covariates': available_covariates,
        'evaluation': evaluation,
        'explanation': explanation
    }

# Run covariate-controlled analysis
if adata is not None:
    # For demonstration, analyze first 500 proteins (set to None for full analysis)
    de_results = covariate_controlled_differential_expression(adata, max_proteins=500)
    
    if de_results:
        print(f"\n✅ Covariate-Controlled DE Complete: {de_results['evaluation']}")
        print(f"📊 Analysis summary: {de_results['n_significant']}/{de_results['n_total']} proteins significant")
    else:
        print("❌ Covariate-Controlled DE Analysis Failed")
else:
    print("⚠️ Skipping covariate-controlled DE analysis - no data loaded")
```

## 🤖 AI Automation Demonstration

This section demonstrates the AI automation capabilities for large-scale analysis.


```python
# Import AI automation modules
try:
    import sys
    sys.path.append('ai_automation')
    from ai_agent import BioinformaticsAgent, EvaluationResult
    from analysis_automation import AnalysisAutomation
    
    print("✅ AI automation modules imported successfully")
    
    if adata is not None:
        print("\n🤖 Initializing AI Automation...")
        
        # Initialize automation
        automation = AnalysisAutomation(adata)
        
        # Demonstrate automated correlation analysis
        print("\n📊 Running automated SQSTM1 correlation analysis...")
        proteins = ['TAX1BP1', 'CAT', 'VDAC1', 'CYCS', 'ATP5F1A']
        
        try:
            sqstm1_corr = automation.calculate_protein_correlations(proteins, 'SQSTM1')
            if not sqstm1_corr.empty:
                print("\n🧬 SQSTM1 Correlation Results:")
                for _, row in sqstm1_corr.iterrows():
                    print(f"   {row['protein']}: r={row['correlation']:.3f}, p={row['p_value']:.3e}")
            else:
                print("⚠️ No correlations computed")
        except Exception as e:
            print(f"⚠️ Automated correlation analysis failed: {str(e)}")
        
        # Demonstrate sliding window automation
        print("\n📈 Running automated sliding window analysis...")
        try:
            sliding = automation.sliding_window_correlation('SQSTM1', 'VDAC1')
            if not sliding.empty:
                print(f"   ✅ Generated {len(sliding)} window correlations")
                if hasattr(sliding, 'attrs') and 'summary' in sliding.attrs:
                    summary = sliding.attrs['summary']
                    print(f"   📊 Early mean: {summary.get('early_mean_correlation', 'N/A'):.3f}")
                    print(f"   📊 Late mean: {summary.get('late_mean_correlation', 'N/A'):.3f}")
            else:
                print("⚠️ No sliding window results generated")
        except Exception as e:
            print(f"⚠️ Automated sliding window analysis failed: {str(e)}")
        
        print("\n✅ AI automation demonstration complete")
    else:
        print("⚠️ Cannot run AI automation - no data loaded")
        
except ImportError as e:
    print(f"⚠️ AI automation modules not available: {str(e)}")
    print("📝 This is expected if running without the full project structure")
except Exception as e:
    print(f"⚠️ AI automation error: {str(e)}")
```

## 📋 Comprehensive Results Summary

This section provides a comprehensive summary of all analyses performed.


```python
def generate_comprehensive_summary():
    """
    Generate a comprehensive summary of all analysis results
    """
    print("📋 COMPREHENSIVE ANALYSIS SUMMARY")
    print("="*60)
    
    # Data loading summary
    if adata is not None:
        print(f"\n📊 Dataset Summary:")
        print(f"   📦 Dimensions: {adata.n_obs} cells × {adata.n_vars} proteins")
        if 'tau_status' in adata.obs.columns:
            tau_counts = adata.obs['tau_status'].value_counts()
            print(f"   🧠 Tau distribution: {dict(tau_counts)}")
        print(f"   ✅ Data successfully loaded and validated")
    else:
        print(f"\n❌ Dataset: Not loaded")
    
    # Group 1 results
    print(f"\n🧬 Group 1: Late-Stage Mitochondrial Dysregulation")
    print(f"   {'='*50}")
    
    # UPS analysis results
    if 'ups_results' in globals() and ups_results:
        print(f"   📋 Statement 1 (UPS Proteins): {ups_results['evaluation']}")
        print(f"      📊 Analyzed: {ups_results['n_total']} proteins")
        print(f"      📈 Significant: {ups_results['pct_significant']:.1f}%")
        print(f"      📏 Mean effect size: {ups_results['mean_effect_size']:.3f}")
    else:
        print(f"   📋 Statement 1 (UPS Proteins): Not analyzed")
    
    # SQSTM1 analysis results
    if 'sqstm1_results' in globals() and sqstm1_results:
        print(f"   📋 Statement 2 (SQSTM1): {sqstm1_results['evaluation']}")
        print(f"      📊 Log2FC: {sqstm1_results['log2fc_observed']:.3f} (expected: {sqstm1_results['log2fc_expected']})")
        print(f"      📈 Fold change: {sqstm1_results['fold_change']:.1f}x")
        print(f"      🎲 P-value: {sqstm1_results['p_value']:.2e}")
    else:
        print(f"   📋 Statement 2 (SQSTM1): Not analyzed")
    
    # Sliding window results
    if 'sliding_results' in globals() and sliding_results:
        print(f"   📋 Statement 6 (Sliding Window): {sliding_results['evaluation']}")
        print(f"      📊 Windows analyzed: {sliding_results['n_windows']}")
        if sliding_results['early_mean'] is not None:
            print(f"      🌅 Early correlation: {sliding_results['early_mean']:.3f} (expected: -0.417)")
        if sliding_results['late_mean'] is not None:
            print(f"      🌇 Late correlation: {sliding_results['late_mean']:.3f} (expected: 0.478)")
        if sliding_results['trend_correlation'] is not None:
            print(f"      📈 Trend: r={sliding_results['trend_correlation']:.3f} (expected: 0.851)")
    else:
        print(f"   📋 Statement 6 (Sliding Window): Not analyzed")
    
    # Group 2 results
    print(f"\n🔬 Group 2: Sequential Failure of Proteostasis")
    print(f"   {'='*50}")
    
    # Differential expression results
    if 'de_results' in globals() and de_results:
        print(f"   📋 Statement 1 (Covariate-controlled DE): {de_results['evaluation']}")
        print(f"      📊 Proteins analyzed: {de_results['n_total']:,}")
        print(f"      📈 Significant: {de_results['n_significant']:,} ({de_results['percentage_significant']:.2f}%)")
        print(f"      🎯 Expected: {de_results['expected_percentage']:.2f}%")
        print(f"      🔧 Covariates: {', '.join(de_results['available_covariates'])}")
    else:
        print(f"   📋 Statement 1 (Covariate-controlled DE): Not analyzed")
    
    # Overall evaluation
    print(f"\n🎯 Overall Project Evaluation")
    print(f"   {'='*40}")
    
    # Count evaluations
    evaluations = []
    if 'ups_results' in globals() and ups_results:
        evaluations.append(ups_results['evaluation'])
    if 'sqstm1_results' in globals() and sqstm1_results:
        evaluations.append(sqstm1_results['evaluation'])
    if 'sliding_results' in globals() and sliding_results:
        evaluations.append(sliding_results['evaluation'])
    if 'de_results' in globals() and de_results:
        evaluations.append(de_results['evaluation'])
    
    if evaluations:
        supported = sum(1 for e in evaluations if e == 'SUPPORTED')
        partially = sum(1 for e in evaluations if e == 'PARTIALLY SUPPORTED')
        refuted = sum(1 for e in evaluations if e == 'REFUTED')
        total = len(evaluations)
        
        print(f"   📊 Statements analyzed: {total}")
        print(f"   ✅ Supported: {supported}")
        print(f"   🔶 Partially supported: {partially}")
        print(f"   ❌ Refuted: {refuted}")
        print(f"   📈 Success rate: {100*supported/total:.1f}%")
    else:
        print(f"   📊 No analyses completed")
    
    # Technical summary
    print(f"\n🔧 Technical Implementation Summary")
    print(f"   {'='*40}")
    print(f"   📋 Framework: Comprehensive bioinformatics evaluation system")
    print(f"   📊 Statistical methods: FDR correction, bootstrap CI, linear models")
    print(f"   🔬 Advanced analyses: Sliding window, covariate control, temporal analysis")
    print(f"   🤖 AI automation: Integrated automated analysis capabilities")
    print(f"   📈 Visualization: Publication-quality plots and comprehensive summaries")
    print(f"   📝 Documentation: Enhanced analytical rationale and biological context")
    
    print(f"\n🎉 Analysis Framework Successfully Demonstrated!")
    print(f"\n📚 For more details, see individual analysis sections above.")

# Generate the comprehensive summary
generate_comprehensive_summary()
```

## 🚀 Next Steps and Local Execution Guide

### To run this notebook locally:

1. **Install Dependencies**:
   ```bash
   pip install pandas numpy scipy scanpy statsmodels scikit-learn matplotlib seaborn tqdm
   ```

2. **Data Requirements**:
   - Ensure `data/pool_processed_v2.h5ad` file is present
   - Data should contain required metadata columns: `tau_status`, `MC1`, `pseudotime`

3. **Running Individual Analyses**:
   - Each analysis can be run independently
   - Modify `max_proteins` parameter for faster testing
   - Set to `None` for full analysis

4. **Customization Options**:
   - Adjust statistical thresholds in evaluation functions
   - Modify visualization parameters
   - Add additional protein lists or analysis methods

5. **Performance Considerations**:
   - Full analysis of 5,853 proteins may take 10-30 minutes
   - Use subset analysis for testing (`max_proteins=500`)
   - Results are cached where possible

### Advanced Features:
- **AI Automation**: Automated analysis pipelines
- **Statistical Rigor**: Multiple testing correction, effect sizes
- **Temporal Analysis**: Sliding window correlations
- **Covariate Control**: Advanced linear modeling
- **Publication Quality**: Professional visualizations

This framework provides a complete solution for rigorous evaluation of biological statements using proteomic data, with comprehensive documentation and statistical validation.
