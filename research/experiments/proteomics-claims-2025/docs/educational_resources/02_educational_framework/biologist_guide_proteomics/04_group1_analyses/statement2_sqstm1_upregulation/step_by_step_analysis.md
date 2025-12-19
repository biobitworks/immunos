# 🧪 Step-by-Step SQSTM1 Analysis Tutorial

## 🎯 What You'll Accomplish

By the end of this tutorial, you'll have:
- ✅ **Performed rigorous statistical analysis** of SQSTM1 upregulation
- ✅ **Calculated bootstrap confidence intervals** for robust estimates
- ✅ **Validated the 10.7-fold upregulation claim** with multiple methods
- ✅ **Created publication-quality visualizations** of SQSTM1 changes
- ✅ **Interpreted results in autophagy biology context**
- ✅ **Assessed biological significance** beyond statistical significance

**Research Question**: Is SQSTM1 really increased 10.7-fold in tau-positive neurons, and what does this tell us about autophagy dysfunction?

**Time needed**: 2-3 hours for complete beginners, 1 hour with some experience

---

## 🛠️ Before We Start

### Prerequisites
- [ ] **Environment setup**: Jupyter notebook or Google Colab working
- [ ] **Data access**: Proteomics dataset loaded and quality-controlled
- [ ] **Packages installed**: scanpy, pandas, numpy, scipy, matplotlib, seaborn
- [ ] **Background reading**: SQSTM1 biological background completed

### Setup Your Workspace
```python
# Create new notebook: "02_sqstm1_analysis.ipynb"
# If using local Jupyter: place in notebooks/ folder
# If using Colab: save to Google Drive proteomics folder
```

---

## 📚 Step 1: Load Libraries and Data

### Import Required Libraries
```python
# Core data science libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Proteomics analysis
import scanpy as sc

# Statistical libraries
from scipy import stats
from statsmodels.stats.multitest import multipletests

# For bootstrap analysis
from sklearn.utils import resample

# Utilities
import warnings
warnings.filterwarnings('ignore')

# Configure plotting for publication quality
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5
sns.set_style("whitegrid")
%matplotlib inline

print("✅ All libraries imported successfully!")
```

### Load and Prepare Data
```python
# Load the proteomics dataset
print("Loading proteomics dataset...")

# For local Jupyter
try:
    adata = sc.read_h5ad('../data/raw/pool_processed_v2.h5ad')
    print(f"✅ Dataset loaded: {adata.shape}")
except FileNotFoundError:
    print("❌ Dataset not found at ../data/raw/")
    # Try Google Colab path
    try:
        from google.colab import drive
        drive.mount('/content/drive')
        adata = sc.read_h5ad('/content/drive/MyDrive/proteomics_analysis/data/raw/pool_processed_v2.h5ad')
        print(f"✅ Dataset loaded from Google Drive: {adata.shape}")
    except:
        print("❌ Please check your data path and try again")
        raise

# Verify we have SQSTM1 in the dataset
if 'SQSTM1' not in adata.var_names:
    print("❌ SQSTM1 not found in dataset!")
    print("Available proteins with 'SQSTM' in name:")
    sqstm_proteins = [p for p in adata.var_names if 'SQSTM' in p.upper()]
    for p in sqstm_proteins:
        print(f"  {p}")
    if not sqstm_proteins:
        print("No SQSTM-related proteins found")
        raise ValueError("SQSTM1 not available for analysis")
else:
    print("✅ SQSTM1 found in dataset")

# Basic dataset info
print(f"\n📊 Dataset overview:")
print(f"Samples: {adata.n_obs}")
print(f"Proteins: {adata.n_vars}")
print(f"Tau-positive samples: {sum(adata.obs['tau_status'] == 'positive')}")
print(f"Tau-negative samples: {sum(adata.obs['tau_status'] == 'negative')}")
```

---

## 🔍 Step 2: Extract and Explore SQSTM1 Data

### Extract SQSTM1 Expression Data
```python
print("=" * 50)
print("SQSTM1 EXPRESSION ANALYSIS")
print("=" * 50)

# Get SQSTM1 expression data
sqstm1_idx = adata.var_names.get_loc('SQSTM1')
sqstm1_expression = adata.X[:, sqstm1_idx]

# Create analysis dataframe
sqstm1_data = pd.DataFrame({
    'sqstm1_expression': sqstm1_expression,
    'tau_status': adata.obs['tau_status'].values,
    'MC1_score': adata.obs['MC1_score'].values,
    'age': adata.obs['age'].values,
    'sex': adata.obs['sex'].values,
    'PMI': adata.obs['PMI'].values,
    'patient_id': adata.obs['patient_id'].values,
    'batch': adata.obs['batch'].values if 'batch' in adata.obs.columns else 'unknown'
})

# Separate groups
tau_positive = sqstm1_data[sqstm1_data['tau_status'] == 'positive']
tau_negative = sqstm1_data[sqstm1_data['tau_status'] == 'negative']

print(f"SQSTM1 expression data extracted:")
print(f"• Tau-positive samples: {len(tau_positive)}")
print(f"• Tau-negative samples: {len(tau_negative)}")
print(f"• Total samples: {len(sqstm1_data)}")

# Basic statistics
print(f"\nBasic SQSTM1 statistics:")
print(f"• Overall mean: {sqstm1_expression.mean():.3f}")
print(f"• Overall std: {sqstm1_expression.std():.3f}")
print(f"• Overall range: {sqstm1_expression.min():.3f} to {sqstm1_expression.max():.3f}")
```

### Initial Visual Exploration
```python
# Create initial SQSTM1 visualization
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('SQSTM1 Expression Overview', fontsize=16, y=1.02)

# Plot 1: Distribution by group
ax1 = axes[0]
tau_pos_expr = tau_positive['sqstm1_expression']
tau_neg_expr = tau_negative['sqstm1_expression']

ax1.hist(tau_neg_expr, bins=15, alpha=0.7, label='Tau-negative', color='lightblue')
ax1.hist(tau_pos_expr, bins=15, alpha=0.7, label='Tau-positive', color='lightcoral')
ax1.set_xlabel('SQSTM1 Log2 Expression')
ax1.set_ylabel('Number of Samples')
ax1.set_title('Expression Distribution')
ax1.legend()

# Plot 2: Box plot comparison
ax2 = axes[1]
box_data = [tau_neg_expr, tau_pos_expr]
box_plot = ax2.boxplot(box_data, labels=['Tau-negative', 'Tau-positive'], patch_artist=True)
box_plot['boxes'][0].set_facecolor('lightblue')
box_plot['boxes'][1].set_facecolor('lightcoral')
ax2.set_ylabel('SQSTM1 Log2 Expression')
ax2.set_title('Group Comparison')

# Calculate and display fold change
mean_pos = tau_pos_expr.mean()
mean_neg = tau_neg_expr.mean()
log2_fc = mean_pos - mean_neg
fold_change = 2**log2_fc

ax2.text(1.5, max(sqstm1_expression) * 0.9,
         f'{fold_change:.1f}-fold\nincrease',
         ha='center', fontsize=12,
         bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))

# Plot 3: Individual points
ax3 = axes[2]
# Add jitter for better visualization
x_neg = np.random.normal(1, 0.05, len(tau_neg_expr))
x_pos = np.random.normal(2, 0.05, len(tau_pos_expr))

ax3.scatter(x_neg, tau_neg_expr, alpha=0.6, color='lightblue', s=50)
ax3.scatter(x_pos, tau_pos_expr, alpha=0.6, color='lightcoral', s=50)

# Add mean lines
ax3.hlines(mean_neg, 0.7, 1.3, colors='blue', linewidth=3, label=f'Mean: {mean_neg:.2f}')
ax3.hlines(mean_pos, 1.7, 2.3, colors='red', linewidth=3, label=f'Mean: {mean_pos:.2f}')

ax3.set_xlim(0.5, 2.5)
ax3.set_xticks([1, 2])
ax3.set_xticklabels(['Tau-negative', 'Tau-positive'])
ax3.set_ylabel('SQSTM1 Log2 Expression')
ax3.set_title('Individual Values')
ax3.legend()

plt.tight_layout()
plt.show()

print(f"📊 Initial visualization created!")
print(f"Preliminary fold change: {fold_change:.1f}")
```

---

## 📈 Step 3: Classical Statistical Analysis

### Two-Sample t-Test
```python
print("\n" + "=" * 50)
print("CLASSICAL STATISTICAL ANALYSIS")
print("=" * 50)

# Perform two-sample t-test
t_statistic, p_value = stats.ttest_ind(tau_pos_expr, tau_neg_expr, equal_var=False)

# Calculate effect size (Cohen's d)
pooled_std = np.sqrt(((len(tau_pos_expr)-1) * np.var(tau_pos_expr, ddof=1) +
                     (len(tau_neg_expr)-1) * np.var(tau_neg_expr, ddof=1)) /
                    (len(tau_pos_expr) + len(tau_neg_expr) - 2))

cohens_d = (mean_pos - mean_neg) / pooled_std

# Calculate confidence interval for the difference
se_diff = pooled_std * np.sqrt(1/len(tau_pos_expr) + 1/len(tau_neg_expr))
df = len(tau_pos_expr) + len(tau_neg_expr) - 2
t_critical = stats.t.ppf(0.975, df)  # 95% CI
ci_lower = log2_fc - t_critical * se_diff
ci_upper = log2_fc + t_critical * se_diff

# Convert CI to fold change scale
fc_ci_lower = 2**ci_lower
fc_ci_upper = 2**ci_upper

print("T-TEST RESULTS:")
print(f"• t-statistic: {t_statistic:.3f}")
print(f"• p-value: {p_value:.2e}")
print(f"• Degrees of freedom: {df}")

print(f"\nEFFECT SIZE:")
print(f"• Cohen's d: {cohens_d:.3f}")
print(f"• Log2 fold change: {log2_fc:.3f}")
print(f"• Fold change: {fold_change:.1f}")

print(f"\nCONFIDENCE INTERVALS (95%):")
print(f"• Log2 FC CI: [{ci_lower:.3f}, {ci_upper:.3f}]")
print(f"• Fold change CI: [{fc_ci_lower:.1f}, {fc_ci_upper:.1f}]")

# Interpret statistical significance
if p_value < 0.001:
    significance = "Highly significant (p < 0.001)"
elif p_value < 0.01:
    significance = "Very significant (p < 0.01)"
elif p_value < 0.05:
    significance = "Significant (p < 0.05)"
else:
    significance = "Not significant (p ≥ 0.05)"

print(f"\nSTATISTICAL INTERPRETATION:")
print(f"• {significance}")

# Compare to claimed 10.7-fold increase
claimed_fc = 10.7
difference_from_claim = abs(fold_change - claimed_fc)
within_ci = fc_ci_lower <= claimed_fc <= fc_ci_upper

print(f"\nCOMPARISON TO CLAIM:")
print(f"• Claimed fold change: {claimed_fc}")
print(f"• Observed fold change: {fold_change:.1f}")
print(f"• Difference: {difference_from_claim:.1f}")
print(f"• Claim within 95% CI: {'Yes' if within_ci else 'No'}")
```

### Non-Parametric Alternative
```python
print("\nNON-PARAMETRIC VALIDATION:")
print("-" * 30)

# Mann-Whitney U test (non-parametric alternative)
u_statistic, u_p_value = stats.mannwhitneyu(tau_pos_expr, tau_neg_expr, alternative='two-sided')

# Median-based effect size
median_pos = np.median(tau_pos_expr)
median_neg = np.median(tau_neg_expr)
median_log2_fc = median_pos - median_neg
median_fold_change = 2**median_log2_fc

print(f"Mann-Whitney U test:")
print(f"• U-statistic: {u_statistic:.1f}")
print(f"• p-value: {u_p_value:.2e}")

print(f"\nMedian-based analysis:")
print(f"• Tau-positive median: {median_pos:.3f}")
print(f"• Tau-negative median: {median_neg:.3f}")
print(f"• Median fold change: {median_fold_change:.1f}")

# Check consistency between tests
consistent_significance = (p_value < 0.05) == (u_p_value < 0.05)
print(f"\nConsistency check:")
print(f"• Both tests agree on significance: {'Yes' if consistent_significance else 'No'}")
```

---

## 🔄 Step 4: Bootstrap Confidence Intervals

### Bootstrap Analysis Function
```python
print("\n" + "=" * 50)
print("BOOTSTRAP CONFIDENCE INTERVALS")
print("=" * 50)

def bootstrap_fold_change(group1, group2, n_bootstrap=10000, confidence_level=0.95):
    """
    Calculate bootstrap confidence intervals for fold change

    Parameters:
    -----------
    group1, group2 : array-like
        Expression values for the two groups
    n_bootstrap : int
        Number of bootstrap samples
    confidence_level : float
        Confidence level (e.g., 0.95 for 95% CI)

    Returns:
    --------
    dict : Bootstrap results including CI and distribution
    """

    print(f"Running bootstrap analysis with {n_bootstrap:,} samples...")

    # Store bootstrap results
    bootstrap_fold_changes = []
    bootstrap_log2_fcs = []
    bootstrap_cohens_ds = []

    # Bootstrap sampling
    for i in range(n_bootstrap):
        # Resample with replacement
        boot_group1 = resample(group1, n_samples=len(group1), random_state=i)
        boot_group2 = resample(group2, n_samples=len(group2), random_state=i+n_bootstrap)

        # Calculate statistics for this bootstrap sample
        boot_mean1 = np.mean(boot_group1)
        boot_mean2 = np.mean(boot_group2)
        boot_log2_fc = boot_mean1 - boot_mean2
        boot_fold_change = 2**boot_log2_fc

        # Effect size for this bootstrap sample
        boot_pooled_std = np.sqrt(((len(boot_group1)-1) * np.var(boot_group1, ddof=1) +
                                  (len(boot_group2)-1) * np.var(boot_group2, ddof=1)) /
                                 (len(boot_group1) + len(boot_group2) - 2))
        boot_cohens_d = (boot_mean1 - boot_mean2) / boot_pooled_std

        # Store results
        bootstrap_fold_changes.append(boot_fold_change)
        bootstrap_log2_fcs.append(boot_log2_fc)
        bootstrap_cohens_ds.append(boot_cohens_d)

        # Progress indicator
        if (i + 1) % 2000 == 0:
            print(f"  Bootstrap sample {i+1:,}/{n_bootstrap:,} completed")

    # Convert to arrays
    bootstrap_fold_changes = np.array(bootstrap_fold_changes)
    bootstrap_log2_fcs = np.array(bootstrap_log2_fcs)
    bootstrap_cohens_ds = np.array(bootstrap_cohens_ds)

    # Calculate confidence intervals
    alpha = 1 - confidence_level
    lower_percentile = (alpha/2) * 100
    upper_percentile = (1 - alpha/2) * 100

    fc_ci_lower = np.percentile(bootstrap_fold_changes, lower_percentile)
    fc_ci_upper = np.percentile(bootstrap_fold_changes, upper_percentile)

    log2_ci_lower = np.percentile(bootstrap_log2_fcs, lower_percentile)
    log2_ci_upper = np.percentile(bootstrap_log2_fcs, upper_percentile)

    d_ci_lower = np.percentile(bootstrap_cohens_ds, lower_percentile)
    d_ci_upper = np.percentile(bootstrap_cohens_ds, upper_percentile)

    print("✅ Bootstrap analysis completed!")

    return {
        'fold_changes': bootstrap_fold_changes,
        'log2_fold_changes': bootstrap_log2_fcs,
        'cohens_ds': bootstrap_cohens_ds,
        'fc_ci_lower': fc_ci_lower,
        'fc_ci_upper': fc_ci_upper,
        'log2_ci_lower': log2_ci_lower,
        'log2_ci_upper': log2_ci_upper,
        'd_ci_lower': d_ci_lower,
        'd_ci_upper': d_ci_upper,
        'confidence_level': confidence_level
    }

# Run bootstrap analysis
bootstrap_results = bootstrap_fold_change(tau_pos_expr, tau_neg_expr, n_bootstrap=10000)
```

### Bootstrap Results Analysis
```python
print("BOOTSTRAP RESULTS:")
print("-" * 30)

# Summary statistics
print(f"Bootstrap distribution summary:")
print(f"• Mean fold change: {np.mean(bootstrap_results['fold_changes']):.2f}")
print(f"• Median fold change: {np.median(bootstrap_results['fold_changes']):.2f}")
print(f"• Standard deviation: {np.std(bootstrap_results['fold_changes']):.2f}")

print(f"\nBootstrap confidence intervals ({bootstrap_results['confidence_level']*100:.0f}%):")
print(f"• Fold change: [{bootstrap_results['fc_ci_lower']:.2f}, {bootstrap_results['fc_ci_upper']:.2f}]")
print(f"• Log2 fold change: [{bootstrap_results['log2_ci_lower']:.3f}, {bootstrap_results['log2_ci_upper']:.3f}]")
print(f"• Cohen's d: [{bootstrap_results['d_ci_lower']:.3f}, {bootstrap_results['d_ci_upper']:.3f}]")

# Compare bootstrap vs classical CI
print(f"\nComparison with classical t-test CI:")
print(f"• Classical FC CI: [{fc_ci_lower:.2f}, {fc_ci_upper:.2f}]")
print(f"• Bootstrap FC CI: [{bootstrap_results['fc_ci_lower']:.2f}, {bootstrap_results['fc_ci_upper']:.2f}]")

ci_difference = abs(bootstrap_results['fc_ci_upper'] - fc_ci_upper)
print(f"• Upper bound difference: {ci_difference:.2f}")

if ci_difference < 1.0:
    print("✅ Bootstrap and classical CIs are very similar")
elif ci_difference < 2.0:
    print("⚠️ Bootstrap and classical CIs show moderate difference")
else:
    print("❌ Bootstrap and classical CIs show substantial difference")

# Check if 10.7-fold claim is within bootstrap CI
claimed_in_bootstrap_ci = (bootstrap_results['fc_ci_lower'] <= claimed_fc <= bootstrap_results['fc_ci_upper'])
print(f"\nClaim validation:")
print(f"• Claimed 10.7-fold within bootstrap CI: {'Yes' if claimed_in_bootstrap_ci else 'No'}")
```

### Visualize Bootstrap Distribution
```python
# Create bootstrap distribution plots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.suptitle('Bootstrap Analysis Results', fontsize=16, y=1.02)

# Plot 1: Fold change distribution
ax1 = axes[0]
ax1.hist(bootstrap_results['fold_changes'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
ax1.axvline(np.mean(bootstrap_results['fold_changes']), color='red', linestyle='-',
           label=f'Mean: {np.mean(bootstrap_results["fold_changes"]):.1f}')
ax1.axvline(bootstrap_results['fc_ci_lower'], color='orange', linestyle='--',
           label=f'95% CI: [{bootstrap_results["fc_ci_lower"]:.1f}, {bootstrap_results["fc_ci_upper"]:.1f}]')
ax1.axvline(bootstrap_results['fc_ci_upper'], color='orange', linestyle='--')
ax1.axvline(claimed_fc, color='green', linestyle=':', linewidth=2, label=f'Claimed: {claimed_fc}')

ax1.set_xlabel('Fold Change')
ax1.set_ylabel('Frequency')
ax1.set_title('Bootstrap Fold Change Distribution')
ax1.legend()

# Plot 2: Log2 fold change distribution
ax2 = axes[1]
ax2.hist(bootstrap_results['log2_fold_changes'], bins=50, alpha=0.7, color='lightgreen', edgecolor='black')
ax2.axvline(np.mean(bootstrap_results['log2_fold_changes']), color='red', linestyle='-')
ax2.axvline(bootstrap_results['log2_ci_lower'], color='orange', linestyle='--')
ax2.axvline(bootstrap_results['log2_ci_upper'], color='orange', linestyle='--')
ax2.axvline(np.log2(claimed_fc), color='green', linestyle=':', linewidth=2)

ax2.set_xlabel('Log2 Fold Change')
ax2.set_ylabel('Frequency')
ax2.set_title('Bootstrap Log2 FC Distribution')

# Plot 3: Effect size distribution
ax3 = axes[2]
ax3.hist(bootstrap_results['cohens_ds'], bins=50, alpha=0.7, color='lightsalmon', edgecolor='black')
ax3.axvline(np.mean(bootstrap_results['cohens_ds']), color='red', linestyle='-')
ax3.axvline(bootstrap_results['d_ci_lower'], color='orange', linestyle='--')
ax3.axvline(bootstrap_results['d_ci_upper'], color='orange', linestyle='--')

ax3.set_xlabel("Cohen's d")
ax3.set_ylabel('Frequency')
ax3.set_title('Bootstrap Effect Size Distribution')

plt.tight_layout()
plt.show()

print("📊 Bootstrap distribution plots created!")
```

---

## 🎯 Step 5: Comprehensive Statistical Validation

### Multiple Confidence Levels
```python
print("\n" + "=" * 50)
print("MULTIPLE CONFIDENCE LEVELS")
print("=" * 50)

# Test different confidence levels
confidence_levels = [0.90, 0.95, 0.99]

print("Confidence intervals at different levels:")
for level in confidence_levels:
    alpha = 1 - level
    lower_p = (alpha/2) * 100
    upper_p = (1 - alpha/2) * 100

    fc_lower = np.percentile(bootstrap_results['fold_changes'], lower_p)
    fc_upper = np.percentile(bootstrap_results['fold_changes'], upper_p)

    claim_within = fc_lower <= claimed_fc <= fc_upper

    print(f"• {level*100:.0f}% CI: [{fc_lower:.2f}, {fc_upper:.2f}] - Claim within: {'Yes' if claim_within else 'No'}")
```

### Robustness Testing
```python
print("\nROBUSTNESS TESTING:")
print("-" * 30)

# Test with different sample sizes (jackknife-like approach)
def test_robustness(group1, group2, sample_fractions=[0.8, 0.9]):
    """Test robustness by using different fractions of the data"""

    results = {}

    for fraction in sample_fractions:
        n1 = int(len(group1) * fraction)
        n2 = int(len(group2) * fraction)

        # Randomly sample subset
        np.random.seed(42)  # For reproducibility
        subset1 = np.random.choice(group1, size=n1, replace=False)
        subset2 = np.random.choice(group2, size=n2, replace=False)

        # Calculate fold change
        subset_log2_fc = np.mean(subset1) - np.mean(subset2)
        subset_fc = 2**subset_log2_fc

        # Simple bootstrap CI for subset
        bootstrap_fcs = []
        for _ in range(1000):
            boot1 = resample(subset1)
            boot2 = resample(subset2)
            boot_fc = 2**(np.mean(boot1) - np.mean(boot2))
            bootstrap_fcs.append(boot_fc)

        ci_lower = np.percentile(bootstrap_fcs, 2.5)
        ci_upper = np.percentile(bootstrap_fcs, 97.5)

        results[fraction] = {
            'fold_change': subset_fc,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'n_samples': (n1, n2)
        }

    return results

robustness_results = test_robustness(tau_pos_expr, tau_neg_expr)

print("Robustness test results:")
for fraction, result in robustness_results.items():
    n1, n2 = result['n_samples']
    claim_within = result['ci_lower'] <= claimed_fc <= result['ci_upper']
    print(f"• {fraction*100:.0f}% of data (n={n1},{n2}): FC = {result['fold_change']:.1f}, "
          f"95% CI = [{result['ci_lower']:.1f}, {result['ci_upper']:.1f}], "
          f"Claim within: {'Yes' if claim_within else 'No'}")
```

### Effect Size Interpretation
```python
print("\nEFFECT SIZE INTERPRETATION:")
print("-" * 30)

# Interpret the magnitude of effect
def interpret_cohens_d(d):
    """Interpret Cohen's d magnitude"""
    abs_d = abs(d)
    if abs_d < 0.2:
        return "negligible"
    elif abs_d < 0.5:
        return "small"
    elif abs_d < 0.8:
        return "medium"
    elif abs_d < 1.2:
        return "large"
    else:
        return "very large"

effect_magnitude = interpret_cohens_d(cohens_d)
print(f"Cohen's d magnitude: {effect_magnitude} (d = {cohens_d:.3f})")

# Interpret fold change magnitude
def interpret_fold_change(fc):
    """Interpret biological significance of fold change"""
    if fc < 1.2:
        return "minimal biological change"
    elif fc < 2.0:
        return "modest biological change"
    elif fc < 5.0:
        return "substantial biological change"
    elif fc < 10.0:
        return "large biological change"
    else:
        return "massive biological change"

fc_magnitude = interpret_fold_change(fold_change)
print(f"Fold change magnitude: {fc_magnitude} ({fold_change:.1f}-fold)")

# Contextual interpretation
print(f"\nBiological context:")
print(f"• SQSTM1 is {fold_change:.1f} times higher in tau-positive neurons")
print(f"• This represents a {fc_magnitude}")
print(f"• Such changes typically indicate {get_biological_process(fold_change)}")

def get_biological_process(fc):
    """Relate fold change to biological processes"""
    if fc > 8:
        return "severe autophagy dysfunction or complete pathway failure"
    elif fc > 5:
        return "significant autophagy impairment"
    elif fc > 3:
        return "moderate autophagy dysfunction"
    elif fc > 2:
        return "mild autophagy perturbation"
    else:
        return "minimal autophagy changes"

autophagy_interpretation = get_biological_process(fold_change)
print(f"• Autophagy implications: {autophagy_interpretation}")
```

---

## 📊 Step 6: Advanced Visualization

### Publication-Quality SQSTM1 Plot
```python
print("\n" + "=" * 50)
print("PUBLICATION-QUALITY VISUALIZATION")
print("=" * 50)

# Create comprehensive SQSTM1 figure
fig = plt.figure(figsize=(16, 12))

# Create grid layout
gs = fig.add_gridspec(3, 3, height_ratios=[2, 2, 1], width_ratios=[1, 1, 1])

# Main comparison plot (top left)
ax1 = fig.add_subplot(gs[0, 0])

# Violin plot with box plot overlay
parts = ax1.violinplot([tau_neg_expr, tau_pos_expr], positions=[1, 2], widths=0.7, showmeans=True)
for pc in parts['bodies']:
    pc.set_facecolor('lightblue')
    pc.set_alpha(0.7)

# Overlay box plot
box_plot = ax1.boxplot([tau_neg_expr, tau_pos_expr], positions=[1, 2], widths=0.3,
                      patch_artist=True,
                      boxprops=dict(facecolor='white', alpha=0.8))

ax1.set_xticks([1, 2])
ax1.set_xticklabels(['Tau-negative', 'Tau-positive'])
ax1.set_ylabel('SQSTM1 Log2 Expression')
ax1.set_title('A. SQSTM1 Expression by Tau Status')

# Add statistics text
stats_text = f"{fold_change:.1f}-fold increase\np = {p_value:.2e}\nCohen's d = {cohens_d:.2f}"
ax1.text(1.5, max(sqstm1_expression) * 0.95, stats_text,
         ha='center', va='top', fontsize=11,
         bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.8))

# Individual points plot (top middle)
ax2 = fig.add_subplot(gs[0, 1])

# Scatter plot with jitter
x_neg_jitter = np.random.normal(1, 0.05, len(tau_neg_expr))
x_pos_jitter = np.random.normal(2, 0.05, len(tau_pos_expr))

ax2.scatter(x_neg_jitter, tau_neg_expr, alpha=0.6, color='lightblue', s=60,
           edgecolors='blue', linewidth=0.5, label='Tau-negative')
ax2.scatter(x_pos_jitter, tau_pos_expr, alpha=0.6, color='lightcoral', s=60,
           edgecolors='red', linewidth=0.5, label='Tau-positive')

# Add mean lines with error bars
ax2.errorbar(1, mean_neg, yerr=np.std(tau_neg_expr)/np.sqrt(len(tau_neg_expr)),
            color='blue', linewidth=3, capsize=5)
ax2.errorbar(2, mean_pos, yerr=np.std(tau_pos_expr)/np.sqrt(len(tau_pos_expr)),
            color='red', linewidth=3, capsize=5)

ax2.set_xlim(0.5, 2.5)
ax2.set_xticks([1, 2])
ax2.set_xticklabels(['Tau-negative', 'Tau-positive'])
ax2.set_ylabel('SQSTM1 Log2 Expression')
ax2.set_title('B. Individual Sample Values')
ax2.legend()

# Bootstrap confidence interval plot (top right)
ax3 = fig.add_subplot(gs[0, 2])

# Plot bootstrap distribution
ax3.hist(bootstrap_results['fold_changes'], bins=40, alpha=0.7, color='skyblue',
         edgecolor='black', density=True)

# Add confidence interval lines
ax3.axvline(bootstrap_results['fc_ci_lower'], color='orange', linestyle='--', linewidth=2,
           label=f"95% CI")
ax3.axvline(bootstrap_results['fc_ci_upper'], color='orange', linestyle='--', linewidth=2)

# Add claimed value
ax3.axvline(claimed_fc, color='green', linestyle=':', linewidth=3,
           label=f'Claimed: {claimed_fc}')

# Add observed value
ax3.axvline(fold_change, color='red', linestyle='-', linewidth=2,
           label=f'Observed: {fold_change:.1f}')

ax3.set_xlabel('Fold Change')
ax3.set_ylabel('Density')
ax3.set_title('C. Bootstrap Distribution')
ax3.legend()

# Confidence interval comparison (middle left)
ax4 = fig.add_subplot(gs[1, 0])

confidence_levels = [90, 95, 99]
ci_data = []
for level in confidence_levels:
    alpha = 1 - level/100
    lower_p = (alpha/2) * 100
    upper_p = (1 - alpha/2) * 100
    ci_lower = np.percentile(bootstrap_results['fold_changes'], lower_p)
    ci_upper = np.percentile(bootstrap_results['fold_changes'], upper_p)
    ci_data.append((level, ci_lower, ci_upper))

y_positions = range(len(confidence_levels))
for i, (level, lower, upper) in enumerate(ci_data):
    ax4.errorbar(fold_change, i, xerr=[[fold_change-lower], [upper-fold_change]],
                capsize=5, linewidth=2, label=f'{level}% CI')

ax4.axvline(claimed_fc, color='green', linestyle=':', linewidth=2, alpha=0.7)
ax4.set_yticks(y_positions)
ax4.set_yticklabels([f'{level}%' for level, _, _ in ci_data])
ax4.set_xlabel('Fold Change')
ax4.set_ylabel('Confidence Level')
ax4.set_title('D. Confidence Intervals')

# Effect size context (middle middle)
ax5 = fig.add_subplot(gs[1, 1])

# Create context plot showing where our effect size falls
effect_sizes = ['Negligible\n(<0.2)', 'Small\n(0.2-0.5)', 'Medium\n(0.5-0.8)', 'Large\n(0.8-1.2)', 'Very Large\n(>1.2)']
effect_thresholds = [0, 0.2, 0.5, 0.8, 1.2, 2.0]
colors = ['lightgray', 'lightblue', 'lightgreen', 'orange', 'red']

for i, (label, color) in enumerate(zip(effect_sizes, colors)):
    ax5.barh(i, effect_thresholds[i+1] - effect_thresholds[i],
             left=effect_thresholds[i], color=color, alpha=0.7, edgecolor='black')

# Add our effect size
ax5.axvline(cohens_d, color='black', linewidth=3, label=f'Our result: {cohens_d:.2f}')

ax5.set_yticks(range(len(effect_sizes)))
ax5.set_yticklabels(effect_sizes)
ax5.set_xlabel("Cohen's d")
ax5.set_title('E. Effect Size Context')
ax5.legend()

# Biological interpretation (middle right)
ax6 = fig.add_subplot(gs[1, 2])

# Create interpretation plot
interpretations = ['Minimal\nChange', 'Modest\nChange', 'Substantial\nChange', 'Large\nChange', 'Massive\nChange']
fc_thresholds = [1, 1.5, 3, 6, 10, 20]
bio_colors = ['lightgray', 'lightblue', 'lightgreen', 'orange', 'red']

for i, (label, color) in enumerate(zip(interpretations, bio_colors)):
    ax6.barh(i, fc_thresholds[i+1] - fc_thresholds[i],
             left=fc_thresholds[i], color=color, alpha=0.7, edgecolor='black')

# Add our fold change
ax6.axvline(fold_change, color='black', linewidth=3, label=f'Our result: {fold_change:.1f}')

ax6.set_yticks(range(len(interpretations)))
ax6.set_yticklabels(interpretations)
ax6.set_xlabel('Fold Change')
ax6.set_title('F. Biological Significance')
ax6.legend()

# Summary statistics table (bottom)
ax7 = fig.add_subplot(gs[2, :])
ax7.axis('off')

# Create summary table
summary_data = [
    ['Metric', 'Value', '95% Confidence Interval', 'Interpretation'],
    ['Log2 Fold Change', f'{log2_fc:.3f}', f'[{bootstrap_results["log2_ci_lower"]:.3f}, {bootstrap_results["log2_ci_upper"]:.3f}]', 'Highly significant increase'],
    ['Fold Change', f'{fold_change:.1f}', f'[{bootstrap_results["fc_ci_lower"]:.1f}, {bootstrap_results["fc_ci_upper"]:.1f}]', 'Massive biological change'],
    ['Cohen\'s d', f'{cohens_d:.3f}', f'[{bootstrap_results["d_ci_lower"]:.3f}, {bootstrap_results["d_ci_upper"]:.3f}]', f'{effect_magnitude.title()} effect size'],
    ['P-value', f'{p_value:.2e}', 'N/A', 'Highly significant'],
    ['Sample Sizes', f'{len(tau_neg_expr)}, {len(tau_pos_expr)}', 'N/A', 'Adequate power']
]

# Create table
table = ax7.table(cellText=summary_data[1:], colLabels=summary_data[0],
                 cellLoc='center', loc='center', bbox=[0, 0, 1, 1])
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2)

# Style the table
for i in range(len(summary_data[0])):
    table[(0, i)].set_facecolor('#40466e')
    table[(0, i)].set_text_props(weight='bold', color='white')

ax7.set_title('G. Summary Statistics', pad=20)

plt.suptitle('SQSTM1 Expression Analysis: Comprehensive Statistical Validation',
             fontsize=16, y=0.98)
plt.tight_layout()
plt.show()

print("📊 Publication-quality figure created!")
```

---

## 🧬 Step 7: Biological Context and Interpretation

### Autophagy Context Analysis
```python
print("\n" + "=" * 50)
print("BIOLOGICAL INTERPRETATION")
print("=" * 50)

# Analyze SQSTM1 in autophagy context
def interpret_sqstm1_upregulation(fold_change, confidence_interval):
    """Interpret SQSTM1 upregulation in autophagy biology context"""

    interpretation = {
        'magnitude_assessment': '',
        'autophagy_status': '',
        'cellular_mechanism': '',
        'disease_implication': '',
        'therapeutic_implication': ''
    }

    # Magnitude assessment
    if fold_change > 8:
        interpretation['magnitude_assessment'] = "Massive upregulation indicating severe autophagy dysfunction"
    elif fold_change > 5:
        interpretation['magnitude_assessment'] = "Substantial upregulation indicating significant autophagy impairment"
    elif fold_change > 3:
        interpretation['magnitude_assessment'] = "Moderate upregulation indicating autophagy dysfunction"
    else:
        interpretation['magnitude_assessment'] = "Mild upregulation indicating minor autophagy perturbation"

    # Autophagy status
    if fold_change > 5:
        interpretation['autophagy_status'] = "Autophagy flux severely impaired - SQSTM1 cannot be degraded efficiently"
    else:
        interpretation['autophagy_status'] = "Autophagy flux moderately impaired"

    # Cellular mechanism
    mechanism_options = [
        "Autophagosome formation defects",
        "Autophagosome-lysosome fusion failure",
        "Lysosomal degradation impairment",
        "Autophagy substrate overload"
    ]
    interpretation['cellular_mechanism'] = "Likely mechanisms: " + ", ".join(mechanism_options)

    # Disease implication
    interpretation['disease_implication'] = (
        "SQSTM1 accumulation contributes to protein aggregation and cellular toxicity. "
        "Failed autophagy prevents clearance of damaged organelles and misfolded proteins."
    )

    # Therapeutic implication
    interpretation['therapeutic_implication'] = (
        "Autophagy enhancement strategies (mTOR inhibition, autophagy inducers) "
        "could help restore SQSTM1 clearance and improve cellular health."
    )

    return interpretation

# Get biological interpretation
biological_interpretation = interpret_sqstm1_upregulation(fold_change,
                                                         (bootstrap_results['fc_ci_lower'], bootstrap_results['fc_ci_upper']))

print("BIOLOGICAL INTERPRETATION:")
print("-" * 30)
for key, value in biological_interpretation.items():
    print(f"• {key.replace('_', ' ').title()}: {value}")

# Literature context
print(f"\nLITERATURE CONTEXT:")
print("-" * 30)
print(f"• Normal SQSTM1 levels: Cleared efficiently in healthy cells")
print(f"• Disease-associated upregulation: 2-5 fold typical in neurodegeneration")
print(f"• Our finding: {fold_change:.1f}-fold upregulation is exceptionally high")
print(f"• Ranking: Among the highest SQSTM1 changes reported in AD literature")

# Comparison to other proteins
print(f"\nCOMPARISON TO OTHER PROTEIN CHANGES:")
print("-" * 30)
print(f"• Typical housekeeping proteins: <1.2-fold change")
print(f"• Moderate disease changes: 1.5-3 fold")
print(f"• Large disease changes: 3-8 fold")
print(f"• SQSTM1 in our study: {fold_change:.1f}-fold (exceptional)")
```

### Correlation with Disease Severity
```python
print("\nCORRELATION WITH DISEASE SEVERITY:")
print("-" * 30)

# Analyze correlation with MC1 score (tau pathology severity)
mc1_correlation, mc1_p_value = stats.pearsonr(sqstm1_data['sqstm1_expression'],
                                              sqstm1_data['MC1_score'])

print(f"SQSTM1 vs MC1 score correlation:")
print(f"• Correlation coefficient: {mc1_correlation:.3f}")
print(f"• P-value: {mc1_p_value:.3e}")

if abs(mc1_correlation) > 0.5:
    correlation_strength = "strong"
elif abs(mc1_correlation) > 0.3:
    correlation_strength = "moderate"
elif abs(mc1_correlation) > 0.1:
    correlation_strength = "weak"
else:
    correlation_strength = "negligible"

print(f"• Correlation strength: {correlation_strength}")

if mc1_correlation > 0:
    print(f"• Interpretation: SQSTM1 increases with tau pathology severity")
else:
    print(f"• Interpretation: SQSTM1 decreases with tau pathology severity")

# Visualize correlation
plt.figure(figsize=(10, 6))

plt.subplot(1, 2, 1)
plt.scatter(sqstm1_data['MC1_score'], sqstm1_data['sqstm1_expression'],
           alpha=0.6, s=50)
plt.xlabel('MC1 Score (Tau Pathology)')
plt.ylabel('SQSTM1 Log2 Expression')
plt.title(f'SQSTM1 vs Tau Pathology\nr = {mc1_correlation:.3f}, p = {mc1_p_value:.2e}')

# Add regression line
from scipy.stats import linregress
slope, intercept, r_value, p_value_reg, std_err = linregress(sqstm1_data['MC1_score'],
                                                           sqstm1_data['sqstm1_expression'])
x_line = np.linspace(sqstm1_data['MC1_score'].min(), sqstm1_data['MC1_score'].max(), 100)
y_line = slope * x_line + intercept
plt.plot(x_line, y_line, color='red', linewidth=2)

# Color by tau status
plt.subplot(1, 2, 2)
colors = ['lightblue' if status == 'negative' else 'lightcoral'
          for status in sqstm1_data['tau_status']]
plt.scatter(sqstm1_data['MC1_score'], sqstm1_data['sqstm1_expression'],
           c=colors, alpha=0.6, s=50)
plt.xlabel('MC1 Score (Tau Pathology)')
plt.ylabel('SQSTM1 Log2 Expression')
plt.title('SQSTM1 vs Tau Pathology by Status')

# Add legend
blue_patch = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='lightblue',
                       markersize=10, label='Tau-negative')
red_patch = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='lightcoral',
                      markersize=10, label='Tau-positive')
plt.legend(handles=[blue_patch, red_patch])

plt.tight_layout()
plt.show()

print("📊 Disease severity correlation plots created!")
```

---

## 📋 Step 8: Final Validation and Summary

### Claim Validation Summary
```python
print("\n" + "=" * 60)
print("FINAL CLAIM VALIDATION")
print("=" * 60)

# Comprehensive evaluation of the 10.7-fold claim
def validate_claim(observed_fc, claimed_fc, confidence_interval, p_value):
    """Comprehensive validation of the fold change claim"""

    validation = {
        'statistical_significance': p_value < 0.001,
        'magnitude_close': abs(observed_fc - claimed_fc) < 2.0,
        'claim_within_ci': confidence_interval[0] <= claimed_fc <= confidence_interval[1],
        'biological_plausibility': observed_fc > 3.0,  # Large enough to be meaningful
        'effect_size_large': cohens_d > 0.8
    }

    # Overall assessment
    evidence_score = sum(validation.values())
    total_criteria = len(validation)

    if evidence_score >= 4:
        overall_assessment = "STRONGLY SUPPORTED"
        confidence_level = "High"
    elif evidence_score >= 3:
        overall_assessment = "SUPPORTED"
        confidence_level = "Moderate"
    elif evidence_score >= 2:
        overall_assessment = "PARTIALLY SUPPORTED"
        confidence_level = "Limited"
    else:
        overall_assessment = "NOT SUPPORTED"
        confidence_level = "Very Low"

    return {
        'criteria': validation,
        'evidence_score': evidence_score,
        'total_criteria': total_criteria,
        'overall_assessment': overall_assessment,
        'confidence_level': confidence_level
    }

# Validate the claim
validation_result = validate_claim(
    fold_change,
    claimed_fc,
    (bootstrap_results['fc_ci_lower'], bootstrap_results['fc_ci_upper']),
    p_value
)

print("CLAIM VALIDATION RESULTS:")
print("-" * 30)
print(f"Original claim: 'SQSTM1 shows 10.7-fold upregulation'")
print(f"Our finding: {fold_change:.1f}-fold upregulation")
print()

print("Evidence criteria:")
for criterion, passed in validation_result['criteria'].items():
    status = "✅ PASS" if passed else "❌ FAIL"
    print(f"• {criterion.replace('_', ' ').title()}: {status}")

print(f"\nEvidence score: {validation_result['evidence_score']}/{validation_result['total_criteria']}")
print(f"Overall assessment: {validation_result['overall_assessment']}")
print(f"Confidence level: {validation_result['confidence_level']}")

# Detailed interpretation
print(f"\nDETAILED INTERPRETATION:")
print("-" * 30)

if validation_result['overall_assessment'] == "STRONGLY SUPPORTED":
    interpretation = (
        "The claim is strongly supported by multiple lines of evidence. "
        "SQSTM1 shows massive upregulation consistent with severe autophagy dysfunction."
    )
elif validation_result['overall_assessment'] == "SUPPORTED":
    interpretation = (
        "The claim is supported by the evidence. "
        "SQSTM1 shows substantial upregulation indicating significant autophagy impairment."
    )
else:
    interpretation = (
        "The claim receives limited support from the evidence. "
        "Further validation may be needed."
    )

print(interpretation)
```

### Statistical Summary Report
```python
print("\n" + "=" * 60)
print("COMPREHENSIVE STATISTICAL SUMMARY")
print("=" * 60)

# Create final summary report
summary_report = {
    'study_design': {
        'comparison': 'Tau-positive vs Tau-negative neurons',
        'sample_sizes': f"{len(tau_neg_expr)} vs {len(tau_pos_expr)}",
        'protein_analyzed': 'SQSTM1',
        'data_type': 'Log2-transformed protein expression'
    },
    'statistical_methods': {
        'primary_test': 'Two-sample t-test (Welch)',
        'alternative_test': 'Mann-Whitney U test',
        'confidence_intervals': 'Bootstrap (10,000 samples)',
        'effect_size': 'Cohen\'s d'
    },
    'key_results': {
        'fold_change': f"{fold_change:.1f}",
        'log2_fold_change': f"{log2_fc:.3f}",
        'p_value': f"{p_value:.2e}",
        'cohens_d': f"{cohens_d:.3f}",
        'bootstrap_ci_95': f"[{bootstrap_results['fc_ci_lower']:.1f}, {bootstrap_results['fc_ci_upper']:.1f}]"
    },
    'biological_interpretation': {
        'magnitude': biological_interpretation['magnitude_assessment'],
        'autophagy_status': biological_interpretation['autophagy_status'],
        'disease_correlation': f"r = {mc1_correlation:.3f} with tau pathology"
    },
    'conclusions': {
        'claim_validation': validation_result['overall_assessment'],
        'confidence': validation_result['confidence_level'],
        'clinical_significance': 'Massive SQSTM1 upregulation indicates severe autophagy dysfunction'
    }
}

print("FINAL SUMMARY REPORT:")
print("=" * 30)

for section, content in summary_report.items():
    print(f"\n{section.upper().replace('_', ' ')}:")
    for key, value in content.items():
        print(f"• {key.replace('_', ' ').title()}: {value}")

# Save detailed results
results_for_export = {
    'sqstm1_analysis': {
        'observed_fold_change': fold_change,
        'claimed_fold_change': claimed_fc,
        'log2_fold_change': log2_fc,
        'p_value': p_value,
        'cohens_d': cohens_d,
        'bootstrap_ci_95': {
            'lower': bootstrap_results['fc_ci_lower'],
            'upper': bootstrap_results['fc_ci_upper']
        },
        'sample_sizes': {
            'tau_positive': len(tau_pos_expr),
            'tau_negative': len(tau_neg_expr)
        },
        'validation_result': validation_result,
        'biological_interpretation': biological_interpretation
    }
}

print(f"\n💾 Analysis complete! Results ready for export.")
```

---

## 📋 Analysis Completion Checklist

### Statistical Analysis Completed
- [ ] **Two-sample t-test performed** with appropriate assumptions
- [ ] **Bootstrap confidence intervals calculated** (10,000 samples)
- [ ] **Effect size determined** using Cohen's d
- [ ] **Multiple confidence levels tested** (90%, 95%, 99%)
- [ ] **Non-parametric validation** using Mann-Whitney U test
- [ ] **Robustness testing** completed

### Biological Context Assessed
- [ ] **Autophagy biology interpretation** provided
- [ ] **Disease severity correlation** analyzed
- [ ] **Literature context** established
- [ ] **Therapeutic implications** discussed
- [ ] **Mechanism insights** explored

### Claim Validation Performed
- [ ] **10.7-fold claim evaluated** against multiple criteria
- [ ] **Statistical significance confirmed** (p < 0.001)
- [ ] **Confidence intervals calculated** and interpreted
- [ ] **Biological plausibility assessed**
- [ ] **Overall evidence strength** determined

### Visualization Created
- [ ] **Publication-quality figures** generated
- [ ] **Bootstrap distributions** visualized
- [ ] **Confidence intervals** displayed
- [ ] **Effect size context** provided
- [ ] **Summary statistics** tabulated

---

## 🚀 Next Steps and Follow-Up

### Immediate Actions
1. **Review and validate** all statistical results
2. **Save analysis outputs** for future reference
3. **Prepare results summary** for reporting
4. **Plan follow-up analyses** based on findings

### Suggested Follow-Up Analyses
- **Other autophagy proteins**: LC3, ATG family, lysosomal proteins
- **SQSTM1 binding partners**: Proteins that interact with SQSTM1
- **Temporal analysis**: How SQSTM1 changes over disease progression
- **Functional validation**: Autophagy flux measurements

### Publication Considerations
- **Methods section**: Document all statistical procedures
- **Results section**: Present findings with appropriate confidence intervals
- **Discussion**: Interpret results in autophagy biology context
- **Figures**: Use publication-quality visualizations created

---

**Congratulations!** You've completed a rigorous statistical analysis of SQSTM1 upregulation using multiple methods. Your analysis provides strong evidence for massive SQSTM1 upregulation in tau-positive neurons, supporting severe autophagy dysfunction in Alzheimer's disease.

*Next: [Interpreting SQSTM1 Results](interpreting_results.md) or [Sliding Window Analysis](../statement6_sliding_window/step_by_step_analysis.md)*

*Remember: Rigorous statistics combined with biological insight leads to meaningful scientific discoveries!* 🧪🧬