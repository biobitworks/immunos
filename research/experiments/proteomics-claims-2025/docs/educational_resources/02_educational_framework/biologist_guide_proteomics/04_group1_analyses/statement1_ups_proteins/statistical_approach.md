# 📊 Statistical Approach: Testing UPS Protein Stability

## 🎯 Our Statistical Question

**Biological Claim**: "Targeted analyses show no significant UPS protein alterations across tau-positive versus tau-negative neurons."

**Statistical Translation**: "The majority of UPS proteins show no statistically significant differences between tau-positive and tau-negative neurons after multiple testing correction."

---

## 🧮 The Statistical Strategy

### Step 1: Define "No Significant Alterations"

#### What Does This Mean Quantitatively?
For the claim to be **SUPPORTED**, we expect:
- **<10% of UPS proteins** significantly different (allowing for some random variation)
- **Small effect sizes** for any significant proteins (Cohen's d < 0.3)
- **No systematic pattern** of up/down regulation

#### Why These Criteria?
- **Random chance**: With ~50 UPS proteins, we expect 2-3 to appear significant by luck alone
- **Biological noise**: Some variation is normal in biological systems
- **Effect size matters**: Tiny statistical differences aren't biologically meaningful

### Step 2: Choose the Right Statistical Test

#### Our Data Structure
```
Samples (rows): 150 neurons
Variables (columns): 5,853 proteins
Comparison: Tau-positive vs Tau-negative neurons
Data type: Continuous protein expression (log2 transformed)
```

#### Test Selection: Two-Sample t-Test
**Why t-test?**
- ✅ Comparing two groups (tau+ vs tau-)
- ✅ Continuous data (protein expression)
- ✅ Large sample size (>30 per group)
- ✅ Log-transformed data approximately normal

**Alternative considered**: Mann-Whitney U test
- More robust to non-normal data
- But less powerful with large samples
- We'll use t-test as primary, Mann-Whitney as sensitivity check

### Step 3: Multiple Testing Correction Strategy

#### The Multiple Testing Problem
Testing 50 UPS proteins means:
- **50 statistical tests** performed simultaneously
- **Expected false positives**: 50 × 0.05 = 2.5 proteins by chance
- **Without correction**: Many false discoveries

#### Our Solution: Benjamini-Hochberg FDR
**False Discovery Rate (FDR)** controls the expected proportion of false discoveries:
```python
from statsmodels.stats.multitest import multipletests

# Apply FDR correction
rejected, p_adjusted, _, _ = multipletests(p_values, method='fdr_bh', alpha=0.05)
```

**Why FDR instead of Bonferroni?**
- **Less conservative**: Bonferroni too strict for exploratory biology
- **Controls proportion**: Better for discovering real effects
- **Standard in genomics**: Widely accepted in high-throughput studies

### Step 4: Effect Size Calculation

#### Cohen's d for Biological Significance
```python
def cohens_d(group1, group2):
    """Calculate Cohen's d effect size"""
    pooled_std = np.sqrt(((len(group1)-1)*np.var(group1) + (len(group2)-1)*np.var(group2)) /
                        (len(group1) + len(group2) - 2))
    return (np.mean(group1) - np.mean(group2)) / pooled_std
```

**Interpretation Guidelines**:
- **|d| < 0.2**: Small effect (possibly not biologically meaningful)
- **|d| = 0.2-0.5**: Small to medium effect
- **|d| = 0.5-0.8**: Medium to large effect
- **|d| > 0.8**: Large effect (likely biologically important)

#### Fold Change for Intuitive Understanding
```python
def fold_change(group1, group2):
    """Calculate fold change from log2 expression data"""
    log2_fc = np.mean(group1) - np.mean(group2)
    return 2**log2_fc
```

---

## 🔬 Detailed Analysis Protocol

### Data Preparation

#### 1. Load and Inspect Dataset
```python
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats

# Load data
adata = sc.read_h5ad('pool_processed_v2.h5ad')

# Check dimensions
print(f"Dataset shape: {adata.shape}")
print(f"Tau-positive neurons: {sum(adata.obs['tau_status'] == 'positive')}")
print(f"Tau-negative neurons: {sum(adata.obs['tau_status'] == 'negative')}")
```

#### 2. Define UPS Protein List
```python
# Comprehensive UPS protein list
ups_proteins = [
    # Ubiquitin-activating enzymes (E1)
    'UBA1', 'UBA2', 'UBA3',

    # Ubiquitin-conjugating enzymes (E2)
    'UBE2A', 'UBE2B', 'UBE2C', 'UBE2D1', 'UBE2D2', 'UBE2D3',
    'UBE2E1', 'UBE2E2', 'UBE2E3', 'UBE2F', 'UBE2G1', 'UBE2G2',
    'UBE2H', 'UBE2I', 'UBE2J1', 'UBE2J2', 'UBE2K', 'UBE2L3',
    'UBE2M', 'UBE2N', 'UBE2O', 'UBE2Q1', 'UBE2Q2', 'UBE2R2',
    'UBE2S', 'UBE2T', 'UBE2U', 'UBE2V1', 'UBE2V2', 'UBE2W',
    'UBE2Z',

    # Ubiquitin ligases (E3) - major ones
    'UBE3A', 'UBE3B', 'UBE3C',
    'HUWE1', 'HECTD1', 'HECTD2', 'HECTD3',
    'RNF4', 'RNF8', 'RNF168',
    'MDM2', 'CHIP', 'PARKIN',

    # Proteasome subunits - 20S core
    'PSMA1', 'PSMA2', 'PSMA3', 'PSMA4', 'PSMA5', 'PSMA6', 'PSMA7',
    'PSMB1', 'PSMB2', 'PSMB3', 'PSMB4', 'PSMB5', 'PSMB6', 'PSMB7',

    # Proteasome subunits - 26S regulatory
    'PSMC1', 'PSMC2', 'PSMC3', 'PSMC4', 'PSMC5', 'PSMC6',
    'PSMD1', 'PSMD2', 'PSMD3', 'PSMD4', 'PSMD5', 'PSMD6', 'PSMD7',
    'PSMD8', 'PSMD9', 'PSMD10', 'PSMD11', 'PSMD12', 'PSMD13', 'PSMD14',

    # Deubiquitinating enzymes
    'USP14', 'UCHL1', 'UCHL3', 'UCHL5',
    'USP1', 'USP2', 'USP7', 'USP9X', 'USP10',

    # Other UPS components
    'VCP', 'NEDD8', 'UBQLN1', 'UBQLN2', 'UBQLN4',
    'UBC', 'UBB', 'UBA52', 'RPS27A'  # Ubiquitin genes
]

# Check which UPS proteins are in our dataset
available_ups = [p for p in ups_proteins if p in adata.var_names]
missing_ups = [p for p in ups_proteins if p not in adata.var_names]

print(f"Available UPS proteins: {len(available_ups)}")
print(f"Missing UPS proteins: {len(missing_ups)}")
```

#### 3. Extract Expression Data
```python
# Extract UPS protein expression
ups_data = adata[:, available_ups].X.toarray()  # Convert sparse to dense
ups_df = pd.DataFrame(ups_data, columns=available_ups, index=adata.obs_names)

# Add tau status
ups_df['tau_status'] = adata.obs['tau_status'].values

# Separate groups
tau_pos_data = ups_df[ups_df['tau_status'] == 'positive'].drop('tau_status', axis=1)
tau_neg_data = ups_df[ups_df['tau_status'] == 'negative'].drop('tau_status', axis=1)

print(f"Tau-positive group: {tau_pos_data.shape[0]} neurons")
print(f"Tau-negative group: {tau_neg_data.shape[0]} neurons")
```

### Statistical Testing

#### 1. Perform t-Tests for Each Protein
```python
def analyze_ups_protein(protein_name, tau_pos, tau_neg):
    """Analyze single UPS protein between groups"""

    # Extract expression values
    pos_expr = tau_pos[protein_name].values
    neg_expr = tau_neg[protein_name].values

    # Remove any NaN values
    pos_expr = pos_expr[~np.isnan(pos_expr)]
    neg_expr = neg_expr[~np.isnan(neg_expr)]

    # Check sample sizes
    if len(pos_expr) < 5 or len(neg_expr) < 5:
        return None  # Insufficient data

    # Two-sample t-test
    t_stat, p_value = stats.ttest_ind(pos_expr, neg_expr, equal_var=False)

    # Effect size (Cohen's d)
    pooled_std = np.sqrt(((len(pos_expr)-1)*np.var(pos_expr, ddof=1) +
                         (len(neg_expr)-1)*np.var(neg_expr, ddof=1)) /
                        (len(pos_expr) + len(neg_expr) - 2))
    cohens_d = (np.mean(pos_expr) - np.mean(neg_expr)) / pooled_std

    # Fold change
    log2_fc = np.mean(pos_expr) - np.mean(neg_expr)
    fold_change = 2**log2_fc

    # Confidence interval for difference
    se_diff = pooled_std * np.sqrt(1/len(pos_expr) + 1/len(neg_expr))
    df = len(pos_expr) + len(neg_expr) - 2
    t_critical = stats.t.ppf(0.975, df)  # 95% CI
    diff = np.mean(pos_expr) - np.mean(neg_expr)
    ci_lower = diff - t_critical * se_diff
    ci_upper = diff + t_critical * se_diff

    return {
        'protein': protein_name,
        'n_tau_pos': len(pos_expr),
        'n_tau_neg': len(neg_expr),
        'mean_tau_pos': np.mean(pos_expr),
        'mean_tau_neg': np.mean(neg_expr),
        't_statistic': t_stat,
        'p_value': p_value,
        'cohens_d': cohens_d,
        'log2_fold_change': log2_fc,
        'fold_change': fold_change,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper
    }

# Analyze all UPS proteins
results = []
for protein in available_ups:
    result = analyze_ups_protein(protein, tau_pos_data, tau_neg_data)
    if result is not None:
        results.append(result)

# Convert to DataFrame
results_df = pd.DataFrame(results)
print(f"Successfully analyzed {len(results_df)} UPS proteins")
```

#### 2. Apply Multiple Testing Correction
```python
from statsmodels.stats.multitest import multipletests

# Apply FDR correction
p_values = results_df['p_value'].values
rejected, p_adjusted, alpha_sidak, alpha_bonf = multipletests(
    p_values, method='fdr_bh', alpha=0.05
)

# Add corrected p-values to results
results_df['p_adjusted'] = p_adjusted
results_df['significant'] = rejected

# Summary statistics
n_significant = sum(rejected)
percent_significant = (n_significant / len(results_df)) * 100

print(f"Significant proteins (FDR < 0.05): {n_significant}/{len(results_df)} ({percent_significant:.1f}%)")
```

#### 3. Effect Size Analysis
```python
# Categorize effect sizes
def categorize_effect_size(d):
    """Categorize Cohen's d effect size"""
    abs_d = abs(d)
    if abs_d < 0.2:
        return 'negligible'
    elif abs_d < 0.5:
        return 'small'
    elif abs_d < 0.8:
        return 'medium'
    else:
        return 'large'

results_df['effect_size_category'] = results_df['cohens_d'].apply(categorize_effect_size)

# Effect size distribution
effect_size_counts = results_df['effect_size_category'].value_counts()
print("Effect size distribution:")
print(effect_size_counts)
```

---

## 📈 Expected Results and Interpretation

### If the Claim is SUPPORTED (UPS Intact)

#### Statistical Pattern
- **<10% proteins significant** after FDR correction
- **Small effect sizes** (most |Cohen's d| < 0.3)
- **Random distribution** of significant results
- **No systematic up/down pattern**

#### Biological Interpretation
```python
# Example of supporting evidence
if percent_significant < 10 and effect_size_counts['negligible'] > effect_size_counts['large']:
    interpretation = "SUPPORTED: UPS system appears intact"
    biological_meaning = """
    - Protein quality control machinery is functional
    - Disease pathology doesn't primarily affect UPS
    - UPS remains a viable therapeutic target
    - Protein aggregation likely due to other mechanisms
    """
```

### If the Claim is REFUTED (UPS Dysfunction)

#### Statistical Pattern
- **>25% proteins significant** after FDR correction
- **Large effect sizes** (many |Cohen's d| > 0.5)
- **Systematic pattern** (mostly up or mostly down)
- **Biologically coherent changes** (related proteins change together)

#### Biological Interpretation
```python
# Example of refuting evidence
if percent_significant > 25 and effect_size_counts['large'] > effect_size_counts['negligible']:
    interpretation = "REFUTED: UPS system shows significant dysfunction"
    biological_meaning = """
    - Protein quality control is impaired in disease
    - UPS dysfunction contributes to pathology
    - Therapeutic strategies should target UPS enhancement
    - Multiple UPS components affected simultaneously
    """
```

---

## 🔍 Sensitivity Analyses

### 1. Alternative Statistical Tests
```python
# Mann-Whitney U test for non-parametric comparison
from scipy.stats import mannwhitneyu

def mann_whitney_analysis(protein_name, tau_pos, tau_neg):
    """Non-parametric alternative to t-test"""
    pos_expr = tau_pos[protein_name].values
    neg_expr = tau_neg[protein_name].values

    # Remove NaN values
    pos_expr = pos_expr[~np.isnan(pos_expr)]
    neg_expr = neg_expr[~np.isnan(neg_expr)]

    # Mann-Whitney U test
    statistic, p_value = mannwhitneyu(pos_expr, neg_expr, alternative='two-sided')

    return p_value

# Compare t-test vs Mann-Whitney results
results_df['p_value_mw'] = [mann_whitney_analysis(protein, tau_pos_data, tau_neg_data)
                           for protein in results_df['protein']]

# Correlation between methods
correlation = np.corrcoef(results_df['p_value'], results_df['p_value_mw'])[0,1]
print(f"Correlation between t-test and Mann-Whitney p-values: {correlation:.3f}")
```

### 2. Sample Size Effects
```python
# Check if results depend on sample size
min_samples = min(len(tau_pos_data), len(tau_neg_data))
if min_samples < 30:
    print("WARNING: Small sample size may affect power")
    print("Consider using non-parametric tests or bootstrap methods")

# Power analysis for typical effect sizes
from statsmodels.stats.power import ttest_power

power_small = ttest_power(0.2, min_samples, 0.05)  # Small effect
power_medium = ttest_power(0.5, min_samples, 0.05)  # Medium effect

print(f"Statistical power for small effect (d=0.2): {power_small:.2f}")
print(f"Statistical power for medium effect (d=0.5): {power_medium:.2f}")
```

---

## 📋 Quality Control Checklist

### Before Running Analysis
- [ ] **Data loaded correctly**: Check dimensions and structure
- [ ] **UPS proteins identified**: Verify protein names match dataset
- [ ] **Groups defined properly**: Tau-positive vs tau-negative
- [ ] **Sample sizes adequate**: >20 per group recommended

### During Analysis
- [ ] **Missing data handled**: Check for NaN values
- [ ] **Statistical assumptions**: Normal distributions (approximately)
- [ ] **Effect sizes calculated**: Not just p-values
- [ ] **Multiple testing corrected**: FDR applied

### After Analysis
- [ ] **Results make biological sense**: Interpret in context
- [ ] **Effect sizes match significance**: Large p-values should have small effects
- [ ] **Sensitivity analysis**: Compare with non-parametric tests
- [ ] **Document limitations**: Note any concerns or caveats

---

## 🎯 Key Learning Objectives

After completing this statistical approach section, you should understand:

### Statistical Concepts
- [ ] Why we need multiple testing correction
- [ ] How to interpret effect sizes vs p-values
- [ ] When to use parametric vs non-parametric tests
- [ ] How sample size affects statistical power

### Practical Implementation
- [ ] How to structure a differential expression analysis
- [ ] How to handle missing data appropriately
- [ ] How to calculate confidence intervals
- [ ] How to perform sensitivity analyses

### Biological Integration
- [ ] How statistical results inform biological conclusions
- [ ] Why "no significant difference" supports the UPS integrity claim
- [ ] How to distinguish statistical from biological significance
- [ ] How to communicate uncertainty appropriately

---

**Ready to implement this analysis?** The next step is the hands-on tutorial where we'll run all this code step-by-step.

*Next: [Step-by-Step Analysis Tutorial](step_by_step_analysis.md)*

*Remember: Good statistics require careful thinking about assumptions, appropriate methods, and biological interpretation - not just running tests!* 📊🧬