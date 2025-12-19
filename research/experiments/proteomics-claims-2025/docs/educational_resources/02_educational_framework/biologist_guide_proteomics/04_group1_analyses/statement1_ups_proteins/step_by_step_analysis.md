# 🧪 Step-by-Step Analysis: UPS Protein Stability

## 🎯 What We're Going to Do

By the end of this tutorial, you'll have:
- ✅ Loaded and explored the proteomic dataset
- ✅ Identified UPS proteins in the data
- ✅ Performed statistical tests comparing tau-positive vs tau-negative neurons
- ✅ Applied multiple testing correction
- ✅ Calculated effect sizes and confidence intervals
- ✅ Interpreted results in biological context

**Time needed**: 2-3 hours for complete beginners, 1 hour for some experience

---

## 🛠️ Before We Start

### Prerequisites Checklist
- [ ] Python and Jupyter installed (see software setup guide)
- [ ] Dataset downloaded (`pool_processed_v2.h5ad`)
- [ ] Required packages installed
- [ ] Understanding of basic UPS biology (see biological background)

### Setup Your Workspace
```bash
# Create analysis directory
mkdir ups_protein_analysis
cd ups_protein_analysis

# Start Jupyter notebook
jupyter notebook
```

---

## 📂 Step 1: Load Required Libraries and Data

### Import Essential Packages
Create a new Jupyter notebook and run this cell:

```python
# Core data science libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Statistical libraries
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.power import ttest_power

# Bioinformatics libraries
import scanpy as sc

# Utilities
import warnings
warnings.filterwarnings('ignore')

# Set up plotting
plt.style.use('default')
sns.set_palette("colorblind")
%matplotlib inline

print("✅ All libraries imported successfully!")
```

**Troubleshooting**: If you get import errors, refer to the software setup guide or run:
```bash
pip install pandas numpy matplotlib seaborn scipy statsmodels scanpy
```

### Load the Dataset
```python
# Load the proteomic dataset
print("Loading dataset...")
adata = sc.read_h5ad('../../../data/pool_processed_v2.h5ad')

print(f"✅ Dataset loaded successfully!")
print(f"📊 Dataset shape: {adata.shape} (neurons × proteins)")
print(f"🧠 Number of neurons: {adata.n_obs}")
print(f"🧬 Number of proteins: {adata.n_vars}")
```

### Explore the Data Structure
```python
# Check what information we have about each neuron
print("📋 Available metadata (first 5 columns):")
print(adata.obs.head())

print("\n🏷️ Available metadata columns:")
print(list(adata.obs.columns))

print("\n🎯 Tau status distribution:")
print(adata.obs['tau_status'].value_counts())
```

**What you should see**:
- Dataset with ~150 neurons and ~5,853 proteins
- Metadata including tau_status, MC1_score, pseudotime
- Two groups: tau-positive and tau-negative neurons

---

## 🔍 Step 2: Identify UPS Proteins in the Dataset

### Define UPS Protein List
```python
# Comprehensive list of UPS (Ubiquitin-Proteasome System) proteins
ups_proteins = [
    # Ubiquitin-activating enzymes (E1)
    'UBA1', 'UBA2', 'UBA3',

    # Ubiquitin-conjugating enzymes (E2)
    'UBE2A', 'UBE2B', 'UBE2C', 'UBE2D1', 'UBE2D2', 'UBE2D3',
    'UBE2E1', 'UBE2E2', 'UBE2E3', 'UBE2F', 'UBE2G1', 'UBE2G2',
    'UBE2H', 'UBE2I', 'UBE2J1', 'UBE2J2', 'UBE2K', 'UBE2L3',
    'UBE2M', 'UBE2N', 'UBE2O', 'UBE2Q1', 'UBE2Q2', 'UBE2R2',
    'UBE2S', 'UBE2T', 'UBE2U', 'UBE2V1', 'UBE2V2', 'UBE2W', 'UBE2Z',

    # Ubiquitin ligases (E3) - major ones
    'UBE3A', 'UBE3B', 'UBE3C',
    'HUWE1', 'HECTD1', 'HECTD2', 'HECTD3',
    'RNF4', 'RNF8', 'RNF168',
    'MDM2', 'PARKIN',

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

print(f"📋 Total UPS proteins in our list: {len(ups_proteins)}")
```

### Check Which UPS Proteins Are Available
```python
# Check which UPS proteins are actually in our dataset
available_ups = []
missing_ups = []

for protein in ups_proteins:
    if protein in adata.var_names:
        available_ups.append(protein)
    else:
        missing_ups.append(protein)

print(f"✅ Available UPS proteins: {len(available_ups)}")
print(f"❌ Missing UPS proteins: {len(missing_ups)}")
print(f"📊 Coverage: {len(available_ups)/len(ups_proteins)*100:.1f}%")

print(f"\n🔍 Available UPS proteins:")
for i, protein in enumerate(available_ups):
    if i % 6 == 0:  # New line every 6 proteins for readability
        print()
    print(f"{protein:8}", end="")

if missing_ups:
    print(f"\n\n❓ Missing UPS proteins (first 10):")
    print(missing_ups[:10])
```

**Understanding the output**:
- You should have 40-60 available UPS proteins (depends on dataset)
- Some proteins might be missing due to low expression or detection limits
- This is normal in proteomic studies

---

## 📊 Step 3: Extract and Prepare UPS Expression Data

### Create UPS-Focused Dataset
```python
# Extract only UPS proteins
print("🔬 Extracting UPS protein expression data...")

# Get expression data for UPS proteins only
ups_adata = adata[:, available_ups].copy()

# Convert sparse matrix to dense for easier manipulation
ups_expression = ups_adata.X.toarray()

# Create a DataFrame for easier analysis
ups_df = pd.DataFrame(
    ups_expression,
    columns=available_ups,
    index=adata.obs_names
)

# Add metadata
ups_df['tau_status'] = adata.obs['tau_status'].values
ups_df['MC1_score'] = adata.obs['MC1_score'].values
ups_df['patient_id'] = adata.obs['patient_id'].values

print(f"✅ UPS dataset created: {ups_df.shape}")
print(f"📊 Columns: {ups_df.shape[1]-3} proteins + 3 metadata columns")
```

### Explore UPS Expression Patterns
```python
# Basic statistics for UPS proteins
protein_cols = [col for col in ups_df.columns if col not in ['tau_status', 'MC1_score', 'patient_id']]

print("📈 UPS protein expression summary:")
print(ups_df[protein_cols].describe().round(2))

# Check for missing values
missing_data = ups_df[protein_cols].isnull().sum()
if missing_data.sum() > 0:
    print(f"\n⚠️ Missing data found in {missing_data.sum()} entries")
    print(missing_data[missing_data > 0])
else:
    print("\n✅ No missing data in UPS proteins")
```

### Visualize Overall UPS Expression
```python
# Create a heatmap of UPS protein expression
plt.figure(figsize=(15, 8))

# Select a subset of proteins for visualization (if too many)
n_proteins_to_show = min(20, len(protein_cols))
proteins_to_plot = protein_cols[:n_proteins_to_show]

# Prepare data for heatmap
heatmap_data = ups_df[proteins_to_plot + ['tau_status']].copy()
heatmap_data = heatmap_data.sort_values('tau_status')  # Group by tau status

# Create heatmap
sns.heatmap(
    heatmap_data[proteins_to_plot].T,  # Transpose so proteins are rows
    cmap='RdBu_r',
    center=0,
    cbar_kws={'label': 'Log2 Expression'},
    xticklabels=False,  # Too many samples to show labels
    yticklabels=True
)

plt.title(f'UPS Protein Expression Heatmap\n(First {n_proteins_to_show} proteins)')
plt.ylabel('UPS Proteins')
plt.xlabel('Neurons (sorted by tau status)')
plt.tight_layout()
plt.show()

print("💡 In this heatmap:")
print("- Each row is a UPS protein")
print("- Each column is a neuron")
print("- Red = higher expression, Blue = lower expression")
print("- Neurons are sorted by tau status")
```

---

## 🔬 Step 4: Statistical Analysis - Comparing Groups

### Separate Groups for Analysis
```python
# Separate tau-positive and tau-negative groups
tau_positive = ups_df[ups_df['tau_status'] == 'positive']
tau_negative = ups_df[ups_df['tau_status'] == 'negative']

print(f"🔴 Tau-positive group: {len(tau_positive)} neurons")
print(f"🔵 Tau-negative group: {len(tau_negative)} neurons")

# Check group balance
total_neurons = len(tau_positive) + len(tau_negative)
pos_percent = len(tau_positive) / total_neurons * 100
neg_percent = len(tau_negative) / total_neurons * 100

print(f"📊 Group distribution: {pos_percent:.1f}% positive, {neg_percent:.1f}% negative")

if min(len(tau_positive), len(tau_negative)) < 20:
    print("⚠️ WARNING: Small group size may reduce statistical power")
```

### Define Analysis Function
```python
def analyze_ups_protein(protein_name, tau_pos_df, tau_neg_df):
    """
    Perform comprehensive statistical analysis for a single UPS protein

    Parameters:
    -----------
    protein_name : str
        Name of the protein to analyze
    tau_pos_df : DataFrame
        Data for tau-positive neurons
    tau_neg_df : DataFrame
        Data for tau-negative neurons

    Returns:
    --------
    dict : Analysis results
    """

    # Extract expression values
    pos_expr = tau_pos_df[protein_name].values
    neg_expr = tau_neg_df[protein_name].values

    # Remove any NaN values (just in case)
    pos_expr = pos_expr[~np.isnan(pos_expr)]
    neg_expr = neg_expr[~np.isnan(neg_expr)]

    # Check if we have enough data
    if len(pos_expr) < 3 or len(neg_expr) < 3:
        print(f"⚠️ Warning: Insufficient data for {protein_name}")
        return None

    # Basic descriptive statistics
    pos_mean = np.mean(pos_expr)
    neg_mean = np.mean(neg_expr)
    pos_std = np.std(pos_expr, ddof=1)
    neg_std = np.std(neg_expr, ddof=1)

    # Two-sample t-test (assuming unequal variances)
    t_stat, p_value = stats.ttest_ind(pos_expr, neg_expr, equal_var=False)

    # Effect size (Cohen's d)
    pooled_std = np.sqrt(((len(pos_expr)-1)*pos_std**2 + (len(neg_expr)-1)*neg_std**2) /
                        (len(pos_expr) + len(neg_expr) - 2))
    cohens_d = (pos_mean - neg_mean) / pooled_std

    # Fold change
    log2_fc = pos_mean - neg_mean  # Data is already log2 transformed
    fold_change = 2**log2_fc

    # 95% Confidence interval for the difference
    se_diff = pooled_std * np.sqrt(1/len(pos_expr) + 1/len(neg_expr))
    df = len(pos_expr) + len(neg_expr) - 2
    t_critical = stats.t.ppf(0.975, df)
    diff = pos_mean - neg_mean
    ci_lower = diff - t_critical * se_diff
    ci_upper = diff + t_critical * se_diff

    return {
        'protein': protein_name,
        'n_tau_pos': len(pos_expr),
        'n_tau_neg': len(neg_expr),
        'mean_tau_pos': pos_mean,
        'std_tau_pos': pos_std,
        'mean_tau_neg': neg_mean,
        'std_tau_neg': neg_std,
        't_statistic': t_stat,
        'p_value': p_value,
        'cohens_d': cohens_d,
        'log2_fold_change': log2_fc,
        'fold_change': fold_change,
        'difference': diff,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'degrees_freedom': df
    }

print("✅ Analysis function defined successfully!")
```

### Run Statistical Tests for All UPS Proteins
```python
# Analyze all available UPS proteins
print("🧮 Running statistical tests for all UPS proteins...")

results = []
failed_proteins = []

for i, protein in enumerate(protein_cols):
    print(f"📊 Analyzing {protein} ({i+1}/{len(protein_cols)})", end="")

    result = analyze_ups_protein(protein, tau_positive, tau_negative)

    if result is not None:
        results.append(result)
        print(" ✅")
    else:
        failed_proteins.append(protein)
        print(" ❌")

print(f"\n✅ Analysis complete!")
print(f"📊 Successfully analyzed: {len(results)} proteins")
print(f"❌ Failed analysis: {len(failed_proteins)} proteins")

if failed_proteins:
    print(f"Failed proteins: {failed_proteins}")
```

### Convert Results to DataFrame
```python
# Convert results to DataFrame for easier manipulation
results_df = pd.DataFrame(results)

# Display first few results
print("📋 First 5 analysis results:")
print(results_df[['protein', 'p_value', 'cohens_d', 'fold_change']].head())

print(f"\n📊 Analysis summary:")
print(f"Mean p-value: {results_df['p_value'].mean():.4f}")
print(f"Median p-value: {results_df['p_value'].median():.4f}")
print(f"Mean |Cohen's d|: {results_df['cohens_d'].abs().mean():.3f}")
print(f"Mean fold change: {results_df['fold_change'].mean():.3f}")
```

---

## 🎯 Step 5: Multiple Testing Correction

### Apply False Discovery Rate (FDR) Correction
```python
# Apply Benjamini-Hochberg FDR correction
print("🔧 Applying multiple testing correction...")

p_values = results_df['p_value'].values
rejected, p_adjusted, alpha_sidak, alpha_bonf = multipletests(
    p_values,
    method='fdr_bh',  # Benjamini-Hochberg FDR
    alpha=0.05
)

# Add corrected p-values to results
results_df['p_adjusted'] = p_adjusted
results_df['significant_fdr'] = rejected

# Summary of significance
n_significant_uncorrected = sum(results_df['p_value'] < 0.05)
n_significant_corrected = sum(rejected)

print(f"📊 Significance summary:")
print(f"Significant before correction (p < 0.05): {n_significant_uncorrected}/{len(results_df)} ({n_significant_uncorrected/len(results_df)*100:.1f}%)")
print(f"Significant after FDR correction: {n_significant_corrected}/{len(results_df)} ({n_significant_corrected/len(results_df)*100:.1f}%)")

# Show significant proteins
if n_significant_corrected > 0:
    significant_proteins = results_df[results_df['significant_fdr']].copy()
    significant_proteins = significant_proteins.sort_values('p_adjusted')

    print(f"\n🎯 Significant UPS proteins (FDR < 0.05):")
    print(significant_proteins[['protein', 'p_value', 'p_adjusted', 'cohens_d', 'fold_change']].to_string(index=False, float_format='%.4f'))
else:
    print("\n✅ No UPS proteins show significant differences after FDR correction")
```

### Visualize P-Value Distribution
```python
# Create p-value histogram
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

# Uncorrected p-values
ax1.hist(results_df['p_value'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
ax1.axvline(0.05, color='red', linestyle='--', label='p = 0.05')
ax1.set_xlabel('Uncorrected P-values')
ax1.set_ylabel('Frequency')
ax1.set_title('Distribution of Uncorrected P-values')
ax1.legend()

# Corrected p-values
ax2.hist(results_df['p_adjusted'], bins=20, alpha=0.7, color='lightcoral', edgecolor='black')
ax2.axvline(0.05, color='red', linestyle='--', label='FDR = 0.05')
ax2.set_xlabel('FDR-Corrected P-values')
ax2.set_ylabel('Frequency')
ax2.set_title('Distribution of FDR-Corrected P-values')
ax2.legend()

plt.tight_layout()
plt.show()

print("💡 Understanding the histograms:")
print("- Left: Original p-values (many may be < 0.05 by chance)")
print("- Right: FDR-corrected p-values (controls false discovery rate)")
print("- If no correction was needed, both plots would look similar")
```

---

## 📏 Step 6: Effect Size Analysis

### Categorize Effect Sizes
```python
def categorize_effect_size(cohens_d):
    """Categorize Cohen's d effect size"""
    abs_d = abs(cohens_d)
    if abs_d < 0.2:
        return 'negligible'
    elif abs_d < 0.5:
        return 'small'
    elif abs_d < 0.8:
        return 'medium'
    else:
        return 'large'

# Add effect size categories
results_df['effect_size_category'] = results_df['cohens_d'].apply(categorize_effect_size)

# Count effect size categories
effect_size_counts = results_df['effect_size_category'].value_counts()
print("📊 Effect size distribution:")
print(effect_size_counts)

# Calculate percentages
effect_size_percentages = (effect_size_counts / len(results_df) * 100).round(1)
print("\n📊 Effect size percentages:")
for category, percentage in effect_size_percentages.items():
    print(f"{category:>12}: {percentage:>5.1f}%")
```

### Visualize Effect Sizes
```python
# Create effect size visualization
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Effect size histogram
ax1.hist(results_df['cohens_d'], bins=20, alpha=0.7, color='green', edgecolor='black')
ax1.axvline(0, color='black', linestyle='-', alpha=0.3)
ax1.axvline(-0.2, color='orange', linestyle='--', alpha=0.7, label='Small effect')
ax1.axvline(0.2, color='orange', linestyle='--', alpha=0.7)
ax1.axvline(-0.5, color='red', linestyle='--', alpha=0.7, label='Medium effect')
ax1.axvline(0.5, color='red', linestyle='--', alpha=0.7)
ax1.set_xlabel("Cohen's d")
ax1.set_ylabel('Frequency')
ax1.set_title('Distribution of Effect Sizes')
ax1.legend()

# Effect size categories pie chart
colors = ['lightblue', 'gold', 'orange', 'red']
ax2.pie(effect_size_counts.values, labels=effect_size_counts.index,
        autopct='%1.1f%%', colors=colors[:len(effect_size_counts)])
ax2.set_title('Effect Size Categories')

plt.tight_layout()
plt.show()
```

### Volcano Plot: P-values vs Effect Sizes
```python
# Create volcano plot
plt.figure(figsize=(10, 6))

# Color points based on significance and effect size
colors = []
for _, row in results_df.iterrows():
    if row['significant_fdr'] and abs(row['cohens_d']) > 0.5:
        colors.append('red')  # Significant and large effect
    elif row['significant_fdr']:
        colors.append('orange')  # Significant but small effect
    elif abs(row['cohens_d']) > 0.5:
        colors.append('blue')  # Large effect but not significant
    else:
        colors.append('gray')  # Neither significant nor large effect

plt.scatter(results_df['cohens_d'], -np.log10(results_df['p_adjusted']),
           c=colors, alpha=0.7, s=50)

# Add significance threshold line
plt.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='FDR = 0.05')

# Add effect size threshold lines
plt.axvline(-0.5, color='blue', linestyle='--', alpha=0.7, label='Medium effect')
plt.axvline(0.5, color='blue', linestyle='--', alpha=0.7)

plt.xlabel("Cohen's d (Effect Size)")
plt.ylabel('-log10(FDR-adjusted p-value)')
plt.title('Volcano Plot: UPS Protein Changes')
plt.legend()

# Add labels for significant proteins
if n_significant_corrected > 0:
    for _, row in results_df[results_df['significant_fdr']].iterrows():
        plt.annotate(row['protein'],
                    (row['cohens_d'], -np.log10(row['p_adjusted'])),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=8, alpha=0.8)

plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

print("🎯 Volcano plot interpretation:")
print("- X-axis: Effect size (Cohen's d)")
print("- Y-axis: Statistical significance (-log10 p-value)")
print("- Red points: Significant AND large effect")
print("- Orange points: Significant but small effect")
print("- Blue points: Large effect but not significant")
print("- Gray points: Neither significant nor large effect")
```

---

## 📋 Step 7: Interpret Results

### Determine Claim Evaluation
```python
# Evaluate the biological claim based on our results
percent_significant = (n_significant_corrected / len(results_df)) * 100
large_effects = sum(results_df['effect_size_category'] == 'large')
medium_effects = sum(results_df['effect_size_category'] == 'medium')

print("🎯 BIOLOGICAL CLAIM EVALUATION")
print("=" * 50)
print(f"Claim: 'No significant UPS protein alterations across tau-positive vs tau-negative neurons'")
print()

print("📊 STATISTICAL EVIDENCE:")
print(f"• Total UPS proteins analyzed: {len(results_df)}")
print(f"• Significant proteins (FDR < 0.05): {n_significant_corrected} ({percent_significant:.1f}%)")
print(f"• Large effect sizes (|d| > 0.8): {large_effects}")
print(f"• Medium effect sizes (|d| > 0.5): {medium_effects}")

print("\n🧬 BIOLOGICAL INTERPRETATION:")

if percent_significant < 10 and large_effects == 0:
    evaluation = "SUPPORTED"
    explanation = """
    ✅ CLAIM IS SUPPORTED

    Evidence:
    • Very few UPS proteins show significant changes (<10%)
    • No proteins show large effect sizes
    • Pattern consistent with intact UPS system

    Biological meaning:
    • UPS system appears functional in late-stage disease
    • Protein quality control machinery is preserved
    • UPS remains a viable therapeutic target
    • Protein aggregation likely due to other mechanisms
    """

elif percent_significant > 25 or large_effects > 3:
    evaluation = "REFUTED"
    explanation = """
    ❌ CLAIM IS REFUTED

    Evidence:
    • Many UPS proteins show significant changes (>25%)
    • Multiple proteins show large effect sizes
    • Pattern suggests UPS dysfunction

    Biological meaning:
    • UPS system is impaired in diseased neurons
    • Protein quality control is compromised
    • UPS dysfunction contributes to pathology
    • Therapeutic strategies should target UPS enhancement
    """

else:
    evaluation = "UNCLEAR"
    explanation = """
    ❓ CLAIM IS UNCLEAR

    Evidence:
    • Moderate number of proteins show changes (10-25%)
    • Mixed pattern of effect sizes
    • Some evidence for both intact and impaired UPS

    Biological meaning:
    • UPS system may be partially compromised
    • Different UPS components may be differentially affected
    • Further investigation needed
    • Consider temporal or subtype-specific analysis
    """

print(explanation)

# Store evaluation result
claim_evaluation = {
    'claim': 'No significant UPS protein alterations',
    'evaluation': evaluation,
    'percent_significant': percent_significant,
    'n_significant': n_significant_corrected,
    'n_total': len(results_df),
    'large_effects': large_effects,
    'evidence_strength': 'strong' if percent_significant < 5 or percent_significant > 30 else 'moderate'
}
```

### Summary Statistics Table
```python
# Create comprehensive summary table
summary_stats = {
    'Metric': [
        'Total UPS proteins analyzed',
        'Proteins with p < 0.05 (uncorrected)',
        'Proteins with FDR < 0.05',
        'Proteins with negligible effects (|d| < 0.2)',
        'Proteins with small effects (0.2 ≤ |d| < 0.5)',
        'Proteins with medium effects (0.5 ≤ |d| < 0.8)',
        'Proteins with large effects (|d| ≥ 0.8)',
        'Mean |Cohen\'s d|',
        'Mean fold change',
        'Proteins upregulated (tau+ > tau-)',
        'Proteins downregulated (tau+ < tau-)'
    ],
    'Value': [
        len(results_df),
        sum(results_df['p_value'] < 0.05),
        n_significant_corrected,
        sum(results_df['effect_size_category'] == 'negligible'),
        sum(results_df['effect_size_category'] == 'small'),
        sum(results_df['effect_size_category'] == 'medium'),
        sum(results_df['effect_size_category'] == 'large'),
        f"{results_df['cohens_d'].abs().mean():.3f}",
        f"{results_df['fold_change'].mean():.3f}",
        sum(results_df['cohens_d'] > 0),
        sum(results_df['cohens_d'] < 0)
    ]
}

summary_df = pd.DataFrame(summary_stats)
print("\n📊 COMPREHENSIVE SUMMARY TABLE:")
print(summary_df.to_string(index=False))
```

### Export Results
```python
# Save detailed results
results_df.to_csv('ups_protein_analysis_results.csv', index=False)

# Save summary
with open('ups_analysis_summary.txt', 'w') as f:
    f.write("UPS Protein Analysis Summary\n")
    f.write("=" * 30 + "\n\n")
    f.write(f"Claim: No significant UPS protein alterations\n")
    f.write(f"Evaluation: {evaluation}\n")
    f.write(f"Evidence strength: {claim_evaluation['evidence_strength']}\n\n")
    f.write(f"Statistical summary:\n")
    f.write(f"- Total proteins: {len(results_df)}\n")
    f.write(f"- Significant (FDR < 0.05): {n_significant_corrected} ({percent_significant:.1f}%)\n")
    f.write(f"- Large effects: {large_effects}\n")
    f.write(f"- Mean |Cohen's d|: {results_df['cohens_d'].abs().mean():.3f}\n")

print("💾 Results saved:")
print("• ups_protein_analysis_results.csv - Detailed results")
print("• ups_analysis_summary.txt - Summary evaluation")
```

---

## 🎯 Step 8: Quality Control and Validation

### Check Statistical Assumptions
```python
# Check normality assumption for a few proteins
print("🔍 QUALITY CONTROL CHECKS")
print("=" * 30)

# Test normality for first 5 proteins
from scipy.stats import shapiro

print("📊 Normality tests (Shapiro-Wilk) for first 5 proteins:")
for protein in protein_cols[:5]:
    pos_data = tau_positive[protein].values
    neg_data = tau_negative[protein].values

    _, p_pos = shapiro(pos_data)
    _, p_neg = shapiro(neg_data)

    print(f"{protein:>8}: tau+ p={p_pos:.3f}, tau- p={p_neg:.3f}", end="")
    if p_pos < 0.05 or p_neg < 0.05:
        print(" ⚠️ (non-normal)")
    else:
        print(" ✅ (normal)")

print("\n💡 Note: With large sample sizes, t-tests are robust to mild non-normality")
```

### Power Analysis
```python
# Calculate statistical power for detecting effects
print("\n⚡ STATISTICAL POWER ANALYSIS:")

min_group_size = min(len(tau_positive), len(tau_negative))

# Power to detect different effect sizes
effect_sizes = [0.2, 0.5, 0.8]
for effect_size in effect_sizes:
    power = ttest_power(effect_size, min_group_size, 0.05)
    print(f"Power to detect Cohen's d = {effect_size}: {power:.2f}")

print(f"\n📊 With {min_group_size} samples per group:")
if min_group_size > 50:
    print("✅ Excellent power to detect medium-large effects")
elif min_group_size > 25:
    print("✅ Good power to detect medium-large effects")
else:
    print("⚠️ Limited power - may miss small-medium effects")
```

### Sensitivity Analysis
```python
# Compare with non-parametric test for a few proteins
print("\n🔄 SENSITIVITY ANALYSIS:")
print("Comparing t-test vs Mann-Whitney U test for significant proteins:")

if n_significant_corrected > 0:
    sig_proteins = results_df[results_df['significant_fdr']]['protein'].values

    for protein in sig_proteins[:min(3, len(sig_proteins))]:  # Test up to 3
        pos_data = tau_positive[protein].values
        neg_data = tau_negative[protein].values

        # T-test (already done)
        t_p = results_df[results_df['protein'] == protein]['p_value'].iloc[0]

        # Mann-Whitney U test
        _, mw_p = stats.mannwhitneyu(pos_data, neg_data, alternative='two-sided')

        print(f"{protein}: t-test p={t_p:.4f}, Mann-Whitney p={mw_p:.4f}")

        if (t_p < 0.05) == (mw_p < 0.05):
            print("  ✅ Consistent results")
        else:
            print("  ⚠️ Different conclusions - investigate further")
else:
    print("No significant proteins to compare")
```

---

## ✅ Step 9: Checklist and Next Steps

### Analysis Completion Checklist
```python
print("📋 ANALYSIS COMPLETION CHECKLIST:")
print("=" * 35)

checklist = [
    ("Data loaded successfully", True),
    ("UPS proteins identified", len(available_ups) > 20),
    ("Groups compared statistically", len(results_df) > 0),
    ("Multiple testing correction applied", 'p_adjusted' in results_df.columns),
    ("Effect sizes calculated", 'cohens_d' in results_df.columns),
    ("Results interpreted biologically", 'evaluation' in locals()),
    ("Quality checks performed", True),
    ("Results exported", True)
]

for task, completed in checklist:
    status = "✅" if completed else "❌"
    print(f"{status} {task}")

all_complete = all(completed for _, completed in checklist)
print(f"\n🎯 Analysis {'complete' if all_complete else 'incomplete'}!")
```

### What You've Learned
```python
print("\n🎓 WHAT YOU'VE ACCOMPLISHED:")
print("• Performed differential expression analysis on UPS proteins")
print("• Applied proper multiple testing correction")
print("• Calculated and interpreted effect sizes")
print("• Made evidence-based biological conclusions")
print("• Conducted quality control checks")
print("• Exported results for further use")

print(f"\n🔬 BIOLOGICAL CONCLUSION:")
print(f"The claim that UPS proteins show no significant alterations is: {evaluation}")
print(f"This conclusion is based on {len(results_df)} proteins with {claim_evaluation['evidence_strength']} evidence.")
```

### Next Steps
```python
print("\n🚀 SUGGESTED NEXT STEPS:")

if evaluation == "SUPPORTED":
    print("1. Investigate why UPS remains intact when other systems fail")
    print("2. Explore UPS as a therapeutic target")
    print("3. Analyze temporal patterns in UPS function")
    print("4. Compare with autophagy system analysis")

elif evaluation == "REFUTED":
    print("1. Identify which UPS components are most affected")
    print("2. Investigate causes of UPS dysfunction")
    print("3. Explore UPS-enhancing therapeutic strategies")
    print("4. Analyze correlation with disease severity")

else:  # UNCLEAR
    print("1. Increase sample size for better power")
    print("2. Perform subgroup analyses")
    print("3. Investigate temporal changes")
    print("4. Consider alternative statistical approaches")

print("\n📚 Continue to next tutorial:")
print("• SQSTM1 upregulation analysis")
print("• Sliding window correlation analysis")
print("• Group 2 proteome-wide analysis")
```

---

## 🎉 Congratulations!

You've completed your first comprehensive proteomic analysis! You now know how to:

✅ **Load and explore** proteomic datasets
✅ **Identify relevant proteins** for biological questions
✅ **Perform statistical comparisons** between groups
✅ **Apply multiple testing correction** properly
✅ **Calculate and interpret effect sizes**
✅ **Make evidence-based biological conclusions**
✅ **Conduct quality control** checks
✅ **Export and document** your results

**This analysis framework can be applied to any proteomic comparison study!**

---

*Next: [Interpreting Results in Biological Context](interpreting_results.md)*

*Remember: The goal isn't just to run statistical tests, but to answer meaningful biological questions with appropriate rigor!* 🧬📊