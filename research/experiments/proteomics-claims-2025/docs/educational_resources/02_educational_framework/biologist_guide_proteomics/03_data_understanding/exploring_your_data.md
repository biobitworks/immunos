# 🔍 Hands-On Data Exploration Tutorial

## 🎯 What You'll Accomplish

By the end of this tutorial, you'll have:
- ✅ **Loaded your proteomics dataset** successfully
- ✅ **Explored sample characteristics** and metadata
- ✅ **Visualized protein expression patterns**
- ✅ **Identified key proteins** for subsequent analyses
- ✅ **Assessed data quality** and potential issues
- ✅ **Created your first publication-quality plots**

**Time needed**: 1-2 hours for complete beginners, 30-45 minutes with some experience

---

## 🛠️ Before We Start

### Prerequisites
- [ ] **Environment setup**: Jupyter notebook or Google Colab working
- [ ] **Data access**: Dataset downloaded and accessible
- [ ] **Packages installed**: scanpy, pandas, numpy, matplotlib, seaborn
- [ ] **Basic understanding**: Read the dataset overview guide

### Setup Your Workspace
```python
# Create new notebook called: "01_data_exploration.ipynb"
# If using local Jupyter: place in notebooks/ folder
# If using Colab: save to Google Drive proteomics folder
```

---

## 📚 Step 1: Load Libraries and Data

### Import Essential Libraries
```python
# Core data science libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Proteomics/single-cell analysis
import scanpy as sc

# Statistical tools
from scipy import stats

# Utilities
import warnings
warnings.filterwarnings('ignore')

# Configure plotting
plt.rcParams['figure.figsize'] = (8, 6)
plt.rcParams['font.size'] = 12
sns.set_style("whitegrid")
%matplotlib inline

print("✅ All libraries imported successfully!")
```

**Troubleshooting**: If you get import errors, refer to the software setup guide.

### Load the Proteomics Dataset

#### For Local Jupyter
```python
# Load the dataset (adjust path as needed)
print("Loading proteomics dataset...")

try:
    adata = sc.read_h5ad('../data/raw/pool_processed_v2.h5ad')
    print(f"✅ Dataset loaded successfully!")
    print(f"📊 Dataset shape: {adata.shape}")
    print(f"🧠 Samples (neurons): {adata.n_obs}")
    print(f"🧬 Proteins: {adata.n_vars}")

except FileNotFoundError:
    print("❌ Dataset not found!")
    print("Check file location and name:")
    print("Expected: ../data/raw/pool_processed_v2.h5ad")
    import os
    print(f"Current directory: {os.getcwd()}")
    print(f"Files in ../data/raw/: {os.listdir('../data/raw/') if os.path.exists('../data/raw/') else 'Directory not found'}")
```

#### For Google Colab
```python
# Mount Google Drive first
from google.colab import drive
drive.mount('/content/drive')

# Load the dataset
print("Loading proteomics dataset from Google Drive...")

try:
    adata = sc.read_h5ad('/content/drive/MyDrive/proteomics_analysis/data/raw/pool_processed_v2.h5ad')
    print(f"✅ Dataset loaded successfully!")
    print(f"📊 Dataset shape: {adata.shape}")
    print(f"🧠 Samples (neurons): {adata.n_obs}")
    print(f"🧬 Proteins: {adata.n_vars}")

except FileNotFoundError:
    print("❌ Dataset not found!")
    print("Check file location in Google Drive")
    import os
    print("Available folders:")
    print(os.listdir('/content/drive/MyDrive/proteomics_analysis/data/'))
```

### Initial Data Inspection
```python
# Basic dataset information
print("=" * 50)
print("DATASET OVERVIEW")
print("=" * 50)

print(f"Data matrix shape: {adata.X.shape}")
print(f"Data matrix type: {type(adata.X)}")
print(f"Memory usage: {adata.X.nbytes / 1024**2:.1f} MB")

# Check if data is sparse or dense
if hasattr(adata.X, 'toarray'):
    print("Data format: Sparse matrix")
    # Convert to dense for easier manipulation
    adata.X = adata.X.toarray()
    print("Converted to dense matrix for analysis")
else:
    print("Data format: Dense matrix")

print(f"Expression value range: {np.min(adata.X):.2f} to {np.max(adata.X):.2f}")
```

---

## 📋 Step 2: Explore Sample Metadata

### Examine Available Metadata
```python
print("=" * 50)
print("SAMPLE METADATA")
print("=" * 50)

print("Available metadata columns:")
for i, col in enumerate(adata.obs.columns, 1):
    print(f"{i:2d}. {col}")

print(f"\nFirst 5 samples metadata:")
print(adata.obs.head())
```

### Analyze Key Variables

#### Tau Status Distribution
```python
# Our primary comparison variable
print("TAU STATUS DISTRIBUTION")
print("-" * 30)

tau_counts = adata.obs['tau_status'].value_counts()
print(tau_counts)

tau_percentages = adata.obs['tau_status'].value_counts(normalize=True) * 100
print(f"\nPercentages:")
for status, pct in tau_percentages.items():
    print(f"{status:>8}: {pct:5.1f}%")

# Check for balanced groups
min_group = tau_counts.min()
max_group = tau_counts.max()
balance_ratio = min_group / max_group

print(f"\nGroup balance:")
print(f"Smallest group: {min_group} samples")
print(f"Largest group: {max_group} samples")
print(f"Balance ratio: {balance_ratio:.2f}")

if balance_ratio > 0.7:
    print("✅ Groups are well balanced")
elif balance_ratio > 0.5:
    print("⚠️ Groups somewhat imbalanced")
else:
    print("❌ Groups highly imbalanced - consider statistical adjustments")
```

#### Patient Demographics
```python
print("\nPATIENT DEMOGRAPHICS")
print("-" * 30)

# Age distribution
print("Age statistics:")
age_stats = adata.obs['age'].describe()
print(age_stats)

# Sex distribution
print("\nSex distribution:")
sex_counts = adata.obs['sex'].value_counts()
print(sex_counts)

# Age by tau status
print("\nAge by tau status:")
age_by_tau = adata.obs.groupby('tau_status')['age'].describe()
print(age_by_tau)

# Test for age differences between groups
tau_pos_age = adata.obs[adata.obs['tau_status'] == 'positive']['age']
tau_neg_age = adata.obs[adata.obs['tau_status'] == 'negative']['age']

t_stat, p_val = stats.ttest_ind(tau_pos_age, tau_neg_age)
print(f"\nAge difference test:")
print(f"Tau-positive mean age: {tau_pos_age.mean():.1f} years")
print(f"Tau-negative mean age: {tau_neg_age.mean():.1f} years")
print(f"T-test p-value: {p_val:.3f}")

if p_val < 0.05:
    print("⚠️ Significant age difference between groups - should control for age in analysis")
else:
    print("✅ No significant age difference between groups")
```

#### Technical Quality Variables
```python
print("\nTECHNICAL QUALITY VARIABLES")
print("-" * 30)

# Post-mortem interval
print("Post-mortem interval (PMI) statistics:")
pmi_stats = adata.obs['PMI'].describe()
print(pmi_stats)

# Batch distribution
print("\nBatch distribution:")
batch_counts = adata.obs['batch'].value_counts().sort_index()
print(batch_counts)

# Check for batch-tau status confounding
print("\nTau status by batch:")
batch_tau_table = pd.crosstab(adata.obs['batch'], adata.obs['tau_status'])
print(batch_tau_table)

# Chi-square test for independence
chi2, p_val, dof, expected = stats.chi2_contingency(batch_tau_table)
print(f"\nBatch-tau independence test:")
print(f"Chi-square p-value: {p_val:.3f}")

if p_val < 0.05:
    print("⚠️ Batch and tau status are not independent - must control for batch effects")
else:
    print("✅ Batch and tau status are independent")
```

### Visualize Sample Characteristics

#### Create Demographic Plots
```python
# Set up subplot layout
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Sample Demographics and Technical Variables', fontsize=16, y=1.02)

# Plot 1: Tau status distribution
ax1 = axes[0, 0]
tau_counts.plot(kind='bar', ax=ax1, color=['lightcoral', 'lightblue'])
ax1.set_title('Tau Status Distribution')
ax1.set_xlabel('Tau Status')
ax1.set_ylabel('Number of Samples')
ax1.tick_params(axis='x', rotation=0)

# Add sample counts on bars
for i, v in enumerate(tau_counts.values):
    ax1.text(i, v + 1, str(v), ha='center', va='bottom')

# Plot 2: Age distribution by tau status
ax2 = axes[0, 1]
tau_pos_age = adata.obs[adata.obs['tau_status'] == 'positive']['age']
tau_neg_age = adata.obs[adata.obs['tau_status'] == 'negative']['age']

ax2.hist([tau_neg_age, tau_pos_age], bins=10, alpha=0.7,
         label=['Tau-negative', 'Tau-positive'], color=['lightblue', 'lightcoral'])
ax2.set_title('Age Distribution by Tau Status')
ax2.set_xlabel('Age (years)')
ax2.set_ylabel('Number of Samples')
ax2.legend()

# Plot 3: PMI distribution
ax3 = axes[1, 0]
adata.obs['PMI'].hist(bins=15, ax=ax3, alpha=0.7, color='lightgreen')
ax3.set_title('Post-Mortem Interval Distribution')
ax3.set_xlabel('PMI (hours)')
ax3.set_ylabel('Number of Samples')

# Plot 4: Batch composition
ax4 = axes[1, 1]
batch_tau_table.plot(kind='bar', stacked=True, ax=ax4,
                     color=['lightblue', 'lightcoral'])
ax4.set_title('Tau Status by Batch')
ax4.set_xlabel('Batch')
ax4.set_ylabel('Number of Samples')
ax4.legend(title='Tau Status')
ax4.tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.show()

print("📊 Demographics plots created successfully!")
```

---

## 🧬 Step 3: Explore Protein Expression

### Basic Expression Statistics
```python
print("=" * 50)
print("PROTEIN EXPRESSION OVERVIEW")
print("=" * 50)

# Overall expression statistics
print("Expression matrix statistics:")
print(f"Total values: {adata.X.size:,}")
print(f"Non-zero values: {np.count_nonzero(adata.X):,}")
print(f"Zero values: {np.sum(adata.X == 0):,}")
print(f"Missing/NaN values: {np.sum(np.isnan(adata.X)):,}")

print(f"\nExpression value distribution:")
print(f"Minimum: {np.min(adata.X):.2f}")
print(f"25th percentile: {np.percentile(adata.X, 25):.2f}")
print(f"Median: {np.median(adata.X):.2f}")
print(f"75th percentile: {np.percentile(adata.X, 75):.2f}")
print(f"Maximum: {np.max(adata.X):.2f}")
print(f"Mean: {np.mean(adata.X):.2f}")
print(f"Standard deviation: {np.std(adata.X):.2f}")
```

### Protein Detection Analysis
```python
# Calculate detection rates for each protein
print("\nPROTEIN DETECTION ANALYSIS")
print("-" * 30)

# Detection rate = fraction of samples where protein is detected (>0)
detection_rates = np.mean(adata.X > 0, axis=0)
adata.var['detection_rate'] = detection_rates

# Calculate mean expression for each protein
mean_expression = np.mean(adata.X, axis=0)
adata.var['mean_expression'] = mean_expression

# Detection rate statistics
print("Detection rate statistics:")
print(f"Mean detection rate: {np.mean(detection_rates):.2f}")
print(f"Median detection rate: {np.median(detection_rates):.2f}")
print(f"Proteins detected in >90% samples: {np.sum(detection_rates > 0.9):,}")
print(f"Proteins detected in >75% samples: {np.sum(detection_rates > 0.75):,}")
print(f"Proteins detected in >50% samples: {np.sum(detection_rates > 0.5):,}")
print(f"Proteins detected in <25% samples: {np.sum(detection_rates < 0.25):,}")

# Highly detected proteins (good for analysis)
well_detected = adata.var[adata.var['detection_rate'] > 0.75]
print(f"\nWell-detected proteins (>75% samples): {len(well_detected):,}")

# Poorly detected proteins (exclude from analysis)
poorly_detected = adata.var[adata.var['detection_rate'] < 0.25]
print(f"Poorly detected proteins (<25% samples): {len(poorly_detected):,}")
```

### Expression Distribution Visualization
```python
# Create expression distribution plots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Protein Expression Patterns', fontsize=16, y=1.02)

# Plot 1: Overall expression distribution
ax1 = axes[0, 0]
ax1.hist(adata.X.flatten(), bins=50, alpha=0.7, color='skyblue', edgecolor='black')
ax1.set_title('Overall Expression Distribution')
ax1.set_xlabel('Log2 Expression Level')
ax1.set_ylabel('Frequency')
ax1.axvline(np.mean(adata.X), color='red', linestyle='--', label=f'Mean: {np.mean(adata.X):.1f}')
ax1.legend()

# Plot 2: Detection rate distribution
ax2 = axes[0, 1]
ax2.hist(detection_rates, bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
ax2.set_title('Protein Detection Rate Distribution')
ax2.set_xlabel('Detection Rate (fraction of samples)')
ax2.set_ylabel('Number of Proteins')
ax2.axvline(0.75, color='red', linestyle='--', label='75% threshold')
ax2.legend()

# Plot 3: Detection rate vs mean expression
ax3 = axes[1, 0]
ax3.scatter(mean_expression, detection_rates, alpha=0.6, s=10)
ax3.set_title('Detection Rate vs Mean Expression')
ax3.set_xlabel('Mean Log2 Expression')
ax3.set_ylabel('Detection Rate')
ax3.axhline(0.75, color='red', linestyle='--', alpha=0.7)
ax3.axvline(np.median(mean_expression), color='orange', linestyle='--', alpha=0.7)

# Plot 4: Sample-wise expression distribution
ax4 = axes[1, 1]
sample_means = np.mean(adata.X, axis=1)
ax4.hist(sample_means, bins=20, alpha=0.7, color='lightsalmon', edgecolor='black')
ax4.set_title('Sample-wise Mean Expression')
ax4.set_xlabel('Mean Log2 Expression per Sample')
ax4.set_ylabel('Number of Samples')

plt.tight_layout()
plt.show()

print("📊 Expression distribution plots created!")
```

---

## 🔍 Step 4: Identify Key Proteins

### Find Our Target Proteins
```python
print("=" * 50)
print("KEY PROTEIN IDENTIFICATION")
print("=" * 50)

# Define proteins of interest for our analyses
key_proteins = {
    'SQSTM1': 'Autophagy receptor (our main focus)',
    'VDAC1': 'Mitochondrial channel (SQSTM1 target)',
    'PSMA1': 'Proteasome subunit (UPS system)',
    'PSMB1': 'Proteasome subunit (UPS system)',
    'PSMC1': 'Proteasome subunit (UPS system)',
    'PSMD1': 'Proteasome subunit (UPS system)',
    'UBA1': 'E1 ubiquitin enzyme (UPS system)',
    'UBE2B': 'E2 ubiquitin enzyme (UPS system)',
    'ACTB': 'Beta-actin (housekeeping control)',
    'GAPDH': 'Glycolysis enzyme (housekeeping control)'
}

# Check which proteins are present in our dataset
present_proteins = []
missing_proteins = []

print("Protein availability check:")
for protein, description in key_proteins.items():
    if protein in adata.var_names:
        present_proteins.append(protein)
        # Get detection rate and mean expression
        idx = adata.var_names.get_loc(protein)
        det_rate = detection_rates[idx]
        mean_expr = mean_expression[idx]
        print(f"✅ {protein:<8}: {description} (detection: {det_rate:.1%}, mean: {mean_expr:.1f})")
    else:
        missing_proteins.append(protein)
        print(f"❌ {protein:<8}: {description} - NOT FOUND")

print(f"\nSummary:")
print(f"Present proteins: {len(present_proteins)}/{len(key_proteins)}")
print(f"Missing proteins: {len(missing_proteins)}")

if missing_proteins:
    print(f"\nMissing proteins: {', '.join(missing_proteins)}")
    print("This is normal - not all proteins are detected in every experiment")
```

### Explore High-Impact Proteins
```python
# Find proteins with highest variance (most interesting biologically)
print("\nTOP VARIABLE PROTEINS")
print("-" * 30)

protein_vars = np.var(adata.X, axis=0)
adata.var['variance'] = protein_vars

# Top 10 most variable proteins
top_variable_indices = np.argsort(protein_vars)[-10:]
top_variable_proteins = adata.var_names[top_variable_indices]

print("Top 10 most variable proteins:")
for i, protein in enumerate(reversed(top_variable_proteins), 1):
    idx = adata.var_names.get_loc(protein)
    variance = protein_vars[idx]
    mean_expr = mean_expression[idx]
    det_rate = detection_rates[idx]
    print(f"{i:2d}. {protein:<12} (var: {variance:.2f}, mean: {mean_expr:.1f}, det: {det_rate:.1%})")
```

### Expression Patterns of Key Proteins
```python
# Visualize expression of key proteins
if len(present_proteins) >= 4:
    # Select 4 key proteins for visualization
    proteins_to_plot = present_proteins[:4]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Key Protein Expression Patterns', fontsize=16, y=1.02)

    for i, protein in enumerate(proteins_to_plot):
        ax = axes[i//2, i%2]

        # Get protein expression
        protein_idx = adata.var_names.get_loc(protein)
        protein_expr = adata.X[:, protein_idx]

        # Split by tau status
        tau_pos_expr = protein_expr[adata.obs['tau_status'] == 'positive']
        tau_neg_expr = protein_expr[adata.obs['tau_status'] == 'negative']

        # Create violin plot
        data_for_plot = [tau_neg_expr, tau_pos_expr]
        parts = ax.violinplot(data_for_plot, positions=[1, 2], widths=0.7)

        # Customize violin plot
        for pc in parts['bodies']:
            pc.set_facecolor('lightblue')
            pc.set_alpha(0.7)

        # Add box plot overlay
        ax.boxplot(data_for_plot, positions=[1, 2], widths=0.3,
                  patch_artist=True,
                  boxprops=dict(facecolor='white', alpha=0.8))

        # Calculate fold change
        mean_pos = np.mean(tau_pos_expr)
        mean_neg = np.mean(tau_neg_expr)
        log2_fc = mean_pos - mean_neg
        fold_change = 2**log2_fc

        ax.set_title(f'{protein}\n({fold_change:.1f}-fold change)')
        ax.set_xticks([1, 2])
        ax.set_xticklabels(['Tau-negative', 'Tau-positive'])
        ax.set_ylabel('Log2 Expression')

        # Add statistics
        t_stat, p_val = stats.ttest_ind(tau_neg_expr, tau_pos_expr)
        if p_val < 0.001:
            p_text = "p < 0.001"
        elif p_val < 0.01:
            p_text = f"p = {p_val:.3f}"
        else:
            p_text = f"p = {p_val:.2f}"

        ax.text(1.5, max(max(tau_pos_expr), max(tau_neg_expr)) * 1.05,
                p_text, ha='center', fontsize=10)

    plt.tight_layout()
    plt.show()

    print("📊 Key protein expression plots created!")
else:
    print("⚠️ Not enough key proteins found for plotting")
```

---

## 📊 Step 5: Sample Quality Assessment

### Assess Sample Quality
```python
print("=" * 50)
print("SAMPLE QUALITY ASSESSMENT")
print("=" * 50)

# Calculate quality metrics for each sample
sample_metrics = pd.DataFrame(index=adata.obs_names)

# Number of detected proteins per sample
sample_metrics['n_proteins_detected'] = np.sum(adata.X > 0, axis=1)

# Mean expression per sample
sample_metrics['mean_expression'] = np.mean(adata.X, axis=1)

# Expression variance per sample (technical noise indicator)
sample_metrics['expression_variance'] = np.var(adata.X, axis=1)

# Add metadata
sample_metrics['tau_status'] = adata.obs['tau_status']
sample_metrics['age'] = adata.obs['age']
sample_metrics['PMI'] = adata.obs['PMI']
sample_metrics['batch'] = adata.obs['batch']

print("Sample quality metrics:")
print(sample_metrics.describe())

# Quality thresholds
min_proteins = 2000  # Minimum proteins detected
max_variance = np.percentile(sample_metrics['expression_variance'], 95)  # 95th percentile

print(f"\nQuality assessment:")
print(f"Samples with <{min_proteins} proteins: {np.sum(sample_metrics['n_proteins_detected'] < min_proteins)}")
print(f"Samples with high variance (>{max_variance:.1f}): {np.sum(sample_metrics['expression_variance'] > max_variance)}")

# Flag potentially problematic samples
problematic_samples = sample_metrics[
    (sample_metrics['n_proteins_detected'] < min_proteins) |
    (sample_metrics['expression_variance'] > max_variance)
]

if len(problematic_samples) > 0:
    print(f"\n⚠️ {len(problematic_samples)} potentially problematic samples identified:")
    print(problematic_samples)
else:
    print("\n✅ All samples pass quality thresholds")
```

### Quality Control Visualizations
```python
# Create quality control plots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Sample Quality Control', fontsize=16, y=1.02)

# Plot 1: Proteins detected vs mean expression
ax1 = axes[0, 0]
colors = ['lightcoral' if status == 'positive' else 'lightblue'
          for status in sample_metrics['tau_status']]
scatter = ax1.scatter(sample_metrics['mean_expression'],
                     sample_metrics['n_proteins_detected'],
                     c=colors, alpha=0.7, s=50)
ax1.set_xlabel('Mean Expression')
ax1.set_ylabel('Proteins Detected')
ax1.set_title('Sample Quality Overview')

# Add threshold lines
ax1.axhline(min_proteins, color='red', linestyle='--', alpha=0.7,
           label=f'Min proteins: {min_proteins}')
ax1.legend()

# Plot 2: Expression variance by tau status
ax2 = axes[0, 1]
tau_neg_var = sample_metrics[sample_metrics['tau_status'] == 'negative']['expression_variance']
tau_pos_var = sample_metrics[sample_metrics['tau_status'] == 'positive']['expression_variance']

ax2.boxplot([tau_neg_var, tau_pos_var], labels=['Tau-negative', 'Tau-positive'])
ax2.set_ylabel('Expression Variance')
ax2.set_title('Technical Noise by Group')

# Plot 3: PMI vs proteins detected
ax3 = axes[1, 0]
ax3.scatter(sample_metrics['PMI'], sample_metrics['n_proteins_detected'],
           c=colors, alpha=0.7, s=50)
ax3.set_xlabel('Post-Mortem Interval (hours)')
ax3.set_ylabel('Proteins Detected')
ax3.set_title('PMI Effect on Detection')

# Add correlation
corr, p_val = stats.pearsonr(sample_metrics['PMI'], sample_metrics['n_proteins_detected'])
ax3.text(0.05, 0.95, f'r = {corr:.3f}\np = {p_val:.3f}',
         transform=ax3.transAxes, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# Plot 4: Batch effects on mean expression
ax4 = axes[1, 1]
batch_groups = [sample_metrics[sample_metrics['batch'] == batch]['mean_expression'].values
                for batch in sorted(sample_metrics['batch'].unique())]
ax4.boxplot(batch_groups, labels=sorted(sample_metrics['batch'].unique()))
ax4.set_xlabel('Batch')
ax4.set_ylabel('Mean Expression')
ax4.set_title('Batch Effects')
ax4.tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.show()

print("📊 Quality control plots created!")
```

---

## 🎯 Step 6: Initial Biological Insights

### Preliminary Analysis: SQSTM1 Expression
```python
print("=" * 50)
print("PRELIMINARY ANALYSIS: SQSTM1")
print("=" * 50)

if 'SQSTM1' in present_proteins:
    # Get SQSTM1 expression data
    sqstm1_idx = adata.var_names.get_loc('SQSTM1')
    sqstm1_expr = adata.X[:, sqstm1_idx]

    # Split by tau status
    tau_pos_sqstm1 = sqstm1_expr[adata.obs['tau_status'] == 'positive']
    tau_neg_sqstm1 = sqstm1_expr[adata.obs['tau_status'] == 'negative']

    # Calculate statistics
    mean_pos = np.mean(tau_pos_sqstm1)
    mean_neg = np.mean(tau_neg_sqstm1)
    log2_fc = mean_pos - mean_neg
    fold_change = 2**log2_fc

    # Statistical test
    t_stat, p_val = stats.ttest_ind(tau_pos_sqstm1, tau_neg_sqstm1)

    print("SQSTM1 expression analysis:")
    print(f"Tau-positive mean: {mean_pos:.2f}")
    print(f"Tau-negative mean: {mean_neg:.2f}")
    print(f"Log2 fold change: {log2_fc:.2f}")
    print(f"Fold change: {fold_change:.1f}")
    print(f"T-statistic: {t_stat:.2f}")
    print(f"P-value: {p_val:.3e}")

    # Compare to expected 10.7-fold change
    expected_fc = 10.7
    if abs(fold_change - expected_fc) < 2:
        print(f"✅ Close to expected {expected_fc}-fold change!")
    else:
        print(f"⚠️ Different from expected {expected_fc}-fold change")

    # Create SQSTM1 visualization
    plt.figure(figsize=(10, 6))

    # Subplot 1: Box plot
    plt.subplot(1, 2, 1)
    data_for_box = [tau_neg_sqstm1, tau_pos_sqstm1]
    box_plot = plt.boxplot(data_for_box, labels=['Tau-negative', 'Tau-positive'],
                          patch_artist=True)
    box_plot['boxes'][0].set_facecolor('lightblue')
    box_plot['boxes'][1].set_facecolor('lightcoral')
    plt.ylabel('SQSTM1 Log2 Expression')
    plt.title('SQSTM1 Expression by Tau Status')

    # Add statistics text
    plt.text(1.5, max(sqstm1_expr) * 0.95,
             f'{fold_change:.1f}-fold change\np = {p_val:.2e}',
             ha='center', fontsize=12,
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Subplot 2: Individual points
    plt.subplot(1, 2, 2)
    x_neg = np.random.normal(1, 0.05, len(tau_neg_sqstm1))
    x_pos = np.random.normal(2, 0.05, len(tau_pos_sqstm1))

    plt.scatter(x_neg, tau_neg_sqstm1, alpha=0.6, color='lightblue', s=50, label='Tau-negative')
    plt.scatter(x_pos, tau_pos_sqstm1, alpha=0.6, color='lightcoral', s=50, label='Tau-positive')

    plt.xlim(0.5, 2.5)
    plt.xticks([1, 2], ['Tau-negative', 'Tau-positive'])
    plt.ylabel('SQSTM1 Log2 Expression')
    plt.title('Individual Sample Values')
    plt.legend()

    plt.tight_layout()
    plt.show()

    print("📊 SQSTM1 analysis plots created!")

else:
    print("❌ SQSTM1 not found in dataset - cannot perform preliminary analysis")
```

### Correlation Analysis Preview
```python
print("\nCORRELATION ANALYSIS PREVIEW")
print("-" * 30)

if 'SQSTM1' in present_proteins and 'VDAC1' in present_proteins:
    # Get expression data
    sqstm1_idx = adata.var_names.get_loc('SQSTM1')
    vdac1_idx = adata.var_names.get_loc('VDAC1')

    sqstm1_expr = adata.X[:, sqstm1_idx]
    vdac1_expr = adata.X[:, vdac1_idx]

    # Overall correlation
    overall_corr, overall_p = stats.pearsonr(sqstm1_expr, vdac1_expr)

    # Correlation by group
    tau_pos_mask = adata.obs['tau_status'] == 'positive'
    tau_neg_mask = adata.obs['tau_status'] == 'negative'

    pos_corr, pos_p = stats.pearsonr(sqstm1_expr[tau_pos_mask], vdac1_expr[tau_pos_mask])
    neg_corr, neg_p = stats.pearsonr(sqstm1_expr[tau_neg_mask], vdac1_expr[tau_neg_mask])

    print("SQSTM1-VDAC1 correlations:")
    print(f"Overall: r = {overall_corr:.3f}, p = {overall_p:.3f}")
    print(f"Tau-positive: r = {pos_corr:.3f}, p = {pos_p:.3f}")
    print(f"Tau-negative: r = {neg_corr:.3f}, p = {neg_p:.3f}")

    # Visualization
    plt.figure(figsize=(12, 4))

    # Overall correlation
    plt.subplot(1, 3, 1)
    plt.scatter(sqstm1_expr, vdac1_expr, alpha=0.6, s=50)
    plt.xlabel('SQSTM1 Expression')
    plt.ylabel('VDAC1 Expression')
    plt.title(f'Overall Correlation\nr = {overall_corr:.3f}')

    # Tau-negative correlation
    plt.subplot(1, 3, 2)
    plt.scatter(sqstm1_expr[tau_neg_mask], vdac1_expr[tau_neg_mask],
               alpha=0.6, s=50, color='lightblue')
    plt.xlabel('SQSTM1 Expression')
    plt.ylabel('VDAC1 Expression')
    plt.title(f'Tau-Negative\nr = {neg_corr:.3f}')

    # Tau-positive correlation
    plt.subplot(1, 3, 3)
    plt.scatter(sqstm1_expr[tau_pos_mask], vdac1_expr[tau_pos_mask],
               alpha=0.6, s=50, color='lightcoral')
    plt.xlabel('SQSTM1 Expression')
    plt.ylabel('VDAC1 Expression')
    plt.title(f'Tau-Positive\nr = {pos_corr:.3f}')

    plt.tight_layout()
    plt.show()

    print("📊 Correlation analysis plots created!")

else:
    print("❌ Both SQSTM1 and VDAC1 needed for correlation analysis")
```

---

## 📋 Step 7: Summary and Next Steps

### Create Analysis Summary
```python
print("=" * 60)
print("DATA EXPLORATION SUMMARY")
print("=" * 60)

# Dataset overview
print("DATASET OVERVIEW:")
print(f"• Samples: {adata.n_obs}")
print(f"• Proteins: {adata.n_vars}")
print(f"• Well-detected proteins (>75%): {np.sum(detection_rates > 0.75):,}")

# Sample characteristics
tau_counts = adata.obs['tau_status'].value_counts()
print(f"\nSAMPLE GROUPS:")
print(f"• Tau-positive: {tau_counts.get('positive', 0)}")
print(f"• Tau-negative: {tau_counts.get('negative', 0)}")

print(f"\nAGE DEMOGRAPHICS:")
print(f"• Age range: {adata.obs['age'].min():.0f}-{adata.obs['age'].max():.0f} years")
print(f"• Mean age: {adata.obs['age'].mean():.1f} years")

# Key proteins status
print(f"\nKEY PROTEINS:")
print(f"• Found: {len(present_proteins)}/{len(key_proteins)}")
if missing_proteins:
    print(f"• Missing: {', '.join(missing_proteins)}")

# Quality assessment
problematic_count = len(problematic_samples) if 'problematic_samples' in locals() else 0
print(f"\nQUALITY ASSESSMENT:")
print(f"• High-quality samples: {adata.n_obs - problematic_count}/{adata.n_obs}")
if problematic_count > 0:
    print(f"• Flagged samples: {problematic_count}")

# SQSTM1 preliminary results
if 'SQSTM1' in present_proteins:
    print(f"\nSQSTM1 PRELIMINARY RESULTS:")
    print(f"• Fold change: {fold_change:.1f}")
    print(f"• P-value: {p_val:.2e}")
    if p_val < 0.001:
        print(f"• Status: Highly significant upregulation ✅")
    elif p_val < 0.05:
        print(f"• Status: Significant upregulation ✅")
    else:
        print(f"• Status: Not significant ❌")

print("\n" + "=" * 60)
```

### Save Your Work
```python
# Save exploration results for future use
exploration_results = {
    'dataset_shape': adata.shape,
    'present_proteins': present_proteins,
    'missing_proteins': missing_proteins,
    'sample_metadata': sample_metrics,
    'protein_detection_rates': detection_rates,
    'protein_mean_expression': mean_expression
}

# For local Jupyter
try:
    import pickle
    with open('../results/exploration_results.pkl', 'wb') as f:
        pickle.dump(exploration_results, f)
    print("💾 Exploration results saved to ../results/exploration_results.pkl")
except:
    print("⚠️ Could not save results - check results/ directory exists")

# For Google Colab
try:
    with open('/content/drive/MyDrive/proteomics_analysis/results/exploration_results.pkl', 'wb') as f:
        pickle.dump(exploration_results, f)
    print("💾 Exploration results saved to Google Drive")
except:
    print("⚠️ Could not save to Google Drive")

# Also save key plots
try:
    plt.figure(figsize=(10, 8))
    # Recreate key summary plot
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Add your most important plots here
    plt.suptitle('Data Exploration Summary', fontsize=16)

    # Save
    plt.tight_layout()
    plt.savefig('../plots/data_exploration_summary.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("📊 Summary plots saved")
except:
    print("⚠️ Could not save plots")
```

### Readiness Check
```python
print("=" * 50)
print("READINESS FOR ANALYSIS")
print("=" * 50)

readiness_checks = {
    'Data loaded successfully': adata is not None and adata.shape[0] > 0,
    'Key proteins present': len(present_proteins) >= 3,
    'Sample groups balanced': tau_counts.min() / tau_counts.max() > 0.5,
    'Quality thresholds met': problematic_count / adata.n_obs < 0.2,
    'SQSTM1 available for analysis': 'SQSTM1' in present_proteins,
    'Metadata variables present': all(col in adata.obs.columns for col in ['tau_status', 'age', 'PMI'])
}

print("Readiness checklist:")
all_ready = True
for check, status in readiness_checks.items():
    symbol = "✅" if status else "❌"
    print(f"{symbol} {check}")
    if not status:
        all_ready = False

print(f"\nOVERALL STATUS: {'✅ READY FOR ANALYSIS!' if all_ready else '⚠️ ADDRESS ISSUES BEFORE PROCEEDING'}")
```

---

## 🎯 Learning Check

After completing this exploration, you should be able to:

### Data Understanding
- [ ] Load and inspect proteomics datasets in H5AD format
- [ ] Understand what log2 expression values represent
- [ ] Interpret detection rates and data quality metrics
- [ ] Identify potential confounding variables

### Visualization Skills
- [ ] Create informative plots of sample characteristics
- [ ] Visualize protein expression patterns by groups
- [ ] Generate quality control visualizations
- [ ] Interpret correlation plots

### Analysis Preparation
- [ ] Assess whether data is suitable for planned analyses
- [ ] Identify which proteins are available for study
- [ ] Understand sample size and power considerations
- [ ] Recognize potential technical issues

### Biological Interpretation
- [ ] Connect statistical patterns to biological mechanisms
- [ ] Distinguish biological from technical variation
- [ ] Interpret preliminary results in disease context
- [ ] Plan follow-up analyses based on initial findings

---

## 🚀 Next Steps

### Immediate Actions
1. **Review any flagged issues** from quality assessment
2. **Plan your analysis strategy** based on available proteins
3. **Choose your first analysis** (UPS proteins, SQSTM1, or temporal)
4. **Organize your files** for reproducible analysis

### Continue Learning
- **Quality Control Deep Dive**: Learn advanced QC methods
- **Group 1 Analyses**: Start with focused protein studies
- **Statistical Methods**: Understand the theory behind the methods
- **Biological Interpretation**: Connect results to disease mechanisms

### Analysis Readiness
You're now ready to begin the focused analyses! Your next choice:
- **UPS Protein Analysis**: Test protein quality control system integrity
- **SQSTM1 Analysis**: Validate the 10.7-fold upregulation claim
- **Temporal Analysis**: Understand disease progression dynamics

---

**Congratulations!** You've successfully explored your proteomics dataset and gained crucial insights into the data structure, quality, and biological patterns. You're now equipped to perform sophisticated statistical analyses with confidence!

*Next: [Quality Control Deep Dive](quality_control_basics.md) or jump to [UPS Protein Analysis](../04_group1_analyses/statement1_ups_proteins/step_by_step_analysis.md)*

*Remember: Good data exploration is the foundation of all successful analyses - the time you've invested here will pay dividends in every subsequent analysis!* 🔍🧬