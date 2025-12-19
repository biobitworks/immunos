# 🔍 Quality Control Basics for Proteomics Data

## 🎯 What You'll Master

By the end of this guide, you'll know how to:
- ✅ **Identify common quality issues** in proteomics datasets
- ✅ **Perform systematic quality control** checks
- ✅ **Decide which samples/proteins to exclude** from analysis
- ✅ **Correct for technical artifacts** when possible
- ✅ **Document quality control decisions** for reproducibility
- ✅ **Assess whether your data is suitable** for planned analyses

**Why This Matters**: Poor quality data leads to unreliable results. Quality control is your first line of defense against false discoveries.

---

## 🧪 Understanding Proteomics Data Quality

### What Can Go Wrong in Proteomics?

#### Sample-Level Issues
1. **Poor sample preservation**: Degraded proteins due to improper storage
2. **Contamination**: Non-neuronal cells mixed with neurons
3. **Technical failures**: Problems during sample processing
4. **Batch effects**: Systematic differences between processing batches
5. **Outlier samples**: Extreme values due to various factors

#### Protein-Level Issues
1. **Low detection sensitivity**: Some proteins below detection limits
2. **Non-specific detection**: False positive protein identifications
3. **Quantification errors**: Inaccurate abundance measurements
4. **Missing values**: Proteins not detected in some samples
5. **Systematic biases**: Consistent errors affecting specific protein types

#### Study Design Issues
1. **Confounding variables**: Age, sex, PMI correlated with disease status
2. **Unbalanced groups**: Unequal sample sizes between conditions
3. **Insufficient power**: Too few samples for reliable statistical inference
4. **Technical replication**: Lack of proper technical controls

---

## 📊 Comprehensive Quality Control Pipeline

### Step 1: Load Data and Setup

```python
# Essential libraries for quality control
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# Load your dataset
print("Loading dataset for quality control...")
adata = sc.read_h5ad('your_dataset_path.h5ad')  # Adjust path as needed

print(f"Dataset dimensions: {adata.shape}")
print(f"Samples: {adata.n_obs}")
print(f"Proteins: {adata.n_vars}")

# Configure plotting
plt.style.use('default')
sns.set_palette("Set2")
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 11
```

### Step 2: Sample-Level Quality Control

#### Basic Sample Statistics
```python
print("=" * 50)
print("SAMPLE-LEVEL QUALITY CONTROL")
print("=" * 50)

# Calculate sample-level metrics
sample_qc = pd.DataFrame(index=adata.obs_names)

# Number of proteins detected per sample
sample_qc['n_proteins'] = np.sum(adata.X > 0, axis=1)

# Total protein abundance per sample
sample_qc['total_abundance'] = np.sum(adata.X, axis=1)

# Mean protein abundance per sample
sample_qc['mean_abundance'] = np.mean(adata.X, axis=1)

# Median protein abundance per sample
sample_qc['median_abundance'] = np.median(adata.X, axis=1)

# Standard deviation (technical noise indicator)
sample_qc['abundance_std'] = np.std(adata.X, axis=1)

# Coefficient of variation
sample_qc['abundance_cv'] = sample_qc['abundance_std'] / sample_qc['mean_abundance']

# Add metadata
for col in adata.obs.columns:
    sample_qc[col] = adata.obs[col]

print("Sample quality metrics summary:")
print(sample_qc[['n_proteins', 'mean_abundance', 'abundance_std', 'abundance_cv']].describe())
```

#### Identify Sample Outliers
```python
# Define quality thresholds
print("\nIDENTIFYING SAMPLE OUTLIERS")
print("-" * 30)

# Thresholds for outlier detection
min_proteins_threshold = np.percentile(sample_qc['n_proteins'], 5)  # 5th percentile
max_proteins_threshold = np.percentile(sample_qc['n_proteins'], 95)  # 95th percentile

mean_abundance_lower = np.percentile(sample_qc['mean_abundance'], 5)
mean_abundance_upper = np.percentile(sample_qc['mean_abundance'], 95)

cv_threshold = np.percentile(sample_qc['abundance_cv'], 95)  # High variability

print(f"Quality thresholds:")
print(f"Proteins detected: {min_proteins_threshold:.0f} - {max_proteins_threshold:.0f}")
print(f"Mean abundance: {mean_abundance_lower:.2f} - {mean_abundance_upper:.2f}")
print(f"CV threshold: <{cv_threshold:.3f}")

# Flag outlier samples
outlier_flags = (
    (sample_qc['n_proteins'] < min_proteins_threshold) |
    (sample_qc['n_proteins'] > max_proteins_threshold) |
    (sample_qc['mean_abundance'] < mean_abundance_lower) |
    (sample_qc['mean_abundance'] > mean_abundance_upper) |
    (sample_qc['abundance_cv'] > cv_threshold)
)

sample_qc['is_outlier'] = outlier_flags
outlier_samples = sample_qc[outlier_flags]

print(f"\nOutlier samples identified: {len(outlier_samples)}/{len(sample_qc)} ({len(outlier_samples)/len(sample_qc)*100:.1f}%)")

if len(outlier_samples) > 0:
    print("\nOutlier samples details:")
    print(outlier_samples[['n_proteins', 'mean_abundance', 'abundance_cv', 'tau_status']].head(10))
```

#### Visualize Sample Quality
```python
# Create sample quality plots
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('Sample Quality Control Assessment', fontsize=16, y=1.02)

# Plot 1: Proteins detected distribution
ax1 = axes[0, 0]
ax1.hist(sample_qc['n_proteins'], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
ax1.axvline(min_proteins_threshold, color='red', linestyle='--', label=f'5th percentile: {min_proteins_threshold:.0f}')
ax1.axvline(max_proteins_threshold, color='red', linestyle='--', label=f'95th percentile: {max_proteins_threshold:.0f}')
ax1.set_xlabel('Number of Proteins Detected')
ax1.set_ylabel('Number of Samples')
ax1.set_title('Proteins Detected per Sample')
ax1.legend()

# Plot 2: Mean abundance distribution
ax2 = axes[0, 1]
ax2.hist(sample_qc['mean_abundance'], bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
ax2.axvline(mean_abundance_lower, color='red', linestyle='--')
ax2.axvline(mean_abundance_upper, color='red', linestyle='--')
ax2.set_xlabel('Mean Protein Abundance')
ax2.set_ylabel('Number of Samples')
ax2.set_title('Mean Abundance per Sample')

# Plot 3: Coefficient of variation
ax3 = axes[0, 2]
ax3.hist(sample_qc['abundance_cv'], bins=30, alpha=0.7, color='lightsalmon', edgecolor='black')
ax3.axvline(cv_threshold, color='red', linestyle='--', label=f'95th percentile: {cv_threshold:.3f}')
ax3.set_xlabel('Coefficient of Variation')
ax3.set_ylabel('Number of Samples')
ax3.set_title('Sample Variability')
ax3.legend()

# Plot 4: Proteins vs mean abundance (outlier detection)
ax4 = axes[1, 0]
normal_samples = sample_qc[~sample_qc['is_outlier']]
outlier_samp = sample_qc[sample_qc['is_outlier']]

ax4.scatter(normal_samples['mean_abundance'], normal_samples['n_proteins'],
           alpha=0.6, color='blue', s=30, label='Normal samples')
if len(outlier_samp) > 0:
    ax4.scatter(outlier_samp['mean_abundance'], outlier_samp['n_proteins'],
               alpha=0.8, color='red', s=50, label='Outlier samples')
ax4.set_xlabel('Mean Abundance')
ax4.set_ylabel('Proteins Detected')
ax4.set_title('Sample Quality Overview')
ax4.legend()

# Plot 5: Quality by tau status
ax5 = axes[1, 1]
tau_groups = sample_qc.groupby('tau_status')['n_proteins'].apply(list).to_dict()
ax5.boxplot(tau_groups.values(), labels=tau_groups.keys())
ax5.set_ylabel('Proteins Detected')
ax5.set_title('Quality by Tau Status')

# Plot 6: Quality by batch
ax6 = axes[1, 2]
if 'batch' in sample_qc.columns:
    batch_groups = sample_qc.groupby('batch')['mean_abundance'].apply(list).to_dict()
    ax6.boxplot(batch_groups.values(), labels=batch_groups.keys())
    ax6.set_ylabel('Mean Abundance')
    ax6.set_title('Quality by Batch')
    ax6.tick_params(axis='x', rotation=45)
else:
    ax6.text(0.5, 0.5, 'No batch info available', ha='center', va='center', transform=ax6.transAxes)
    ax6.set_title('Batch Information')

plt.tight_layout()
plt.show()

print("📊 Sample quality plots created!")
```

### Step 3: Protein-Level Quality Control

#### Basic Protein Statistics
```python
print("\n" + "=" * 50)
print("PROTEIN-LEVEL QUALITY CONTROL")
print("=" * 50)

# Calculate protein-level metrics
protein_qc = pd.DataFrame(index=adata.var_names)

# Detection rate (fraction of samples where protein is detected)
protein_qc['detection_rate'] = np.mean(adata.X > 0, axis=0)

# Mean expression across all samples
protein_qc['mean_expression'] = np.mean(adata.X, axis=0)

# Standard deviation across samples
protein_qc['expression_std'] = np.std(adata.X, axis=0)

# Coefficient of variation
protein_qc['expression_cv'] = protein_qc['expression_std'] / protein_qc['mean_expression']

# Number of zero values
protein_qc['n_zeros'] = np.sum(adata.X == 0, axis=0)

# Fraction of zeros
protein_qc['zero_fraction'] = protein_qc['n_zeros'] / adata.n_obs

print("Protein quality metrics summary:")
print(protein_qc[['detection_rate', 'mean_expression', 'expression_cv']].describe())
```

#### Identify Low-Quality Proteins
```python
print("\nIDENTIFYING LOW-QUALITY PROTEINS")
print("-" * 30)

# Quality thresholds for proteins
min_detection_rate = 0.25  # Protein must be detected in at least 25% of samples
max_zero_fraction = 0.75   # Less than 75% zeros
min_mean_expression = np.percentile(protein_qc['mean_expression'], 5)  # 5th percentile
max_cv = np.percentile(protein_qc['expression_cv'], 95)  # 95th percentile CV

print(f"Protein quality thresholds:")
print(f"Minimum detection rate: {min_detection_rate}")
print(f"Maximum zero fraction: {max_zero_fraction}")
print(f"Minimum mean expression: {min_mean_expression:.2f}")
print(f"Maximum CV: {max_cv:.2f}")

# Flag low-quality proteins
low_quality_flags = (
    (protein_qc['detection_rate'] < min_detection_rate) |
    (protein_qc['zero_fraction'] > max_zero_fraction) |
    (protein_qc['mean_expression'] < min_mean_expression) |
    (protein_qc['expression_cv'] > max_cv)
)

protein_qc['is_low_quality'] = low_quality_flags
low_quality_proteins = protein_qc[low_quality_flags]

print(f"\nLow-quality proteins: {len(low_quality_proteins)}/{len(protein_qc)} ({len(low_quality_proteins)/len(protein_qc)*100:.1f}%)")

# High-quality proteins for analysis
high_quality_proteins = protein_qc[~low_quality_flags]
print(f"High-quality proteins: {len(high_quality_proteins)}/{len(protein_qc)} ({len(high_quality_proteins)/len(protein_qc)*100:.1f}%)")

if len(low_quality_proteins) > 0:
    print(f"\nExamples of low-quality proteins:")
    print(low_quality_proteins[['detection_rate', 'mean_expression', 'expression_cv']].head())
```

#### Visualize Protein Quality
```python
# Create protein quality plots
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('Protein Quality Control Assessment', fontsize=16, y=1.02)

# Plot 1: Detection rate distribution
ax1 = axes[0, 0]
ax1.hist(protein_qc['detection_rate'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
ax1.axvline(min_detection_rate, color='red', linestyle='--', label=f'Threshold: {min_detection_rate}')
ax1.set_xlabel('Detection Rate')
ax1.set_ylabel('Number of Proteins')
ax1.set_title('Protein Detection Rates')
ax1.legend()

# Plot 2: Mean expression distribution
ax2 = axes[0, 1]
ax2.hist(protein_qc['mean_expression'], bins=50, alpha=0.7, color='lightgreen', edgecolor='black')
ax2.axvline(min_mean_expression, color='red', linestyle='--')
ax2.set_xlabel('Mean Expression (log2)')
ax2.set_ylabel('Number of Proteins')
ax2.set_title('Mean Expression Distribution')

# Plot 3: Coefficient of variation
ax3 = axes[0, 2]
# Remove infinite and very high CV values for visualization
cv_for_plot = protein_qc['expression_cv'].replace([np.inf, -np.inf], np.nan).dropna()
cv_for_plot = cv_for_plot[cv_for_plot < np.percentile(cv_for_plot, 99)]

ax3.hist(cv_for_plot, bins=50, alpha=0.7, color='lightsalmon', edgecolor='black')
ax3.set_xlabel('Coefficient of Variation')
ax3.set_ylabel('Number of Proteins')
ax3.set_title('Expression Variability')

# Plot 4: Detection rate vs mean expression
ax4 = axes[1, 0]
high_qual = protein_qc[~protein_qc['is_low_quality']]
low_qual = protein_qc[protein_qc['is_low_quality']]

ax4.scatter(high_qual['mean_expression'], high_qual['detection_rate'],
           alpha=0.6, color='blue', s=10, label='High quality')
if len(low_qual) > 0:
    ax4.scatter(low_qual['mean_expression'], low_qual['detection_rate'],
               alpha=0.6, color='red', s=10, label='Low quality')

ax4.axhline(min_detection_rate, color='red', linestyle='--', alpha=0.7)
ax4.axvline(min_mean_expression, color='red', linestyle='--', alpha=0.7)
ax4.set_xlabel('Mean Expression')
ax4.set_ylabel('Detection Rate')
ax4.set_title('Protein Quality Overview')
ax4.legend()

# Plot 5: Quality protein counts
ax5 = axes[1, 1]
quality_counts = protein_qc['is_low_quality'].value_counts()
labels = ['High Quality', 'Low Quality']
ax5.pie([quality_counts[False], quality_counts[True]], labels=labels,
        autopct='%1.1f%%', colors=['lightblue', 'lightcoral'])
ax5.set_title('Protein Quality Distribution')

# Plot 6: Expression range by quality
ax6 = axes[1, 2]
if len(low_qual) > 0:
    ax6.boxplot([high_qual['mean_expression'], low_qual['mean_expression']],
               labels=['High Quality', 'Low Quality'])
    ax6.set_ylabel('Mean Expression')
    ax6.set_title('Expression by Quality')
else:
    ax6.text(0.5, 0.5, 'All proteins pass quality filters',
             ha='center', va='center', transform=ax6.transAxes)

plt.tight_layout()
plt.show()

print("📊 Protein quality plots created!")
```

### Step 4: Technical Artifacts and Batch Effects

#### Assess Batch Effects
```python
print("\n" + "=" * 50)
print("BATCH EFFECTS ASSESSMENT")
print("=" * 50)

if 'batch' in adata.obs.columns:
    # Calculate batch effect strength
    def calculate_batch_effect(expression_data, batch_labels):
        """Calculate variance explained by batch effects"""
        from sklearn.preprocessing import LabelEncoder

        # Encode batch labels
        le = LabelEncoder()
        batch_encoded = le.fit_transform(batch_labels)

        # Calculate total variance and batch variance for each protein
        total_var = np.var(expression_data, axis=0)

        batch_effects = []
        for i in range(expression_data.shape[1]):
            protein_expr = expression_data[:, i]
            # One-way ANOVA to test batch effects
            batch_groups = [protein_expr[batch_labels == batch]
                           for batch in np.unique(batch_labels)]
            try:
                f_stat, p_val = stats.f_oneway(*batch_groups)
                batch_effects.append(p_val)
            except:
                batch_effects.append(1.0)  # No effect if calculation fails

        return np.array(batch_effects)

    batch_p_values = calculate_batch_effect(adata.X, adata.obs['batch'])

    # Proteins significantly affected by batch
    batch_affected = np.sum(batch_p_values < 0.05)
    batch_percentage = (batch_affected / len(batch_p_values)) * 100

    print(f"Proteins significantly affected by batch: {batch_affected}/{len(batch_p_values)} ({batch_percentage:.1f}%)")

    if batch_percentage > 20:
        print("⚠️ Substantial batch effects detected - consider batch correction")
    elif batch_percentage > 10:
        print("⚠️ Moderate batch effects detected - monitor in analysis")
    else:
        print("✅ Minimal batch effects detected")

    # Visualize batch effects
    plt.figure(figsize=(12, 8))

    # Plot 1: PCA colored by batch
    plt.subplot(2, 2, 1)
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(StandardScaler().fit_transform(adata.X))

    batches = adata.obs['batch'].unique()
    colors = plt.cm.Set3(np.linspace(0, 1, len(batches)))

    for batch, color in zip(batches, colors):
        mask = adata.obs['batch'] == batch
        plt.scatter(X_pca[mask, 0], X_pca[mask, 1],
                   label=f'Batch {batch}', alpha=0.7, color=color)

    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    plt.title('PCA by Batch')
    plt.legend()

    # Plot 2: Batch effect p-values
    plt.subplot(2, 2, 2)
    plt.hist(-np.log10(batch_p_values), bins=50, alpha=0.7, color='orange')
    plt.axvline(-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
    plt.xlabel('-log10(p-value)')
    plt.ylabel('Number of Proteins')
    plt.title('Batch Effect Strength')
    plt.legend()

    # Plot 3: Mean expression by batch
    plt.subplot(2, 2, 3)
    batch_means = adata.obs.groupby('batch').apply(
        lambda x: np.mean(adata.X[x.index])
    )
    plt.bar(range(len(batch_means)), batch_means.values)
    plt.xlabel('Batch')
    plt.ylabel('Mean Expression')
    plt.title('Average Expression by Batch')
    plt.xticks(range(len(batch_means)), batch_means.index)

    # Plot 4: Sample clustering by batch
    plt.subplot(2, 2, 4)
    # Hierarchical clustering would go here - simplified version:
    plt.scatter(X_pca[:, 0], X_pca[:, 1],
               c=[hash(str(batch)) for batch in adata.obs['batch']],
               alpha=0.7, cmap='tab10')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('Sample Clustering')

    plt.tight_layout()
    plt.show()

else:
    print("No batch information available for assessment")
```

#### Check for Other Technical Artifacts

```python
print("\nOTHER TECHNICAL ARTIFACTS")
print("-" * 30)

# PMI effects
if 'PMI' in adata.obs.columns:
    pmi_correlations = []
    for i in range(adata.X.shape[1]):
        corr, p_val = stats.pearsonr(adata.obs['PMI'], adata.X[:, i])
        pmi_correlations.append(abs(corr))

    strong_pmi_effects = np.sum(np.array(pmi_correlations) > 0.3)
    print(f"Proteins strongly correlated with PMI (|r| > 0.3): {strong_pmi_effects}/{adata.n_vars}")

    if strong_pmi_effects > adata.n_vars * 0.1:
        print("⚠️ PMI effects detected - consider controlling for PMI in analysis")
    else:
        print("✅ Minimal PMI effects")

# Age effects
if 'age' in adata.obs.columns:
    age_correlations = []
    for i in range(adata.X.shape[1]):
        corr, p_val = stats.pearsonr(adata.obs['age'], adata.X[:, i])
        age_correlations.append(abs(corr))

    strong_age_effects = np.sum(np.array(age_correlations) > 0.3)
    print(f"Proteins strongly correlated with age (|r| > 0.3): {strong_age_effects}/{adata.n_vars}")

    if strong_age_effects > adata.n_vars * 0.1:
        print("⚠️ Age effects detected - consider controlling for age in analysis")
    else:
        print("✅ Minimal age effects")

# Sex effects
if 'sex' in adata.obs.columns:
    sex_effects = []
    for i in range(adata.X.shape[1]):
        male_expr = adata.X[adata.obs['sex'] == 'M', i]
        female_expr = adata.X[adata.obs['sex'] == 'F', i]
        if len(male_expr) > 0 and len(female_expr) > 0:
            t_stat, p_val = stats.ttest_ind(male_expr, female_expr)
            sex_effects.append(p_val)
        else:
            sex_effects.append(1.0)

    sex_affected = np.sum(np.array(sex_effects) < 0.05)
    print(f"Proteins with significant sex differences (p < 0.05): {sex_affected}/{adata.n_vars}")

    if sex_affected > adata.n_vars * 0.1:
        print("⚠️ Sex effects detected - consider controlling for sex in analysis")
    else:
        print("✅ Minimal sex effects")
```

### Step 5: Final Quality Assessment

#### Overall Data Quality Score
```python
print("\n" + "=" * 50)
print("OVERALL DATA QUALITY ASSESSMENT")
print("=" * 50)

# Calculate overall quality metrics
quality_metrics = {
    'total_samples': adata.n_obs,
    'total_proteins': adata.n_vars,
    'outlier_samples': len(outlier_samples),
    'outlier_percentage': len(outlier_samples) / adata.n_obs * 100,
    'low_quality_proteins': len(low_quality_proteins),
    'low_quality_percentage': len(low_quality_proteins) / adata.n_vars * 100,
    'high_quality_proteins': len(high_quality_proteins),
    'usable_proteins_percentage': len(high_quality_proteins) / adata.n_vars * 100,
    'mean_detection_rate': np.mean(protein_qc['detection_rate']),
    'mean_proteins_per_sample': np.mean(sample_qc['n_proteins']),
}

print("QUALITY SUMMARY:")
print(f"📊 Total samples: {quality_metrics['total_samples']}")
print(f"📊 Total proteins: {quality_metrics['total_proteins']}")
print(f"🚨 Outlier samples: {quality_metrics['outlier_samples']} ({quality_metrics['outlier_percentage']:.1f}%)")
print(f"🚨 Low-quality proteins: {quality_metrics['low_quality_proteins']} ({quality_metrics['low_quality_percentage']:.1f}%)")
print(f"✅ High-quality proteins: {quality_metrics['high_quality_proteins']} ({quality_metrics['usable_proteins_percentage']:.1f}%)")
print(f"📈 Mean detection rate: {quality_metrics['mean_detection_rate']:.1%}")
print(f"📈 Mean proteins per sample: {quality_metrics['mean_proteins_per_sample']:.0f}")

# Overall quality score (0-100)
sample_quality_score = max(0, 100 - quality_metrics['outlier_percentage'] * 2)
protein_quality_score = quality_metrics['usable_proteins_percentage']
detection_quality_score = quality_metrics['mean_detection_rate'] * 100

overall_quality_score = np.mean([sample_quality_score, protein_quality_score, detection_quality_score])

print(f"\n🎯 OVERALL QUALITY SCORE: {overall_quality_score:.0f}/100")

if overall_quality_score >= 85:
    quality_status = "EXCELLENT"
    quality_color = "🟢"
elif overall_quality_score >= 70:
    quality_status = "GOOD"
    quality_color = "🟡"
elif overall_quality_score >= 50:
    quality_status = "ACCEPTABLE"
    quality_color = "🟠"
else:
    quality_status = "POOR"
    quality_color = "🔴"

print(f"{quality_color} Data quality: {quality_status}")
```

#### Recommendations Based on Quality Assessment
```python
print("\n" + "=" * 50)
print("QUALITY-BASED RECOMMENDATIONS")
print("=" * 50)

recommendations = []

# Sample recommendations
if quality_metrics['outlier_percentage'] > 10:
    recommendations.append("🔴 HIGH: Remove outlier samples before analysis")
elif quality_metrics['outlier_percentage'] > 5:
    recommendations.append("🟡 MEDIUM: Consider removing outlier samples")
else:
    recommendations.append("🟢 LOW: Sample quality is acceptable")

# Protein recommendations
if quality_metrics['low_quality_percentage'] > 50:
    recommendations.append("🔴 HIGH: Filter out low-quality proteins (>50% are poor quality)")
elif quality_metrics['low_quality_percentage'] > 25:
    recommendations.append("🟡 MEDIUM: Apply protein quality filters")
else:
    recommendations.append("🟢 LOW: Most proteins are high quality")

# Statistical power recommendations
usable_samples = adata.n_obs - len(outlier_samples)
if usable_samples < 20:
    recommendations.append("🔴 HIGH: Very low sample size - results may be unreliable")
elif usable_samples < 50:
    recommendations.append("🟡 MEDIUM: Small sample size - use conservative statistical methods")
else:
    recommendations.append("🟢 LOW: Sample size adequate for analysis")

# Covariate recommendations
if 'batch_percentage' in locals() and batch_percentage > 20:
    recommendations.append("🔴 HIGH: Strong batch effects - batch correction required")
elif 'batch_percentage' in locals() and batch_percentage > 10:
    recommendations.append("🟡 MEDIUM: Consider batch correction")

print("RECOMMENDATIONS:")
for i, rec in enumerate(recommendations, 1):
    print(f"{i}. {rec}")

# Analysis readiness
analysis_ready = (
    overall_quality_score >= 50 and
    quality_metrics['outlier_percentage'] < 20 and
    quality_metrics['usable_proteins_percentage'] > 60 and
    usable_samples >= 20
)

print(f"\n🎯 ANALYSIS READINESS: {'✅ READY' if analysis_ready else '❌ NOT READY'}")

if not analysis_ready:
    print("Address quality issues before proceeding with statistical analysis")
```

### Step 6: Apply Quality Filters

#### Filter Data Based on QC Results
```python
print("\n" + "=" * 50)
print("APPLYING QUALITY FILTERS")
print("=" * 50)

# Create filtered dataset
print("Creating filtered dataset...")

# Remove outlier samples
if len(outlier_samples) > 0 and len(outlier_samples) < adata.n_obs * 0.2:  # Don't remove >20%
    print(f"Removing {len(outlier_samples)} outlier samples")
    good_samples = ~sample_qc['is_outlier']
    adata_filtered = adata[good_samples, :].copy()
else:
    print("Keeping all samples (outlier removal would be too aggressive)")
    adata_filtered = adata.copy()

# Remove low-quality proteins
if len(low_quality_proteins) > 0:
    print(f"Removing {len(low_quality_proteins)} low-quality proteins")
    good_proteins = ~protein_qc['is_low_quality']
    adata_filtered = adata_filtered[:, good_proteins].copy()
else:
    print("Keeping all proteins")

print(f"\nFiltered dataset dimensions: {adata_filtered.shape}")
print(f"Original dataset dimensions: {adata.shape}")
print(f"Samples retained: {adata_filtered.n_obs}/{adata.n_obs} ({adata_filtered.n_obs/adata.n_obs*100:.1f}%)")
print(f"Proteins retained: {adata_filtered.n_vars}/{adata.n_vars} ({adata_filtered.n_vars/adata.n_vars*100:.1f}%)")

# Update metadata
adata_filtered.obs['passed_qc'] = True
adata_filtered.var['passed_qc'] = True

# Add QC metrics to metadata
for col in ['n_proteins', 'mean_abundance', 'abundance_cv']:
    if col in sample_qc.columns:
        matching_samples = adata_filtered.obs_names.intersection(sample_qc.index)
        adata_filtered.obs.loc[matching_samples, col] = sample_qc.loc[matching_samples, col]

for col in ['detection_rate', 'mean_expression', 'expression_cv']:
    if col in protein_qc.columns:
        matching_proteins = adata_filtered.var_names.intersection(protein_qc.index)
        adata_filtered.var.loc[matching_proteins, col] = protein_qc.loc[matching_proteins, col]
```

#### Save Quality Control Results
```python
print("\nSAVING QUALITY CONTROL RESULTS")
print("-" * 30)

# Save QC report
qc_report = {
    'original_dimensions': adata.shape,
    'filtered_dimensions': adata_filtered.shape,
    'quality_metrics': quality_metrics,
    'overall_quality_score': overall_quality_score,
    'recommendations': recommendations,
    'outlier_samples': outlier_samples.index.tolist() if len(outlier_samples) > 0 else [],
    'low_quality_proteins': low_quality_proteins.index.tolist() if len(low_quality_proteins) > 0 else [],
    'analysis_ready': analysis_ready
}

# Save detailed QC data
sample_qc.to_csv('sample_quality_control.csv')
protein_qc.to_csv('protein_quality_control.csv')

# Save QC report
import json
with open('quality_control_report.json', 'w') as f:
    json.dump(qc_report, f, indent=2)

# Save filtered dataset
adata_filtered.write('dataset_quality_filtered.h5ad')

print("✅ Quality control files saved:")
print("  📄 sample_quality_control.csv")
print("  📄 protein_quality_control.csv")
print("  📄 quality_control_report.json")
print("  📊 dataset_quality_filtered.h5ad")
```

---

## 📋 Quality Control Checklist

### Before Analysis
- [ ] **Sample quality assessed**: Outliers identified and handled appropriately
- [ ] **Protein quality assessed**: Low-detection proteins flagged or removed
- [ ] **Technical artifacts evaluated**: Batch effects, PMI effects, age/sex effects
- [ ] **Overall quality score calculated**: Data quality meets analysis standards
- [ ] **Filtering decisions documented**: Clear rationale for inclusions/exclusions

### Quality Thresholds Met
- [ ] **Sample outliers**: <20% of samples flagged as outliers
- [ ] **Protein detection**: >60% of proteins pass quality filters
- [ ] **Detection rate**: Mean protein detection rate >50%
- [ ] **Batch effects**: <20% of proteins significantly affected by batch
- [ ] **Sample size**: Adequate power for planned statistical tests

### Documentation Complete
- [ ] **QC metrics saved**: Detailed quality metrics for all samples and proteins
- [ ] **Decisions recorded**: Documentation of filtering criteria and rationale
- [ ] **Filtered data saved**: Clean dataset ready for analysis
- [ ] **Quality report generated**: Summary of data quality and recommendations

---

## 🚀 Next Steps After Quality Control

### Immediate Actions
1. **Review QC results** and ensure you understand all quality issues
2. **Validate filtering decisions** with biological knowledge
3. **Plan analysis strategy** based on data quality and limitations
4. **Choose appropriate statistical methods** for your data characteristics

### Analysis Considerations
- **Adjust statistical methods** based on sample size and data quality
- **Include relevant covariates** to control for technical effects
- **Use appropriate multiple testing correction** for number of proteins tested
- **Consider effect size thresholds** appropriate for your data quality

### Documentation for Publication
- **Methods section**: Document QC procedures and filtering criteria
- **Supplementary materials**: Include QC plots and detailed metrics
- **Results interpretation**: Consider quality limitations when interpreting findings
- **Reproducibility**: Ensure others can replicate your QC procedures

---

**Congratulations!** You've completed a comprehensive quality control assessment of your proteomics data. Your dataset is now ready for rigorous statistical analysis with confidence in the data quality.

*Next: Begin your focused analyses with [UPS Protein Analysis](../04_group1_analyses/statement1_ups_proteins/step_by_step_analysis.md) or [SQSTM1 Analysis](../04_group1_analyses/statement2_sqstm1_upregulation/step_by_step_analysis.md)*

*Remember: Time spent on quality control is never wasted - it's the foundation of reliable scientific discoveries!* 🔍🧬