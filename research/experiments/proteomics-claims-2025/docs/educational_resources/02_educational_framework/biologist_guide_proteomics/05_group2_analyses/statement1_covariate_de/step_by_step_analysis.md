# 📊 Step-by-Step: Covariate-Adjusted Differential Expression Analysis

## 🎯 Learning Objectives

By the end of this tutorial, you'll be able to:
- ✅ **Perform differential expression** with covariate adjustment
- ✅ **Build and interpret linear models** for proteomics data
- ✅ **Diagnose and handle batch effects** in your dataset
- ✅ **Compare different adjustment strategies** and their impacts
- ✅ **Validate your findings** with sensitivity analyses

---

## 🚀 Complete Analysis Workflow

### Overview of Steps
1. **Data Preparation** - Load and explore your data
2. **Covariate Assessment** - Identify important factors
3. **Model Building** - Create adjusted analyses
4. **Results Interpretation** - Understand the output
5. **Validation** - Confirm robustness

---

## 📂 Step 1: Data Preparation and Exploration

### 1.1 Load Your Data and Libraries

```python
# Import necessary libraries
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.api import OLS
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import fdrcorrection
import warnings
warnings.filterwarnings('ignore')

# Set up plotting style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

print("Libraries loaded successfully! ✓")
```

### 1.2 Load the Proteomics Data

```python
# For Google Colab users:
from google.colab import drive
drive.mount('/content/drive')
data_path = '/content/drive/MyDrive/proteomics_data/alzheimer_data.h5ad'

# For local users:
# data_path = 'data/alzheimer_data.h5ad'

# Load the data
adata = sc.read_h5ad(data_path)

print(f"Data loaded successfully! ✓")
print(f"Shape: {adata.shape[0]} samples × {adata.shape[1]} proteins")
print(f"\nAvailable metadata columns:")
for col in adata.obs.columns:
    print(f"  - {col}")
```

### 1.3 Examine Your Covariates

```python
# Create a summary of all covariates
def summarize_covariates(adata):
    """
    Create a comprehensive summary of all sample metadata
    """
    summary_data = []

    for col in adata.obs.columns:
        col_data = adata.obs[col]

        # Determine data type and create appropriate summary
        if col_data.dtype in ['float64', 'float32', 'int64', 'int32']:
            summary_data.append({
                'Variable': col,
                'Type': 'Continuous',
                'Missing': col_data.isna().sum(),
                'Mean ± SD': f"{col_data.mean():.2f} ± {col_data.std():.2f}",
                'Range': f"{col_data.min():.2f} - {col_data.max():.2f}"
            })
        else:
            value_counts = col_data.value_counts()
            summary_data.append({
                'Variable': col,
                'Type': 'Categorical',
                'Missing': col_data.isna().sum(),
                'Categories': len(value_counts),
                'Distribution': ', '.join([f"{v}:{c}" for v, c in value_counts.head(3).items()])
            })

    summary_df = pd.DataFrame(summary_data)
    return summary_df

# Display covariate summary
covariate_summary = summarize_covariates(adata)
print("\n📊 Covariate Summary:")
print(covariate_summary.to_string())
```

### 1.4 Check Group Balance

```python
# Check if groups are balanced for key covariates
def check_group_balance(adata, group_col='tau_status'):
    """
    Assess whether groups are balanced for important covariates
    """
    print(f"\n🔍 Checking balance between {group_col} groups:\n")

    groups = adata.obs[group_col].unique()

    # Check continuous variables
    continuous_vars = ['age', 'pmi', 'rin', 'pH']
    for var in continuous_vars:
        if var in adata.obs.columns:
            print(f"\n{var.upper()}:")
            for group in groups:
                group_data = adata.obs[adata.obs[group_col] == group][var]
                print(f"  {group}: {group_data.mean():.2f} ± {group_data.std():.2f}")

            # Perform t-test
            if len(groups) == 2:
                g1_data = adata.obs[adata.obs[group_col] == groups[0]][var]
                g2_data = adata.obs[adata.obs[group_col] == groups[1]][var]
                t_stat, p_val = stats.ttest_ind(g1_data.dropna(), g2_data.dropna())
                print(f"  Difference: p = {p_val:.4f} {'⚠️ IMBALANCED' if p_val < 0.05 else '✓ Balanced'}")

    # Check categorical variables
    categorical_vars = ['sex', 'brain_region', 'batch']
    for var in categorical_vars:
        if var in adata.obs.columns:
            print(f"\n{var.upper()}:")
            contingency = pd.crosstab(adata.obs[var], adata.obs[group_col])
            print(contingency)
            chi2, p_val, _, _ = stats.chi2_contingency(contingency)
            print(f"  Chi-square: p = {p_val:.4f} {'⚠️ IMBALANCED' if p_val < 0.05 else '✓ Balanced'}")

check_group_balance(adata)
```

---

## 🔍 Step 2: Exploratory Covariate Analysis

### 2.1 Visualize Covariate Relationships

```python
# Create comprehensive covariate visualization
def visualize_covariates(adata):
    """
    Create multi-panel figure showing covariate relationships
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Covariate Relationships and Distributions', fontsize=16)

    # Age distribution by group
    ax = axes[0, 0]
    for group in adata.obs['tau_status'].unique():
        group_data = adata.obs[adata.obs['tau_status'] == group]['age']
        ax.hist(group_data, alpha=0.6, label=group, bins=20)
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('Count')
    ax.set_title('Age Distribution by Tau Status')
    ax.legend()

    # PMI distribution
    ax = axes[0, 1]
    ax.hist(adata.obs['pmi'], bins=30, color='coral', edgecolor='black')
    ax.set_xlabel('Post-Mortem Interval (hours)')
    ax.set_ylabel('Count')
    ax.set_title('PMI Distribution')

    # Sex balance
    ax = axes[0, 2]
    sex_counts = pd.crosstab(adata.obs['sex'], adata.obs['tau_status'])
    sex_counts.plot(kind='bar', ax=ax, rot=0)
    ax.set_xlabel('Sex')
    ax.set_ylabel('Count')
    ax.set_title('Sex Distribution by Group')
    ax.legend(title='Tau Status')

    # Batch effects visualization
    ax = axes[1, 0]
    if 'batch' in adata.obs.columns:
        batch_counts = pd.crosstab(adata.obs['batch'], adata.obs['tau_status'])
        batch_counts.plot(kind='bar', ax=ax, stacked=True)
        ax.set_xlabel('Batch')
        ax.set_ylabel('Count')
        ax.set_title('Batch Distribution')
        ax.legend(title='Tau Status')

    # Age vs PMI correlation
    ax = axes[1, 1]
    ax.scatter(adata.obs['age'], adata.obs['pmi'], alpha=0.6)
    ax.set_xlabel('Age (years)')
    ax.set_ylabel('PMI (hours)')
    ax.set_title('Age vs PMI Correlation')

    # Correlation matrix
    ax = axes[1, 2]
    numeric_cols = ['age', 'pmi', 'rin', 'pH']
    numeric_data = adata.obs[numeric_cols].dropna()
    corr_matrix = numeric_data.corr()
    im = ax.imshow(corr_matrix, cmap='coolwarm', vmin=-1, vmax=1)
    ax.set_xticks(range(len(numeric_cols)))
    ax.set_yticks(range(len(numeric_cols)))
    ax.set_xticklabels(numeric_cols, rotation=45)
    ax.set_yticklabels(numeric_cols)
    ax.set_title('Covariate Correlations')
    plt.colorbar(im, ax=ax)

    # Add correlation values
    for i in range(len(numeric_cols)):
        for j in range(len(numeric_cols)):
            text = ax.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}',
                          ha="center", va="center", color="black", fontsize=8)

    plt.tight_layout()
    plt.show()

visualize_covariates(adata)
```

### 2.2 Detect Batch Effects

```python
# Comprehensive batch effect detection
def detect_batch_effects(adata, n_top_proteins=50):
    """
    Multiple methods to detect batch effects
    """
    print("🔍 Detecting Batch Effects...\n")

    # Method 1: PCA visualization
    print("1. PCA-based detection:")

    # Perform PCA
    sc.pp.scale(adata)
    sc.pp.pca(adata)

    # Plot PCA colored by batch
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Color by batch
    ax = axes[0]
    for batch in adata.obs['batch'].unique():
        batch_mask = adata.obs['batch'] == batch
        ax.scatter(adata.obsm['X_pca'][batch_mask, 0],
                  adata.obsm['X_pca'][batch_mask, 1],
                  label=f'Batch {batch}', alpha=0.6)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title('PCA colored by Batch')
    ax.legend()

    # Color by biological group
    ax = axes[1]
    for group in adata.obs['tau_status'].unique():
        group_mask = adata.obs['tau_status'] == group
        ax.scatter(adata.obsm['X_pca'][group_mask, 0],
                  adata.obsm['X_pca'][group_mask, 1],
                  label=group, alpha=0.6)
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title('PCA colored by Tau Status')
    ax.legend()

    # Method 2: ANOVA for batch effects
    print("\n2. ANOVA-based detection (top proteins):")

    # Test top variable proteins
    protein_vars = adata.var['variance'].nlargest(n_top_proteins).index
    batch_pvals = []

    for protein in protein_vars:
        protein_data = pd.DataFrame({
            'expression': adata[:, protein].X.flatten(),
            'batch': adata.obs['batch'].values
        })

        f_stat, p_val = stats.f_oneway(*[group['expression'].values
                                         for name, group in protein_data.groupby('batch')])
        batch_pvals.append(p_val)

    significant_batch_effects = sum([p < 0.05 for p in batch_pvals])
    print(f"  Proteins with significant batch effects: {significant_batch_effects}/{n_top_proteins}")
    print(f"  Percentage: {100*significant_batch_effects/n_top_proteins:.1f}%")

    if significant_batch_effects/n_top_proteins > 0.3:
        print("  ⚠️ Strong batch effects detected! Adjustment recommended.")
    else:
        print("  ✓ Minimal batch effects detected.")

    # Method 3: Variance decomposition
    ax = axes[2]
    ax.hist(batch_pvals, bins=20, edgecolor='black')
    ax.axvline(x=0.05, color='red', linestyle='--', label='p=0.05')
    ax.set_xlabel('Batch Effect P-values')
    ax.set_ylabel('Number of Proteins')
    ax.set_title('Distribution of Batch Effect P-values')
    ax.legend()

    plt.tight_layout()
    plt.show()

    return batch_pvals

batch_pvals = detect_batch_effects(adata)
```

---

## 🏗️ Step 3: Building Linear Models with Covariates

### 3.1 Simple vs Adjusted Models

```python
def compare_models(adata, protein_name, covariates=[]):
    """
    Compare simple and adjusted differential expression models
    """
    print(f"\n📊 Analyzing: {protein_name}")
    print("=" * 50)

    # Extract protein data
    if protein_name not in adata.var_names:
        print(f"❌ Protein {protein_name} not found!")
        return None

    protein_expr = adata[:, protein_name].X.flatten()

    # Create dataframe for modeling
    model_data = pd.DataFrame({
        'expression': protein_expr,
        'tau_status': adata.obs['tau_status'].values,
        'age': adata.obs['age'].values,
        'sex': adata.obs['sex'].values,
        'pmi': adata.obs['pmi'].values,
        'batch': adata.obs['batch'].values
    })

    # Convert tau_status to numeric (0/1)
    model_data['tau_binary'] = (model_data['tau_status'] == 'tau_positive').astype(int)

    results = {}

    # Model 1: Simple (unadjusted)
    print("\n🔹 Model 1: Simple (Unadjusted)")
    print("  Formula: expression ~ tau_status")

    simple_model = ols('expression ~ tau_binary', data=model_data).fit()
    tau_effect_simple = simple_model.params['tau_binary']
    tau_pval_simple = simple_model.pvalues['tau_binary']

    print(f"  Tau Effect: {tau_effect_simple:.3f}")
    print(f"  P-value: {tau_pval_simple:.4f}")
    print(f"  R-squared: {simple_model.rsquared:.3f}")

    results['simple'] = {
        'effect': tau_effect_simple,
        'pvalue': tau_pval_simple,
        'rsquared': simple_model.rsquared
    }

    # Model 2: Age-adjusted
    print("\n🔹 Model 2: Age-Adjusted")
    print("  Formula: expression ~ tau_status + age")

    age_model = ols('expression ~ tau_binary + age', data=model_data).fit()
    tau_effect_age = age_model.params['tau_binary']
    tau_pval_age = age_model.pvalues['tau_binary']

    print(f"  Tau Effect: {tau_effect_age:.3f}")
    print(f"  P-value: {tau_pval_age:.4f}")
    print(f"  R-squared: {age_model.rsquared:.3f}")
    print(f"  Age coefficient: {age_model.params['age']:.4f} (p={age_model.pvalues['age']:.4f})")

    results['age_adjusted'] = {
        'effect': tau_effect_age,
        'pvalue': tau_pval_age,
        'rsquared': age_model.rsquared
    }

    # Model 3: Fully adjusted
    print("\n🔹 Model 3: Fully Adjusted")
    print("  Formula: expression ~ tau_status + age + sex + pmi + batch")

    # Convert categorical variables to dummy variables
    model_data_encoded = pd.get_dummies(model_data, columns=['sex', 'batch'], drop_first=True)

    formula = 'expression ~ tau_binary + age + pmi'
    sex_cols = [col for col in model_data_encoded.columns if col.startswith('sex_')]
    batch_cols = [col for col in model_data_encoded.columns if col.startswith('batch_')]

    if sex_cols:
        formula += ' + ' + ' + '.join(sex_cols)
    if batch_cols:
        formula += ' + ' + ' + '.join(batch_cols)

    full_model = ols(formula, data=model_data_encoded).fit()
    tau_effect_full = full_model.params['tau_binary']
    tau_pval_full = full_model.pvalues['tau_binary']

    print(f"  Tau Effect: {tau_effect_full:.3f}")
    print(f"  P-value: {tau_pval_full:.4f}")
    print(f"  R-squared: {full_model.rsquared:.3f}")

    results['full_adjusted'] = {
        'effect': tau_effect_full,
        'pvalue': tau_pval_full,
        'rsquared': full_model.rsquared
    }

    # Summary comparison
    print("\n📈 Model Comparison Summary:")
    print("-" * 50)
    print(f"{'Model':<20} {'Tau Effect':<12} {'P-value':<12} {'R²':<8}")
    print("-" * 50)
    for model_name, res in results.items():
        sig = "***" if res['pvalue'] < 0.001 else "**" if res['pvalue'] < 0.01 else "*" if res['pvalue'] < 0.05 else ""
        print(f"{model_name:<20} {res['effect']:>11.3f} {res['pvalue']:>11.4f}{sig:<3} {res['rsquared']:>7.3f}")

    # Visualization
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f'{protein_name}: Model Comparison', fontsize=14)

    # Effect sizes
    ax = axes[0]
    models = list(results.keys())
    effects = [results[m]['effect'] for m in models]
    colors = ['red' if results[m]['pvalue'] < 0.05 else 'gray' for m in models]
    ax.bar(range(len(models)), effects, color=colors)
    ax.set_xticks(range(len(models)))
    ax.set_xticklabels(models, rotation=45)
    ax.set_ylabel('Tau Effect Size')
    ax.set_title('Effect Size Comparison')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)

    # P-values
    ax = axes[1]
    pvals = [results[m]['pvalue'] for m in models]
    ax.bar(range(len(models)), [-np.log10(p) for p in pvals], color=colors)
    ax.set_xticks(range(len(models)))
    ax.set_xticklabels(models, rotation=45)
    ax.set_ylabel('-log10(p-value)')
    ax.set_title('Statistical Significance')
    ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
    ax.legend()

    # R-squared
    ax = axes[2]
    rsquared = [results[m]['rsquared'] for m in models]
    ax.bar(range(len(models)), rsquared, color='steelblue')
    ax.set_xticks(range(len(models)))
    ax.set_xticklabels(models, rotation=45)
    ax.set_ylabel('R-squared')
    ax.set_title('Model Fit')
    ax.set_ylim(0, max(rsquared) * 1.2)

    plt.tight_layout()
    plt.show()

    return results

# Example analysis with a key protein
results = compare_models(adata, 'SQSTM1')
```

### 3.2 Proteome-Wide Adjusted Analysis

```python
def proteome_wide_adjusted_de(adata, formula='expression ~ tau_binary + age + sex + pmi + batch',
                              n_proteins=None, show_progress=True):
    """
    Perform adjusted differential expression for all proteins
    """
    if n_proteins is None:
        n_proteins = adata.shape[1]

    proteins = adata.var_names[:n_proteins]

    results = []
    failed_proteins = []

    print(f"🚀 Analyzing {n_proteins} proteins with covariates...")
    print(f"Formula: {formula}")
    print("-" * 50)

    for i, protein in enumerate(proteins):
        if show_progress and i % 100 == 0:
            print(f"Progress: {i}/{n_proteins} proteins analyzed...")

        try:
            # Extract protein expression
            protein_expr = adata[:, protein].X.flatten()

            # Create model data
            model_data = pd.DataFrame({
                'expression': protein_expr,
                'tau_binary': (adata.obs['tau_status'] == 'tau_positive').astype(int),
                'age': adata.obs['age'].values,
                'sex': adata.obs['sex'].values,
                'pmi': adata.obs['pmi'].values,
                'batch': adata.obs['batch'].values
            })

            # Convert categorical variables
            model_data_encoded = pd.get_dummies(model_data, columns=['sex', 'batch'], drop_first=True)

            # Build formula dynamically
            actual_formula = 'expression ~ tau_binary + age + pmi'
            sex_cols = [col for col in model_data_encoded.columns if col.startswith('sex_')]
            batch_cols = [col for col in model_data_encoded.columns if col.startswith('batch_')]

            if sex_cols:
                actual_formula += ' + ' + ' + '.join(sex_cols)
            if batch_cols:
                actual_formula += ' + ' + ' + '.join(batch_cols)

            # Fit model
            model = ols(actual_formula, data=model_data_encoded).fit()

            # Extract results
            tau_effect = model.params['tau_binary']
            tau_se = model.bse['tau_binary']
            tau_pval = model.pvalues['tau_binary']

            # Calculate confidence interval
            ci_lower = tau_effect - 1.96 * tau_se
            ci_upper = tau_effect + 1.96 * tau_se

            # Calculate fold change
            mean_tau_pos = model_data[model_data['tau_binary'] == 1]['expression'].mean()
            mean_tau_neg = model_data[model_data['tau_binary'] == 0]['expression'].mean()
            fold_change = mean_tau_pos / mean_tau_neg if mean_tau_neg != 0 else np.nan

            results.append({
                'protein': protein,
                'coefficient': tau_effect,
                'std_error': tau_se,
                'pvalue': tau_pval,
                'ci_lower': ci_lower,
                'ci_upper': ci_upper,
                'fold_change': fold_change,
                'rsquared': model.rsquared,
                'mean_tau_pos': mean_tau_pos,
                'mean_tau_neg': mean_tau_neg
            })

        except Exception as e:
            failed_proteins.append(protein)
            if len(failed_proteins) <= 5:
                print(f"  Warning: Failed to analyze {protein}: {str(e)}")

    # Create results dataframe
    results_df = pd.DataFrame(results)

    # Apply FDR correction
    results_df['pvalue_adj'] = fdrcorrection(results_df['pvalue'])[1]

    # Add significance categories
    results_df['significant'] = results_df['pvalue_adj'] < 0.05
    results_df['log2fc'] = np.log2(results_df['fold_change'])

    # Sort by p-value
    results_df = results_df.sort_values('pvalue')

    print(f"\n✅ Analysis complete!")
    print(f"  Successfully analyzed: {len(results)}/{n_proteins} proteins")
    if failed_proteins:
        print(f"  Failed proteins: {len(failed_proteins)}")

    # Summary statistics
    print(f"\n📊 Results Summary:")
    print(f"  Significant proteins (FDR < 0.05): {results_df['significant'].sum()}")
    print(f"  Upregulated (FC > 1.5, FDR < 0.05): {((results_df['fold_change'] > 1.5) & results_df['significant']).sum()}")
    print(f"  Downregulated (FC < 0.67, FDR < 0.05): {((results_df['fold_change'] < 0.67) & results_df['significant']).sum()}")

    return results_df

# Run proteome-wide analysis
de_results = proteome_wide_adjusted_de(adata, n_proteins=100)  # Start with 100 for testing
```

---

## 📊 Step 4: Visualizing and Interpreting Results

### 4.1 Volcano Plot with Adjustments

```python
def plot_volcano_adjusted(results_df, title='Adjusted Differential Expression'):
    """
    Create volcano plot for adjusted DE results
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Define colors
    colors = []
    for _, row in results_df.iterrows():
        if row['pvalue_adj'] < 0.05:
            if row['log2fc'] > 0.5:
                colors.append('red')  # Upregulated
            elif row['log2fc'] < -0.5:
                colors.append('blue')  # Downregulated
            else:
                colors.append('orange')  # Significant but small effect
        else:
            colors.append('gray')  # Not significant

    # Create scatter plot
    scatter = ax.scatter(results_df['log2fc'],
                        -np.log10(results_df['pvalue_adj']),
                        c=colors, alpha=0.6, s=30)

    # Add significance thresholds
    ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='FDR = 0.05')
    ax.axvline(x=0.5, color='gray', linestyle='--', alpha=0.3)
    ax.axvline(x=-0.5, color='gray', linestyle='--', alpha=0.3)

    # Label top proteins
    top_proteins = results_df.nsmallest(10, 'pvalue_adj')
    for _, row in top_proteins.iterrows():
        ax.annotate(row['protein'],
                   xy=(row['log2fc'], -np.log10(row['pvalue_adj'])),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=8, alpha=0.7)

    ax.set_xlabel('Log2 Fold Change', fontsize=12)
    ax.set_ylabel('-Log10 Adjusted P-value', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')

    # Add legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8, label='Upregulated'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8, label='Downregulated'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='orange', markersize=8, label='Small effect'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=8, label='Not significant')
    ]
    ax.legend(handles=legend_elements, loc='upper left')

    # Add statistics text
    n_up = ((results_df['log2fc'] > 0.5) & (results_df['pvalue_adj'] < 0.05)).sum()
    n_down = ((results_df['log2fc'] < -0.5) & (results_df['pvalue_adj'] < 0.05)).sum()
    ax.text(0.02, 0.98, f'Upregulated: {n_up}\nDownregulated: {n_down}',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.show()

plot_volcano_adjusted(de_results)
```

### 4.2 Compare Adjusted vs Unadjusted

```python
def compare_adjusted_unadjusted(adata, n_proteins=100):
    """
    Direct comparison of adjusted vs unadjusted results
    """
    print("🔄 Comparing adjusted vs unadjusted analyses...")

    # Run both analyses
    results_comparison = []

    for protein in adata.var_names[:n_proteins]:
        try:
            protein_expr = adata[:, protein].X.flatten()

            # Create model data
            model_data = pd.DataFrame({
                'expression': protein_expr,
                'tau_binary': (adata.obs['tau_status'] == 'tau_positive').astype(int),
                'age': adata.obs['age'].values,
                'sex': adata.obs['sex'].values,
                'pmi': adata.obs['pmi'].values,
                'batch': adata.obs['batch'].values
            })

            # Unadjusted model
            unadj_model = ols('expression ~ tau_binary', data=model_data).fit()
            unadj_coef = unadj_model.params['tau_binary']
            unadj_pval = unadj_model.pvalues['tau_binary']

            # Adjusted model
            model_data_encoded = pd.get_dummies(model_data, columns=['sex', 'batch'], drop_first=True)
            formula = 'expression ~ tau_binary + age + pmi'
            sex_cols = [col for col in model_data_encoded.columns if col.startswith('sex_')]
            batch_cols = [col for col in model_data_encoded.columns if col.startswith('batch_')]
            if sex_cols:
                formula += ' + ' + ' + '.join(sex_cols)
            if batch_cols:
                formula += ' + ' + ' + '.join(batch_cols)

            adj_model = ols(formula, data=model_data_encoded).fit()
            adj_coef = adj_model.params['tau_binary']
            adj_pval = adj_model.pvalues['tau_binary']

            results_comparison.append({
                'protein': protein,
                'unadj_coef': unadj_coef,
                'unadj_pval': unadj_pval,
                'adj_coef': adj_coef,
                'adj_pval': adj_pval,
                'coef_change': adj_coef - unadj_coef,
                'pval_ratio': adj_pval / unadj_pval if unadj_pval > 0 else np.inf
            })

        except:
            continue

    comp_df = pd.DataFrame(results_comparison)

    # Apply FDR correction
    comp_df['unadj_pval_adj'] = fdrcorrection(comp_df['unadj_pval'])[1]
    comp_df['adj_pval_adj'] = fdrcorrection(comp_df['adj_pval'])[1]

    # Categorize changes
    comp_df['category'] = 'No change'
    comp_df.loc[(comp_df['unadj_pval_adj'] < 0.05) & (comp_df['adj_pval_adj'] >= 0.05), 'category'] = 'Lost significance'
    comp_df.loc[(comp_df['unadj_pval_adj'] >= 0.05) & (comp_df['adj_pval_adj'] < 0.05), 'category'] = 'Gained significance'
    comp_df.loc[(comp_df['unadj_pval_adj'] < 0.05) & (comp_df['adj_pval_adj'] < 0.05), 'category'] = 'Remained significant'

    # Visualization
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Impact of Covariate Adjustment', fontsize=14, fontweight='bold')

    # Coefficient comparison
    ax = axes[0, 0]
    ax.scatter(comp_df['unadj_coef'], comp_df['adj_coef'], alpha=0.5)
    ax.plot([comp_df['unadj_coef'].min(), comp_df['unadj_coef'].max()],
            [comp_df['unadj_coef'].min(), comp_df['unadj_coef'].max()],
            'r--', label='y=x')
    ax.set_xlabel('Unadjusted Coefficient')
    ax.set_ylabel('Adjusted Coefficient')
    ax.set_title('Effect Size Comparison')
    ax.legend()

    # P-value comparison
    ax = axes[0, 1]
    ax.scatter(-np.log10(comp_df['unadj_pval']), -np.log10(comp_df['adj_pval']), alpha=0.5)
    ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.3)
    ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.3)
    ax.set_xlabel('-log10(Unadjusted P-value)')
    ax.set_ylabel('-log10(Adjusted P-value)')
    ax.set_title('Statistical Significance Comparison')

    # Category distribution
    ax = axes[1, 0]
    category_counts = comp_df['category'].value_counts()
    colors = {'No change': 'gray', 'Lost significance': 'red',
              'Gained significance': 'green', 'Remained significant': 'blue'}
    ax.pie(category_counts.values, labels=category_counts.index, autopct='%1.1f%%',
           colors=[colors[cat] for cat in category_counts.index])
    ax.set_title('Impact on Significance')

    # Effect size changes
    ax = axes[1, 1]
    ax.hist(comp_df['coef_change'], bins=30, edgecolor='black')
    ax.axvline(x=0, color='red', linestyle='--')
    ax.set_xlabel('Change in Coefficient (Adj - Unadj)')
    ax.set_ylabel('Number of Proteins')
    ax.set_title('Distribution of Effect Size Changes')

    plt.tight_layout()
    plt.show()

    # Print summary
    print("\n📊 Adjustment Impact Summary:")
    print("-" * 40)
    for category, count in category_counts.items():
        print(f"{category}: {count} proteins ({100*count/len(comp_df):.1f}%)")

    print(f"\nMean coefficient change: {comp_df['coef_change'].mean():.3f}")
    print(f"Median p-value ratio: {comp_df['pval_ratio'].median():.2f}")

    return comp_df

comparison_results = compare_adjusted_unadjusted(adata)
```

---

## 🔍 Step 5: Sensitivity and Validation Analyses

### 5.1 Leave-One-Covariate-Out Analysis

```python
def sensitivity_analysis(adata, protein_name):
    """
    Test sensitivity to individual covariates
    """
    print(f"\n🔬 Sensitivity Analysis for {protein_name}")
    print("=" * 50)

    if protein_name not in adata.var_names:
        print(f"❌ Protein {protein_name} not found!")
        return None

    protein_expr = adata[:, protein_name].X.flatten()

    # Create full model data
    model_data = pd.DataFrame({
        'expression': protein_expr,
        'tau_binary': (adata.obs['tau_status'] == 'tau_positive').astype(int),
        'age': adata.obs['age'].values,
        'sex': adata.obs['sex'].values,
        'pmi': adata.obs['pmi'].values,
        'batch': adata.obs['batch'].values
    })

    model_data_encoded = pd.get_dummies(model_data, columns=['sex', 'batch'], drop_first=True)

    # Full model
    full_formula = 'expression ~ tau_binary + age + pmi'
    sex_cols = [col for col in model_data_encoded.columns if col.startswith('sex_')]
    batch_cols = [col for col in model_data_encoded.columns if col.startswith('batch_')]
    if sex_cols:
        full_formula += ' + ' + ' + '.join(sex_cols)
    if batch_cols:
        full_formula += ' + ' + ' + '.join(batch_cols)

    full_model = ols(full_formula, data=model_data_encoded).fit()
    full_coef = full_model.params['tau_binary']
    full_pval = full_model.pvalues['tau_binary']

    print(f"Full Model: coef={full_coef:.3f}, p={full_pval:.4f}")
    print("\nLeave-one-out results:")

    sensitivity_results = []

    # Test each covariate
    covariates_to_test = ['age', 'pmi'] + sex_cols + batch_cols

    for covariate in covariates_to_test:
        # Create formula without this covariate
        reduced_formula = full_formula.replace(f' + {covariate}', '').replace(f'{covariate} + ', '')

        reduced_model = ols(reduced_formula, data=model_data_encoded).fit()
        reduced_coef = reduced_model.params['tau_binary']
        reduced_pval = reduced_model.pvalues['tau_binary']

        change_pct = 100 * (reduced_coef - full_coef) / abs(full_coef) if full_coef != 0 else 0

        sensitivity_results.append({
            'excluded': covariate,
            'coefficient': reduced_coef,
            'pvalue': reduced_pval,
            'change_pct': change_pct
        })

        impact = "HIGH" if abs(change_pct) > 20 else "MEDIUM" if abs(change_pct) > 10 else "LOW"
        print(f"  Without {covariate:15}: coef={reduced_coef:.3f}, p={reduced_pval:.4f}, "
              f"change={change_pct:+.1f}% [{impact}]")

    # Visualization
    sens_df = pd.DataFrame(sensitivity_results)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f'Sensitivity Analysis: {protein_name}', fontsize=14)

    # Coefficient changes
    ax = axes[0]
    colors = ['red' if abs(c) > 20 else 'orange' if abs(c) > 10 else 'green'
              for c in sens_df['change_pct']]
    ax.barh(range(len(sens_df)), sens_df['change_pct'], color=colors)
    ax.set_yticks(range(len(sens_df)))
    ax.set_yticklabels(sens_df['excluded'])
    ax.set_xlabel('% Change in Coefficient')
    ax.set_title('Impact of Removing Each Covariate')
    ax.axvline(x=0, color='black', linewidth=0.5)

    # P-value changes
    ax = axes[1]
    ax.barh(range(len(sens_df)), -np.log10(sens_df['pvalue']))
    ax.set_yticks(range(len(sens_df)))
    ax.set_yticklabels(sens_df['excluded'])
    ax.set_xlabel('-log10(p-value)')
    ax.set_title('P-values Without Each Covariate')
    ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.show()

    return sens_df

# Run sensitivity analysis
sensitivity_results = sensitivity_analysis(adata, 'SQSTM1')
```

### 5.2 Stratified Analysis

```python
def stratified_analysis(adata, protein_name):
    """
    Analyze effects within subgroups
    """
    print(f"\n📊 Stratified Analysis for {protein_name}")
    print("=" * 50)

    if protein_name not in adata.var_names:
        print(f"❌ Protein {protein_name} not found!")
        return None

    protein_expr = adata[:, protein_name].X.flatten()

    # Prepare data
    analysis_data = pd.DataFrame({
        'expression': protein_expr,
        'tau_status': adata.obs['tau_status'].values,
        'age': adata.obs['age'].values,
        'sex': adata.obs['sex'].values,
        'age_group': pd.cut(adata.obs['age'], bins=[0, 65, 75, 100],
                            labels=['<65', '65-75', '>75'])
    })

    stratified_results = []

    # Stratify by sex
    print("\n🔹 Stratified by Sex:")
    for sex in analysis_data['sex'].unique():
        subset = analysis_data[analysis_data['sex'] == sex]

        tau_pos = subset[subset['tau_status'] == 'tau_positive']['expression']
        tau_neg = subset[subset['tau_status'] == 'tau_negative']['expression']

        if len(tau_pos) > 0 and len(tau_neg) > 0:
            t_stat, p_val = stats.ttest_ind(tau_pos, tau_neg)
            effect_size = (tau_pos.mean() - tau_neg.mean()) / subset['expression'].std()

            print(f"  {sex}: Effect={effect_size:.3f}, p={p_val:.4f}, "
                  f"n_tau+={len(tau_pos)}, n_tau-={len(tau_neg)}")

            stratified_results.append({
                'stratum': f'Sex: {sex}',
                'effect_size': effect_size,
                'pvalue': p_val,
                'n_tau_pos': len(tau_pos),
                'n_tau_neg': len(tau_neg)
            })

    # Stratify by age group
    print("\n🔹 Stratified by Age Group:")
    for age_group in analysis_data['age_group'].unique():
        subset = analysis_data[analysis_data['age_group'] == age_group]

        tau_pos = subset[subset['tau_status'] == 'tau_positive']['expression']
        tau_neg = subset[subset['tau_status'] == 'tau_negative']['expression']

        if len(tau_pos) > 0 and len(tau_neg) > 0:
            t_stat, p_val = stats.ttest_ind(tau_pos, tau_neg)
            effect_size = (tau_pos.mean() - tau_neg.mean()) / subset['expression'].std()

            print(f"  {age_group}: Effect={effect_size:.3f}, p={p_val:.4f}, "
                  f"n_tau+={len(tau_pos)}, n_tau-={len(tau_neg)}")

            stratified_results.append({
                'stratum': f'Age: {age_group}',
                'effect_size': effect_size,
                'pvalue': p_val,
                'n_tau_pos': len(tau_pos),
                'n_tau_neg': len(tau_neg)
            })

    # Visualization
    strat_df = pd.DataFrame(stratified_results)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(f'Stratified Analysis: {protein_name}', fontsize=14)

    # Forest plot
    ax = axes[0]
    y_pos = range(len(strat_df))
    ax.errorbar(strat_df['effect_size'], y_pos,
                xerr=1.96 * strat_df['effect_size'].std() / np.sqrt(len(strat_df)),
                fmt='o', capsize=5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(strat_df['stratum'])
    ax.set_xlabel('Effect Size')
    ax.set_title('Forest Plot of Stratified Effects')
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5)

    # Sample sizes
    ax = axes[1]
    x = np.arange(len(strat_df))
    width = 0.35
    ax.bar(x - width/2, strat_df['n_tau_pos'], width, label='Tau+', color='coral')
    ax.bar(x + width/2, strat_df['n_tau_neg'], width, label='Tau-', color='skyblue')
    ax.set_xticks(x)
    ax.set_xticklabels(strat_df['stratum'], rotation=45, ha='right')
    ax.set_ylabel('Sample Size')
    ax.set_title('Sample Sizes by Stratum')
    ax.legend()

    plt.tight_layout()
    plt.show()

    return strat_df

# Run stratified analysis
stratified_results = stratified_analysis(adata, 'SQSTM1')
```

---

## 📈 Step 6: Advanced Topics and Best Practices

### 6.1 Interaction Effects

```python
def test_interaction_effects(adata, protein_name):
    """
    Test for interaction between tau status and covariates
    """
    print(f"\n🔄 Testing Interaction Effects for {protein_name}")
    print("=" * 50)

    protein_expr = adata[:, protein_name].X.flatten()

    model_data = pd.DataFrame({
        'expression': protein_expr,
        'tau_binary': (adata.obs['tau_status'] == 'tau_positive').astype(int),
        'age': adata.obs['age'].values,
        'sex_binary': (adata.obs['sex'] == 'M').astype(int)
    })

    # Model with interaction
    interaction_model = ols('expression ~ tau_binary * age + sex_binary',
                           data=model_data).fit()

    # Model without interaction
    no_interaction_model = ols('expression ~ tau_binary + age + sex_binary',
                              data=model_data).fit()

    # Compare models
    interaction_pval = interaction_model.pvalues['tau_binary:age']

    print(f"\n🔹 Tau × Age Interaction:")
    print(f"  Coefficient: {interaction_model.params['tau_binary:age']:.4f}")
    print(f"  P-value: {interaction_pval:.4f}")

    if interaction_pval < 0.05:
        print("  ⚠️ SIGNIFICANT interaction detected!")
        print("  Interpretation: The effect of tau pathology depends on age")
    else:
        print("  ✓ No significant interaction")

    # Visualization
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot by tau status
    for tau_status in [0, 1]:
        subset = model_data[model_data['tau_binary'] == tau_status]
        ax.scatter(subset['age'], subset['expression'],
                  label=f'Tau {"+" if tau_status else "-"}',
                  alpha=0.5)

        # Add regression line
        z = np.polyfit(subset['age'], subset['expression'], 1)
        p = np.poly1d(z)
        age_range = np.linspace(subset['age'].min(), subset['age'].max(), 100)
        ax.plot(age_range, p(age_range), "--", alpha=0.8)

    ax.set_xlabel('Age (years)')
    ax.set_ylabel(f'{protein_name} Expression')
    ax.set_title(f'Age × Tau Status Interaction\n(p = {interaction_pval:.4f})')
    ax.legend()

    plt.tight_layout()
    plt.show()

    return interaction_model.summary()

# Test interactions
interaction_summary = test_interaction_effects(adata, 'SQSTM1')
```

### 6.2 Diagnostic Checks

```python
def model_diagnostics(adata, protein_name):
    """
    Comprehensive model diagnostics
    """
    print(f"\n🔍 Model Diagnostics for {protein_name}")
    print("=" * 50)

    protein_expr = adata[:, protein_name].X.flatten()

    model_data = pd.DataFrame({
        'expression': protein_expr,
        'tau_binary': (adata.obs['tau_status'] == 'tau_positive').astype(int),
        'age': adata.obs['age'].values,
        'pmi': adata.obs['pmi'].values
    })

    # Fit model
    model = ols('expression ~ tau_binary + age + pmi', data=model_data).fit()

    # Get residuals and fitted values
    residuals = model.resid
    fitted = model.fittedvalues

    # Diagnostic plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f'Model Diagnostics: {protein_name}', fontsize=14)

    # 1. Residuals vs Fitted
    ax = axes[0, 0]
    ax.scatter(fitted, residuals, alpha=0.5)
    ax.axhline(y=0, color='red', linestyle='--')
    ax.set_xlabel('Fitted Values')
    ax.set_ylabel('Residuals')
    ax.set_title('Residuals vs Fitted')

    # 2. Q-Q plot
    ax = axes[0, 1]
    stats.probplot(residuals, dist="norm", plot=ax)
    ax.set_title('Q-Q Plot')

    # 3. Scale-Location
    ax = axes[1, 0]
    standardized_residuals = residuals / residuals.std()
    ax.scatter(fitted, np.sqrt(np.abs(standardized_residuals)), alpha=0.5)
    ax.set_xlabel('Fitted Values')
    ax.set_ylabel('√|Standardized Residuals|')
    ax.set_title('Scale-Location Plot')

    # 4. Residual histogram
    ax = axes[1, 1]
    ax.hist(residuals, bins=30, edgecolor='black', density=True)

    # Overlay normal distribution
    mu, std = residuals.mean(), residuals.std()
    xmin, xmax = ax.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    p = stats.norm.pdf(x, mu, std)
    ax.plot(x, p, 'r-', linewidth=2, label='Normal')

    ax.set_xlabel('Residuals')
    ax.set_ylabel('Density')
    ax.set_title('Residual Distribution')
    ax.legend()

    plt.tight_layout()
    plt.show()

    # Statistical tests
    print("\n📊 Statistical Tests:")

    # Shapiro-Wilk test for normality
    shapiro_stat, shapiro_p = stats.shapiro(residuals[:min(5000, len(residuals))])
    print(f"  Shapiro-Wilk test: p = {shapiro_p:.4f}")
    if shapiro_p < 0.05:
        print("    ⚠️ Residuals may not be normally distributed")
    else:
        print("    ✓ Residuals appear normally distributed")

    # Breusch-Pagan test for homoscedasticity
    from statsmodels.stats.diagnostic import het_breuschpagan
    bp_stat, bp_p, _, _ = het_breuschpagan(residuals, model.model.exog)
    print(f"  Breusch-Pagan test: p = {bp_p:.4f}")
    if bp_p < 0.05:
        print("    ⚠️ Heteroscedasticity detected")
    else:
        print("    ✓ Homoscedasticity assumption appears satisfied")

    # VIF for multicollinearity
    from statsmodels.stats.outliers_influence import variance_inflation_factor
    vif_data = pd.DataFrame()
    vif_data["Variable"] = model.model.exog_names[1:]  # Exclude intercept
    vif_data["VIF"] = [variance_inflation_factor(model.model.exog, i)
                       for i in range(1, model.model.exog.shape[1])]
    print(f"\n  Variance Inflation Factors:")
    for _, row in vif_data.iterrows():
        vif_warning = "⚠️" if row['VIF'] > 10 else "✓"
        print(f"    {row['Variable']}: {row['VIF']:.2f} {vif_warning}")

    return model

# Run diagnostics
model_with_diagnostics = model_diagnostics(adata, 'SQSTM1')
```

---

## 💡 Tips and Best Practices

### Do's and Don'ts

```python
# Create a reference guide
print("""
✅ DO's:
--------
1. ALWAYS check covariate balance between groups
2. ALWAYS perform diagnostic checks on your models
3. ALWAYS validate findings with sensitivity analyses
4. ALWAYS report both adjusted and unadjusted results
5. ALWAYS consider biological plausibility

❌ DON'Ts:
----------
1. DON'T adjust for mediators (part of causal pathway)
2. DON'T include highly correlated covariates (VIF > 10)
3. DON'T ignore batch effects in proteomics data
4. DON'T over-adjust (too many covariates, small sample)
5. DON'T interpret without considering effect sizes

📊 COVARIATE SELECTION GUIDE:
-----------------------------
ALWAYS ADJUST FOR:
- Batch/technical factors
- Age (in aging/disease studies)
- Sex (if mixed population)

CONSIDER ADJUSTING FOR:
- PMI (if variable)
- Brain region (if multiple)
- RIN (if RNA/protein quality varies)

RARELY ADJUST FOR:
- Disease severity markers
- Symptoms
- Downstream effects

🔍 WARNING SIGNS:
-----------------
- R² > 0.9 (overfitting)
- VIF > 10 (multicollinearity)
- Residuals not normal (model assumptions violated)
- Huge changes after adjustment (confounding or error)
""")
```

---

## 🎯 Summary and Next Steps

### What You've Learned

1. **How to identify and assess covariates** in your proteomics data
2. **Building linear models** with appropriate adjustments
3. **Comparing adjusted vs unadjusted** analyses
4. **Performing sensitivity analyses** to validate findings
5. **Diagnostic checking** to ensure model validity

### Your Analysis Checklist

```python
analysis_checklist = """
☐ 1. Load data and examine metadata
☐ 2. Check group balance for covariates
☐ 3. Detect batch effects (PCA, ANOVA)
☐ 4. Select appropriate covariates
☐ 5. Build adjusted models
☐ 6. Compare to unadjusted results
☐ 7. Perform sensitivity analyses
☐ 8. Check model diagnostics
☐ 9. Test for interactions if relevant
☐ 10. Validate in independent data
"""
print(analysis_checklist)
```

### Common Questions

**Q: How many covariates should I include?**
A: Rule of thumb: 1 covariate per 10-15 samples. With 150 samples, limit to ~10-15 covariates maximum.

**Q: What if my groups are imbalanced for a covariate?**
A: This makes adjustment MORE important, not less. The imbalance could drive false positives.

**Q: Should I use propensity score matching instead?**
A: Both approaches are valid. Linear models are simpler and more interpretable for most analyses.

**Q: What if adjustment changes everything?**
A: This suggests strong confounding. Your "discoveries" might have been driven by covariates, not biology.

---

## 📚 Further Resources

- [Linear Models in R Tutorial](https://www.datacamp.com/tutorial/linear-models-r)
- [Combat for Batch Correction](https://bioconductor.org/packages/ComBat/)
- [limma User Guide](https://bioconductor.org/packages/limma/)
- [Statistical Methods for DE](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106)

---

**Congratulations! You've completed the covariate-adjusted DE analysis tutorial. Your results are now more robust and reliable, accounting for confounding factors that could mislead your biological interpretations.**

*Next: [Statistical Methods Deep Dive](statistical_methods.md)*