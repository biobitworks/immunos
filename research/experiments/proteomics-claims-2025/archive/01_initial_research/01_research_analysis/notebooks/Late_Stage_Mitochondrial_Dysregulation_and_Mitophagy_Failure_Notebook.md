# Late-Stage Mitochondrial Dysregulation and Mitophagy Failure: Evaluation Notebook

This is a working notebook for performing Kosmos data analysis evaluation. You'll find instructions for how to use this notebook below.

If you are not familiar with Colab Notebooks, please visit the welcome notebook at: https://colab.research.google.com/notebooks/intro.ipynb. In short, they function similarly to standard Jupyter Notebooks, but are more shareable for our purpose.

## Instructions - IMPORTANT, READ CAREFULLY

This notebook has been created for your specific use, and you should work directly in it. DO NOT download the notebook and work locally and reupload.

All cells of the notebook should be run in the place and left as-is so outputs can be inspected.

All code should be written in Python as well as use terminal commands. If you need to conduct single-cell analysis, use the scanpy package.

Please be highly verbose using markdown cells to outline the high-level steps taken including rationale behind decisions made for using various tools or carrying out specific steps. You should use inline comments diligently to explain each code block's purpose.

It is also helpful if you are liberal in your use of separate code cells rather than including too much code together in a single cell.

You may import external packages at the top or load them throughout the code. Feel free to do so in whichever way is most natural to you.

## File handling

You'll find any necessary input files in the same folder you found this notebook in. They should be uploaded to the notebook environment by clicking on the folder icon on the left. These are stored in a /content/ directory. They can be programmatically accessed simply with the filename.

Output files you write are also written to the /content/ directory. You can view this directory at any time by clicking the folder icon.

As always, if you have any questions do not hesitate to reach out to Jon via Slack or email (jon@futurehouse.org)

---

# Research Analysis Overview

## Biological Context: Mitochondrial Dysfunction and Mitophagy Failure in Neurodegeneration

Mitochondrial dysfunction is a hallmark of neurodegenerative diseases, particularly in late-stage disease progression. The failure of mitophagy (selective autophagy of mitochondria) leads to the accumulation of damaged mitochondria, contributing to neuronal death.

### Key Molecular Players:
1. **SQSTM1/p62**: Autophagy receptor that can target mitochondria for mitophagy
2. **VDAC1**: Voltage-dependent anion channel, mitochondrial outer membrane protein
3. **CYCS**: Cytochrome c, critical for both respiration and apoptosis
4. **UPS Proteins**: Ubiquitin-proteasome system components
5. **Autophagy Machinery**: BECN1, CTSD, ATG12, ULK1, CTSL

### Research Hypothesis:
In late-stage neurodegeneration, mitophagy becomes increasingly dysregulated, leading to a shift from protective to pathological processes. This manifests as altered protein expression patterns and correlation changes along disease progression (pseudotime).

### Analytical Framework:
This notebook will rigorously evaluate 8 specific biological claims about mitochondrial dysregulation and mitophagy failure using:
- Differential expression analysis
- Correlation analysis (global and temporal)
- Sliding window analysis along pseudotime
- Biphasic pattern detection
- Statistical validation with appropriate multiple testing correction

# Data loading

Load any local or remote data in this section.

For each data file or set of related files, include a brief description as a comment above or at the end of the loading call.


```python
# Install scanpy via pip. Remove hashtag if needed.
# Scanpy is essential for reading and manipulating h5ad format proteomics data
# It provides the infrastructure for single-cell/proteomics analysis
!pip install scanpy
```


```python
# We need to load the package to read and load the data file.
# For python, scanpy is a package needed to conduct single cell analysis.
# So, we import the scanpy package as sc for convenience.
import scanpy as sc
```


```python
# Import essential data manipulation and analysis libraries
import pandas as pd  # Data manipulation and analysis
import numpy as np   # Numerical computing and array operations
import scipy.stats as stats  # Statistical functions and tests
from scipy.signal import find_peaks  # For peak detection in sliding window analysis
```


```python
# Import advanced statistical libraries for comprehensive analysis
import statsmodels.api as sm  # Statistical modeling
from statsmodels.stats.multitest import multipletests  # Multiple testing correction
from sklearn.linear_model import LinearRegression  # For regression analysis
from sklearn.preprocessing import StandardScaler  # Data standardization
```


```python
# Import visualization libraries for publication-quality plots
import matplotlib.pyplot as plt  # Basic plotting functionality
import seaborn as sns           # Statistical data visualization
from matplotlib.patches import Rectangle  # For custom plot elements
import matplotlib.patches as mpatches     # For legend customization

# Configure plotting parameters for publication quality
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['legend.frameon'] = True
sns.set_style("whitegrid")  # Clean, professional appearance
```


```python
# We can use the scanpy.read_h5ad() function to load the data.
# We need to alter the function to sc.read_h5ad() due to the import.
# The data file was uploaded to the /content/ directory per instructions.
adata = sc.read_h5ad('/content/pool_processed_v2.h5ad')
```


```python
# We want to inspect the data to verify.
# Printing the data provides the contents with dimensions.
print(adata)
print(f"\nDataset successfully loaded:")
print(f"- Samples (neuronal pools): {adata.n_obs}")
print(f"- Proteins quantified: {adata.n_vars}")
print(f"- Data matrix type: {type(adata.X)}")
```

# Data exploration

Include data exploration here, including reading the structure of the data files such as column names, sample ID formats, etc.

## Understanding Dataset Structure for Mitochondrial Analysis

Before analyzing mitochondrial dysfunction and mitophagy failure, we need to understand:
- Sample metadata structure and key variables
- Protein annotation format
- Expression data characteristics
- Presence of target proteins for our analysis


```python
# AnnData is a core package within scanpy.
# AnnData is a Python package that handles annotated data matrices.

# We can access the shape of the data, which provides the dimensions.
# (n_observations, n_variables)
print(f"Data shape: {adata.shape}")
print(f"Observations (samples): {adata.n_obs}")
print(f"Variables (proteins): {adata.n_vars}")
```


```python
# We can verify that the observations contain information about each sample.
print("=== SAMPLE METADATA OVERVIEW ===")
print("Available columns in .obs:")
for col in adata.obs.columns:
    dtype = adata.obs[col].dtype
    unique_vals = adata.obs[col].nunique()
    print(f"  - {col}: {dtype} ({unique_vals} unique values)")

print("\n=== SAMPLE METADATA PREVIEW ===")
print(adata.obs.head())
```


```python
# Examine key variables critical for mitochondrial analysis
print("=== KEY VARIABLES FOR MITOCHONDRIAL ANALYSIS ===")

# Tau status distribution - critical for differential expression
print("Tau Status Distribution:")
tau_counts = adata.obs['TauStatus'].value_counts()
print(tau_counts)
print(f"Tau-positive: {tau_counts.get('positive', 0)} samples")
print(f"Tau-negative: {tau_counts.get('negative', 0)} samples")

# Pseudotime distribution - critical for temporal analysis
print("\nPseudotime Statistics:")
print(adata.obs['pseudotime'].describe())

# MC1 scores - important for biphasic analysis
print("\nMC1 Score Statistics:")
print(adata.obs['MC1'].describe())
```


```python
# We can confirm that the variables contain information about each protein.
print("=== PROTEIN ANNOTATION OVERVIEW ===")
print("Available columns in .var:")
for col in adata.var.columns:
    dtype = adata.var[col].dtype
    non_null = adata.var[col].notna().sum()
    print(f"  - {col}: {dtype} ({non_null}/{len(adata.var)} non-null)")

print("\n=== PROTEIN ANNOTATION PREVIEW ===")
print(adata.var.head())
```


```python
# We can view the values of the matrix.
# We can view the expression using the pandas DataFrame.
# We can see the values are the Log2-transformed protein expression levels.
print("=== EXPRESSION DATA CHARACTERISTICS ===")

# Convert to DataFrame for analysis
expr_df = adata.to_df()

print(f"Expression data shape: {expr_df.shape}")
print(f"Data type: {expr_df.dtypes.iloc[0]}")
print(f"\nExpression statistics:")
print(f"  Range: {expr_df.values.min():.3f} to {expr_df.values.max():.3f}")
print(f"  Mean: {expr_df.values.mean():.3f}")
print(f"  Std: {expr_df.values.std():.3f}")
print(f"  Missing values: {expr_df.isnull().sum().sum()}")

print(f"\n=== SAMPLE EXPRESSION VALUES ===")
print("First 5 samples x 5 proteins:")
print(expr_df.iloc[:5, :5])
```


```python
# Identify key proteins for mitochondrial dysfunction analysis
print("=== IDENTIFYING TARGET PROTEINS ===")

# Define target proteins for our analysis based on the claims
target_proteins = {
    'SQSTM1': ['SQSTM1', 'sequestosome', 'p62'],  # Autophagy receptor
    'VDAC1': ['VDAC1', 'voltage-dependent anion'],  # Mitochondrial channel
    'CYCS': ['CYCS', 'cytochrome c'],  # Cytochrome c
    'BECN1': ['BECN1', 'beclin'],  # Autophagy protein
    'CTSD': ['CTSD', 'cathepsin D'],  # Lysosomal protease
    'ATG12': ['ATG12', 'autophagy related 12'],  # Autophagy protein
    'ULK1': ['ULK1', 'unc-51 like autophagy'],  # Autophagy initiation
    'CTSL': ['CTSL', 'cathepsin L'],  # Lysosomal protease
    'TAX1BP1': ['TAX1BP1', 'Tax1 binding protein'],  # Autophagy receptor
    'CAT': ['CAT', 'catalase'],  # Antioxidant enzyme
    'PRDX1': ['PRDX1', 'peroxiredoxin'],  # Antioxidant enzyme
    'KEAP1': ['KEAP1', 'kelch like ECH'],  # Oxidative stress regulator
    'TFRC': ['TFRC', 'transferrin receptor']  # Iron transport
}

# Search for each target protein in the dataset
found_proteins = {}
missing_proteins = []

for protein_name, search_terms in target_proteins.items():
    # Search in gene names and descriptions
    matches = []
    
    for term in search_terms:
        # Case-insensitive search in GeneName column
        gene_matches = adata.var[adata.var['GeneName'].str.contains(term, case=False, na=False)]
        # Case-insensitive search in Description column
        desc_matches = adata.var[adata.var['Description'].str.contains(term, case=False, na=False)]
        
        matches.extend(gene_matches.index.tolist())
        matches.extend(desc_matches.index.tolist())
    
    # Remove duplicates
    matches = list(set(matches))
    
    if matches:
        found_proteins[protein_name] = matches
        print(f"{protein_name}: Found {len(matches)} match(es)")
        for match in matches[:3]:  # Show first 3 matches
            gene = adata.var.loc[match, 'GeneName']
            desc = adata.var.loc[match, 'Description'][:50] + '...' if len(adata.var.loc[match, 'Description']) > 50 else adata.var.loc[match, 'Description']
            print(f"  - {match}: {gene} | {desc}")
        if len(matches) > 3:
            print(f"  ... and {len(matches)-3} more")
    else:
        missing_proteins.append(protein_name)
        print(f"{protein_name}: NOT FOUND")

print(f"\n=== SUMMARY ===")
print(f"Proteins found: {len(found_proteins)}/{len(target_proteins)}")
print(f"Missing proteins: {missing_proteins}")
```

# Claim 1: Targeted analyses show no significant UPS protein alterations across tau-positive versus tau-negative neurons.

## Biological Context

The Ubiquitin-Proteasome System (UPS) is the primary pathway for degrading misfolded and damaged proteins in cells. In neurodegenerative diseases, UPS function is often impaired, leading to protein accumulation. However, this claim suggests that UPS proteins themselves are not significantly altered between tau-positive and tau-negative neurons.

### UPS Components to Analyze:
- **Proteasome subunits**: PSMA1-7, PSMB1-7, PSMD1-14
- **E1 enzymes**: UBA1, UBA6
- **E2 enzymes**: UBE2 family proteins
- **E3 ligases**: Various RING, HECT, and RBR domain proteins
- **Deubiquitinases**: USP family proteins

## Analytical Approach

1. **Identify UPS proteins** in the dataset using systematic annotation
2. **Perform differential expression analysis** comparing tau+ vs tau- neurons
3. **Apply multiple testing correction** to control false discovery rate
4. **Assess effect sizes** to determine biological significance
5. **Validate results** with appropriate statistical tests

## Analysis code

Run your analysis code in this section. Important: include detailed markdown explaining your analytical steps and rationale for choosing certain modelling and data processing approaches and similar decisions. We also ask that you comment your code itself thoroughly.


```python
# Identify UPS (Ubiquitin-Proteasome System) proteins in the dataset
print("=== IDENTIFYING UPS PROTEINS ===")

# Define comprehensive UPS protein categories and search terms
ups_categories = {
    'proteasome_20s': ['PSMA', 'PSMB'],  # 20S proteasome subunits
    'proteasome_19s': ['PSMD', 'PSMC'],  # 19S regulatory particle
    'e1_enzymes': ['UBA1', 'UBA6'],  # Ubiquitin-activating enzymes
    'e2_enzymes': ['UBE2', 'UBC'],  # Ubiquitin-conjugating enzymes
    'e3_ligases': ['UBR', 'RING', 'HECT', 'RBR', 'TRIM', 'MDM2', 'HUWE1'],  # E3 ligases
    'deubiquitinases': ['USP', 'UCH', 'OTU', 'JAMM'],  # Deubiquitinases
    'ubiquitin': ['UBB', 'UBC', 'UBA52', 'RPS27A'],  # Ubiquitin proteins
    'proteasome_assembly': ['POMP', 'PAC', 'PSMF1', 'PSMG'],  # Assembly factors
    'proteasome_activators': ['PSME', 'PA28', 'PA200'],  # Proteasome activators
}

# Additional UPS-related terms for broader search
ups_general_terms = [
    'proteasome', 'ubiquitin', 'protease', 'degradation',
    'deubiquitinating', 'ubiquitin ligase', 'proteasome subunit'
]

# Search for UPS proteins in the dataset
ups_proteins_found = {}
all_ups_indices = set()

print("Searching for UPS proteins by category...")

for category, terms in ups_categories.items():
    category_matches = set()
    
    for term in terms:
        # Search in gene names (primary search)
        gene_matches = adata.var[adata.var['GeneName'].str.contains(term, case=False, na=False)].index
        category_matches.update(gene_matches)
        
        # Search in descriptions (secondary search)
        desc_matches = adata.var[adata.var['Description'].str.contains(term, case=False, na=False)].index
        category_matches.update(desc_matches)
    
    if category_matches:
        ups_proteins_found[category] = list(category_matches)
        all_ups_indices.update(category_matches)
        print(f"  {category}: {len(category_matches)} proteins")
    else:
        ups_proteins_found[category] = []
        print(f"  {category}: 0 proteins")

# Additional search using general terms
print("\nSearching using general UPS terms...")
general_matches = set()

for term in ups_general_terms:
    # Search in descriptions for general terms
    matches = adata.var[adata.var['Description'].str.contains(term, case=False, na=False)].index
    general_matches.update(matches)

# Remove proteins already found in specific categories
additional_ups = general_matches - all_ups_indices
all_ups_indices.update(additional_ups)

if additional_ups:
    ups_proteins_found['general_ups'] = list(additional_ups)
    print(f"  Additional UPS proteins: {len(additional_ups)}")

print(f"\n=== UPS PROTEIN IDENTIFICATION SUMMARY ===")
print(f"Total UPS proteins identified: {len(all_ups_indices)}")
print(f"Categories with proteins: {sum(1 for cat in ups_proteins_found.values() if cat)}")

# Convert to list for easier handling
ups_protein_indices = list(all_ups_indices)

# Display sample of identified UPS proteins
if ups_protein_indices:
    print(f"\nSample of identified UPS proteins:")
    sample_indices = ups_protein_indices[:10]  # Show first 10
    for idx in sample_indices:
        gene = adata.var.loc[idx, 'GeneName']
        desc = adata.var.loc[idx, 'Description'][:60] + '...' if len(adata.var.loc[idx, 'Description']) > 60 else adata.var.loc[idx, 'Description']
        print(f"  {idx}: {gene} | {desc}")
    
    if len(ups_protein_indices) > 10:
        print(f"  ... and {len(ups_protein_indices)-10} more UPS proteins")
else:
    print("WARNING: No UPS proteins identified in the dataset!")
```


```python
# Perform differential expression analysis for UPS proteins
print("=== UPS DIFFERENTIAL EXPRESSION ANALYSIS ===")

if len(ups_protein_indices) == 0:
    print("ERROR: No UPS proteins found for analysis")
else:
    # Extract UPS protein expression data
    ups_expr_df = adata.to_df().iloc[:, ups_protein_indices]
    
    # Get tau status for grouping
    tau_status = adata.obs['TauStatus']
    
    print(f"Analyzing {len(ups_protein_indices)} UPS proteins")
    print(f"Sample sizes: Tau+ = {(tau_status == 'positive').sum()}, Tau- = {(tau_status == 'negative').sum()}")
    
    # Perform statistical tests for each UPS protein
    ups_de_results = []
    
    for i, protein_idx in enumerate(ups_protein_indices):
        # Get expression values for current protein
        protein_expr = ups_expr_df.iloc[:, i]
        
        # Split by tau status
        tau_pos_expr = protein_expr[tau_status == 'positive']
        tau_neg_expr = protein_expr[tau_status == 'negative']
        
        # Perform Welch's t-test (unequal variances)
        # This is more robust than standard t-test when sample sizes differ
        t_stat, p_value = stats.ttest_ind(tau_pos_expr, tau_neg_expr, equal_var=False)
        
        # Calculate effect size (Cohen's d)
        pooled_std = np.sqrt(((len(tau_pos_expr)-1)*tau_pos_expr.var() + 
                             (len(tau_neg_expr)-1)*tau_neg_expr.var()) / 
                            (len(tau_pos_expr) + len(tau_neg_expr) - 2))
        
        if pooled_std > 0:  # Avoid division by zero
            cohens_d = (tau_pos_expr.mean() - tau_neg_expr.mean()) / pooled_std
        else:
            cohens_d = 0
        
        # Calculate log2 fold change
        log2_fc = tau_pos_expr.mean() - tau_neg_expr.mean()
        
        # Get protein information
        gene_name = adata.var.loc[protein_idx, 'GeneName']
        
        # Store results
        ups_de_results.append({
            'protein_index': protein_idx,
            'gene_name': gene_name,
            'log2_fold_change': log2_fc,
            'p_value': p_value,
            't_statistic': t_stat,
            'cohens_d': cohens_d,
            'tau_pos_mean': tau_pos_expr.mean(),
            'tau_neg_mean': tau_neg_expr.mean(),
            'tau_pos_std': tau_pos_expr.std(),
            'tau_neg_std': tau_neg_expr.std()
        })
    
    # Convert results to DataFrame for easier manipulation
    ups_results_df = pd.DataFrame(ups_de_results)
    
    print(f"\nCompleted differential expression analysis for {len(ups_results_df)} UPS proteins")
    print(f"Raw p-value range: {ups_results_df['p_value'].min():.2e} to {ups_results_df['p_value'].max():.2e}")
    print(f"Log2 fold change range: {ups_results_df['log2_fold_change'].min():.3f} to {ups_results_df['log2_fold_change'].max():.3f}")
```


```python
# Apply multiple testing correction to UPS protein analysis
print("=== MULTIPLE TESTING CORRECTION ===")

if 'ups_results_df' in locals() and len(ups_results_df) > 0:
    # Apply Benjamini-Hochberg FDR correction
    # This controls the expected proportion of false discoveries
    rejected, corrected_pvals, alpha_sidak, alpha_bonf = multipletests(
        ups_results_df['p_value'], 
        alpha=0.05, 
        method='fdr_bh'
    )
    
    # Add corrected results to DataFrame
    ups_results_df['fdr_corrected_pvalue'] = corrected_pvals
    ups_results_df['significant_fdr'] = rejected
    
    # Also apply Bonferroni correction for comparison
    bonf_rejected, bonf_corrected, _, _ = multipletests(
        ups_results_df['p_value'], 
        alpha=0.05, 
        method='bonferroni'
    )
    
    ups_results_df['bonferroni_corrected_pvalue'] = bonf_corrected
    ups_results_df['significant_bonferroni'] = bonf_rejected
    
    # Calculate additional metrics
    ups_results_df['abs_log2_fc'] = np.abs(ups_results_df['log2_fold_change'])
    ups_results_df['abs_cohens_d'] = np.abs(ups_results_df['cohens_d'])
    
    print(f"Multiple testing correction results:")
    print(f"  Total UPS proteins tested: {len(ups_results_df)}")
    print(f"  Significant (raw p < 0.05): {(ups_results_df['p_value'] < 0.05).sum()}")
    print(f"  Significant (FDR < 0.05): {rejected.sum()}")
    print(f"  Significant (Bonferroni < 0.05): {bonf_rejected.sum()}")
    
    # Calculate percentages
    total_ups = len(ups_results_df)
    raw_sig_pct = (ups_results_df['p_value'] < 0.05).sum() / total_ups * 100
    fdr_sig_pct = rejected.sum() / total_ups * 100
    bonf_sig_pct = bonf_rejected.sum() / total_ups * 100
    
    print(f"\nPercentages:")
    print(f"  Raw significance: {raw_sig_pct:.1f}%")
    print(f"  FDR significance: {fdr_sig_pct:.1f}%")
    print(f"  Bonferroni significance: {bonf_sig_pct:.1f}%")
    
    # Analyze effect sizes
    print(f"\nEffect size analysis:")
    print(f"  Mean |log2FC|: {ups_results_df['abs_log2_fc'].mean():.3f}")
    print(f"  Mean |Cohen's d|: {ups_results_df['abs_cohens_d'].mean():.3f}")
    print(f"  Large effects (|log2FC| > 1.0): {(ups_results_df['abs_log2_fc'] > 1.0).sum()}")
    print(f"  Large effects (|Cohen's d| > 0.8): {(ups_results_df['abs_cohens_d'] > 0.8).sum()}")
    
else:
    print("ERROR: No UPS differential expression results available for correction")
```


```python
# Detailed analysis of UPS protein results
print("=== DETAILED UPS PROTEIN RESULTS ===")

if 'ups_results_df' in locals() and len(ups_results_df) > 0:
    # Sort by FDR-corrected p-value
    ups_sorted = ups_results_df.sort_values('fdr_corrected_pvalue')
    
    print(f"Top 10 UPS proteins by statistical significance (FDR-corrected):")
    print("Rank | Gene | Log2FC | Raw p-val | FDR p-val | Cohen's d | Significant")
    print("-" * 80)
    
    for i, (idx, row) in enumerate(ups_sorted.head(10).iterrows()):
        rank = i + 1
        gene = row['gene_name'][:10]  # Truncate long names
        log2fc = row['log2_fold_change']
        raw_p = row['p_value']
        fdr_p = row['fdr_corrected_pvalue']
        cohens_d = row['cohens_d']
        sig_marker = "***" if row['significant_fdr'] else ""
        
        print(f"{rank:2d}   | {gene:10s} | {log2fc:6.3f} | {raw_p:8.2e} | {fdr_p:8.2e} | {cohens_d:7.3f} | {sig_marker}")
    
    # Show any significant proteins
    significant_ups = ups_results_df[ups_results_df['significant_fdr']]
    
    if len(significant_ups) > 0:
        print(f"\n=== SIGNIFICANT UPS PROTEINS (FDR < 0.05) ===")
        print(f"Found {len(significant_ups)} significant UPS proteins:")
        
        for idx, row in significant_ups.iterrows():
            gene = row['gene_name']
            log2fc = row['log2_fold_change']
            fdr_p = row['fdr_corrected_pvalue']
            cohens_d = row['cohens_d']
            direction = "UP" if log2fc > 0 else "DOWN"
            
            print(f"  {gene}: {direction} {abs(log2fc):.3f} log2FC, FDR p-val = {fdr_p:.2e}, Cohen's d = {cohens_d:.3f}")
            
            # Add protein description
            protein_desc = adata.var.loc[row['protein_index'], 'Description']
            print(f"    Description: {protein_desc[:100]}...")
    else:
        print(f"\n=== NO SIGNIFICANT UPS PROTEINS FOUND ===")
        print(f"Zero UPS proteins show significant alterations (FDR < 0.05)")
        print(f"This supports the claim that UPS proteins are not significantly altered.")
    
    # Summary statistics for claim validation
    print(f"\n=== CLAIM VALIDATION STATISTICS ===")
    print(f"Total UPS proteins analyzed: {len(ups_results_df)}")
    print(f"Significant alterations (FDR < 0.05): {(ups_results_df['significant_fdr']).sum()}")
    print(f"Percentage significant: {(ups_results_df['significant_fdr']).sum() / len(ups_results_df) * 100:.1f}%")
    
    # Effect size summary
    mean_abs_fc = ups_results_df['abs_log2_fc'].mean()
    mean_abs_d = ups_results_df['abs_cohens_d'].mean()
    
    print(f"\nEffect size summary:")
    print(f"  Mean absolute log2 fold change: {mean_abs_fc:.3f}")
    print(f"  Mean absolute Cohen's d: {mean_abs_d:.3f}")
    
    # Biological interpretation of effect sizes
    if mean_abs_d < 0.2:
        effect_interpretation = "negligible"
    elif mean_abs_d < 0.5:
        effect_interpretation = "small"
    elif mean_abs_d < 0.8:
        effect_interpretation = "medium"
    else:
        effect_interpretation = "large"
    
    print(f"  Overall effect size interpretation: {effect_interpretation}")
    
else:
    print("ERROR: No UPS results available for detailed analysis")
```

## Results code

Put results in this section. This includes generating any final plots or calculating any final values, as well as verbose markdown explaining any conclusions based on the interpretation of the results.


```python
# Create comprehensive visualization for UPS protein analysis
print("=== GENERATING UPS ANALYSIS VISUALIZATION ===")

if 'ups_results_df' in locals() and len(ups_results_df) > 0:
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Volcano plot of UPS proteins
    x = ups_results_df['log2_fold_change']
    y = -np.log10(ups_results_df['fdr_corrected_pvalue'])
    
    # Color points by significance
    colors = ['red' if sig else 'gray' for sig in ups_results_df['significant_fdr']]
    
    axes[0,0].scatter(x, y, c=colors, alpha=0.7, s=50)
    axes[0,0].axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.7, label='FDR = 0.05')
    axes[0,0].axvline(x=0, color='black', linestyle='-', alpha=0.3)
    axes[0,0].set_xlabel('Log2 Fold Change (Tau+ vs Tau-)')
    axes[0,0].set_ylabel('-Log10(FDR-corrected p-value)')
    axes[0,0].set_title(f'UPS Proteins Volcano Plot\n({len(ups_results_df)} proteins analyzed)')
    axes[0,0].legend()
    axes[0,0].grid(True, alpha=0.3)
    
    # Add text annotation
    sig_count = (ups_results_df['significant_fdr']).sum()
    axes[0,0].text(0.02, 0.98, f'Significant: {sig_count}/{len(ups_results_df)}', 
                  transform=axes[0,0].transAxes, verticalalignment='top',
                  bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Plot 2: Distribution of p-values
    axes[0,1].hist(ups_results_df['p_value'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    axes[0,1].axvline(x=0.05, color='red', linestyle='--', label='p = 0.05')
    axes[0,1].set_xlabel('Raw p-value')
    axes[0,1].set_ylabel('Frequency')
    axes[0,1].set_title('Distribution of Raw p-values')
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)
    
    # Plot 3: Distribution of fold changes
    axes[1,0].hist(ups_results_df['log2_fold_change'], bins=20, alpha=0.7, color='lightgreen', edgecolor='black')
    axes[1,0].axvline(x=0, color='red', linestyle='--', label='No change')
    axes[1,0].set_xlabel('Log2 Fold Change')
    axes[1,0].set_ylabel('Frequency')
    axes[1,0].set_title('Distribution of Log2 Fold Changes')
    axes[1,0].legend()
    axes[1,0].grid(True, alpha=0.3)
    
    # Plot 4: Effect sizes (Cohen's d)
    axes[1,1].hist(ups_results_df['cohens_d'], bins=20, alpha=0.7, color='orange', edgecolor='black')
    axes[1,1].axvline(x=0, color='red', linestyle='--', label='No effect')
    # Add effect size interpretation lines
    axes[1,1].axvline(x=0.2, color='blue', linestyle=':', alpha=0.5, label='Small effect')
    axes[1,1].axvline(x=-0.2, color='blue', linestyle=':', alpha=0.5)
    axes[1,1].axvline(x=0.8, color='purple', linestyle=':', alpha=0.5, label='Large effect')
    axes[1,1].axvline(x=-0.8, color='purple', linestyle=':', alpha=0.5)
    axes[1,1].set_xlabel("Cohen's d")
    axes[1,1].set_ylabel('Frequency')
    axes[1,1].set_title('Distribution of Effect Sizes')
    axes[1,1].legend()
    axes[1,1].grid(True, alpha=0.3)
    
    plt.suptitle('UPS Protein Differential Expression Analysis', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.show()
    
    print("UPS protein analysis visualization completed")
    
else:
    print("Cannot generate visualization - no UPS results available")
```


```python
# Final assessment of Claim 1
print("=== CLAIM 1 COMPREHENSIVE ASSESSMENT ===")
print("\nClaim: 'Targeted analyses show no significant UPS protein alterations")
print("across tau-positive versus tau-negative neurons.'")

if 'ups_results_df' in locals() and len(ups_results_df) > 0:
    # Key statistics for claim validation
    total_ups_proteins = len(ups_results_df)
    significant_proteins = (ups_results_df['significant_fdr']).sum()
    percentage_significant = significant_proteins / total_ups_proteins * 100
    
    print(f"\n=== QUANTITATIVE RESULTS ===")
    print(f"Total UPS proteins analyzed: {total_ups_proteins}")
    print(f"Significant alterations (FDR < 0.05): {significant_proteins}")
    print(f"Percentage significant: {percentage_significant:.1f}%")
    print(f"Percentage non-significant: {100 - percentage_significant:.1f}%")
    
    # Statistical power assessment
    mean_abs_effect = ups_results_df['abs_cohens_d'].mean()
    median_pvalue = ups_results_df['p_value'].median()
    
    print(f"\n=== EFFECT SIZE ASSESSMENT ===")
    print(f"Mean absolute effect size (Cohen's d): {mean_abs_effect:.3f}")
    print(f"Median raw p-value: {median_pvalue:.3f}")
    
    # Categorize effect sizes
    negligible_effects = (ups_results_df['abs_cohens_d'] < 0.2).sum()
    small_effects = ((ups_results_df['abs_cohens_d'] >= 0.2) & (ups_results_df['abs_cohens_d'] < 0.5)).sum()
    medium_effects = ((ups_results_df['abs_cohens_d'] >= 0.5) & (ups_results_df['abs_cohens_d'] < 0.8)).sum()
    large_effects = (ups_results_df['abs_cohens_d'] >= 0.8).sum()
    
    print(f"\nEffect size distribution:")
    print(f"  Negligible (|d| < 0.2): {negligible_effects} ({negligible_effects/total_ups_proteins*100:.1f}%)")
    print(f"  Small (0.2 ≤ |d| < 0.5): {small_effects} ({small_effects/total_ups_proteins*100:.1f}%)")
    print(f"  Medium (0.5 ≤ |d| < 0.8): {medium_effects} ({medium_effects/total_ups_proteins*100:.1f}%)")
    print(f"  Large (|d| ≥ 0.8): {large_effects} ({large_effects/total_ups_proteins*100:.1f}%)")
    
    # Direction of changes analysis
    upregulated = (ups_results_df['log2_fold_change'] > 0).sum()
    downregulated = (ups_results_df['log2_fold_change'] < 0).sum()
    unchanged = (ups_results_df['log2_fold_change'] == 0).sum()
    
    print(f"\n=== DIRECTION OF CHANGES ===")
    print(f"Upregulated (log2FC > 0): {upregulated} ({upregulated/total_ups_proteins*100:.1f}%)")
    print(f"Downregulated (log2FC < 0): {downregulated} ({downregulated/total_ups_proteins*100:.1f}%)")
    print(f"Unchanged (log2FC = 0): {unchanged} ({unchanged/total_ups_proteins*100:.1f}%)")
    
    # Biological interpretation
    print(f"\n=== BIOLOGICAL INTERPRETATION ===")
    
    if significant_proteins == 0:
        biological_conclusion = "STRONG SUPPORT"
        interpretation = "No UPS proteins show significant alterations, strongly supporting the claim."
    elif percentage_significant < 5:
        biological_conclusion = "SUPPORT"
        interpretation = f"Very few UPS proteins ({percentage_significant:.1f}%) show significant alterations, generally supporting the claim."
    elif percentage_significant < 10:
        biological_conclusion = "WEAK SUPPORT"
        interpretation = f"A small proportion ({percentage_significant:.1f}%) of UPS proteins show alterations, providing weak support for the claim."
    else:
        biological_conclusion = "NO SUPPORT"
        interpretation = f"A substantial proportion ({percentage_significant:.1f}%) of UPS proteins show significant alterations, contradicting the claim."
    
    print(f"Biological conclusion: {biological_conclusion}")
    print(f"Interpretation: {interpretation}")
    
    # Technical considerations
    print(f"\n=== TECHNICAL CONSIDERATIONS ===")
    print(f"Statistical method: Welch's t-test with FDR correction")
    print(f"Sample size: Tau+ = {(adata.obs['TauStatus'] == 'positive').sum()}, Tau- = {(adata.obs['TauStatus'] == 'negative').sum()}")
    print(f"Multiple testing: Benjamini-Hochberg FDR control")
    print(f"Effect size metric: Cohen's d")
    
    # List any significant proteins for transparency
    if significant_proteins > 0:
        print(f"\n=== SIGNIFICANT UPS PROTEINS IDENTIFIED ===")
        significant_ups = ups_results_df[ups_results_df['significant_fdr']]
        for idx, row in significant_ups.iterrows():
            gene = row['gene_name']
            log2fc = row['log2_fold_change']
            fdr_p = row['fdr_corrected_pvalue']
            direction = "upregulated" if log2fc > 0 else "downregulated"
            print(f"  {gene}: {direction} (log2FC = {log2fc:.3f}, FDR p-val = {fdr_p:.2e})")
    
    # Final validation conclusion
    print(f"\n=== FINAL CLAIM VALIDATION ===")
    print(f"Status: {biological_conclusion}")
    print(f"\nEvidence summary:")
    print(f"• {total_ups_proteins} UPS proteins analyzed using rigorous statistical methods")
    print(f"• {significant_proteins} proteins ({percentage_significant:.1f}%) show significant alterations")
    print(f"• Mean effect size is {mean_abs_effect:.3f} (Cohen's d), indicating {effect_interpretation} effects")
    print(f"• Results are consistent with maintained UPS function in tau pathology")
    
    if percentage_significant < 5:
        print(f"\nConclusion: The claim is SUPPORTED by the data. UPS proteins show")
        print(f"minimal significant alterations between tau-positive and tau-negative neurons.")
    else:
        print(f"\nConclusion: The claim requires QUALIFICATION. While most UPS proteins")
        print(f"are unchanged, {significant_proteins} proteins do show significant alterations.")
        
else:
    print("\n=== ANALYSIS LIMITATION ===")
    print("Could not complete UPS protein analysis due to protein identification issues.")
    print("This may indicate that UPS proteins are either:")
    print("1. Not present in the dataset (technical limitation)")
    print("2. Named differently than expected (annotation issue)")
    print("3. Below detection threshold (biological reality)")
    print("\nRecommendation: Verify protein annotation methods and detection sensitivity.")
```
