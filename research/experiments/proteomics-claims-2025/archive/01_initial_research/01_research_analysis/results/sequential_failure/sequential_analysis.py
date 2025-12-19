#!/usr/bin/env python
# coding: utf-8

# # Sequential Failure of Proteostasis Mechanisms: Evaluation Notebook
# 
# This is a working notebook for performing Kosmos data analysis evaluation. You'll find instructions for how to use this notebook below.
# 
# If you are not familiar with Colab Notebooks, please visit the welcome notebook at: https://colab.research.google.com/notebooks/intro.ipynb. In short, they function similarly to standard Jupyter Notebooks, but are more shareable for our purpose.

# ## Instructions - IMPORTANT, READ CAREFULLY
# 
# This notebook has been created for your specific use, and you should work directly in it. DO NOT download the notebook and work locally and reupload.
# 
# All cells of the notebook should be run in the place and left as-is so outputs can be inspected.
# 
# All code should be written in Python as well as use terminal commands. If you need to conduct single-cell analysis, use the scanpy package.
# 
# Please be highly verbose using markdown cells to outline the high-level steps taken including rationale behind decisions made for using various tools or carrying out specific steps. You should use inline comments diligently to explain each code block's purpose.
# 
# It is also helpful if you are liberal in your use of separate code cells rather than including too much code together in a single cell.
# 
# You may import external packages at the top or load them throughout the code. Feel free to do so in whichever way is most natural to you.

# ## File handling
# 
# You'll find any necessary input files in the same folder you found this notebook in. They should be uploaded to the notebook environment by clicking on the folder icon on the left. These are stored in a /content/ directory. They can be programmatically accessed simply with the filename.
# 
# Output files you write are also written to the /content/ directory. You can view this directory at any time by clicking the folder icon.

# In[ ]:


# Install scanpy for single-cell/proteomics analysis
# Scanpy is the standard package for analyzing high-dimensional biological data
# We need this to read and manipulate the h5ad format proteomics data
get_ipython().system('pip install scanpy')


# As always, if you have any questions do not hesitate to reach out to Jon via Slack or email (jon@futurehouse.org)
# 
# ---
# 
# # Research Analysis Overview
# 
# ## Biological Context: Sequential Proteostasis Failure in Neurodegeneration
# 
# Protein homeostasis (proteostasis) is critical for neuronal survival and function. In neurodegenerative diseases like Alzheimer's, multiple protein quality control systems fail sequentially, leading to toxic protein accumulation and neuronal death.
# 
# ### Key Proteostasis Systems:
# 1. **Ubiquitin-Proteasome System (UPS)**: Primary degradation pathway for misfolded proteins
# 2. **Autophagy-Lysosome Pathway**: Bulk degradation system, including mitophagy
# 3. **V-ATPase Complex**: Critical for lysosomal acidification and function
# 4. **Molecular Chaperones**: Protein folding assistance
# 
# ### Research Hypothesis:
# We hypothesize that proteostasis systems fail in a sequential, ordered manner during neurodegeneration, with specific timing and molecular signatures that can be detected through proteomics analysis.
# 
# ### Analytical Approach:
# This notebook will rigorously evaluate 8 specific biological claims about sequential proteostasis failure using statistical analysis of proteomics data from tau-positive and tau-negative neurons.

# # Data loading
# 
# Load any local or remote data in this section.
# 
# For each data file or set of related files, include a brief description as a comment above or at the end of the loading call.

# In[ ]:


# Import core scientific computing libraries
import scanpy as sc  # Single-cell analysis package (works for proteomics data)
import pandas as pd  # Data manipulation and analysis
import numpy as np   # Numerical computing


# In[ ]:


# Import statistical analysis libraries
import scipy.stats as stats  # Statistical tests and distributions
import statsmodels.api as sm  # Advanced statistical modeling
from statsmodels.stats.multitest import multipletests  # Multiple testing correction
from sklearn.linear_model import LinearRegression  # Linear regression analysis
from sklearn.preprocessing import StandardScaler   # Data standardization


# In[ ]:


# Import visualization libraries
import matplotlib.pyplot as plt  # Basic plotting
import seaborn as sns           # Statistical data visualization
from matplotlib.patches import Rectangle  # For creating custom plot elements

# Set publication-quality plotting parameters
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 1.5
sns.set_style("whitegrid")


# In[ ]:


# Load the main proteomics dataset
# This is an h5ad file containing log2-transformed protein expression data
# from tau-positive and tau-negative neurons in Alzheimer's disease
adata = sc.read_h5ad('/content/pool_processed_v2.h5ad')

# Display basic information about the loaded dataset
print(f"Dataset loaded successfully:")
print(f"- Samples (neurons): {adata.n_obs}")
print(f"- Proteins quantified: {adata.n_vars}")
print(f"- Data format: {type(adata.X)}")


# # Data exploration
# 
# Include data exploration here, including reading the structure of the data files such as column names, sample ID formats, etc.
# 
# ## Understanding the Dataset Structure
# 
# Before proceeding with analysis, we need to thoroughly understand our data structure, including:
# - Sample metadata and covariates
# - Protein annotation information
# - Data distribution and quality metrics
# - Missing values and outliers

# In[ ]:


# Examine the overall structure of our AnnData object
# AnnData is the standard format for storing annotated data matrices
print("=== DATASET OVERVIEW ===")
print(adata)
print(f"\nData matrix shape: {adata.shape}")
print(f"Data type: {adata.X.dtype}")
print(f"Is sparse: {hasattr(adata.X, 'todense')}")


# In[ ]:


# Examine sample metadata (observations)
# This contains information about each neuronal sample
print("=== SAMPLE METADATA COLUMNS ===")
for col in adata.obs.columns:
    print(f"- {col}: {adata.obs[col].dtype}")
    
print("\n=== SAMPLE METADATA PREVIEW ===")
print(adata.obs.head())


# In[ ]:


# Examine protein metadata (variables)
# This contains annotation information for each protein
print("=== PROTEIN METADATA COLUMNS ===")
for col in adata.var.columns:
    print(f"- {col}: {adata.var[col].dtype}")
    
print("\n=== PROTEIN METADATA PREVIEW ===")
print(adata.var.head())


# In[ ]:


# Examine key variables for our analysis
print("=== TAU STATUS DISTRIBUTION ===")
print(adata.obs['TauStatus'].value_counts())

print("\n=== PATIENT ID DISTRIBUTION ===")
print(adata.obs['PatientID'].value_counts())

print("\n=== COVARIATE SUMMARY STATISTICS ===")
print("Age at death:")
print(adata.obs['Age at death'].describe())

print("\nPMI hours:")
print(adata.obs['PMI hours'].describe())

print("\nMC1 scores:")
print(adata.obs['MC1'].describe())

print("\nPseudotime:")
print(adata.obs['pseudotime'].describe())


# In[ ]:


# Examine the expression data distribution
# Convert to DataFrame for easier manipulation
expr_df = adata.to_df()

print("=== EXPRESSION DATA OVERVIEW ===")
print(f"Expression range: {expr_df.values.min():.3f} to {expr_df.values.max():.3f}")
print(f"Mean expression: {expr_df.values.mean():.3f}")
print(f"Standard deviation: {expr_df.values.std():.3f}")
print(f"Missing values: {expr_df.isnull().sum().sum()}")

# Display sample of expression values
print("\n=== EXPRESSION DATA SAMPLE ===")
print(expr_df.iloc[:5, :5])  # First 5 samples and 5 proteins


# In[ ]:


# Perform data quality checks
print("=== DATA QUALITY ASSESSMENT ===")

# Check for missing values in metadata
print("Missing values in sample metadata:")
missing_obs = adata.obs.isnull().sum()
for col, count in missing_obs.items():
    if count > 0:
        print(f"  {col}: {count} missing")
    else:
        print(f"  {col}: Complete")

# Check for missing values in protein metadata
print("\nMissing values in protein metadata:")
missing_var = adata.var.isnull().sum()
for col, count in missing_var.items():
    if count > 0:
        print(f"  {col}: {count} missing")
    else:
        print(f"  {col}: Complete")

# Check for extreme outliers in expression data
print("\nExpression data outlier analysis:")
Q1 = np.percentile(expr_df.values, 25)
Q3 = np.percentile(expr_df.values, 75)
IQR = Q3 - Q1
lower_bound = Q1 - 1.5 * IQR
upper_bound = Q3 + 1.5 * IQR

outliers = np.sum((expr_df.values < lower_bound) | (expr_df.values > upper_bound))
total_values = expr_df.size
print(f"  Outliers detected: {outliers} / {total_values} ({outliers/total_values*100:.2f}%)")


# # Claim 1: A covariate-controlled differential expression analysis (age, PMI, and PatientID) across 5,853 proteins (BH-FDR) identified 2,115 proteins (36.14%) significantly altered between tau-positive and tau-negative neurons.
# 
# ## Analytical Approach
# 
# This claim requires a rigorous differential expression analysis controlling for potential confounding variables:
# 
# 1. **Covariate Control**: We must include age, post-mortem interval (PMI), and patient ID as covariates
# 2. **Statistical Model**: Linear regression with tau status as the primary predictor
# 3. **Multiple Testing**: Benjamini-Hochberg False Discovery Rate (BH-FDR) correction
# 4. **Effect Size**: Calculate fold changes and confidence intervals
# 5. **Validation**: Ensure results are biologically meaningful and statistically robust
# 
# ### Rationale for Statistical Approach:
# - **Linear regression**: Appropriate for continuous expression data with covariates
# - **BH-FDR correction**: Controls false discovery rate while maintaining power
# - **Covariate inclusion**: Critical to avoid confounding by patient characteristics

# ## Analysis code
# 
# Run your analysis code in this section. Important: include detailed markdown explaining your analytical steps and rationale for choosing certain modelling and data processing approaches and similar decisions. We also ask that you comment your code itself thoroughly.

# In[ ]:


# Prepare data for differential expression analysis
# We need to set up our predictor variables and ensure proper data types

# Create a copy of the AnnData object to avoid modifying the original
adata_de = adata.copy()

# Convert tau status to binary encoding for regression
# tau-positive = 1, tau-negative = 0
adata_de.obs['tau_binary'] = (adata_de.obs['TauStatus'] == 'positive').astype(int)

# Create patient ID dummy variables for covariate control
# This accounts for patient-specific effects
patient_dummies = pd.get_dummies(adata_de.obs['PatientID'], prefix='Patient')

# Prepare the design matrix for regression
design_df = pd.DataFrame({
    'tau_positive': adata_de.obs['tau_binary'],
    'age': adata_de.obs['Age at death'],
    'pmi': adata_de.obs['PMI hours']
})

# Add patient dummy variables (excluding one to avoid multicollinearity)
# This is the standard approach for categorical variables in regression
design_df = pd.concat([design_df, patient_dummies.iloc[:, :-1]], axis=1)

print(f"Design matrix shape: {design_df.shape}")
print(f"Included covariates: {list(design_df.columns)}")
print(f"\nSample size by tau status:")
print(adata_de.obs['TauStatus'].value_counts())


# In[ ]:


# Define function for covariate-controlled differential expression
def perform_covariate_de_analysis(expression_matrix, design_matrix, target_column='tau_positive'):
    """
    Perform covariate-controlled differential expression analysis
    
    Parameters:
    - expression_matrix: DataFrame with samples as rows, proteins as columns
    - design_matrix: DataFrame with covariates and target variable
    - target_column: Name of the column representing the condition of interest
    
    Returns:
    - results_df: DataFrame with statistical results for each protein
    """
    
    # Initialize results storage
    results = []
    
    # Add constant term for regression intercept
    X = sm.add_constant(design_matrix)
    
    print(f"Running regression analysis for {expression_matrix.shape[1]} proteins...")
    print(f"Including covariates: {list(design_matrix.columns)}")
    
    # Iterate through each protein
    for i, protein in enumerate(expression_matrix.columns):
        if i % 1000 == 0:  # Progress indicator
            print(f"  Processed {i} proteins...")
        
        try:
            # Get expression values for current protein
            y = expression_matrix[protein]
            
            # Fit linear regression model
            model = sm.OLS(y, X).fit()
            
            # Extract results for tau status coefficient
            tau_coef = model.params[target_column]
            tau_pvalue = model.pvalues[target_column]
            tau_stderr = model.bse[target_column]
            
            # Calculate confidence interval
            ci_lower, ci_upper = model.conf_int().loc[target_column]
            
            # Store results
            results.append({
                'protein': protein,
                'log2_fold_change': tau_coef,
                'p_value': tau_pvalue,
                'std_error': tau_stderr,
                'ci_lower': ci_lower,
                'ci_upper': ci_upper,
                'r_squared': model.rsquared
            })
            
        except Exception as e:
            # Handle any errors in regression fitting
            print(f"  Warning: Error fitting model for protein {protein}: {e}")
            results.append({
                'protein': protein,
                'log2_fold_change': np.nan,
                'p_value': np.nan,
                'std_error': np.nan,
                'ci_lower': np.nan,
                'ci_upper': np.nan,
                'r_squared': np.nan
            })
    
    print(f"Completed regression analysis for all {len(results)} proteins.")
    return pd.DataFrame(results)


# In[ ]:


# Run the covariate-controlled differential expression analysis
print("=== COVARIATE-CONTROLLED DIFFERENTIAL EXPRESSION ANALYSIS ===")

# Get expression data as DataFrame
expr_df = adata_de.to_df()

# Ensure design matrix aligns with expression data
assert expr_df.shape[0] == design_df.shape[0], "Sample sizes must match"

# Run the analysis
de_results = perform_covariate_de_analysis(expr_df, design_df)

# Remove any proteins with failed fits
de_results_clean = de_results.dropna(subset=['p_value'])

print(f"\nSuccessful analyses: {len(de_results_clean)} / {len(de_results)} proteins")
print(f"Failed analyses: {len(de_results) - len(de_results_clean)} proteins")


# In[ ]:


# Apply multiple testing correction (Benjamini-Hochberg FDR)
print("=== MULTIPLE TESTING CORRECTION ===")

# Apply Benjamini-Hochberg FDR correction
# This controls the expected proportion of false discoveries
rejected, corrected_pvals, alpha_sidak, alpha_bonf = multipletests(
    de_results_clean['p_value'], 
    alpha=0.05, 
    method='fdr_bh'
)

# Add corrected p-values and significance flags to results
de_results_clean['fdr_corrected_pvalue'] = corrected_pvals
de_results_clean['significant_fdr'] = rejected

# Calculate additional statistics
de_results_clean['abs_log2_fc'] = np.abs(de_results_clean['log2_fold_change'])
de_results_clean['neg_log10_pval'] = -np.log10(de_results_clean['fdr_corrected_pvalue'])

print(f"Total proteins analyzed: {len(de_results_clean)}")
print(f"Proteins significant (FDR < 0.05): {rejected.sum()}")
print(f"Percentage significant: {rejected.sum() / len(de_results_clean) * 100:.2f}%")

# Verify we have the expected number of proteins from the claim
expected_total = 5853
expected_significant = 2115
expected_percentage = 36.14

print(f"\n=== COMPARISON WITH CLAIM ===")
print(f"Expected total proteins: {expected_total}")
print(f"Observed total proteins: {len(de_results_clean)}")
print(f"Expected significant: {expected_significant} ({expected_percentage}%)")
print(f"Observed significant: {rejected.sum()} ({rejected.sum() / len(de_results_clean) * 100:.2f}%)")


# In[ ]:


# Calculate effect sizes and summary statistics for significant proteins
print("=== EFFECT SIZE ANALYSIS ===")

# Focus on significant proteins
significant_proteins = de_results_clean[de_results_clean['significant_fdr']]

print(f"Analysis of {len(significant_proteins)} significant proteins:")
print(f"\nLog2 fold change distribution:")
print(significant_proteins['log2_fold_change'].describe())

print(f"\nAbsolute log2 fold change distribution:")
print(significant_proteins['abs_log2_fc'].describe())

# Count upregulated vs downregulated
upregulated = significant_proteins[significant_proteins['log2_fold_change'] > 0]
downregulated = significant_proteins[significant_proteins['log2_fold_change'] < 0]

print(f"\nDirection of changes:")
print(f"Upregulated (log2FC > 0): {len(upregulated)} proteins ({len(upregulated)/len(significant_proteins)*100:.1f}%)")
print(f"Downregulated (log2FC < 0): {len(downregulated)} proteins ({len(downregulated)/len(significant_proteins)*100:.1f}%)")

# Identify proteins with large effect sizes
large_effects = significant_proteins[significant_proteins['abs_log2_fc'] > 1.0]
print(f"\nProteins with large effects (|log2FC| > 1.0): {len(large_effects)}")

if len(large_effects) > 0:
    print("\nTop 10 proteins by absolute fold change:")
    top_proteins = significant_proteins.nlargest(10, 'abs_log2_fc')
    for idx, row in top_proteins.iterrows():
        print(f"  {row['protein']}: log2FC = {row['log2_fold_change']:.3f}, FDR p-val = {row['fdr_corrected_pvalue']:.2e}")


# ## Results code
# 
# Put results in this section. This includes generating any final plots or calculating any final values, as well as verbose markdown explaining any conclusions based on the interpretation of the results.

# In[ ]:


# Create volcano plot to visualize differential expression results
print("=== GENERATING VOLCANO PLOT ===")

plt.figure(figsize=(12, 8))

# Plot all proteins
plt.scatter(de_results_clean['log2_fold_change'], 
           de_results_clean['neg_log10_pval'],
           c='lightgray', alpha=0.6, s=20, label='Non-significant')

# Highlight significant proteins
significant_mask = de_results_clean['significant_fdr']
plt.scatter(de_results_clean.loc[significant_mask, 'log2_fold_change'], 
           de_results_clean.loc[significant_mask, 'neg_log10_pval'],
           c='red', alpha=0.7, s=25, label='Significant (FDR < 0.05)')

# Add significance threshold line
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.7, 
           label='FDR = 0.05 threshold')

# Add fold change threshold lines
plt.axvline(x=1, color='blue', linestyle=':', alpha=0.5, label='|log2FC| = 1')
plt.axvline(x=-1, color='blue', linestyle=':', alpha=0.5)

plt.xlabel('Log2 Fold Change (Tau+ vs Tau-)', fontsize=14)
plt.ylabel('-Log10(FDR-corrected p-value)', fontsize=14)
plt.title('Volcano Plot: Covariate-Controlled Differential Expression\nTau-Positive vs Tau-Negative Neurons', 
         fontsize=16, pad=20)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True, alpha=0.3)

# Add text box with summary statistics
summary_text = f"Total proteins: {len(de_results_clean)}\n" + \
               f"Significant: {significant_mask.sum()} ({significant_mask.sum()/len(de_results_clean)*100:.1f}%)\n" + \
               f"Upregulated: {len(upregulated)}\n" + \
               f"Downregulated: {len(downregulated)}"

plt.text(0.02, 0.98, summary_text, transform=plt.gca().transAxes, 
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.show()

print(f"Volcano plot generated showing {len(de_results_clean)} proteins")
print(f"Red points represent {significant_mask.sum()} significantly altered proteins (FDR < 0.05)")


# In[ ]:


# Create distribution plots for fold changes and p-values
print("=== GENERATING DISTRIBUTION PLOTS ===")

fig, axes = plt.subplots(2, 2, figsize=(15, 10))

# Plot 1: Distribution of log2 fold changes
axes[0,0].hist(de_results_clean['log2_fold_change'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
axes[0,0].axvline(x=0, color='red', linestyle='--', label='No change')
axes[0,0].set_xlabel('Log2 Fold Change')
axes[0,0].set_ylabel('Frequency')
axes[0,0].set_title('Distribution of Log2 Fold Changes')
axes[0,0].legend()
axes[0,0].grid(True, alpha=0.3)

# Plot 2: Distribution of raw p-values
axes[0,1].hist(de_results_clean['p_value'], bins=50, alpha=0.7, color='lightgreen', edgecolor='black')
axes[0,1].axvline(x=0.05, color='red', linestyle='--', label='p = 0.05')
axes[0,1].set_xlabel('Raw p-value')
axes[0,1].set_ylabel('Frequency')
axes[0,1].set_title('Distribution of Raw p-values')
axes[0,1].legend()
axes[0,1].grid(True, alpha=0.3)

# Plot 3: Distribution of FDR-corrected p-values
axes[1,0].hist(de_results_clean['fdr_corrected_pvalue'], bins=50, alpha=0.7, color='orange', edgecolor='black')
axes[1,0].axvline(x=0.05, color='red', linestyle='--', label='FDR = 0.05')
axes[1,0].set_xlabel('FDR-corrected p-value')
axes[1,0].set_ylabel('Frequency')
axes[1,0].set_title('Distribution of FDR-corrected p-values')
axes[1,0].legend()
axes[1,0].grid(True, alpha=0.3)

# Plot 4: Fold change by significance
sig_fc = de_results_clean[de_results_clean['significant_fdr']]['log2_fold_change']
nonsig_fc = de_results_clean[~de_results_clean['significant_fdr']]['log2_fold_change']

axes[1,1].hist([nonsig_fc, sig_fc], bins=30, alpha=0.7, 
              label=['Non-significant', 'Significant (FDR < 0.05)'],
              color=['lightgray', 'red'])
axes[1,1].set_xlabel('Log2 Fold Change')
axes[1,1].set_ylabel('Frequency')
axes[1,1].set_title('Fold Changes by Significance Status')
axes[1,1].legend()
axes[1,1].grid(True, alpha=0.3)

plt.suptitle('Differential Expression Analysis: Statistical Distributions', fontsize=16, y=1.02)
plt.tight_layout()
plt.show()

print("Distribution plots generated showing statistical properties of the analysis")


# In[ ]:


# Final validation and claim assessment
print("=== CLAIM 1 VALIDATION SUMMARY ===")
print("\nClaim: 'A covariate-controlled differential expression analysis (age, PMI, and PatientID)")
print("across 5,853 proteins (BH-FDR) identified 2,115 proteins (36.14%) significantly")
print("altered between tau-positive and tau-negative neurons.'")

# Calculate exact percentages and compare to claim
total_proteins_analyzed = len(de_results_clean)
significant_proteins_found = significant_mask.sum()
percentage_significant = significant_proteins_found / total_proteins_analyzed * 100

print(f"\n=== RESULTS COMPARISON ===")
print(f"Expected total proteins: 5,853")
print(f"Actual total proteins: {total_proteins_analyzed}")
print(f"Expected significant proteins: 2,115 (36.14%)")
print(f"Actual significant proteins: {significant_proteins_found} ({percentage_significant:.2f}%)")

# Statistical assessment
print(f"\n=== STATISTICAL ASSESSMENT ===")
print(f"Analysis method: Linear regression with covariates (age, PMI, PatientID)")
print(f"Multiple testing correction: Benjamini-Hochberg FDR")
print(f"Significance threshold: FDR < 0.05")
print(f"Covariates included: {list(design_df.columns)}")

# Biological interpretation
print(f"\n=== BIOLOGICAL INTERPRETATION ===")
print(f"• {percentage_significant:.1f}% of the proteome shows significant alteration")
print(f"• {len(upregulated)} proteins upregulated in tau-positive neurons")
print(f"• {len(downregulated)} proteins downregulated in tau-positive neurons")
print(f"• Large effect sizes (|log2FC| > 1): {len(large_effects)} proteins")

# Claim validation conclusion
percentage_diff = abs(percentage_significant - 36.14)
if percentage_diff < 5.0:  # Within 5% of claimed value
    validation_status = "SUPPORTED"
elif percentage_diff < 10.0:  # Within 10% of claimed value
    validation_status = "PARTIALLY SUPPORTED"
else:
    validation_status = "NOT SUPPORTED"

print(f"\n=== CLAIM VALIDATION ===")
print(f"Status: {validation_status}")
print(f"Reason: Observed {percentage_significant:.2f}% vs claimed 36.14% (difference: {percentage_diff:.2f}%)")
print(f"\nNote: Results depend on exact dataset and analysis parameters used.")
print(f"The analysis methodology is appropriate and follows best practices for")
print(f"covariate-controlled differential expression analysis.")


# # Claim 2: The autophagy receptor SQSTM1 was the most upregulated protein (+3.41 log2) in the covariate-controlled differential expression analysis.
# 
# ## Biological Context
# 
# SQSTM1 (Sequestosome-1, also known as p62) is a critical autophagy receptor protein that:
# - Recognizes ubiquitinated proteins for autophagic degradation
# - Links the ubiquitin-proteasome system with autophagy
# - Accumulates when autophagy is impaired
# - Is commonly upregulated in neurodegenerative diseases
# 
# ## Analytical Approach
# 
# To validate this claim, we need to:
# 1. Identify SQSTM1 in our differential expression results
# 2. Confirm it has the highest log2 fold change among all proteins
# 3. Verify the specific fold change value (+3.41 log2)
# 4. Assess the statistical significance of this change
# 5. Compare with other highly upregulated proteins

# ## Analysis code
# 
# Run your analysis code in this section. Important: include detailed markdown explaining your analytical steps and rationale for choosing certain modelling and data processing approaches and similar decisions. We also ask that you comment your code itself thoroughly.

# In[ ]:


# Identify SQSTM1 in the protein annotations and differential expression results
print("=== IDENTIFYING SQSTM1 IN THE DATASET ===")

# Search for SQSTM1 in the protein metadata
# We need to check both gene names and protein descriptions
sqstm1_mask = (adata.var['GeneName'].str.contains('SQSTM1', na=False, case=False)) | \
              (adata.var['Description'].str.contains('SQSTM1|sequestosome|p62', na=False, case=False))

sqstm1_proteins = adata.var[sqstm1_mask]

print(f"Found {len(sqstm1_proteins)} protein(s) matching SQSTM1:")
for idx, row in sqstm1_proteins.iterrows():
    print(f"  Index: {idx}")
    print(f"  Gene: {row['GeneName']}")
    print(f"  UniProt: {row['UniprotID']}")
    print(f"  Description: {row['Description']}")
    print()

# If multiple matches, identify the canonical SQSTM1
if len(sqstm1_proteins) > 1:
    print("Multiple SQSTM1-related proteins found. Selecting canonical SQSTM1...")
    # Look for exact gene name match
    exact_match = sqstm1_proteins[sqstm1_proteins['GeneName'] == 'SQSTM1']
    if len(exact_match) > 0:
        sqstm1_index = exact_match.index[0]
        print(f"Selected protein at index {sqstm1_index} as canonical SQSTM1")
    else:
        sqstm1_index = sqstm1_proteins.index[0]
        print(f"No exact match found, using first match at index {sqstm1_index}")
elif len(sqstm1_proteins) == 1:
    sqstm1_index = sqstm1_proteins.index[0]
    print(f"Single SQSTM1 protein found at index {sqstm1_index}")
else:
    print("ERROR: No SQSTM1 protein found in the dataset!")
    sqstm1_index = None


# In[ ]:


# Extract SQSTM1 results from differential expression analysis
if sqstm1_index is not None:
    print("=== SQSTM1 DIFFERENTIAL EXPRESSION RESULTS ===")
    
    # Get the protein name/identifier from the results
    sqstm1_protein_name = adata.var.loc[sqstm1_index, 'GeneName']
    
    # Find SQSTM1 in the differential expression results
    # The protein names in results should match the column names from expression data
    sqstm1_result = de_results_clean[de_results_clean['protein'] == sqstm1_index]
    
    if len(sqstm1_result) == 0:
        # Try alternative protein identifiers
        print(f"Direct index match failed, searching by gene name: {sqstm1_protein_name}")
        # This might happen if protein names are stored differently
        # Let's check what protein identifiers look like in the results
        print(f"Sample protein names in results: {list(de_results_clean['protein'].head())}")
        
        # Try to match by checking if gene name appears in protein identifier
        potential_matches = de_results_clean[de_results_clean['protein'].astype(str).str.contains('SQSTM1', na=False)]
        
        if len(potential_matches) > 0:
            sqstm1_result = potential_matches.iloc[0:1]  # Take first match
            print(f"Found potential match: {sqstm1_result['protein'].iloc[0]}")
        else:
            print("Could not find SQSTM1 in differential expression results")
            sqstm1_result = None
    
    if sqstm1_result is not None and len(sqstm1_result) > 0:
        sqstm1_data = sqstm1_result.iloc[0]  # Get first (and likely only) result
        
        print(f"SQSTM1 Results:")
        print(f"  Protein identifier: {sqstm1_data['protein']}")
        print(f"  Log2 fold change: {sqstm1_data['log2_fold_change']:.3f}")
        print(f"  Raw p-value: {sqstm1_data['p_value']:.2e}")
        print(f"  FDR-corrected p-value: {sqstm1_data['fdr_corrected_pvalue']:.2e}")
        print(f"  95% CI: [{sqstm1_data['ci_lower']:.3f}, {sqstm1_data['ci_upper']:.3f}]")
        print(f"  Significant (FDR < 0.05): {sqstm1_data['significant_fdr']}")
        
        # Compare to claimed value
        claimed_log2fc = 3.41
        observed_log2fc = sqstm1_data['log2_fold_change']
        difference = abs(observed_log2fc - claimed_log2fc)
        
        print(f"\n  Claimed log2 FC: {claimed_log2fc}")
        print(f"  Observed log2 FC: {observed_log2fc:.3f}")
        print(f"  Difference: {difference:.3f}")
else:
    print("Cannot proceed with SQSTM1 analysis - protein not found")
    sqstm1_data = None


# In[ ]:


# Analyze ranking of SQSTM1 among all proteins
print("=== PROTEIN RANKING ANALYSIS ===")

# Sort all proteins by log2 fold change (descending)
ranked_proteins = de_results_clean.sort_values('log2_fold_change', ascending=False).reset_index(drop=True)

# Find SQSTM1's rank
if sqstm1_data is not None:
    sqstm1_protein_id = sqstm1_data['protein']
    sqstm1_rank = ranked_proteins[ranked_proteins['protein'] == sqstm1_protein_id].index
    
    if len(sqstm1_rank) > 0:
        sqstm1_rank = sqstm1_rank[0] + 1  # Convert to 1-based ranking
        
        print(f"SQSTM1 ranking: #{sqstm1_rank} out of {len(ranked_proteins)} proteins")
        
        # Show top 10 most upregulated proteins
        print(f"\nTop 10 most upregulated proteins:")
        for i, (idx, row) in enumerate(ranked_proteins.head(10).iterrows()):
            rank = i + 1
            is_sqstm1 = "*** SQSTM1 ***" if row['protein'] == sqstm1_protein_id else ""
            print(f"  {rank:2d}. {row['protein']}: {row['log2_fold_change']:.3f} (FDR p-val: {row['fdr_corrected_pvalue']:.2e}) {is_sqstm1}")
        
        # Validate the claim about being "most upregulated"
        is_most_upregulated = sqstm1_rank == 1
        
        print(f"\n=== CLAIM VALIDATION: MOST UPREGULATED ===")
        print(f"Is SQSTM1 the most upregulated protein? {is_most_upregulated}")
        
        if not is_most_upregulated:
            most_upregulated = ranked_proteins.iloc[0]
            print(f"Most upregulated protein is: {most_upregulated['protein']}")
            print(f"  Log2 FC: {most_upregulated['log2_fold_change']:.3f}")
            print(f"  FDR p-value: {most_upregulated['fdr_corrected_pvalue']:.2e}")
    else:
        print("ERROR: Could not find SQSTM1 in ranked results")
else:
    print("Cannot perform ranking analysis - SQSTM1 data not available")


# In[ ]:


# Validate SQSTM1 expression changes using raw data
print("=== SQSTM1 EXPRESSION VALIDATION ===")

if sqstm1_index is not None:
    # Extract SQSTM1 expression values
    sqstm1_expr = adata.to_df().iloc[:, sqstm1_index]
    
    # Group by tau status
    tau_positive_expr = sqstm1_expr[adata.obs['TauStatus'] == 'positive']
    tau_negative_expr = sqstm1_expr[adata.obs['TauStatus'] == 'negative']
    
    print(f"SQSTM1 Expression Summary:")
    print(f"\nTau-positive neurons (n={len(tau_positive_expr)}):")
    print(f"  Mean: {tau_positive_expr.mean():.3f}")
    print(f"  Std:  {tau_positive_expr.std():.3f}")
    print(f"  Range: {tau_positive_expr.min():.3f} - {tau_positive_expr.max():.3f}")
    
    print(f"\nTau-negative neurons (n={len(tau_negative_expr)}):")
    print(f"  Mean: {tau_negative_expr.mean():.3f}")
    print(f"  Std:  {tau_negative_expr.std():.3f}")
    print(f"  Range: {tau_negative_expr.min():.3f} - {tau_negative_expr.max():.3f}")
    
    # Calculate simple fold change for comparison
    simple_log2fc = tau_positive_expr.mean() - tau_negative_expr.mean()
    print(f"\nSimple log2 fold change: {simple_log2fc:.3f}")
    print(f"Covariate-adjusted log2 FC: {sqstm1_data['log2_fold_change']:.3f}" if sqstm1_data is not None else "N/A")
    
    # Perform independent t-test for validation
    t_stat, t_pval = stats.ttest_ind(tau_positive_expr, tau_negative_expr)
    print(f"\nIndependent t-test validation:")
    print(f"  t-statistic: {t_stat:.3f}")
    print(f"  p-value: {t_pval:.2e}")
    
    # Calculate effect size (Cohen's d)
    pooled_std = np.sqrt(((len(tau_positive_expr)-1)*tau_positive_expr.var() + 
                         (len(tau_negative_expr)-1)*tau_negative_expr.var()) / 
                        (len(tau_positive_expr) + len(tau_negative_expr) - 2))
    cohens_d = (tau_positive_expr.mean() - tau_negative_expr.mean()) / pooled_std
    print(f"  Cohen's d: {cohens_d:.3f}")
    
    # Interpretation of effect size
    if abs(cohens_d) < 0.2:
        effect_interpretation = "negligible"
    elif abs(cohens_d) < 0.5:
        effect_interpretation = "small"
    elif abs(cohens_d) < 0.8:
        effect_interpretation = "medium"
    else:
        effect_interpretation = "large"
    
    print(f"  Effect size interpretation: {effect_interpretation}")
else:
    print("Cannot perform expression validation - SQSTM1 not found")


# ## Results code
# 
# Put results in this section. This includes generating any final plots or calculating any final values, as well as verbose markdown explaining any conclusions based on the interpretation of the results.

# In[ ]:


# Create comprehensive visualization for SQSTM1 analysis
print("=== GENERATING SQSTM1 VISUALIZATION ===")

if sqstm1_index is not None and sqstm1_data is not None:
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: SQSTM1 expression by tau status
    sqstm1_expr = adata.to_df().iloc[:, sqstm1_index]
    tau_status = adata.obs['TauStatus']
    
    # Box plot
    tau_pos_data = sqstm1_expr[tau_status == 'positive']
    tau_neg_data = sqstm1_expr[tau_status == 'negative']
    
    box_data = [tau_neg_data, tau_pos_data]
    bp = axes[0,0].boxplot(box_data, labels=['Tau-negative', 'Tau-positive'], patch_artist=True)
    bp['boxes'][0].set_facecolor('lightblue')
    bp['boxes'][1].set_facecolor('lightcoral')
    
    axes[0,0].set_ylabel('SQSTM1 Expression (log2)')
    axes[0,0].set_title(f'SQSTM1 Expression by Tau Status\nLog2FC = {sqstm1_data["log2_fold_change"]:.3f}')
    axes[0,0].grid(True, alpha=0.3)
    
    # Add statistical annotation
    y_max = max(sqstm1_expr) * 1.05
    axes[0,0].plot([1, 2], [y_max, y_max], 'k-', alpha=0.5)
    axes[0,0].text(1.5, y_max * 1.02, f'p = {sqstm1_data["fdr_corrected_pvalue"]:.2e}', 
                  ha='center', va='bottom')
    
    # Plot 2: Ranking comparison (top 20 proteins)
    top20 = ranked_proteins.head(20)
    colors = ['red' if protein == sqstm1_data['protein'] else 'skyblue' for protein in top20['protein']]
    
    bars = axes[0,1].bar(range(len(top20)), top20['log2_fold_change'], color=colors)
    axes[0,1].set_xlabel('Protein Rank')
    axes[0,1].set_ylabel('Log2 Fold Change')
    axes[0,1].set_title('Top 20 Most Upregulated Proteins\n(SQSTM1 highlighted in red)')
    axes[0,1].set_xticks(range(0, 20, 2))
    axes[0,1].set_xticklabels(range(1, 21, 2))
    axes[0,1].grid(True, alpha=0.3)
    
    # Add horizontal line at claimed value
    axes[0,1].axhline(y=3.41, color='green', linestyle='--', alpha=0.7, 
                     label='Claimed value (3.41)')
    axes[0,1].legend()
    
    # Plot 3: SQSTM1 expression vs covariates
    # Scatter plot with age
    age = adata.obs['Age at death']
    colors_scatter = ['red' if status == 'positive' else 'blue' for status in tau_status]
    
    scatter = axes[1,0].scatter(age, sqstm1_expr, c=colors_scatter, alpha=0.7, s=50)
    axes[1,0].set_xlabel('Age at Death')
    axes[1,0].set_ylabel('SQSTM1 Expression (log2)')
    axes[1,0].set_title('SQSTM1 Expression vs Age\n(Red=Tau+, Blue=Tau-)')
    axes[1,0].grid(True, alpha=0.3)
    
    # Add trend line
    z = np.polyfit(age, sqstm1_expr, 1)
    p = np.poly1d(z)
    axes[1,0].plot(age, p(age), "k--", alpha=0.5, label=f'Trend (r={np.corrcoef(age, sqstm1_expr)[0,1]:.3f})')
    axes[1,0].legend()
    
    # Plot 4: Confidence interval visualization
    # Show SQSTM1 fold change with confidence interval
    fc_value = sqstm1_data['log2_fold_change']
    ci_lower = sqstm1_data['ci_lower']
    ci_upper = sqstm1_data['ci_upper']
    
    axes[1,1].errorbar([1], [fc_value], yerr=[[fc_value - ci_lower], [ci_upper - fc_value]], 
                      fmt='ro', markersize=10, capsize=5, capthick=2, label='Observed')
    axes[1,1].axhline(y=3.41, color='green', linestyle='--', linewidth=2, label='Claimed (3.41)')
    axes[1,1].axhline(y=0, color='black', linestyle='-', alpha=0.3, label='No change')
    
    axes[1,1].set_xlim(0.5, 1.5)
    axes[1,1].set_ylabel('Log2 Fold Change')
    axes[1,1].set_title('SQSTM1 Fold Change\nwith 95% Confidence Interval')
    axes[1,1].set_xticks([1])
    axes[1,1].set_xticklabels(['SQSTM1'])
    axes[1,1].legend()
    axes[1,1].grid(True, alpha=0.3)
    
    plt.suptitle('SQSTM1 Comprehensive Analysis', fontsize=16, y=1.02)
    plt.tight_layout()
    plt.show()
    
    print("SQSTM1 visualization completed")
else:
    print("Cannot generate visualization - SQSTM1 data not available")


# In[ ]:


# Final assessment of Claim 2
print("=== CLAIM 2 COMPREHENSIVE ASSESSMENT ===")
print("\nClaim: 'The autophagy receptor SQSTM1 was the most upregulated protein")
print("(+3.41 log2) in the covariate-controlled differential expression analysis.'")

if sqstm1_data is not None:
    # Extract key values
    observed_log2fc = sqstm1_data['log2_fold_change']
    claimed_log2fc = 3.41
    fdr_pvalue = sqstm1_data['fdr_corrected_pvalue']
    is_significant = sqstm1_data['significant_fdr']
    
    print(f"\n=== QUANTITATIVE ASSESSMENT ===")
    print(f"Claimed log2 fold change: +{claimed_log2fc}")
    print(f"Observed log2 fold change: {observed_log2fc:+.3f}")
    print(f"Absolute difference: {abs(observed_log2fc - claimed_log2fc):.3f}")
    print(f"Relative difference: {abs(observed_log2fc - claimed_log2fc)/claimed_log2fc*100:.1f}%")
    
    print(f"\nStatistical significance: {is_significant}")
    print(f"FDR-corrected p-value: {fdr_pvalue:.2e}")
    print(f"95% CI: [{sqstm1_data['ci_lower']:.3f}, {sqstm1_data['ci_upper']:.3f}]")
    
    # Check if claimed value falls within confidence interval
    ci_contains_claim = sqstm1_data['ci_lower'] <= claimed_log2fc <= sqstm1_data['ci_upper']
    print(f"\nClaimed value within 95% CI: {ci_contains_claim}")
    
    # Ranking assessment
    sqstm1_protein_id = sqstm1_data['protein']
    sqstm1_rank = ranked_proteins[ranked_proteins['protein'] == sqstm1_protein_id].index[0] + 1
    
    print(f"\n=== RANKING ASSESSMENT ===")
    print(f"SQSTM1 rank among all proteins: #{sqstm1_rank}")
    print(f"Is most upregulated protein: {sqstm1_rank == 1}")
    
    if sqstm1_rank > 1:
        most_upregulated = ranked_proteins.iloc[0]
        print(f"\nActual most upregulated protein:")
        print(f"  Protein: {most_upregulated['protein']}")
        print(f"  Log2FC: {most_upregulated['log2_fold_change']:.3f}")
        print(f"  FDR p-value: {most_upregulated['fdr_corrected_pvalue']:.2e}")
    
    # Biological assessment
    print(f"\n=== BIOLOGICAL ASSESSMENT ===")
    print(f"SQSTM1 is significantly upregulated: {is_significant}")
    print(f"Magnitude of upregulation: {2**observed_log2fc:.1f}-fold linear increase")
    print(f"Effect size (Cohen's d): {cohens_d:.2f} ({effect_interpretation})")
    
    # Overall claim validation
    print(f"\n=== CLAIM VALIDATION SUMMARY ===")
    
    # Assess each component of the claim
    fold_change_accurate = abs(observed_log2fc - claimed_log2fc) < 0.5  # Within 0.5 log2 units
    statistically_significant = is_significant
    is_most_upregulated = sqstm1_rank == 1
    
    print(f"1. Fold change accuracy (+3.41 log2): {fold_change_accurate}")
    print(f"   • Observed: {observed_log2fc:.3f}, Claimed: {claimed_log2fc}")
    print(f"   • Difference: {abs(observed_log2fc - claimed_log2fc):.3f} log2 units")
    
    print(f"\n2. Statistical significance: {statistically_significant}")
    print(f"   • FDR-corrected p-value: {fdr_pvalue:.2e}")
    
    print(f"\n3. Most upregulated status: {is_most_upregulated}")
    print(f"   • Rank: #{sqstm1_rank} out of {len(ranked_proteins)} proteins")
    
    # Final validation conclusion
    if all([fold_change_accurate, statistically_significant, is_most_upregulated]):
        validation_status = "FULLY SUPPORTED"
    elif fold_change_accurate and statistically_significant:
        validation_status = "PARTIALLY SUPPORTED"
    elif statistically_significant:
        validation_status = "WEAKLY SUPPORTED"
    else:
        validation_status = "NOT SUPPORTED"
    
    print(f"\n=== FINAL VALIDATION ===")
    print(f"Status: {validation_status}")
    
    if validation_status == "FULLY SUPPORTED":
        print("All aspects of the claim are confirmed by the analysis.")
    elif validation_status == "PARTIALLY SUPPORTED":
        print("SQSTM1 shows significant upregulation with correct magnitude,")
        print("but may not be the single most upregulated protein.")
    elif validation_status == "WEAKLY SUPPORTED":
        print("SQSTM1 is significantly upregulated but with different magnitude")
        print("and/or ranking than claimed.")
    else:
        print("The claim is not supported by the statistical analysis.")
        
    print(f"\nNote: Results depend on the exact dataset and analysis parameters.")
    print(f"SQSTM1 upregulation is biologically consistent with autophagy dysfunction")
    print(f"in neurodegeneration, regardless of its exact ranking.")
    
else:
    print("\n=== ANALYSIS LIMITATION ===")
    print("Could not complete SQSTM1 analysis due to protein identification issues.")
    print("This may be due to differences in protein naming conventions or")
    print("data processing pipelines between the original study and this analysis.")

