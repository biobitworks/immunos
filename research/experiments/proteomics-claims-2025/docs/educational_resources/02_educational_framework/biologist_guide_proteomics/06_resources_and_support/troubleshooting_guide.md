# 🔧 Troubleshooting Guide: Solutions for Common Problems

## 🎯 What This Guide Covers

This comprehensive troubleshooting guide helps you solve:
- ✅ **Python and package installation issues**
- ✅ **Data loading and file format problems**
- ✅ **Statistical analysis errors and warnings**
- ✅ **Visualization and plotting issues**
- ✅ **Memory and performance problems**
- ✅ **Results interpretation challenges**

---

## 🐍 Python and Package Issues

### Installation Problems

#### Problem: "Python not found" or "command not found"
```bash
# Error messages:
# 'python' is not recognized as an internal or external command
# bash: python: command not found
# No module named 'pandas'
```

**Solutions:**

**Option 1: Check Python Installation**
```bash
# Test if Python is installed
python --version
python3 --version

# If neither works, Python isn't installed or not in PATH
# Install Python from: https://www.python.org/downloads/
```

**Option 2: Fix PATH Issues**
```bash
# On Windows - Add Python to PATH
# 1. Find Python installation (usually C:\Python39\ or similar)
# 2. Add to system PATH environment variable
# 3. Restart command prompt

# On Mac/Linux - Add to shell profile
echo 'export PATH="/usr/local/bin/python3:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

**Option 3: Use Conda/Anaconda**
```bash
# Download and install Anaconda
# https://www.anaconda.com/products/distribution

# Create new environment
conda create -n proteomics python=3.9
conda activate proteomics
```

#### Problem: Package Installation Fails
```bash
# Common error messages:
# ERROR: Could not install packages due to an EnvironmentError
# Permission denied
# No matching distribution found
```

**Solutions:**

**Permission Issues:**
```bash
# Use --user flag (recommended)
pip install --user pandas numpy scipy

# Or use virtual environment
python -m venv proteomics_env
source proteomics_env/bin/activate  # On Windows: proteomics_env\Scripts\activate
pip install pandas numpy scipy
```

**Network/Firewall Issues:**
```bash
# Use different index
pip install --index-url https://pypi.org/simple/ pandas

# Download and install offline
# 1. Download .whl files from https://pypi.org/
# 2. pip install downloaded_file.whl
```

**Version Conflicts:**
```bash
# Update pip first
python -m pip install --upgrade pip

# Install specific versions
pip install pandas==1.3.3 numpy==1.21.2

# Use conda for complex dependencies
conda install pandas numpy scipy matplotlib seaborn
```

### Package Import Errors

#### Problem: "ModuleNotFoundError" after installation
```python
# Error: ModuleNotFoundError: No module named 'scanpy'
import scanpy as sc
```

**Solutions:**

**Check Installation:**
```bash
# Verify package is installed
pip list | grep scanpy
pip show scanpy

# If not installed
pip install scanpy
```

**Virtual Environment Issues:**
```python
# Make sure you're in the right environment
import sys
print(sys.executable)
print(sys.path)

# Should show your intended Python installation
```

**Jupyter Kernel Issues:**
```bash
# Install ipykernel in your environment
pip install ipykernel

# Add environment to Jupyter
python -m ipykernel install --user --name=proteomics_env

# In Jupyter, select the right kernel:
# Kernel → Change Kernel → proteomics_env
```

#### Problem: Version Compatibility Issues
```python
# Error: ImportError: cannot import name 'soft_unicode' from 'markupsafe'
# Error: AttributeError: module 'pandas' has no attribute 'DataFrame'
```

**Solutions:**

**Create Clean Environment:**
```bash
# Remove problematic environment
conda env remove -n old_env

# Create new environment with specific versions
conda create -n proteomics_new python=3.9
conda activate proteomics_new

# Install compatible versions
conda install pandas=1.3.3 numpy=1.21.2 scipy=1.7.1
conda install matplotlib=3.4.3 seaborn=0.11.2
pip install scanpy==1.8.1
```

**Check Version Compatibility:**
```python
# Check installed versions
import pandas as pd
import numpy as np
import scipy
print(f"Pandas: {pd.__version__}")
print(f"NumPy: {np.__version__}")
print(f"SciPy: {scipy.__version__}")

# Compare with requirements in tutorials
```

---

## 📁 Data Loading Issues

### File Not Found Errors

#### Problem: "FileNotFoundError" when loading data
```python
# Error: FileNotFoundError: [Errno 2] No such file or directory: 'pool_processed_v2.h5ad'
adata = sc.read_h5ad('pool_processed_v2.h5ad')
```

**Solutions:**

**Check File Location:**
```python
import os

# Check current directory
print("Current directory:", os.getcwd())

# List files in current directory
print("Files available:", os.listdir('.'))

# Check if file exists
if os.path.exists('pool_processed_v2.h5ad'):
    print("File found!")
else:
    print("File not found - check location")
```

**Fix File Path:**
```python
# Use absolute path
adata = sc.read_h5ad('/full/path/to/pool_processed_v2.h5ad')

# Or relative path from current location
adata = sc.read_h5ad('../data/pool_processed_v2.h5ad')

# Or change working directory
os.chdir('/path/to/data/folder')
adata = sc.read_h5ad('pool_processed_v2.h5ad')
```

**Google Colab File Access:**
```python
# Mount Google Drive
from google.colab import drive
drive.mount('/content/drive')

# Access file in Drive
adata = sc.read_h5ad('/content/drive/MyDrive/proteomics_data/pool_processed_v2.h5ad')

# Or upload file directly
from google.colab import files
uploaded = files.upload()
# Then select and upload your .h5ad file
```

### File Format Issues

#### Problem: "UnicodeDecodeError" or corrupted file
```python
# Error: UnicodeDecodeError: 'utf-8' codec can't decode byte
# Error: OSError: Unable to open file (truncated file)
```

**Solutions:**

**Check File Integrity:**
```python
import os

# Check file size
file_size = os.path.getsize('pool_processed_v2.h5ad')
print(f"File size: {file_size / (1024**2):.1f} MB")

# If file is unexpectedly small (< 10 MB), it may be corrupted
```

**Re-download Data:**
```bash
# If file is corrupted, re-download from source
# Check download completion
# Verify file hash if provided
```

**Alternative Loading Methods:**
```python
# Try different readers
try:
    adata = sc.read_h5ad('pool_processed_v2.h5ad')
except Exception as e:
    print(f"Error loading with scanpy: {e}")

    # Try with h5py directly
    import h5py
    with h5py.File('pool_processed_v2.h5ad', 'r') as f:
        print("File structure:", list(f.keys()))
```

#### Problem: Wrong file format or extension
```python
# Error: ValueError: Unknown file format
# Error: File format not supported
```

**Solutions:**

**Check File Extension:**
```python
import os

# Check actual file type
filename = 'pool_processed_v2.h5ad'
_, extension = os.path.splitext(filename)
print(f"File extension: {extension}")

# For different formats:
if extension == '.csv':
    data = pd.read_csv(filename)
elif extension == '.xlsx':
    data = pd.read_excel(filename)
elif extension == '.h5ad':
    data = sc.read_h5ad(filename)
elif extension == '.h5':
    data = pd.read_hdf(filename)
```

**Convert Between Formats:**
```python
# If you have CSV instead of H5AD
df = pd.read_csv('proteomics_data.csv', index_col=0)

# Create AnnData object
import anndata as ad
adata = ad.AnnData(X=df.values)
adata.obs_names = df.index
adata.var_names = df.columns

# Save as H5AD
adata.write('converted_data.h5ad')
```

### Memory Issues with Large Files

#### Problem: "MemoryError" when loading data
```python
# Error: MemoryError: Unable to allocate array
# Error: Out of memory
```

**Solutions:**

**Check Available Memory:**
```python
import psutil
import os

# Check system memory
memory = psutil.virtual_memory()
print(f"Available memory: {memory.available / (1024**3):.1f} GB")
print(f"Memory usage: {memory.percent}%")

# Check file size
file_size = os.path.getsize('pool_processed_v2.h5ad')
print(f"File size: {file_size / (1024**3):.1f} GB")
```

**Optimize Memory Usage:**
```python
# Load data in chunks or subsets
# For CSV files
chunk_iter = pd.read_csv('large_file.csv', chunksize=1000)
for chunk in chunk_iter:
    # Process each chunk
    process_chunk(chunk)

# For H5AD files, load metadata first
adata = sc.read_h5ad('pool_processed_v2.h5ad', backed='r')  # Read-only mode
print(adata)  # Check size before loading into memory

# Load subset of data
adata_subset = adata[:100, :1000].copy()  # First 100 samples, 1000 proteins
```

**Increase Virtual Memory:**
```bash
# On Windows - Increase page file size
# System Properties → Advanced → Performance Settings → Virtual Memory

# On Mac/Linux - Check swap space
free -h
# Add swap if needed (requires admin rights)
```

---

## 📊 Statistical Analysis Errors

### Common Statistical Warnings

#### Problem: "RuntimeWarning: invalid value encountered"
```python
# Warning: RuntimeWarning: invalid value encountered in divide
# Warning: RuntimeWarning: invalid value encountered in log
```

**Solutions:**

**Check for Problematic Values:**
```python
import numpy as np

# Check for NaN, inf, or zero values
data = your_expression_data
print(f"NaN values: {np.isnan(data).sum()}")
print(f"Infinite values: {np.isinf(data).sum()}")
print(f"Zero values: {(data == 0).sum()}")
print(f"Negative values: {(data < 0).sum()}")
```

**Handle Missing/Invalid Data:**
```python
# Remove or replace problematic values
data_clean = data.copy()

# Replace zeros with small positive number before log transformation
data_clean[data_clean == 0] = 1e-6

# Remove infinite values
data_clean = data_clean[~np.isinf(data_clean)]

# Handle NaN values
data_clean = data_clean.dropna()  # Remove
# OR
data_clean = data_clean.fillna(data_clean.median())  # Impute
```

**Safe Mathematical Operations:**
```python
# Safe logarithm
def safe_log2(x):
    return np.log2(np.maximum(x, 1e-10))

# Safe division
def safe_divide(a, b):
    return np.divide(a, b, out=np.zeros_like(a), where=(b != 0))

# Use these functions in your analysis
log2_fc = safe_log2(group1.mean()) - safe_log2(group2.mean())
```

#### Problem: "Warning: p-value clamping" or "p-value overflow"
```python
# Warning: p-value was clamped to 1.0
# Error: OverflowError in statistical calculation
```

**Solutions:**

**Check Sample Sizes:**
```python
# Verify adequate sample sizes
group1_size = len(group1)
group2_size = len(group2)

print(f"Group 1 size: {group1_size}")
print(f"Group 2 size: {group2_size}")

# Need at least 3 samples per group for t-test
if min(group1_size, group2_size) < 3:
    print("Warning: Sample size too small for t-test")
    # Use non-parametric test instead
    from scipy.stats import mannwhitneyu
    stat, pval = mannwhitneyu(group1, group2)
```

**Handle Extreme Values:**
```python
# Check for extreme values
print(f"Group 1 range: {group1.min():.3f} to {group1.max():.3f}")
print(f"Group 2 range: {group2.min():.3f} to {group2.max():.3f}")

# Remove outliers if necessary
def remove_outliers(data, z_threshold=3):
    z_scores = np.abs((data - data.mean()) / data.std())
    return data[z_scores < z_threshold]

group1_clean = remove_outliers(group1)
group2_clean = remove_outliers(group2)
```

**Use Robust Statistical Methods:**
```python
# When distributions are problematic
from scipy.stats import mannwhitneyu, wilcoxon

# Non-parametric alternative to t-test
stat, pval = mannwhitneyu(group1, group2, alternative='two-sided')

# For paired data
stat, pval = wilcoxon(group1, group2, alternative='two-sided')
```

### Multiple Testing Correction Issues

#### Problem: "All p-values become non-significant after correction"
```python
# After FDR correction, no proteins are significant
from statsmodels.stats.multitest import multipletests
rejected, pvals_corrected, _, _ = multipletests(pvals, method='fdr_bh')
print(f"Significant after correction: {sum(rejected)}")  # Returns 0
```

**Solutions:**

**Check Original P-values:**
```python
# Examine p-value distribution
import matplotlib.pyplot as plt

plt.hist(pvals, bins=50, alpha=0.7)
plt.xlabel('P-values')
plt.ylabel('Frequency')
plt.title('P-value Distribution')
plt.show()

print(f"Raw p-values < 0.05: {sum(np.array(pvals) < 0.05)}")
print(f"Raw p-values < 0.01: {sum(np.array(pvals) < 0.01)}")
print(f"Minimum p-value: {min(pvals)}")
```

**Adjust Significance Threshold:**
```python
# Try different correction methods
methods = ['fdr_bh', 'fdr_by', 'bonferroni', 'holm']

for method in methods:
    rejected, pvals_corr, _, _ = multipletests(pvals, method=method)
    print(f"{method}: {sum(rejected)} significant")

# Use less conservative method if appropriate
rejected, pvals_fdr, _, _ = multipletests(pvals, method='fdr_bh', alpha=0.1)
```

**Consider Effect Sizes:**
```python
# Focus on biologically meaningful effects
# Don't rely only on p-values
effect_sizes = calculate_cohens_d(group1, group2)

# Combined criteria
significant_stats = (np.array(pvals) < 0.05)
large_effects = (np.abs(effect_sizes) > 0.5)
meaningful = significant_stats & large_effects

print(f"Statistically significant: {sum(significant_stats)}")
print(f"Large effect sizes: {sum(large_effects)}")
print(f"Both criteria: {sum(meaningful)}")
```

---

## 📈 Visualization Problems

### Matplotlib/Seaborn Issues

#### Problem: "Figures not displaying" in Jupyter
```python
# Plots don't appear in notebook
plt.plot([1, 2, 3, 4])
plt.show()  # Nothing appears
```

**Solutions:**

**Enable Inline Plotting:**
```python
# In Jupyter notebook
%matplotlib inline
import matplotlib.pyplot as plt

# Or for interactive plots
%matplotlib widget
```

**Check Backend:**
```python
import matplotlib
print(f"Backend: {matplotlib.get_backend()}")

# Change backend if needed
matplotlib.use('Agg')  # For file output
# OR
matplotlib.use('TkAgg')  # For display
```

**Force Display:**
```python
# Explicitly show plots
plt.figure(figsize=(10, 6))
plt.plot([1, 2, 3, 4])
plt.show()

# Or save to file
plt.savefig('test_plot.png', dpi=300, bbox_inches='tight')
plt.close()  # Close figure to save memory
```

#### Problem: "Figure is too small" or "Text is unreadable"
```python
# Plot appears but is tiny or text is too small
```

**Solutions:**

**Adjust Figure Size:**
```python
# Set figure size globally
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 14

# Or for individual plots
fig, ax = plt.subplots(figsize=(15, 10))

# For subplots
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
```

**Fix Text Size:**
```python
# Adjust font sizes
plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12
})

# Or manually for each plot
plt.xlabel('X-axis', fontsize=14)
plt.ylabel('Y-axis', fontsize=14)
plt.title('Title', fontsize=16)
```

**High-DPI Displays:**
```python
# For high-resolution displays
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 300

# In Jupyter
%config InlineBackend.figure_format = 'retina'
```

#### Problem: "Memory error" with many plots
```python
# MemoryError when creating multiple plots
# Figures accumulate in memory
```

**Solutions:**

**Close Figures:**
```python
# Close individual figures
plt.figure()
plt.plot([1, 2, 3])
plt.savefig('plot1.png')
plt.close()  # Important!

# Close all figures
plt.close('all')
```

**Use Context Managers:**
```python
# Automatic cleanup
with plt.figure() as fig:
    plt.plot([1, 2, 3])
    plt.savefig('plot.png')
# Figure automatically closed

# Or with subplots
fig, ax = plt.subplots()
try:
    ax.plot([1, 2, 3])
    plt.savefig('plot.png')
finally:
    plt.close(fig)
```

**Generate Plots in Batches:**
```python
# For many plots
def create_plot_batch(data_list, batch_size=10):
    for i in range(0, len(data_list), batch_size):
        batch = data_list[i:i+batch_size]

        for j, data in enumerate(batch):
            plt.figure()
            plt.plot(data)
            plt.savefig(f'plot_{i+j}.png')
            plt.close()

        # Clear memory between batches
        plt.close('all')
```

### Seaborn-Specific Issues

#### Problem: "Seaborn plots look different than expected"
```python
# Plots don't match tutorial examples
# Colors or styles are different
```

**Solutions:**

**Set Seaborn Style:**
```python
import seaborn as sns

# Set consistent style
sns.set_style("whitegrid")
sns.set_palette("husl")

# Or specific style
sns.set_theme(style="darkgrid", palette="muted")

# Reset to defaults if needed
sns.reset_defaults()
```

**Version Compatibility:**
```python
# Check seaborn version
print(f"Seaborn version: {sns.__version__}")

# Some function names changed between versions
# Old: sns.distplot() → New: sns.histplot()
# Old: sns.boxplot() → Still works
```

**Explicit Parameters:**
```python
# Be explicit about parameters
sns.boxplot(data=df, x='group', y='value',
           palette='Set1', width=0.6)

# Rather than relying on defaults
sns.boxplot(data=df, x='group', y='value')
```

---

## 💾 Memory and Performance Issues

### Large Dataset Problems

#### Problem: "Analysis is extremely slow"
```python
# Code runs but takes hours
# System becomes unresponsive
```

**Solutions:**

**Profile Your Code:**
```python
import time

# Time critical sections
start_time = time.time()
result = your_slow_function(data)
end_time = time.time()
print(f"Function took {end_time - start_time:.2f} seconds")

# Use line profiler for detailed analysis
%load_ext line_profiler
%lprun -f your_function your_function(data)
```

**Optimize Data Types:**
```python
# Check memory usage
import pandas as pd

print(df.info(memory_usage='deep'))

# Optimize data types
df['integer_column'] = df['integer_column'].astype('int32')  # Instead of int64
df['category_column'] = df['category_column'].astype('category')
df['float_column'] = df['float_column'].astype('float32')  # Instead of float64
```

**Vectorize Operations:**
```python
# Slow: Loop-based
result = []
for i in range(len(data)):
    result.append(data[i] * 2 + 1)

# Fast: Vectorized
result = data * 2 + 1

# Use pandas/numpy operations
df['new_col'] = df['col1'] * df['col2']  # Fast
# Instead of df.apply(lambda x: x['col1'] * x['col2'], axis=1)  # Slow
```

**Process in Chunks:**
```python
def process_large_dataset(df, chunk_size=1000):
    results = []

    for i in range(0, len(df), chunk_size):
        chunk = df.iloc[i:i+chunk_size]
        chunk_result = process_chunk(chunk)
        results.append(chunk_result)

        # Optional: Progress reporting
        if i % (chunk_size * 10) == 0:
            print(f"Processed {i}/{len(df)} rows")

    return pd.concat(results, ignore_index=True)
```

#### Problem: "Jupyter kernel keeps dying"
```python
# Kernel restarts unexpectedly
# "The kernel appears to have died" message
```

**Solutions:**

**Increase Memory Limits:**
```bash
# Start Jupyter with more memory
jupyter notebook --NotebookApp.max_buffer_size=1000000000

# Or in config file
c.NotebookApp.max_buffer_size = 1000000000
```

**Monitor Memory Usage:**
```python
import psutil
import os

def check_memory():
    process = psutil.Process(os.getpid())
    memory_mb = process.memory_info().rss / 1024 / 1024
    print(f"Memory usage: {memory_mb:.1f} MB")
    return memory_mb

# Check regularly
check_memory()

# Set memory alerts
def memory_warning(threshold_mb=8000):
    if check_memory() > threshold_mb:
        print(f"WARNING: Memory usage above {threshold_mb} MB")
        print("Consider restarting kernel or optimizing code")
```

**Clear Variables:**
```python
# Clear large variables when done
del large_dataframe
del analysis_results

# Garbage collection
import gc
gc.collect()

# Reset namespace (nuclear option)
%reset -f
```

---

## 🤔 Results Interpretation Issues

### Statistical Results Don't Make Sense

#### Problem: "All proteins are significant" or "Nothing is significant"
```python
# Unrealistic results
# Either everything or nothing passes significance tests
```

**Solutions:**

**Check Data Distribution:**
```python
# Examine your data
print("Data summary:")
print(df.describe())

# Check for outliers
Q1 = df.quantile(0.25)
Q3 = df.quantile(0.75)
IQR = Q3 - Q1
outliers = ((df < (Q1 - 1.5 * IQR)) | (df > (Q3 + 1.5 * IQR))).sum()
print(f"Outliers per column: {outliers}")

# Visualize distributions
import matplotlib.pyplot as plt
df.hist(bins=50, figsize=(15, 10))
plt.tight_layout()
plt.show()
```

**Verify Sample Labels:**
```python
# Make sure groups are labeled correctly
print("Group distribution:")
print(metadata['group'].value_counts())

# Check for mislabeling
sample_check = pd.DataFrame({
    'sample_id': metadata.index,
    'group': metadata['group'],
    'mean_expression': df.mean(axis=1)
})

# Plot by group
sns.boxplot(data=sample_check, x='group', y='mean_expression')
plt.title('Mean Expression by Group')
plt.show()
```

**Validate Statistical Tests:**
```python
# Test with known controls
# Use housekeeping proteins that shouldn't change
housekeeping = ['ACTB', 'GAPDH', 'TUBB']
for protein in housekeeping:
    if protein in df.columns:
        group1 = df[df['group'] == 'case'][protein]
        group2 = df[df['group'] == 'control'][protein]

        stat, pval = ttest_ind(group1, group2)
        print(f"{protein}: p-value = {pval:.3f}")

# Housekeeping proteins should have p > 0.05
```

#### Problem: "Results contradict literature"
```python
# Your findings don't match published studies
# Effect directions are opposite
```

**Solutions:**

**Check Analysis Direction:**
```python
# Verify which group is "up" vs "down"
group1_mean = df[df['group'] == 'case']['protein_X'].mean()
group2_mean = df[df['group'] == 'control']['protein_X'].mean()

print(f"Case group mean: {group1_mean:.3f}")
print(f"Control group mean: {group2_mean:.3f}")
print(f"Fold change: {group1_mean / group2_mean:.3f}")

# Make sure you're interpreting direction correctly
if group1_mean > group2_mean:
    print("Protein X is higher in cases")
else:
    print("Protein X is lower in cases")
```

**Consider Study Differences:**
```python
# Compare study characteristics
print("Your study characteristics:")
print(f"- Sample type: {sample_type}")
print(f"- Sample size: {sample_size}")
print(f"- Population: {population}")
print(f"- Method: {method}")

# Literature differences might explain discrepancies:
# - Different sample types (tissue vs blood vs CSF)
# - Different disease stages
# - Different populations
# - Different measurement methods
```

**Validate Key Findings:**
```python
# Focus on most robust results
significant_proteins = results[results['fdr_corrected_p'] < 0.05]
large_effects = significant_proteins[abs(significant_proteins['cohens_d']) > 0.8]

print(f"Most robust findings ({len(large_effects)} proteins):")
for idx, row in large_effects.iterrows():
    print(f"- {row['protein']}: FC={row['fold_change']:.2f}, d={row['cohens_d']:.2f}")

# Literature search these specific proteins
```

### Pathway Analysis Problems

#### Problem: "No pathways are enriched"
```python
# Pathway analysis returns no significant results
# All FDR values > 0.05
```

**Solutions:**

**Check Gene List Quality:**
```python
# Verify your gene list
gene_list = ['SQSTM1', 'VDAC1', 'ACTB']  # Your significant genes

print(f"Number of genes: {len(gene_list)}")
print(f"Unique genes: {len(set(gene_list))}")

# Check for valid gene symbols
import requests

def check_gene_symbols(genes):
    # Simple validation (you could use more sophisticated APIs)
    valid_genes = []
    invalid_genes = []

    for gene in genes:
        if gene.isalpha() and gene.isupper() and len(gene) > 1:
            valid_genes.append(gene)
        else:
            invalid_genes.append(gene)

    return valid_genes, invalid_genes

valid, invalid = check_gene_symbols(gene_list)
print(f"Valid gene symbols: {len(valid)}")
print(f"Invalid gene symbols: {invalid}")
```

**Adjust Analysis Parameters:**
```python
# Try different databases
databases = ['GO_Biological_Process', 'KEGG_2021', 'Reactome_2022']

# Try less stringent criteria
p_thresholds = [0.05, 0.1, 0.2]
min_genes = [3, 5, 10]

for db in databases:
    for p_thresh in p_thresholds:
        for min_gene in min_genes:
            # Run pathway analysis with these parameters
            print(f"Testing {db}, p<{p_thresh}, min_genes={min_gene}")
            # Your pathway analysis code here
```

**Use Alternative Approaches:**
```python
# Try different pathway tools
# 1. DAVID (web-based)
# 2. Enrichr (web-based)
# 3. g:Profiler (web-based)
# 4. clusterProfiler (R package)

# Or manual pathway testing
pathway_genes = {
    'Autophagy': ['SQSTM1', 'LC3B', 'ATG5', 'ATG7'],
    'Proteasome': ['PSMA1', 'PSMB1', 'PSMD1'],
    # Add more pathways
}

from scipy.stats import hypergeom

def manual_pathway_test(gene_list, pathway_genes, background_size):
    overlap = len(set(gene_list) & set(pathway_genes))
    pathway_size = len(pathway_genes)
    list_size = len(gene_list)

    pval = hypergeom.sf(overlap-1, background_size, pathway_size, list_size)
    return overlap, pval

# Test each pathway
for pathway, genes in pathway_genes.items():
    overlap, pval = manual_pathway_test(gene_list, genes, 20000)
    print(f"{pathway}: {overlap} genes, p={pval:.3f}")
```

---

## 🆘 Getting Additional Help

### Documentation and Resources

#### Official Documentation
```python
# Package documentation
import pandas as pd
help(pd.read_csv)  # Function help

# Online documentation
# Pandas: https://pandas.pydata.org/docs/
# NumPy: https://numpy.org/doc/
# SciPy: https://docs.scipy.org/
# Matplotlib: https://matplotlib.org/stable/
# Seaborn: https://seaborn.pydata.org/
# Scanpy: https://scanpy.readthedocs.io/
```

#### Community Forums
```python
# Where to get help:
"""
STACK OVERFLOW:
- https://stackoverflow.com/
- Tag your questions with relevant packages
- Include minimal reproducible examples

BIOSTARS:
- https://www.biostars.org/
- Bioinformatics-specific questions
- Proteomics and genomics community

GITHUB ISSUES:
- Package-specific issues
- Bug reports and feature requests
- Check existing issues first

REDDIT COMMUNITIES:
- r/bioinformatics
- r/learnpython
- r/statistics
"""
```

### Creating Good Help Requests

#### How to Ask for Help Effectively
```python
# Include this information:
"""
1. CLEAR PROBLEM DESCRIPTION
   - What you're trying to do
   - What you expected
   - What actually happened

2. REPRODUCIBLE EXAMPLE
   - Minimal code that shows the problem
   - Sample data (if possible)
   - Full error message

3. SYSTEM INFORMATION
   - Operating system
   - Python version
   - Package versions
   - Environment (Jupyter, command line, etc.)

4. WHAT YOU'VE TRIED
   - Solutions you've attempted
   - Error messages received
   - Any partial success
"""
```

#### Example Good Help Request
```python
# Title: "Unable to load H5AD file - FileNotFoundError in scanpy"

# Description:
"""
I'm trying to load a proteomics dataset following the tutorial, but getting a FileNotFoundError.

## What I'm trying to do:
Load an H5AD file containing proteomics data for analysis.

## Code that fails:
```python
import scanpy as sc
adata = sc.read_h5ad('pool_processed_v2.h5ad')
```

## Error message:
FileNotFoundError: [Errno 2] No such file or directory: 'pool_processed_v2.h5ad'

## What I've tried:
1. Verified file exists with os.path.exists() - returns True
2. Tried absolute path - same error
3. File size is 45 MB, seems reasonable

## System info:
- Windows 10
- Python 3.9.7
- Scanpy 1.8.1
- Running in Jupyter notebook

## Sample file structure:
```
/project/
  /data/
    pool_processed_v2.h5ad
  /notebooks/
    analysis.ipynb (current location)
```

Any help would be appreciated!
"""
```

### Emergency Troubleshooting Checklist

#### When Everything Goes Wrong
```python
# Emergency reset procedure:
"""
1. RESTART KERNEL
   - Jupyter: Kernel → Restart & Clear Output
   - Clear all variables and imports

2. CHECK BASIC FUNCTIONALITY
   - import pandas as pd
   - print("Hello world")
   - pd.DataFrame({'a': [1, 2, 3]})

3. RELOAD DATA
   - Start with minimal data loading
   - Check file existence and permissions
   - Try alternative loading methods

4. SIMPLIFY ANALYSIS
   - Remove complex operations
   - Test with subset of data
   - Add debugging prints

5. ENVIRONMENT RESET
   - Deactivate/reactivate virtual environment
   - Reinstall problematic packages
   - Create fresh environment if needed

6. SYSTEM RESTART
   - Close all applications
   - Restart computer
   - Fresh start with minimal programs
"""
```

---

## 🎯 Prevention is Better Than Cure

### Best Practices to Avoid Problems

#### Code Organization
```python
# Write defensive code
def safe_analysis(data):
    """Example of defensive programming"""

    # Input validation
    if data is None:
        raise ValueError("Data cannot be None")

    if len(data) == 0:
        raise ValueError("Data cannot be empty")

    # Check data types
    if not isinstance(data, pd.DataFrame):
        try:
            data = pd.DataFrame(data)
        except:
            raise TypeError("Data must be convertible to DataFrame")

    # Check for missing values
    if data.isnull().any().any():
        print("Warning: Missing values detected")
        data = data.dropna()

    # Proceed with analysis
    return analysis_result
```

#### Version Control
```python
# Track your environment
# Create requirements.txt
"""
pandas==1.3.3
numpy==1.21.2
scipy==1.7.1
matplotlib==3.4.3
seaborn==0.11.2
scanpy==1.8.1
"""

# Save environment
pip freeze > requirements.txt

# Recreate environment
pip install -r requirements.txt
```

#### Documentation Habits
```python
# Document your analysis
"""
# Analysis Log
Date: 2024-01-15
Objective: Differential expression analysis of tau-positive vs tau-negative neurons

## Data:
- File: pool_processed_v2.h5ad
- Samples: 150 (75 tau+, 75 tau-)
- Proteins: 5,853

## Methods:
- Statistical test: t-test with FDR correction
- Effect size: Cohen's d
- Significance threshold: FDR < 0.05

## Results:
- Significant proteins: 423 (7.2%)
- Top upregulated: SQSTM1 (FC=10.7)
- Top downregulated: VDAC1 (FC=0.3)

## Issues encountered:
- Initial memory error with full dataset
- Solved by processing in chunks

## Next steps:
- Pathway analysis
- Validation of top hits
"""
```

---

## 🚀 You've Got This!

### Remember: Every Expert Was Once a Beginner

Troubleshooting is a normal part of data analysis. Even experienced researchers encounter these problems regularly. The key is to:

1. **Stay calm** and work through problems systematically
2. **Read error messages carefully** - they often contain the solution
3. **Search online** - someone else has likely faced the same issue
4. **Ask for help** when you're stuck - the community is helpful
5. **Document solutions** for future reference

### Building Troubleshooting Skills

Each problem you solve makes you a better analyst:
- **Technical skills** improve with practice
- **Problem-solving ability** develops over time
- **Confidence** grows with each success
- **Research skills** become more efficient

### When to Take a Break

Sometimes the best solution is to step away:
- Take a coffee break and come back fresh
- Sleep on difficult problems
- Explain the problem to someone else (rubber duck debugging)
- Try a different approach entirely

---

**You have all the tools and knowledge needed to overcome any technical challenge. Trust the process, stay persistent, and remember that every problem has a solution!** 🔧✨

*Remember: The best programmers are not those who never encounter problems, but those who solve them efficiently!*