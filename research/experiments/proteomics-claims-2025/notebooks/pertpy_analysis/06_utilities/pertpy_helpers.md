# PertPy Helper Functions for DGE Analysis

## Overview
Utility functions for differential expression analysis using PertPy. These functions provide a robust framework for running PyDESeq2 with fallback methods, creating visualizations, and evaluating scientific claims.

## Main Functions

### 1. `run_pydeseq2_safe()`
Run PyDESeq2 with automatic fallback to traditional methods if it fails.

```python
def run_pydeseq2_safe(adata, design="~tau_status", contrast_column="tau_status",
                      baseline="negative", compare="positive")
```

**Parameters:**
- `adata`: AnnData object with expression data
- `design`: DESeq2 design formula
- `contrast_column`: Column for contrast comparison
- `baseline`: Reference group
- `compare`: Comparison group

**Returns:**
- DataFrame with log2FC, pvalue, padj columns

**Features:**
- Automatic handling of count vs log-transformed data
- Fallback to t-test and Mann-Whitney U if PyDESeq2 fails
- Cohen's d effect size calculation
- FDR correction using Benjamini-Hochberg

---

### 2. `find_proteins_in_dataset()`
Find proteins in dataset with fuzzy matching capabilities.

```python
def find_proteins_in_dataset(adata, protein_list, fuzzy_match=True)
```

**Parameters:**
- `adata`: AnnData object
- `protein_list`: List of protein names to find
- `fuzzy_match`: Enable partial string matching

**Returns:**
- Dictionary with 'found' and 'missing' protein lists

**Use Case:**
Essential for handling different protein naming conventions across datasets.

---

### 3. `create_volcano_plot()`
Generate publication-quality volcano plots from DGE results.

```python
def create_volcano_plot(results_df, title="Volcano Plot",
                       fc_threshold=0.5, pval_threshold=0.05,
                       label_top=10, save_path=None)
```

**Parameters:**
- `results_df`: DataFrame with DGE results
- `title`: Plot title
- `fc_threshold`: Fold change threshold for significance
- `pval_threshold`: P-value threshold
- `label_top`: Number of top proteins to label
- `save_path`: Optional path to save figure

**Returns:**
- matplotlib figure and axes objects

**Features:**
- Color coding by significance and effect size
- Automatic labeling of top hits
- Customizable thresholds
- High-resolution output (300 DPI)

---

### 4. `analyze_protein_set()`
Comprehensive analysis of a specific protein set.

```python
def analyze_protein_set(adata, protein_set, set_name="Protein Set",
                       design="~tau_status", **kwargs)
```

**Parameters:**
- `adata`: Full dataset
- `protein_set`: List of proteins to analyze
- `set_name`: Name for reporting
- `design`: DESeq2 formula
- `**kwargs`: Additional arguments for PyDESeq2

**Returns:**
- Dictionary with:
  - Number of proteins found/missing
  - Significance statistics
  - Mean log2 fold change
  - Full results DataFrame

---

### 5. `temporal_correlation_analysis()`
Analyze temporal correlations for proteins across pseudotime or disease progression.

```python
def temporal_correlation_analysis(adata, protein_list, time_column='pseudotime')
```

**Parameters:**
- `adata`: AnnData object
- `protein_list`: Proteins to analyze
- `time_column`: Column with temporal information

**Returns:**
- DataFrame with:
  - Spearman correlation coefficients
  - P-values with FDR correction
  - Linear regression slopes
  - Direction of change

**Applications:**
- Disease progression analysis
- Temporal dynamics studies
- Identifying early/late response proteins

---

### 6. `evaluate_claim()`
Systematically evaluate scientific claims against observed data.

```python
def evaluate_claim(observed_value, claimed_value, value_type='fold_change',
                  tolerance=0.2)
```

**Parameters:**
- `observed_value`: Value from your analysis
- `claimed_value`: Published/claimed value
- `value_type`: Type of measurement
- `tolerance`: Acceptable deviation (0.2 = 20%)

**Returns:**
- Tuple of (verdict, explanation)

**Verdicts:**
- **SUPPORTED**: Within tolerance
- **PARTIALLY SUPPORTED**: Within 2x tolerance
- **REFUTED**: Outside acceptable range
- **UNSURE**: Missing data

---

## Usage Examples

### Example 1: Basic DGE Analysis
```python
import scanpy as sc
from pertpy_helpers import run_pydeseq2_safe, create_volcano_plot

# Load data
adata = sc.read_h5ad('data.h5ad')

# Run DGE
results = run_pydeseq2_safe(adata, design="~condition")

# Create volcano plot
fig, ax = create_volcano_plot(results, title="My Analysis")
```

### Example 2: Protein Set Analysis
```python
from pertpy_helpers import analyze_protein_set

# Define UPS proteins
ups_proteins = ['PSMA1', 'PSMB5', 'UCHL1', 'USP14']

# Analyze
summary = analyze_protein_set(adata, ups_proteins, "UPS Proteins")
print(f"Significant: {summary['n_significant']}/{summary['n_proteins_found']}")
```

### Example 3: Temporal Analysis
```python
from pertpy_helpers import temporal_correlation_analysis

# Analyze mitochondrial proteins over pseudotime
mito_proteins = ['COX4I1', 'ATP5A1', 'VDAC1']
temporal_results = temporal_correlation_analysis(adata, mito_proteins)

# Get significantly changing proteins
significant = temporal_results[temporal_results['significant']]
```

### Example 4: Claim Evaluation
```python
from pertpy_helpers import evaluate_claim

# Test a fold change claim
observed_fc = 1.32
claimed_fc = 3.413
verdict, explanation = evaluate_claim(observed_fc, claimed_fc, 'fold_change')
print(f"Verdict: {verdict}")
print(f"Reason: {explanation}")
```

---

## Error Handling

The functions include robust error handling:

1. **PyDESeq2 Failures**: Automatic fallback to traditional statistics
2. **Missing Proteins**: Clear reporting of found vs missing
3. **NaN Values**: Proper handling in all calculations
4. **Small Sample Sizes**: Warnings when n < 10

---

## Statistical Methods

### Multiple Testing Correction
- Benjamini-Hochberg FDR applied consistently
- Default threshold: FDR < 0.05

### Effect Size Metrics
- Cohen's d for standardized effect size
- Log2 fold change for expression differences

### Correlation Methods
- Spearman correlation (non-parametric)
- Linear regression for trend analysis

---

## Best Practices

1. **Always check protein coverage**:
   ```python
   protein_info = find_proteins_in_dataset(adata, my_proteins)
   print(f"Found: {len(protein_info['found'])}")
   print(f"Missing: {protein_info['missing']}")
   ```

2. **Use appropriate design formulas**:
   - Simple: `"~condition"`
   - With covariates: `"~condition + age + sex"`
   - With interaction: `"~condition * time"`

3. **Save results systematically**:
   ```python
   results.to_csv('results/dge_analysis.csv')
   fig.savefig('figures/volcano_plot.png', dpi=300)
   ```

4. **Document parameters**:
   - Record all thresholds used
   - Note any proteins excluded
   - Save sessionInfo equivalent

---

## Dependencies

Required packages:
- pertpy >= 0.5.0
- scanpy >= 1.9.0
- pandas >= 1.5.0
- numpy >= 1.23.0
- scipy >= 1.9.0
- statsmodels >= 0.13.0
- matplotlib >= 3.5.0
- seaborn >= 0.12.0
- scikit-learn >= 1.1.0

---

## Citation

If using these utilities in publication, cite:
- PertPy paper (when available)
- PyDESeq2 implementation
- Original DESeq2 method (Love et al., 2014)