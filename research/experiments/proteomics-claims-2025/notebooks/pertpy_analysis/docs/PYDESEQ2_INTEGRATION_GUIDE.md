# 🧬 PyDESeq2 Integration Guide for PertPy DGE Analysis

## 📊 Current Approach vs PyDESeq2

### Current Implementation (Simplified T-tests)
```python
# Current approach in your notebooks
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

# Simple differential expression
tau_pos_expr = adata.X[tau_positive_mask, :]
tau_neg_expr = adata.X[tau_negative_mask, :]

# T-test for each protein
t_stat, p_val = stats.ttest_ind(tau_pos_expr, tau_neg_expr)

# FDR correction
p_adjusted = fdrcorrection(p_values)[1]
```

**Pros:**
- ✅ Simple and fast
- ✅ Works with any data type
- ✅ Minimal dependencies
- ✅ Easy to understand

**Cons:**
- ❌ Assumes normal distribution
- ❌ Doesn't account for count-based nature of proteomics data
- ❌ No dispersion modeling
- ❌ Less sophisticated normalization

### PyDESeq2 Approach (Advanced DEG Analysis)
```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# Create DESeq2 dataset
dds = DeseqDataSet(
    adata=adata,
    design_factors="tau_status",
    ref_level=["tau_status", "negative"],
    refit_cooks=True
)

# Run DESeq2 pipeline
dds.deseq2()

# Statistical inference
inference = DefaultInference(n_cpus=8)
dds.calculate_qvalues(inference)

# Get results
stat_res = DeseqStats(dds, inference=inference)
stat_res.summary()
```

**Pros:**
- ✅ Gold standard for differential expression
- ✅ Models count overdispersion
- ✅ Sophisticated normalization (size factors)
- ✅ Log fold change shrinkage
- ✅ Better handling of low-count proteins
- ✅ Robust outlier detection

**Cons:**
- ❌ More complex setup
- ❌ Requires count data (not intensities)
- ❌ Longer computation time
- ❌ Steeper learning curve

---

## 🔄 How PyDESeq2 Fits Your Analysis

### 1. **Enhanced Statistical Power**
PyDESeq2 would improve your analyses by:
- **Better handling of biological variability** through dispersion modeling
- **More accurate p-values** for proteins with low expression
- **Reduced false positives** through independent filtering
- **Shrinkage estimators** for more reliable fold changes

### 2. **Integration Points**

#### Option A: Full PyDESeq2 Replacement
Replace current t-test approach entirely:

```python
# Modified notebook section
def run_pydeseq2_analysis(adata, protein_list):
    """Run PyDESeq2 differential expression analysis"""

    # Subset to proteins of interest
    adata_subset = adata[:, protein_list].copy()

    # Initialize PyDESeq2
    dds = DeseqDataSet(
        adata=adata_subset,
        design_factors="tau_status",
        ref_level=["tau_status", "negative"]
    )

    # Run analysis
    dds.deseq2()

    # Get results
    stat_res = DeseqStats(dds)
    results_df = stat_res.summary()

    return results_df
```

#### Option B: Hybrid Approach
Use PyDESeq2 for discovery, t-tests for quick validation:

```python
# Discovery phase with PyDESeq2
significant_proteins = run_pydeseq2_discovery(adata)

# Quick validation with t-tests
for protein in significant_proteins:
    quick_ttest_validation(adata, protein)
```

#### Option C: Comparative Analysis
Run both methods and compare:

```python
# Run both analyses
ttest_results = run_simple_ttest(adata, proteins)
deseq2_results = run_pydeseq2_analysis(adata, proteins)

# Compare results
concordance = compare_methods(ttest_results, deseq2_results)
```

---

## 📝 Recommended Integration Strategy

### For Your Current Notebooks:

1. **Keep simplified version as default** (current implementation)
   - Maintains accessibility for users
   - Works with any data type
   - Fast execution

2. **Add PyDESeq2 as advanced option**:

```python
# In each notebook, add advanced analysis section
## Advanced Analysis with PyDESeq2 (Optional)

```python
# For users with count-based proteomics data
USE_PYDESEQ2 = False  # Set to True for advanced analysis

if USE_PYDESEQ2:
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        print("🔬 Running advanced PyDESeq2 analysis...")

        # PyDESeq2 pipeline
        dds = DeseqDataSet(
            adata=adata_proteins,
            design_factors="tau_status",
            ref_level=["tau_status", "negative"]
        )

        dds.deseq2()
        stat_res = DeseqStats(dds)
        results = stat_res.summary()

        # Extract key metrics
        significant = results[results['padj'] < 0.05]
        print(f"PyDESeq2 found {len(significant)} significant proteins")

    except ImportError:
        print("⚠️ PyDESeq2 not installed, using simple t-tests")
        USE_PYDESEQ2 = False
```

---

## 🎯 Specific Benefits for Your Claims

### Claims That Would Benefit Most from PyDESeq2:

1. **"SQSTM1/p62 is upregulated"** - Better fold change estimation
2. **"Mitochondrial complexes I-V decreased"** - Handle low-abundance subunits
3. **"Temporal dynamics"** - Model time as continuous covariate
4. **"V-ATPase dysregulation"** - Account for subunit correlations
5. **"Retromer complex"** - Multi-protein complex analysis

### PyDESeq2 Advanced Features for Your Analysis:

```python
# 1. Time-series analysis for temporal claims
dds = DeseqDataSet(
    adata=adata,
    design_factors=["tau_status", "pseudotime"],
    continuous_factors=["pseudotime"]
)

# 2. Interaction effects for complex claims
dds = DeseqDataSet(
    adata=adata,
    design_factors="tau_status + treatment + tau_status:treatment"
)

# 3. Batch correction
dds = DeseqDataSet(
    adata=adata,
    design_factors="batch + tau_status"
)
```

---

## 💻 Implementation Example

### Complete PyDESeq2 Integration for Claim 2:

```python
# claim2_sqstm1_upregulation_colab.md - Enhanced Version

## Advanced PyDESeq2 Analysis

```python
# Install PyDESeq2
!pip install -q pydeseq2

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pandas as pd

# Prepare data for PyDESeq2
# Ensure counts are integers (required by DESeq2)
if not np.issubdtype(adata.X.dtype, np.integer):
    print("Converting to counts...")
    adata.X = np.round(adata.X).astype(int)

# Focus on SQSTM1 and related proteins
autophagy_proteins = ['SQSTM1', 'NBR1', 'OPTN', 'NDP52', 'TAX1BP1']
available = [p for p in autophagy_proteins if p in adata.var_names]

# Create DESeq2 dataset
dds = DeseqDataSet(
    adata=adata[:, available],
    design_factors="tau_status",
    ref_level=["tau_status", "negative"],
    refit_cooks=True,
    n_cpus=4  # Use parallel processing
)

# Run DESeq2
print("Running DESeq2 analysis...")
dds.deseq2()

# Get statistical results
stat_res = DeseqStats(dds)
stat_res.summary()

# Extract results
results_df = stat_res.results_df
results_df = results_df.sort_values('padj')

# Check SQSTM1 specifically
if 'SQSTM1' in results_df.index:
    sqstm1_result = results_df.loc['SQSTM1']

    print(f"\n🎯 SQSTM1/p62 PyDESeq2 Results:")
    print(f"  Base mean: {sqstm1_result['baseMean']:.2f}")
    print(f"  Log2 fold change: {sqstm1_result['log2FoldChange']:.3f}")
    print(f"  Adjusted p-value: {sqstm1_result['padj']:.4f}")

    # Verdict based on PyDESeq2
    if sqstm1_result['padj'] < 0.05 and sqstm1_result['log2FoldChange'] > 0.5:
        verdict_deseq2 = "✅ STRONGLY SUPPORTED"
    elif sqstm1_result['padj'] < 0.05 and sqstm1_result['log2FoldChange'] > 0:
        verdict_deseq2 = "✅ SUPPORTED"
    else:
        verdict_deseq2 = "❌ NOT SUPPORTED"

    print(f"\n  PyDESeq2 Verdict: {verdict_deseq2}")

# Volcano plot with PyDESeq2 results
plt.figure(figsize=(10, 6))
plt.scatter(results_df['log2FoldChange'],
           -np.log10(results_df['padj']),
           alpha=0.6)

# Highlight SQSTM1
if 'SQSTM1' in results_df.index:
    sqstm1_data = results_df.loc['SQSTM1']
    plt.scatter(sqstm1_data['log2FoldChange'],
               -np.log10(sqstm1_data['padj']),
               color='red', s=100, label='SQSTM1/p62')

plt.xlabel('Log2 Fold Change (PyDESeq2)')
plt.ylabel('-Log10(adjusted p-value)')
plt.title('PyDESeq2 Differential Expression Analysis')
plt.axhline(y=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5)
plt.axvline(x=0, color='gray', linestyle='-', alpha=0.5)
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

---

## 📊 When to Use Each Method

### Use Simple T-tests When:
- Quick exploratory analysis needed
- Data is not count-based (e.g., normalized intensities)
- Sample size is small (n < 3 per group)
- Teaching/demonstration purposes
- Computational resources limited

### Use PyDESeq2 When:
- Publication-quality results required
- Working with count-based proteomics/RNA-seq
- Complex experimental designs (multiple factors)
- Time-series or dose-response analysis
- Need robust handling of outliers
- Want log fold change shrinkage

---

## 🚀 Migration Path

### Phase 1: Add Optional PyDESeq2 (Current)
- Keep existing t-test approach
- Add PyDESeq2 as optional advanced analysis
- Document when to use each method

### Phase 2: Parallel Analysis
- Run both methods by default
- Compare and report concordance
- Build confidence in PyDESeq2 results

### Phase 3: PyDESeq2 Primary (Future)
- Make PyDESeq2 the default for appropriate data
- Keep t-tests as fallback/validation
- Full integration with PertPy framework

---

## 📚 Resources

- [PyDESeq2 Documentation](https://pydeseq2.readthedocs.io/)
- [Original DESeq2 Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
- [PertPy Integration](https://pertpy.readthedocs.io/)

---

## 💡 Recommendation

**For your current analysis:** The simplified t-test approach is appropriate for initial exploration and demonstration. However, consider adding PyDESeq2 as an optional advanced analysis section in each notebook for users who want more sophisticated statistical modeling.

**Key advantage of PyDESeq2 for your work:** It would provide more reliable p-values and fold changes, especially for low-abundance proteins and complex experimental designs, making your biological conclusions more robust.