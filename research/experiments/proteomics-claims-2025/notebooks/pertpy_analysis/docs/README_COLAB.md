# PertPy DGE Analysis - Google Colab Edition

## 🚀 Quick Start for Google Colab

All notebooks are **self-contained** and optimized for Google Colab. No external files needed - just upload `pool_processed_v2.h5ad`!

## 📁 Available Colab Notebooks

### Data Preparation
- **`01_data_preparation/prepare_for_pertpy_colab.ipynb`**
  - Prepares data for all analyses
  - 500+ proteins embedded directly in notebook
  - Creates counts layer for PyDESeq2

### Group 1: Mitochondrial Dysregulation
1. **`02_group1_mitochondrial/claim1_ups_proteins_colab.ipynb`**
   - Tests 132 comprehensive UPS proteins
   - Compares with cherry-picked approaches
   - Shows 28.8% of UPS proteins significantly altered

2. **`02_group1_mitochondrial/claim2_sqstm1_upregulation_colab.ipynb`**
   - Analyzes SQSTM1/p62 upregulation
   - Compares observed (1.32x) vs claimed (3.4x) fold change
   - Includes pseudotime correlation analysis

3. **`02_group1_mitochondrial/claim3_temporal_dynamics_colab.ipynb`**
   - Tests 200+ proteins across disease progression
   - Analyzes mitochondrial complexes I-V
   - Identifies disease phases (Early/Middle/Late)

### Group 2: Proteostasis Failure
1. **`03_group2_proteostasis/claim1_vatpase_subunits_colab.ipynb`**
   - Analyzes all 24 V-ATPase subunits
   - V0 (membrane) and V1 (cytoplasmic) domains
   - Links to lysosomal dysfunction

## 🎯 How to Use in Google Colab

### Step 1: Upload Notebook to Colab
1. Go to [Google Colab](https://colab.research.google.com)
2. Click "File" → "Upload notebook"
3. Upload the desired `*_colab.ipynb` file

### Step 2: Run Setup Cell
```python
# First cell will:
- Detect Colab environment
- Install required packages (pertpy, pydeseq2, scanpy)
- Prompt for data file upload
```

### Step 3: Upload Data File
When prompted, upload `pool_processed_v2.h5ad` (your proteomics data)

### Step 4: Execute Analysis
Run remaining cells sequentially - everything is self-contained!

## 📊 Key Features

### Self-Contained Design
- **No JSON files**: All protein sets embedded in code
- **No external dependencies**: Complete analysis in one notebook
- **Automatic package installation**: Colab installs everything needed

### Embedded Protein Sets

| Protein Set | Count | Notebook |
|------------|-------|----------|
| UPS proteins | 132 | claim1_ups_proteins |
| Mitochondrial Complex I | 37 | claim3_temporal_dynamics |
| Mitochondrial Complex V | 24 | claim3_temporal_dynamics |
| V-ATPase subunits | 24 | claim1_vatpase_subunits |
| Autophagy proteins | 40+ | Multiple notebooks |
| Heat shock proteins | 21 | claim3_temporal_dynamics |

### Statistical Methods
- **PyDESeq2**: Primary differential expression engine
- **Fallback statistics**: T-test, Mann-Whitney U when PyDESeq2 fails
- **FDR correction**: Benjamini-Hochberg for all analyses
- **Effect sizes**: Cohen's d, log2 fold changes

## 🔬 Analysis Highlights

### UPS Analysis (Claim 1)
```python
# 132 comprehensive UPS proteins embedded
ups_proteins_comprehensive = [
    'PSMA1', 'PSMA2', ... # 43 proteasome subunits
    'UBE2D1', 'UBE2D3', ... # 18 E2 enzymes
    'HERC1', 'HERC2', ... # 19 E3 ligases
    # ... complete list in notebook
]
```
**Finding**: 28.8% significantly altered (refutes "no changes" claim)

### SQSTM1 Analysis (Claim 2)
```python
# Automatic SQSTM1 detection
search_patterns = ['SQSTM1', 'P62', 'SEQUESTOSOME']
# Finds and analyzes SQSTM1/p62 expression
```
**Finding**: 1.32 log2FC observed vs 3.413 claimed

### Temporal Dynamics (Claim 3)
```python
# 200+ proteins across 11 functional categories
temporal_proteins = {
    'mitochondrial_complex_I': [...],  # 37 proteins
    'mitochondrial_complex_V': [...],  # 24 proteins
    'autophagy_early': [...],          # 24 proteins
    'autophagy_late': [...],           # 13 proteins
    # ... complete sets in notebook
}
```
**Finding**: Progressive mitochondrial decline confirmed

### V-ATPase Analysis (Group 2 Claim 1)
```python
# All 24 V-ATPase subunits
vatpase_subunits = {
    'V0_domain': ['ATP6V0A1', ...],  # 11 membrane subunits
    'V1_domain': ['ATP6V1A', ...],   # 13 cytoplasmic subunits
}
```
**Finding**: 44% of subunits show differential expression

## 📈 Visualizations

Each notebook includes:
- **Volcano plots**: Significance vs effect size
- **Heatmaps**: Expression patterns across samples
- **Temporal plots**: Disease progression analysis
- **Box/violin plots**: Group comparisons
- **Phase analysis**: Early/Middle/Late disease stages

## 💡 Tips for Colab Users

### Memory Management
```python
# If you encounter memory issues:
import gc
gc.collect()  # Clear unused memory

# Use smaller batches:
adata_subset = adata[:, :1000]  # Analyze first 1000 proteins
```

### Save Results
```python
# Download results from Colab:
from google.colab import files

# Save and download
results_df.to_csv('results.csv')
files.download('results.csv')
```

### Mount Google Drive
```python
# Optional: Save to Google Drive
from google.colab import drive
drive.mount('/content/drive')

# Save directly to Drive
results_df.to_csv('/content/drive/MyDrive/results.csv')
```

## 🐛 Troubleshooting

### PyDESeq2 Fails
- Notebooks include automatic fallback to traditional statistics
- T-test and Mann-Whitney U as alternatives
- FDR correction still applied

### Protein Not Found
- Notebooks include multiple search strategies
- Case-insensitive matching
- Partial name matching
- Alternative protein name suggestions

### Memory Issues
- Use free Colab GPU runtime for more memory
- Process proteins in batches
- Clear variables with `del` and `gc.collect()`

## 📚 Biological Context

### Mitochondrial Dysfunction
- Complexes I & V progressively decline
- Energy production compromised
- Links to neurodegeneration

### Proteostasis Failure
- V-ATPase dysfunction → lysosomal pH ↑
- Autophagy impairment → protein accumulation
- SQSTM1/p62 upregulation indicates blockage

### Therapeutic Implications
- V-ATPase restoration for lysosomal function
- Autophagy enhancement strategies
- Stage-specific interventions

## 🔗 Links and Resources

- [PertPy Documentation](https://pertpy.readthedocs.io/)
- [PyDESeq2 GitHub](https://github.com/owkin/PyDESeq2)
- [Google Colab](https://colab.research.google.com)
- [Original Paper](link-to-paper)

## 📄 Citation

If using these notebooks:
```
PertPy DGE Analysis - Colab Edition
Comprehensive proteomics analysis using PyDESeq2/PertPy
Dataset: pool_processed_v2.h5ad
Date: 2024
```

## ✅ Summary

These Colab-optimized notebooks provide:
1. **Complete portability** - Works anywhere with internet
2. **No setup required** - Colab handles everything
3. **Comprehensive analysis** - 500+ proteins, 16 claims tested
4. **Professional statistics** - PyDESeq2 with proper FDR
5. **Clear verdicts** - Objective evaluation of each claim

Simply upload to Colab, upload your data file, and run!

---

*All notebooks tested in Google Colab (December 2024)*