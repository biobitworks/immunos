# PertPy DGE Analysis - Complete Notebook Index

## 📚 Available Formats

Each analysis is available in **two formats** for maximum flexibility:
- **`.ipynb`** - Jupyter notebook for direct execution in Colab
- **`.md`** - Markdown with complete code for easy copy-paste

## 🗂️ Complete File List

### 📁 01_data_preparation/
| File | Description | Proteins |
|------|-------------|----------|
| `prepare_for_pertpy_colab.ipynb` | Colab notebook for data prep | 500+ embedded |
| **`prepare_for_pertpy_colab.md`** | **Markdown version with all code** | **Complete copy-paste ready** |

### 📁 02_group1_mitochondrial/
| File | Description | Key Features |
|------|-------------|--------------|
| `claim1_ups_proteins_colab.ipynb` | UPS protein analysis | 132 proteins embedded |
| **`claim1_ups_proteins_colab.md`** | **Markdown with all code** | **Complete UPS list included** |
| `claim2_sqstm1_upregulation_colab.ipynb` | SQSTM1/p62 analysis | Auto-detection code |
| **`claim2_sqstm1_upregulation_colab.md`** | **Markdown with all code** | **Pseudotime regression** |
| `claim3_temporal_dynamics_colab.ipynb` | Temporal progression | 200+ proteins |
| **`claim3_temporal_dynamics_colab.md`** | **Markdown with all code** | **Phase detection analysis** |

### 📁 03_group2_proteostasis/
| File | Description | Key Features |
|------|-------------|--------------|
| `claim1_vatpase_subunits_colab.ipynb` | V-ATPase analysis | 24 subunits embedded |
| **`claim1_vatpase_subunits_colab.md`** | **Markdown with all code** | **V0/V1 domain analysis** |

## 🚀 Quick Start Guide

### Option 1: Use Jupyter Notebooks in Colab
1. Upload any `.ipynb` file to [Google Colab](https://colab.research.google.com)
2. Run first cell to install packages
3. Upload `pool_processed_v2.h5ad` when prompted
4. Execute remaining cells

### Option 2: Copy from Markdown Files
1. Open the `.md` version of desired analysis
2. Copy code cells as needed
3. Paste into your environment (Colab, Jupyter, IDE)
4. All protein lists are embedded in the code

## 📋 Complete Protein Coverage

### UPS Proteins (132 total)
```python
# Embedded in claim1_ups_proteins_colab.md
- Proteasome subunits: 43
- E1 enzymes: 7
- E2 enzymes: 18
- E3 ligases: 19
- Deubiquitinases: 28
- Regulators: 9
- Modifiers: 8
```

### Mitochondrial Complexes
```python
# Embedded in claim3_temporal_dynamics_colab.md
- Complex I: 37 proteins
- Complex II: 4 proteins
- Complex III: 10 proteins
- Complex IV: 16 proteins
- Complex V: 24 proteins
```

### V-ATPase Subunits (24 total)
```python
# Embedded in claim1_vatpase_subunits_colab.md
- V0 domain: 11 subunits
- V1 domain: 13 subunits
```

### Autophagy & Proteostasis
```python
# Embedded across multiple notebooks
- Core autophagy: 24 proteins
- Autophagy receptors: 6 proteins
- LC3/GABARAP family: 7 proteins
- Lysosomal markers: 14 proteins
- Heat shock proteins: 20 proteins
```

## 💡 Key Advantages of Markdown Files

1. **Complete Code Visibility**: See all code without running
2. **Easy Copy-Paste**: Copy specific sections as needed
3. **Version Control Friendly**: Track changes in git
4. **Platform Independent**: View in any text editor
5. **Search Friendly**: Find proteins/functions quickly
6. **Documentation**: Code and explanations together

## 📊 Analysis Results Summary

| Claim | Verdict | Key Finding | Notebook |
|-------|---------|-------------|----------|
| UPS preservation | **REFUTED** | 28.8% proteins altered | claim1_ups_proteins |
| SQSTM1 upregulation | **PARTIALLY SUPPORTED** | 1.32x vs 3.4x claimed | claim2_sqstm1 |
| Temporal dynamics | **SUPPORTED** | Progressive decline | claim3_temporal |
| V-ATPase differential | **SUPPORTED** | 44% subunits changed | claim1_vatpase |

## 🔧 Technical Details

### All Notebooks Include:
- ✅ Automatic Colab detection
- ✅ Package installation code
- ✅ File upload handling
- ✅ Data standardization
- ✅ PyDESeq2 with fallback
- ✅ FDR correction
- ✅ Comprehensive visualizations
- ✅ Statistical summaries

### No External Dependencies:
- ❌ No JSON files needed
- ❌ No protein list files
- ❌ No configuration files
- ✅ Just upload `pool_processed_v2.h5ad`

## 📝 Example: Using Markdown Files

```python
# Example: Extract UPS protein list from markdown

# 1. Open claim1_ups_proteins_colab.md
# 2. Find the protein definition section
# 3. Copy the entire ups_proteins_comprehensive list
# 4. Paste into your analysis

ups_proteins_comprehensive = [
    # ... complete 132 protein list from markdown ...
]
```

## 🎯 Recommended Workflow

1. **Start with data prep**: `prepare_for_pertpy_colab.md`
2. **Choose your analysis**: Select relevant claim notebook
3. **Copy code sections**: Use markdown for easy extraction
4. **Customize as needed**: Modify protein lists or parameters
5. **Run in your environment**: Colab, local Jupyter, or IDE

## 📚 Additional Resources

- [README_COLAB.md](README_COLAB.md) - Detailed Colab instructions
- [all_claims_summary.md](../04_comprehensive_analysis/all_claims_summary.md) - Results overview
- [Statistical methods reference](../06_utilities/statistical_methods.md)

---

**All notebooks are self-contained with embedded protein sets - no external files required!**

*Last updated: December 2024*