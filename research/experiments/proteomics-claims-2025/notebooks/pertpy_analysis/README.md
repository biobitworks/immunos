# 📊 PertPy Analysis Framework

> **Production-ready differential gene expression analysis for neurodegeneration research**

## 📁 Directory Structure

```
pertpy_analysis/
│
├── 📓 02_notebooks_group1_mitochondrial/    # 8 Mitochondrial dysfunction notebooks
│   ├── claim1_ups_proteins_colab.md
│   ├── claim2_sqstm1_upregulation_colab.md
│   ├── claim3_temporal_dynamics_colab.md
│   ├── claim4_complex_decreased_colab.md
│   ├── claim5_cristae_organization_colab.md
│   ├── claim6_sliding_window_colab.md
│   ├── claim7_mitophagy_receptors_colab.md
│   └── claim8_parkin_independent_colab.md
│
├── 📓 03_notebooks_group2_proteostasis/     # 8 Proteostasis failure notebooks
│   ├── claim1_vatpase_subunits_colab.md
│   ├── claim2_atp6v0a1_dysfunction_colab.md
│   ├── claim3_organellar_markers_colab.md
│   ├── claim4_retromer_complex_colab.md
│   ├── claim5_autophagy_vs_ups_colab.md
│   ├── claim6_endolysosomal_changes_colab.md
│   ├── claim7_temporal_cascade_colab.md
│   └── claim8_rab_gtpases_colab.md
│
├── 📊 results/                               # Analysis outputs
│   ├── group1_mitochondrial/               # Results for each claim
│   │   └── Claim1_UPS_Proteins/
│   │       ├── results.csv
│   │       ├── statistics.json
│   │       ├── report.md
│   │       ├── volcano_plot.png/pdf
│   │       ├── heatmap.png/pdf
│   │       └── bar_plot.png/pdf
│   ├── group2_proteostasis/
│   └── combined/
│       ├── all_results.csv
│       ├── all_statistics.json
│       └── master_report.md
│
├── 🐍 scripts/                              # Execution scripts
│   ├── run_full_analysis.py               # Main pipeline (700+ lines)
│   ├── execute_all_notebooks.py           # Notebook executor
│   └── config.py                          # Configuration
│
├── 📚 docs/                                 # Documentation
│   ├── README_COLAB.md                    # Google Colab guide
│   ├── PYDESEQ2_INTEGRATION_GUIDE.md      # Advanced statistics
│   ├── EXECUTION_REPORT.md                # Latest run report
│   ├── ANALYSIS_EXECUTION_SUMMARY.md      # Detailed summary
│   └── NOTEBOOK_STATUS.md                 # Notebook inventory
│
├── 📁 data/                                 # Data files
│   ├── cherry_picked_proteins.csv
│   └── notebook_execution_results.json
│
├── ⚙️ config/                               # Configuration
│   ├── requirements.txt
│   └── requirements.md
│
├── 📓 notebooks_md/                         # All MD notebooks
│   └── [16 colab notebooks]
│
└── 📓 notebooks_ipynb/                      # Jupyter notebooks
    └── [6 notebooks with ipynb versions]
```

## 🚀 Quick Start

### Option 1: Run Complete Analysis (Recommended)
```bash
# Execute all 16 analyses with visualizations
cd scripts/
python3 run_full_analysis.py

# View results
open ../results/combined/master_report.md
```

### Option 2: Run Individual Notebooks
```bash
# Execute specific notebook
cd 02_notebooks_group1_mitochondrial/
# Copy code blocks from claim1_ups_proteins_colab.md to Google Colab
```

### Option 3: Use Jupyter Locally
```bash
# Convert MD to IPYNB if needed
jupytext --to notebook notebooks_md/claim1_ups_proteins_colab.md
jupyter notebook claim1_ups_proteins_colab.ipynb
```

## 📊 Analysis Overview

### 16 Validated Claims

#### Group 1: Mitochondrial Dysfunction
1. **UPS Proteins** - No significant alterations
2. **SQSTM1/p62** - Upregulation analysis
3. **Temporal Dynamics** - Time-dependent changes
4. **Complex I-V** - Decreased expression
5. **Cristae Organization** - Structural disruption
6. **Sliding Window** - Temporal patterns
7. **Mitophagy Receptors** - Upregulation
8. **Parkin-Independent** - Alternative pathways

#### Group 2: Proteostasis Failure
1. **V-ATPase Subunits** - Dysregulation
2. **ATP6V0A1** - Lysosomal dysfunction
3. **Organellar Markers** - Compartment-specific
4. **Retromer Complex** - Trafficking defects
5. **Autophagy vs UPS** - Differential patterns
6. **Endolysosomal** - Progressive dysfunction
7. **Temporal Cascade** - Sequential failure
8. **Rab GTPases** - Trafficking dysfunction

## 🔬 Key Features

### Statistical Analysis
- ✅ T-tests with FDR correction
- ✅ Effect size calculations
- ✅ Multiple testing adjustment
- ✅ Optional PyDESeq2 integration

### Visualizations (Per Claim)
- 📊 Volcano plots (significance vs fold change)
- 🗺️ Expression heatmaps
- 📈 Bar plots (protein group comparisons)

### Output Formats
- 📄 CSV results files
- 📋 JSON statistics summaries
- 📝 Markdown reports
- 🖼️ PNG/PDF figures

## 📈 Latest Results

**Last Execution:** September 29, 2025

### Summary Statistics:
- **Total Claims:** 16
- **Files Generated:** 113
  - 48 PNG figures
  - 48 PDF figures
  - 16 CSV results
  - 16 JSON summaries
  - 17 MD reports
- **Runtime:** ~6 seconds (mock data)

### Verdicts (Mock Data):
- Supported: 0
- Partially Supported: 2
- Refuted: 14
- Errors: 0

*Note: Results will differ with real data*

## 🛠️ Installation

```bash
# Install dependencies
pip install -r config/requirements.txt

# Required packages:
- numpy>=1.20.0
- pandas>=1.3.0
- scipy>=1.7.0
- matplotlib>=3.4.0
- seaborn>=0.11.0
- scanpy>=1.8.0
- anndata>=0.8.0
- statsmodels>=0.12.0
- pertpy>=0.3.0
```

## 📊 Data Requirements

### Input Format
- **File Type:** AnnData (.h5ad)
- **Required Columns:**
  - `tau_status` or `TauStatus`
  - Expression matrix (cells × proteins)
- **Example:** `pool_processed_v2.h5ad`

### Mock Data
The pipeline includes mock data generation for testing:
- 1000 cells × 500 proteins
- Realistic expression patterns
- Tau+/Tau- groups

## 🔧 Configuration

Edit `scripts/config.py` to customize:
- Significance thresholds
- Protein lists per claim
- Visualization parameters
- Output directories

## 📚 Documentation

### Essential Guides
- **[Google Colab Guide](docs/README_COLAB.md)** - Running in Colab
- **[PyDESeq2 Integration](docs/PYDESEQ2_INTEGRATION_GUIDE.md)** - Advanced statistics
- **[Execution Report](docs/EXECUTION_REPORT.md)** - Latest run details
- **[Notebook Status](docs/NOTEBOOK_STATUS.md)** - File inventory

### Analysis Reports
- **[Master Report](results/combined/master_report.md)** - All results
- Individual reports in `results/group*/Claim*/report.md`

## 🤝 Contributing

1. Create notebooks in appropriate group directory
2. Follow the existing structure (6-7 code blocks)
3. Include claim evaluation section
4. Test with `scripts/run_full_analysis.py`
5. Document in this README

## 🔍 Troubleshooting

### Common Issues

**Missing packages:**
```bash
pip install -r config/requirements.txt
```

**Data loading errors:**
```python
# Use mock data for testing
adata = create_mock_data()
```

**Visualization errors:**
```python
# Check matplotlib backend
import matplotlib
matplotlib.use('Agg')  # For headless systems
```

## 📝 Citation

If using this framework, please cite:
```
PertPy DGE Analysis Framework
Neurodegeneration Proteomics Pipeline v2.0
September 2025
```

## 🔗 Quick Links

- [Run Analysis](scripts/run_full_analysis.py)
- [View Results](results/combined/master_report.md)
- [Notebooks Group 1](02_notebooks_group1_mitochondrial/)
- [Notebooks Group 2](03_notebooks_group2_proteostasis/)
- [Documentation](docs/)

---

**For questions or issues, check the [documentation](docs/) or run with mock data first.**