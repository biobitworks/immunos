# 📑 PertPy Analysis - File Index

## 🚀 Start Here

### For Analysis
1. **[Run Full Analysis](scripts/run_full_analysis.py)** - Execute all 16 claims
2. **[View Latest Results](results/combined/master_report.md)** - Summary of last run
3. **[Check Individual Results](results/)** - Detailed per-claim outputs

### For Development
1. **[Main README](README.md)** - Complete documentation
2. **[Notebooks Group 1](02_notebooks_group1_mitochondrial/)** - Mitochondrial claims
3. **[Notebooks Group 2](03_notebooks_group2_proteostasis/)** - Proteostasis claims

## 📊 Analysis Notebooks (16 Claims)

### Group 1: Mitochondrial Dysfunction
| # | Claim | Notebook | Results |
|---|-------|----------|---------|
| 1 | UPS Proteins | [claim1_ups_proteins_colab.md](02_notebooks_group1_mitochondrial/claim1_ups_proteins_colab.md) | [Results](results/group1_mitochondrial/Claim1_UPS_Proteins/) |
| 2 | SQSTM1 Upregulation | [claim2_sqstm1_upregulation_colab.md](02_notebooks_group1_mitochondrial/claim2_sqstm1_upregulation_colab.md) | [Results](results/group1_mitochondrial/Claim2_SQSTM1_Upregulation/) |
| 3 | Temporal Dynamics | [claim3_temporal_dynamics_colab.md](02_notebooks_group1_mitochondrial/claim3_temporal_dynamics_colab.md) | [Results](results/group1_mitochondrial/Claim3_Temporal_Dynamics/) |
| 4 | Complex Decreased | [claim4_complex_decreased_colab.md](02_notebooks_group1_mitochondrial/claim4_complex_decreased_colab.md) | [Results](results/group1_mitochondrial/Claim4_Complex_Decreased/) |
| 5 | Cristae Organization | [claim5_cristae_organization_colab.md](02_notebooks_group1_mitochondrial/claim5_cristae_organization_colab.md) | [Results](results/group1_mitochondrial/Claim5_Cristae_Organization/) |
| 6 | Sliding Window | [claim6_sliding_window_colab.md](02_notebooks_group1_mitochondrial/claim6_sliding_window_colab.md) | [Results](results/group1_mitochondrial/Claim6_Sliding_Window/) |
| 7 | Mitophagy Receptors | [claim7_mitophagy_receptors_colab.md](02_notebooks_group1_mitochondrial/claim7_mitophagy_receptors_colab.md) | [Results](results/group1_mitochondrial/Claim7_Mitophagy_Receptors/) |
| 8 | Parkin-Independent | [claim8_parkin_independent_colab.md](02_notebooks_group1_mitochondrial/claim8_parkin_independent_colab.md) | [Results](results/group1_mitochondrial/Claim8_Parkin_Independent/) |

### Group 2: Proteostasis Failure
| # | Claim | Notebook | Results |
|---|-------|----------|---------|
| 1 | V-ATPase Subunits | [claim1_vatpase_subunits_colab.md](03_notebooks_group2_proteostasis/claim1_vatpase_subunits_colab.md) | [Results](results/group2_proteostasis/Claim1_VATPase_Subunits/) |
| 2 | ATP6V0A1 Dysfunction | [claim2_atp6v0a1_dysfunction_colab.md](03_notebooks_group2_proteostasis/claim2_atp6v0a1_dysfunction_colab.md) | [Results](results/group2_proteostasis/Claim2_ATP6V0A1_Dysfunction/) |
| 3 | Organellar Markers | [claim3_organellar_markers_colab.md](03_notebooks_group2_proteostasis/claim3_organellar_markers_colab.md) | [Results](results/group2_proteostasis/Claim3_Organellar_Markers/) |
| 4 | Retromer Complex | [claim4_retromer_complex_colab.md](03_notebooks_group2_proteostasis/claim4_retromer_complex_colab.md) | [Results](results/group2_proteostasis/Claim4_Retromer_Complex/) |
| 5 | Autophagy vs UPS | [claim5_autophagy_vs_ups_colab.md](03_notebooks_group2_proteostasis/claim5_autophagy_vs_ups_colab.md) | [Results](results/group2_proteostasis/Claim5_Autophagy_vs_UPS/) |
| 6 | Endolysosomal Changes | [claim6_endolysosomal_changes_colab.md](03_notebooks_group2_proteostasis/claim6_endolysosomal_changes_colab.md) | [Results](results/group2_proteostasis/Claim6_Endolysosomal_Changes/) |
| 7 | Temporal Cascade | [claim7_temporal_cascade_colab.md](03_notebooks_group2_proteostasis/claim7_temporal_cascade_colab.md) | [Results](results/group2_proteostasis/Claim7_Temporal_Cascade/) |
| 8 | Rab GTPases | [claim8_rab_gtpases_colab.md](03_notebooks_group2_proteostasis/claim8_rab_gtpases_colab.md) | [Results](results/group2_proteostasis/Claim8_Rab_GTPases/) |

## 🐍 Scripts

| Script | Purpose | Usage |
|--------|---------|-------|
| [run_full_analysis.py](scripts/run_full_analysis.py) | Main pipeline (700+ lines) | `python3 run_full_analysis.py` |
| [execute_all_notebooks.py](scripts/execute_all_notebooks.py) | Notebook executor | `python3 execute_all_notebooks.py` |
| [config.py](scripts/config.py) | Configuration | Import in scripts |

## 📚 Documentation

| Document | Description |
|----------|-------------|
| [README.md](README.md) | Main documentation |
| [README_COLAB.md](docs/README_COLAB.md) | Google Colab guide |
| [PYDESEQ2_INTEGRATION_GUIDE.md](docs/PYDESEQ2_INTEGRATION_GUIDE.md) | Advanced statistics guide |
| [EXECUTION_REPORT.md](docs/EXECUTION_REPORT.md) | Latest execution report |
| [ANALYSIS_EXECUTION_SUMMARY.md](docs/ANALYSIS_EXECUTION_SUMMARY.md) | Detailed execution summary |
| [NOTEBOOK_STATUS.md](docs/NOTEBOOK_STATUS.md) | Notebook inventory |

## 📊 Results Structure

```
results/
├── group1_mitochondrial/
│   ├── Claim1_UPS_Proteins/
│   │   ├── report.md           # Human-readable report
│   │   ├── results.csv         # Raw results
│   │   ├── statistics.json     # Statistical summary
│   │   ├── volcano_plot.png    # Differential expression
│   │   ├── heatmap.png         # Expression patterns
│   │   └── bar_plot.png        # Group comparisons
│   └── ... (7 more claims)
├── group2_proteostasis/
│   └── ... (8 claims)
└── combined/
    ├── master_report.md         # Overall summary
    ├── all_results.csv          # Combined data
    └── all_statistics.json      # All statistics
```

## 📈 Quick Statistics

- **Total Notebooks:** 16
- **Total Scripts:** 4
- **Documentation Files:** 10+
- **Results Generated per Run:** 113 files
  - 48 PNG figures
  - 48 PDF figures
  - 16 CSV results
  - 16 JSON summaries
  - 17 MD reports

## 🔧 Requirements

```bash
# Install all dependencies
pip install -r config/requirements.txt

# Core packages:
numpy pandas scipy matplotlib seaborn
scanpy anndata statsmodels pertpy
```

## 🚦 Status Indicators

| Component | Status |
|-----------|--------|
| Notebooks | ✅ All 16 validated |
| Scripts | ✅ Fully operational |
| Documentation | ✅ Complete |
| Mock Data | ✅ Working |
| Real Data | ⏳ Ready when available |
| Visualizations | ✅ Generating correctly |
| Reports | ✅ Auto-generated |

## 🎯 Next Steps

1. **Run with real data:** Replace mock data with `pool_processed_v2.h5ad`
2. **Review results:** Check `results/combined/master_report.md`
3. **Customize:** Modify protein lists in `scripts/run_full_analysis.py`
4. **Extend:** Add new claims following existing patterns

---

**Navigation:** [README](README.md) | [Run Analysis](scripts/run_full_analysis.py) | [View Results](results/combined/master_report.md) | [Documentation](docs/)