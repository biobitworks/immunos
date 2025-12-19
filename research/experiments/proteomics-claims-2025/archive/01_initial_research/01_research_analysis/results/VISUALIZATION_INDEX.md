# Visualization Index

## Overview
Complete index of all visualizations generated for the proteomics analysis project.

## Quick Access by Category

### 📊 Summary Dashboards
- [[figures/master_analysis_dashboard|Master Analysis Dashboard]] - 4-panel overview
- [[paper_replications/summary_dashboard|Paper Replication Dashboard]] - 9-panel comprehensive view
- [[figures/evaluation_summary|Evaluation Summary]] - Pie charts of claim outcomes

### 🔬 Sequential Failure Visualizations

#### Core Analysis Plots
- [[figures/claim_1_plot|V-ATPase Differential Expression]] - Volcano plot and heatmap
- [[figures/claim_2_plot|ATP6V0A1 Stage Analysis]] - Disease progression patterns

#### Paper Replications
- [[paper_replications/figure3_volcano_histogram|Figure 3: Proteome Remodeling]] - 2,115 DE proteins
- [[paper_replications/figure4_vatpase_mc1|Figure 4: V-ATPase Biphasic Pattern]] - MC1 analysis
- [[paper_replications/figure5_vatpase_pseudotime|Figure 5: Pseudotime Breakpoint]] - Sequential failure
- [[paper_replications/figure6_11_vatpase_mc1_segmented|Figure 6: Critical MC1 Threshold]] - MC1 = 2.831

### 🧬 Mitochondrial Dysregulation Visualizations

#### Core Analysis Plots
- [[figures/mito_claim_1_plot|UPS Protein Analysis]] - Differential expression
- [[figures/mito_claim_2_plot|SQSTM1 Upregulation]] - Validation plots

#### Paper Replications
- [[paper_replications/figure7_autophagy_dysregulation|Figure 7: Autophagy Dysregulation]] - Pathway analysis
- [[paper_replications/figure8_sqstm1_vdac1_correlation|Figure 8: SQSTM1-VDAC1 Coupling]] - Dynamic correlation
- [[paper_replications/figure9_cycs_expression|Figure 9: Cytochrome C Decline]] - Biphasic pattern
- [[paper_replications/figure10_coordinated_decline|Figure 10: Coordinated Failure]] - Mito-lyso decline

## Visualization Statistics

### File Counts
- **Core Analysis Figures**: 6
- **Paper Replications**: 9
- **Total Visualizations**: 15

### Key Findings Visualized

#### Critical Thresholds
| System | Breakpoint | Visualization |
|--------|------------|---------------|
| Proteasome | Pseudotime 0.372 | [[paper_replications/figure5_vatpase_pseudotime|Figure 5]] |
| V-ATPase (time) | Pseudotime 0.654 | [[paper_replications/figure5_vatpase_pseudotime|Figure 5]] |
| V-ATPase (MC1) | MC1 2.831 | [[paper_replications/figure6_11_vatpase_mc1_segmented|Figure 6]] |

#### Key Protein Changes
| Protein | Change | Visualization |
|---------|--------|---------------|
| SQSTM1 | +1.32 fold | [[figures/mito_claim_2_plot|SQSTM1 Plot]] |
| V-ATPase | Biphasic | [[figures/claim_1_plot|V-ATPase Plot]] |
| CYCS | -2.58 Cohen's d | [[paper_replications/figure9_cycs_expression|Figure 9]] |

## Navigation

### By Analysis Type
- **Differential Expression**: [[paper_replications/figure3_volcano_histogram|Volcano Plots]]
- **Temporal Analysis**: [[paper_replications/figure5_vatpase_pseudotime|Pseudotime Plots]]
- **Correlation Analysis**: [[paper_replications/figure8_sqstm1_vdac1_correlation|Running Correlations]]
- **Segmented Regression**: [[paper_replications/figure6_11_vatpase_mc1_segmented|Breakpoint Analysis]]

### By Biological System
- **Proteostasis**: Claims 1-8 (Sequential Failure)
- **Mitochondria**: Claims 1-8 (Mitochondrial Dysregulation)
- **Autophagy**: [[paper_replications/figure7_autophagy_dysregulation|Figure 7]]
- **Lysosomes**: [[paper_replications/figure4_vatpase_mc1|V-ATPase Figures]]

## File Locations

```
01_research_analysis/results/
├── figures/                    # Core analysis plots
│   ├── claim_1_plot.png
│   ├── claim_2_plot.png
│   ├── evaluation_summary.png
│   ├── master_analysis_dashboard.png
│   ├── mito_claim_1_plot.png
│   └── mito_claim_2_plot.png
└── paper_replications/         # Paper figure replications
    ├── figure3_volcano_histogram.png
    ├── figure4_vatpase_mc1.png
    ├── figure5_vatpase_pseudotime.png
    ├── figure6_11_vatpase_mc1_segmented.png
    ├── figure7_autophagy_dysregulation.png
    ├── figure8_sqstm1_vdac1_correlation.png
    ├── figure9_cycs_expression.png
    ├── figure10_coordinated_decline.png
    └── summary_dashboard.png
```

## Quick Links

- [[../README|Results README]]
- [[../../INDEX|Project Index]]
- [[../PAPER_FIGURE_REPLICATION_SUMMARY|Replication Summary]]
- [[../sequential_failure/README|Sequential Failure Results]]
- [[../mitochondrial_dysregulation/README|Mitochondrial Results]]

---

*Last updated: 2024-09-28*
*Total visualizations: 15 high-quality figures*
