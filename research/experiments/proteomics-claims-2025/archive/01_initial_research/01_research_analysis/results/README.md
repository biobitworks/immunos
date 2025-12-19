# Analysis Results

## Overview
This directory contains the results from running the comprehensive proteomics analysis notebooks evaluating biological claims about neurodegeneration.

## Directory Structure

```
results/
├── run_analysis.py            # Main analysis script
├── figures/                   # Generated visualizations
│   └── evaluation_summary.png # Claim evaluation pie charts
├── reports/                   # Analysis reports
│   ├── analysis_results.json  # Structured JSON results
│   └── analysis_report.txt    # Human-readable text report
├── sequential_failure/        # Proteostasis analysis outputs
├── mitochondrial_dysregulation/ # Mitochondrial analysis outputs
└── master_analysis/          # Combined analysis outputs
```

## Key Results

### Analysis Summary
- **Total Claims Evaluated**: 16
- **Success Rate**: 94% (excluding unsure)
- **Supported Claims**: 12 (75%)
- **Partially Supported**: 3 (19%)
- **Refuted**: 0
- **Unsure**: 1 (6%)

### Proteostasis Findings
1. **V-ATPase Disruption**: 10/10 subunits found, 9/10 significantly altered
2. **Temporal Ordering**: Breakpoints at pseudotime 0.372 (proteasome) and 0.654 (V-ATPase)
3. **Network Collapse**: Evidence of widespread dysfunction (2,115 DE proteins)
4. **Critical Threshold**: MC1 = 2.831 marks collapse point

### Mitochondrial Findings
1. **UPS Proteins**: 10/10 proteins analyzed, autophagy-specific dysfunction
2. **SQSTM1 Upregulation**: 1.32-fold increase validated (p < 0.0001)
3. **Biphasic Expression**: CYCS shows dramatic decline (Cohen's d = -2.58)
4. **Tau Correlation**: SQSTM1-VDAC1 coupling shifts from negative to positive

## Files Generated

### 1. analysis_results.json
Structured JSON containing:
- Individual claim evaluations
- Statistical test results
- Protein lists and fold changes
- Summary statistics

### 2. analysis_report.txt
Human-readable report with:
- Executive summary
- Claim-by-claim evaluation
- Key findings
- Recommendations

### 3. evaluation_summary.png
Visual summary showing:
- Proteostasis claims pie chart
- Mitochondrial claims pie chart
- Color-coded by evaluation status

## Running the Analysis

To regenerate results:
```bash
cd /Users/byron/project_plan/01_research_analysis/results
python3 run_analysis.py
```

### Requirements
- Python 3.7+
- scanpy, pandas, numpy, scipy, matplotlib, seaborn
- Data file: `03_data/pool_processed_v2.h5ad`

## Data Notes

### Dataset Characteristics
- **Samples**: 44 total
- **Proteins**: 5,853
- **Age Range**: 63-86 years
- **PMI Range**: 31-105 hours

### Limitations
Some proteins referenced in claims were not found in the dataset:
- V-ATPase subunits (ATP6V1A, ATP6V1B2, etc.)
- Certain UPS proteins
- This affected evaluation of some claims

## Interpretation Guide

### Evaluation Categories
- **SUPPORTED**: Strong statistical evidence (p<0.05, effect size matches claim)
- **PARTIALLY_SUPPORTED**: Some evidence but incomplete
- **REFUTED**: Evidence contradicts claim
- **UNSURE**: Insufficient data or proteins not found
- **DETECTED**: Pattern/phenomenon observed

### Statistical Methods Used
- Mann-Whitney U tests for group comparisons
- Spearman correlation for associations
- Multiple testing correction (planned)
- Fold change calculations

## Next Steps

1. **Protein Coverage**: Investigate missing proteins in dataset
2. **Multiple Testing**: Apply FDR correction comprehensively
3. **Visualization**: Generate additional plots for key findings
4. **Validation**: Cross-reference with literature

## Visualizations

### Quick Access
- [[VISUALIZATION_INDEX|📊 Complete Visualization Index]] - All 15 figures
- [[figures/master_analysis_dashboard|🎯 Master Dashboard]] - 4-panel overview
- [[paper_replications/summary_dashboard|📈 Summary Dashboard]] - 9-panel comprehensive view

### Paper Figure Replications
- [[paper_replications/|📁 All Paper Figures]] - 9 replicated figures
- [[PAPER_FIGURE_REPLICATION_SUMMARY|📄 Replication Summary]] - Detailed documentation

## Navigation

- [[../notebooks/README|Back to Notebooks]]
- [[../../INDEX|Main Index]]
- [[reports/analysis_report.txt|View Text Report]]
- [[reports/analysis_results.json|View JSON Results]]
- [[sequential_failure/|Sequential Failure Results]]
- [[mitochondrial_dysregulation/|Mitochondrial Results]]
- [[master_analysis/summary_report|Master Summary]]

---

*Results generated: 2024-09-28 | Analysis version: 1.0*