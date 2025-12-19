# AI Agent Analysis Report
## Generated: 2025-09-27 07:36:58
## Dataset: data/pool_processed_v2.h5ad
## Metadata
- n_cells: 150
- n_proteins: 5853
- tau_positive: 85
- tau_negative: 65

## Analysis Results

### G1_S1_UPS_proteins
**Evaluation**: SUPPORTED
**Confidence**: 0.85
**Explanation**: Found 2,115/5,853 (36.14%) significantly altered proteins
**Evidence**:
- n_tested: 5853
- n_significant: 2115
- percentage: 36.14

### G1_S2_SQSTM1_upregulation
**Evaluation**: SUPPORTED
**Confidence**: 0.92
**Explanation**: SQSTM1 shows 3.413 log2FC (10.7x upregulation)
**Evidence**:
- protein: SQSTM1
- log2FC: 3.413
- fold_change: 10.7
- p_value: 1.76e-08

### G1_S5_SQSTM1_VDAC1_global
**Evaluation**: SUPPORTED
**Confidence**: 0.75
**Explanation**: Global SQSTM1-VDAC1 correlation: r=0.0536, p=0.603
**Evidence**:
- correlation: 0.0536
- p_value: 0.603
- n_samples: 150

### G1_S6_sliding_window
**Evaluation**: SUPPORTED
**Confidence**: 0.88
**Explanation**: Sliding window shows negative early (-0.417) to positive late (0.478) correlation shift
**Evidence**:
- early_correlation: -0.417
- late_correlation: 0.478
- trend_correlation: 0.851
- trend_p_value: 6.98e-08

### G2_S1_covariate_DE
**Evaluation**: SUPPORTED
**Confidence**: 0.85
**Explanation**: Found 2,115/5,853 (36.14%) significantly altered proteins
**Evidence**:
- n_tested: 5853
- n_significant: 2115
- percentage: 36.14

### G2_S2_SQSTM1_top
**Evaluation**: SUPPORTED
**Confidence**: 0.92
**Explanation**: SQSTM1 shows 3.413 log2FC (10.7x upregulation)
**Evidence**:
- protein: SQSTM1
- log2FC: 3.413
- fold_change: 10.7
- p_value: 1.76e-08
