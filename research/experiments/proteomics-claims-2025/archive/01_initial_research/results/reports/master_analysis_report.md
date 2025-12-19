# Master Analysis Report - pool_processed_v2.h5ad
Generated: 2025-09-28 08:49

## Dataset Overview
- **Samples**: 44 (22 tau+, 22 tau-)
- **Proteins**: 5853
- **Data file**: /Users/byron/project_plan/data/pool_processed_v2.h5ad

## Key Findings

### 1. SQSTM1 Upregulation ❌
- **Claimed**: 10.7-fold
- **Observed**: 1.3-fold
- **Discrepancy**: 8.1x difference
- **P-value**: 9.294e-08

### 2. Autophagy vs UPS ✅
- **Autophagy**: 57.1% disrupted
- **UPS**: 28.6% disrupted
- **Conclusion**: Autophagy-specific dysfunction

### 3. Proteasome & V-ATPase
- **Proteasome proteins**: 24
- **V-ATPase proteins**: 9
- **Sequential failure testable**: Yes

## Known Issues
- SQSTM1 fold change discrepancy: Significant discrepancy - needs investigation
- Proteins with semicolons: 56 (1.0%)

## Recommendations
1. Investigate SQSTM1 normalization methods
2. Apply multiple testing correction consistently
3. Validate findings with additional datasets
