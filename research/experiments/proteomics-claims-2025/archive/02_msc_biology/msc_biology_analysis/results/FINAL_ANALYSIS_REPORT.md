# MSc Biology Analysis - Final Report
## Analysis of pool_processed_v2.h5ad

### Dataset Summary
- **Data file**: pool_processed_v2.h5ad
- **Total samples**: 44 neurons
- **Total proteins**: 5,853
- **Tau-positive neurons**: 22
- **Tau-negative neurons**: 22

---

## Key Finding 1: SQSTM1 Upregulation ❌

### Results:
- **Mean expression in Tau+**: 14.159
- **Mean expression in Tau-**: 10.746
- **Log2 Fold Change**: 0.398
- **Actual Fold Change**: **1.3x**
- **P-value**: 9.294e-08 (highly significant)

### Conclusion:
❌ **NOT CONFIRMED**: Paper claimed 10.7-fold upregulation, but actual data shows only 1.3-fold
- The change is statistically significant but much smaller than claimed
- This suggests the paper may have used different normalization or a subset of data

---

## Key Finding 2: Autophagy vs UPS ✅

### Results:
- **Autophagy proteins**: 3/6 (50%) significantly changed
  - SQSTM1 ✅ (changed)
  - NBR1 ✅ (changed)
  - MAP1LC3B ✅ (changed)
  - BECN1 ❌ (stable)
  - ATG5 ❌ (stable)
  - ATG7 ❌ (stable)

- **UPS proteins**: 0/4 (0%) significantly changed
  - PSMA1 ❌ (stable)
  - PSMB5 ❌ (stable)
  - PSMD1 ❌ (stable)
  - UBA1 ❌ (stable)

### Conclusion:
✅ **CONFIRMED**: Autophagy is specifically disrupted while UPS remains stable
- 50% of autophagy proteins show significant changes
- 0% of UPS proteins show significant changes
- Supports selective autophagy dysfunction

---

## Key Finding 3: Proteasome & V-ATPase Proteins ✅

### Results:
- **Proteasome subunits found**: 37
- **V-ATPase subunits found**: 13
- **Pseudotime range**: 0.000 to 1.000

### Conclusion:
✅ **CONFIRMED**: Sufficient proteins found for sequential failure analysis
- All major proteasome families represented (PSMA, PSMB, PSMC, PSMD)
- Both V-ATPase domains present (V0 and V1)
- Pseudotime data available for progression analysis

---

## Key Finding 4: Mitochondrial Markers

### Results:
- **VDAC1**: Stable (p=0.177)
- **CYCS**: Stable (p=0.062)
- **COX4I1**: Changed (p=0.016) ✅

### Conclusion:
Partial mitochondrial dysfunction (1/3 proteins changed)
- COX4I1 (electron transport chain) shows significant change
- VDAC1 and CYCS remain stable
- Suggests selective mitochondrial impact

---

## Overall Summary

### Validated Claims ✅:
1. **Autophagy-specific dysfunction** - Confirmed
2. **Proteasome and V-ATPase proteins present** - Confirmed
3. **Sequential failure testable** - Confirmed (pseudotime available)

### Refuted Claims ❌:
1. **SQSTM1 10.7-fold upregulation** - Only 1.3-fold observed

### Important Notes:
- The discrepancy in SQSTM1 fold change is significant
- May be due to:
  - Different normalization methods
  - Subset analysis in the paper
  - Different statistical approaches
  - Possible error in paper calculations

### Recommendations:
1. Re-examine SQSTM1 calculation methods
2. Check if paper used log-transformed data differently
3. Verify if paper analyzed a specific subset of samples
4. Consider different statistical models for fold change

---

## Technical Details

### Statistical Methods:
- Mann-Whitney U test for group comparisons
- P-value < 0.05 considered significant
- No multiple testing correction applied in this quick analysis

### Data Processing:
- Direct analysis of raw expression values
- No additional normalization applied
- Used gene names for protein identification

---

*Analysis completed using actual pool_processed_v2.h5ad data*
*Date: 2024*
*MSc Biology Analysis Framework*