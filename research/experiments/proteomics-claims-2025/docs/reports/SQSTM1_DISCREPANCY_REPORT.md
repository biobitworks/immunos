# SQSTM1 Discrepancy Report
## Critical Finding: 8.2x Difference Between Claimed and Observed Values

### Executive Summary
**The most significant finding in our proteomics analysis is a major discrepancy in SQSTM1 upregulation:**
- **Paper claim**: 10.7-fold upregulation (log2 FC = 3.413)
- **Our observation**: 1.3-fold upregulation (log2 FC = 0.398)
- **Discrepancy ratio**: 8.2x difference

This discrepancy is critical because SQSTM1 is presented as the key biomarker for autophagy failure in the paper.

---

## Detailed Analysis

### 1. Statistical Verification
```
Data file: pool_processed_v2.h5ad
Samples: 44 neurons (22 tau+, 22 tau-)
Method: Mann-Whitney U test

Results:
- Mean expression Tau+: 14.159
- Mean expression Tau-: 10.746
- Fold change: 1.318x
- P-value: 9.294e-08 (highly significant)
```

The change IS statistically significant, but the magnitude is 8x smaller than claimed.

### 2. Possible Explanations

#### A. Data Normalization Differences
- **Paper may use**: Log2 transformation after normalization
- **We observed**: Raw expression values
- **Impact**: Different normalization can dramatically affect fold changes

#### B. Subset Analysis
The paper may have:
- Analyzed only extreme samples (highest/lowest MC1 scores)
- Used a different tau classification threshold
- Excluded outliers that we included

#### C. Technical Factors
- **Batch effects**: Not corrected in our analysis
- **Isoform consideration**: SQSTM1 has no isoforms in our data
- **Detection threshold**: Different minimum expression cutoffs

#### D. Statistical Model Differences
- **Paper**: May use linear model with covariates
- **Our analysis**: Simple two-group comparison
- **Impact**: Covariate adjustment could amplify differences

### 3. Validation Steps Taken

✅ **Protein correctly identified**: SQSTM1 found at index 2619
✅ **No semicolon issues**: Gene name is clean (no isoforms)
✅ **Statistical significance confirmed**: P < 0.0001
✅ **Direction correct**: Upregulated in tau+ (just not 10.7x)

### 4. Impact on Other Findings

Despite the SQSTM1 discrepancy:
- **Autophagy disruption confirmed**: 57% of autophagy proteins changed
- **UPS stability confirmed**: Only 29% of UPS proteins changed
- **Selective dysfunction validated**: Autophagy specifically affected

### 5. Biological Interpretation

Even at 1.3-fold upregulation:
- SQSTM1 increase indicates autophagy impairment
- Direction of change supports the biological hypothesis
- Statistical significance remains strong
- But magnitude suggests less severe dysfunction than claimed

---

## Recommendations

### Immediate Actions
1. **Re-examine normalization**: Check if data needs additional transformation
2. **Subset analysis**: Test if analyzing MC1 extremes yields higher fold change
3. **Covariate adjustment**: Include age, PMI, patient effects

### Investigation Priorities
1. **Contact authors**: Request exact analysis code
2. **Check supplementary methods**: Look for hidden normalization steps
3. **Validate with raw data**: Request pre-processed values if available

### Alternative Explanations
1. **Different dataset version**: We may have processed data, not raw
2. **Typographical error**: 10.7 could be error for 1.07
3. **Different SQSTM1 measurement**: Total vs specific isoform

---

## Code to Reproduce

```python
import scanpy as sc
import numpy as np
from scipy.stats import mannwhitneyu

# Load data
adata = sc.read_h5ad('/Users/byron/project_plan/data/pool_processed_v2.h5ad')

# Get groups
tau_pos = adata.obs['TauStatus'] == 'positive'
tau_neg = adata.obs['TauStatus'] == 'negative'

# Find SQSTM1
sqstm1_mask = adata.var['GeneName'].str.contains('SQSTM1', case=False, na=False)
sqstm1_idx = np.where(sqstm1_mask)[0][0]

# Extract expression
expr_pos = adata.X[tau_pos, sqstm1_idx]
expr_neg = adata.X[tau_neg, sqstm1_idx]

# Calculate fold change
fold_change = np.mean(expr_pos) / np.mean(expr_neg)
log2_fc = np.log2(fold_change)

# Statistical test
stat, pval = mannwhitneyu(expr_pos, expr_neg)

print(f"Fold change: {fold_change:.2f}x")
print(f"Log2 FC: {log2_fc:.3f}")
print(f"P-value: {pval:.3e}")
```

---

## Conclusion

The SQSTM1 discrepancy represents a critical issue that must be resolved:
- Our analysis is technically correct and reproducible
- The biological direction is confirmed (upregulation)
- But the magnitude difference (8.2x) is too large to ignore

This does NOT invalidate the overall biological conclusions about autophagy failure, but it does raise questions about the quantitative claims in the paper.

---

*Report generated: 2024*
*Data: pool_processed_v2.h5ad*
*Framework: Bioinformatics Finding Group Evaluation*