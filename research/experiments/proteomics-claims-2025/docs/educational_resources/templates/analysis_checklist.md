# Analysis Checklist & Key Functions

## Initial Setup
- [ ] Load data with scanpy: `adata = sc.read_h5ad('pool_processed_v2.h5ad')`
- [ ] Check data dimensions: `adata.shape`
- [ ] Review metadata columns: `adata.obs.columns`
- [ ] Check protein names: `adata.var_names`
- [ ] Verify tau status groups exist
- [ ] Check MC1 and pseudotime columns

## Key Analysis Functions

### Differential Expression (Tau+ vs Tau-)
```python
# Split by tau status
tau_pos = adata[adata.obs['tau_status'] == 'positive']
tau_neg = adata[adata.obs['tau_status'] == 'negative']

# Calculate log2 fold change
gene_pos = tau_pos[:, gene_name].X.mean()
gene_neg = tau_neg[:, gene_name].X.mean()
log2fc = np.log2(gene_pos / gene_neg)

# T-test
t_stat, p_val = ttest_ind(tau_pos[:, gene_name].X.flatten(),
                          tau_neg[:, gene_name].X.flatten())

# FDR correction
from statsmodels.stats.multitest import multipletests
rejected, p_adjusted, _, _ = multipletests(p_values, method='fdr_bh')
```

### Correlation Analysis
```python
# Simple correlation
from scipy.stats import pearsonr, spearmanr
corr, p_val = pearsonr(x, y)

# Bonferroni correction
bonferroni_alpha = 0.05 / n_tests
```

### Sliding Window Analysis
```python
# Sort by pseudotime
sorted_idx = np.argsort(adata.obs['pseudotime'])
sorted_data = adata[sorted_idx]

# Sliding window
window_size = 20
correlations = []
positions = []

for i in range(len(sorted_data) - window_size + 1):
    window = sorted_data[i:i+window_size]
    x = window[:, 'GENE1'].X.flatten()
    y = window[:, 'GENE2'].X.flatten()
    corr, _ = pearsonr(x, y)
    correlations.append(corr)
    positions.append(window.obs['pseudotime'].mean())
```

### Cohen's d Calculation
```python
def cohens_d(group1, group2):
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    return (np.mean(group1) - np.mean(group2)) / pooled_std
```

## Statement-Specific Checks

### Statement 1 (UPS proteins)
- [ ] Identify UPS proteins in dataset
- [ ] Compare tau+ vs tau- for each
- [ ] Check significance

### Statement 2 (SQSTM1)
- [ ] Calculate log2FC for SQSTM1
- [ ] Check FDR value
- [ ] Correlate with pseudotime (β coefficient)

### Statement 3 (Autophagy proteins)
- [ ] Check SQSTM1, BECN1, CTSD (up?)
- [ ] Check ATG12, ULK1, CTSL (down?)
- [ ] Verify UPS proteins unchanged

### Statement 4 (10 proteins)
- [ ] Test all 10 correlations with SQSTM1
- [ ] Test all 10 correlations with pseudotime
- [ ] Apply Bonferroni (α = 0.0025 for 20 tests)

### Statement 5 (Global correlation)
- [ ] Calculate r for SQSTM1-VDAC1
- [ ] Check p-value

### Statement 6 (Running correlation)
- [ ] Window size = 20
- [ ] Calculate mean r for pseudotime < 0.33
- [ ] Calculate mean r for pseudotime > 0.67
- [ ] Test trend significance

### Statement 7 (CYCS pattern)
- [ ] Group by MC1 < 2.5 vs ≥ 3.0
- [ ] Calculate means and SDs
- [ ] T-test and Cohen's d

### Statement 8 (CYCS-VATPase)
- [ ] Identify V-ATPase subunits
- [ ] Calculate correlation

## Important Notes
- Data is already log2 transformed
- Use robust statistics for cofactors
- Match exact numbers not required, focus on direction/significance
- Document any proteins not found in dataset