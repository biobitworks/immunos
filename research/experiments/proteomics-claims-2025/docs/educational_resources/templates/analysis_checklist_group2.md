# Analysis Checklist - Finding Group 2

## Initial Setup
- [ ] Load data: `adata = sc.read_h5ad('pool_processed_v2.h5ad')`
- [ ] Verify 5,853 proteins present
- [ ] Check covariates exist: age, PMI, PatientID
- [ ] Verify tau status groups
- [ ] Check MC1 and pseudotime columns

## Key Techniques for This Group

### Covariate-Controlled Differential Expression
```python
# Using statsmodels for covariate control
import statsmodels.api as sm
from statsmodels.formula.api import ols

# Example approach
for protein in adata.var_names:
    df = pd.DataFrame({
        'expression': adata[:, protein].X.flatten(),
        'tau_status': adata.obs['tau_status'],
        'age': adata.obs['age'],
        'PMI': adata.obs['PMI'],
        'PatientID': adata.obs['PatientID']
    })

    # Linear model with covariates
    model = ols('expression ~ tau_status + age + PMI + C(PatientID)', data=df)
    results = model.fit()

    # Extract coefficient and p-value for tau_status
```

### Segmented Regression (Breakpoint Analysis)
```python
from scipy.optimize import curve_fit

def piecewise_linear(x, x0, y0, b1, b2):
    """Piecewise linear function with breakpoint at x0"""
    return np.piecewise(x, [x < x0],
                        [lambda x: y0 + b1*(x-x0),
                         lambda x: y0 + b2*(x-x0)])

# Fit to data
popt, _ = curve_fit(piecewise_linear, x_data, y_data, p0=[2.8, 15, -0.003, -0.4])
breakpoint = popt[0]
slope1, slope2 = popt[2], popt[3]
```

### V-ATPase Score Calculation
```python
# Identify V-ATPase subunits
vatpase_genes = ['ATP6V1A', 'ATP6V1B2', 'ATP6V1C1', 'ATP6V1H',
                 'ATP6V0A1', 'ATP6V0D1', 'ATP6V1E1', 'ATP6V1F']

# Filter to available genes
available_vatpase = [g for g in vatpase_genes if g in adata.var_names]

# Calculate mean expression (already log2)
vatpase_score = np.mean(adata[:, available_vatpase].X, axis=1)
```

### Proteasome Score Calculation
```python
# Proteasome subunits typically start with PSM
proteasome_genes = [g for g in adata.var_names if 'PSM' in g]

# Or specific subunits
proteasome_genes = ['PSMA1', 'PSMA2', 'PSMB1', 'PSMB2', 'PSMC1', 'PSMD1']

proteasome_score = np.mean(adata[:, proteasome_genes].X, axis=1)
```

## Statement-Specific Checks

### Statement 1 (DE with covariates)
- [ ] Control for: age, PMI, PatientID
- [ ] Apply BH-FDR correction
- [ ] Count proteins with FDR < 0.05
- [ ] Should find ~2,115 proteins (36.14%)

### Statement 2 (SQSTM1)
- [ ] From DE results, identify top upregulated
- [ ] Check SQSTM1 log2FC (~3.41)
- [ ] Verify it's #1 upregulated

### Statement 3 (Collagens)
- [ ] Check COL1A1, COL1A2, COL6A2
- [ ] Verify log2FC < -4.0
- [ ] Check if largest decreases

### Statement 4 (V-ATPase vs MC1)
- [ ] Filter to tau-positive only
- [ ] Check ATP6V1A, ATP6V1B2, ATP6V1C1, ATP6V1H
- [ ] Calculate correlations with MC1
- [ ] Mean should be ~-0.723

### Statement 5 (V-ATPase segmented)
- [ ] Calculate V-ATPase Score
- [ ] Segmented regression vs MC1
- [ ] Find breakpoint (~2.831)
- [ ] Check slopes before/after

### Statement 6 (Biphasic behavior)
- [ ] V-ATPase vs pseudotime
- [ ] Proteasome vs pseudotime
- [ ] Find both breakpoints
- [ ] V-ATPase at ~0.654
- [ ] Compare difference (0.282 units)

### Statement 7 (Temporal order)
- [ ] Compare breakpoint timings
- [ ] V-ATPase should be later (~0.65)
- [ ] Confirms lysosomal after proteasomal

## Important Notes
- Use covariate control for DE analysis
- FDR correction essential (Benjamini-Hochberg)
- Segmented regression for breakpoint detection
- Focus on tau-positive cells for Statement 4
- Calculate composite scores (mean of subunits)
- Document any missing proteins