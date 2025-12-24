---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/docs/educational_resources/templates/analysis_template_group2.py
relative: research/experiments/proteomics-claims-2025/docs/educational_resources/templates/analysis_template_group2.py
generated_at: 2025-12-23 10:28
---

```python
"""
Analysis Template for Finding Group 2 Evaluation
Focus: Sequential Failure of Proteostasis Mechanisms
Dataset: pool_processed_v2.h5ad
"""

import pandas as pd
import numpy as np
import scanpy as sc
import scipy.stats as stats
from scipy.stats import pearsonr, spearmanr, ttest_ind
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import seaborn as sns

# For segmented regression
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

# Configuration
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Load the dataset
adata = sc.read_h5ad('data/pool_processed_v2.h5ad')

# Initial data exploration
print("Dataset shape:", adata.shape)
print(f"Number of cells: {adata.n_obs}")
print(f"Number of proteins: {adata.n_vars}")
print("\nMetadata columns:", adata.obs.columns.tolist())
print("\nFirst few samples:")
print(adata.obs.head())

# ==========================================
# STATEMENT 1: Covariate-controlled DE Analysis
# ==========================================
def covariate_controlled_de_analysis():
    """
    Perform differential expression with covariates: age, PMI, PatientID
    Compare tau-positive vs tau-negative neurons
    """

    # Prepare data for analysis
    tau_pos = adata[adata.obs['tau_status'] == 'positive']
    tau_neg = adata[adata.obs['tau_status'] == 'negative']

    # Need to control for: age, PMI, PatientID
    # Your covariate-controlled analysis here

    # Count significantly altered proteins
    # Apply BH-FDR correction

    return de_results

# ==========================================
# STATEMENT 2: SQSTM1 Analysis
# ==========================================
def analyze_sqstm1_upregulation():
    """
    Check if SQSTM1 is the most upregulated protein
    Should show +3.41 log2 fold change
    """

    # From the covariate-controlled DE results
    # Identify top upregulated proteins

    return top_upregulated

# ==========================================
# STATEMENT 3: Collagen Analysis
# ==========================================
def analyze_collagen_decreases():
    """
    Check collagen proteins (COL1A1, COL1A2, COL6A2)
    Should show > -4.0 log2 fold changes
    """

    collagens = ['COL1A1', 'COL1A2', 'COL6A2']

    # Extract fold changes from DE analysis

    return collagen_results

# ==========================================
# STATEMENT 4: V-ATPase-MC1 Correlations
# ==========================================
def analyze_vatpase_mc1_correlations():
    """
    Within tau-positive cells, correlate V-ATPase subunits with MC1
    Subunits: ATP6V1A, ATP6V1B2, ATP6V1C1, ATP6V1H
    """

    # Filter to tau-positive cells only
    tau_pos = adata[adata.obs['tau_status'] == 'positive']

    vatpase_subunits = ['ATP6V1A', 'ATP6V1B2', 'ATP6V1C1', 'ATP6V1H']

    correlations = {}
    for subunit in vatpase_subunits:
        if subunit in tau_pos.var_names:
            expr = tau_pos[:, subunit].X.flatten()
            mc1 = tau_pos.obs['MC1']
            corr, p_val = pearsonr(expr, mc1)
            correlations[subunit] = {'r': corr, 'p': p_val}
            print(f"{subunit}: r={corr:.3f}, p={p_val:.3e}")

    # Calculate mean correlation
    mean_r = np.mean([c['r'] for c in correlations.values()])
    print(f"\nMean correlation: {mean_r:.3f}")

    return correlations

# ==========================================
# STATEMENT 5: V-ATPase Score Segmented Analysis
# ==========================================
def analyze_vatpase_score_segmented():
    """
    Calculate V-ATPase Score and perform segmented regression
    Find breakpoint at MC1 = 2.831
    """

    # Identify all V-ATPase subunits in data
    vatpase_genes = [gene for gene in adata.var_names if gene.startswith('ATP6V')]

    # Calculate V-ATPase Score (mean of log2 expression)
    vatpase_expr = adata[:, vatpase_genes].X
    vatpase_score = np.mean(vatpase_expr, axis=1)

    # Segmented regression against MC1
    mc1_values = adata.obs['MC1'].values

    # Define piecewise linear function
    def piecewise_linear(x, x0, b1, b2):
        """Piecewise linear with breakpoint at x0"""
        return np.piecewise(x, [x < x0],
                           [lambda x: b1*x,
                            lambda x: b1*x0 + b2*(x-x0)])

    # Fit segmented regression
    # Your segmented regression code here

    return breakpoint, slopes

# ==========================================
# STATEMENT 6: V-ATPase Biphasic Behavior
# ==========================================
def analyze_vatpase_biphasic():
    """
    Analyze V-ATPase biphasic behavior along pseudotime
    Compare with proteasome breakpoint
    """

    # Calculate V-ATPase Score
    vatpase_genes = [gene for gene in adata.var_names if gene.startswith('ATP6V')]
    vatpase_score = np.mean(adata[:, vatpase_genes].X, axis=1)

    # Calculate Proteasome Score (identify proteasome subunits)
    proteasome_genes = [gene for gene in adata.var_names if 'PSM' in gene]
    proteasome_score = np.mean(adata[:, proteasome_genes].X, axis=1)

    pseudotime = adata.obs['pseudotime'].values

    # Segmented regression for both systems
    # Your biphasic analysis here

    return vatpase_breakpoint, proteasome_breakpoint

# ==========================================
# STATEMENT 7: Temporal Order Analysis
# ==========================================
def analyze_temporal_order():
    """
    Compare breakpoints of V-ATPase vs proteasome
    V-ATPase should be at ~0.65 pseudotime
    """

    # Use results from Statement 6
    # Compare timing of breakpoints

    return timing_comparison

# ==========================================
# Helper Functions
# ==========================================
def calculate_cohens_d(group1, group2):
    """Calculate Cohen's d effect size"""
    n1, n2 = len(group1), len(group2)
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
    return (np.mean(group1) - np.mean(group2)) / pooled_std

def perform_fdr_correction(p_values, method='fdr_bh'):
    """Perform FDR correction on p-values"""
    rejected, p_adjusted, _, _ = multipletests(p_values, method=method)
    return rejected, p_adjusted

# Run all analyses
if __name__ == "__main__":
    print("="*50)
    print("Starting Finding Group 2 Analysis")
    print("Focus: Sequential Failure of Proteostasis")
    print("="*50)

    # Run each analysis function
    # Uncomment as you implement each one

    # covariate_controlled_de_analysis()
    # analyze_sqstm1_upregulation()
    # analyze_collagen_decreases()
    # analyze_vatpase_mc1_correlations()
    # analyze_vatpase_score_segmented()
    # analyze_vatpase_biphasic()
    # analyze_temporal_order()
```
