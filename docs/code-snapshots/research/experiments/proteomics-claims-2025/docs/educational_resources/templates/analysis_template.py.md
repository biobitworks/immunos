---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/docs/educational_resources/templates/analysis_template.py
relative: research/experiments/proteomics-claims-2025/docs/educational_resources/templates/analysis_template.py
generated_at: 2025-12-23 10:28
---

```python
"""
Analysis Template for Finding Group Evaluation
Dataset: pool_processed_v2.h5ad
"""

import pandas as pd
import numpy as np
import scanpy as sc
import scipy.stats as stats
from scipy.stats import pearsonr, spearmanr, ttest_ind
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

# Configuration
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Load the dataset
adata = sc.read_h5ad('data/pool_processed_v2.h5ad')

# Initial data exploration
print("Dataset shape:", adata.shape)
print("\nMetadata columns:", adata.obs.columns.tolist())
print("\nFirst few samples:")
print(adata.obs.head())

# Check tau status distribution
print("\nTau status distribution:")
print(adata.obs['tau_status'].value_counts() if 'tau_status' in adata.obs.columns else "tau_status not found")

# Check MC1 score statistics
if 'MC1' in adata.obs.columns:
    print("\nMC1 statistics:")
    print(adata.obs['MC1'].describe())

# Check pseudotime
if 'pseudotime' in adata.obs.columns:
    print("\nPseudotime statistics:")
    print(adata.obs['pseudotime'].describe())

# ==========================================
# STATEMENT 1: UPS Proteins Analysis
# ==========================================
def analyze_ups_proteins():
    """Analyze UPS protein alterations between tau-positive and tau-negative neurons"""

    # List of UPS proteins to check (you'll need to identify these in the data)
    ups_proteins = []  # Add UPS protein names here

    # Split by tau status
    tau_pos = adata[adata.obs['tau_status'] == 'positive']
    tau_neg = adata[adata.obs['tau_status'] == 'negative']

    # Perform differential expression or comparison
    # Your analysis here

    return results

# ==========================================
# STATEMENT 2: SQSTM1 Analysis
# ==========================================
def analyze_sqstm1():
    """Analyze SQSTM1 upregulation and correlation with pseudotime"""

    # Check if SQSTM1 exists in the data
    if 'SQSTM1' in adata.var_names:
        sqstm1_expr = adata[:, 'SQSTM1'].X

        # Calculate log2FC between tau positive and negative
        # Your analysis here

        # Correlation with pseudotime
        if 'pseudotime' in adata.obs.columns:
            corr, p_val = pearsonr(sqstm1_expr.flatten(), adata.obs['pseudotime'])
            print(f"SQSTM1 vs pseudotime: r={corr:.3f}, p={p_val:.3e}")

    return results

# ==========================================
# STATEMENT 3: Autophagy proteins
# ==========================================
def analyze_autophagy_proteins():
    """Analyze autophagy and UPS proteins relative to SQSTM1"""

    proteins_up = ['BECN1', 'CTSD']
    proteins_down = ['ATG12', 'ULK1', 'CTSL']

    # Your analysis here

    return results

# ==========================================
# STATEMENT 4: Protein correlations
# ==========================================
def analyze_protein_correlations():
    """Analyze correlations with SQSTM1 and pseudotime"""

    proteins = ['TAX1BP1', 'CAT', 'VDAC1', 'CYCS', 'ATP5F1A',
                'UQCRC2', 'COX4I1', 'PRDX1', 'KEAP1', 'TFRC']

    # Bonferroni correction threshold
    bonferroni_alpha = 0.05 / len(proteins)  # 0.005 for 10 proteins

    # Your correlation analysis here

    return results

# ==========================================
# STATEMENT 5: SQSTM1-VDAC1 global correlation
# ==========================================
def analyze_sqstm1_vdac1_global():
    """Analyze global correlation between SQSTM1 and VDAC1"""

    if 'SQSTM1' in adata.var_names and 'VDAC1' in adata.var_names:
        sqstm1 = adata[:, 'SQSTM1'].X.flatten()
        vdac1 = adata[:, 'VDAC1'].X.flatten()

        corr, p_val = pearsonr(sqstm1, vdac1)
        print(f"Global SQSTM1-VDAC1: r={corr:.4f}, p={p_val:.3f}")

    return corr, p_val

# ==========================================
# STATEMENT 6: Running correlation analysis
# ==========================================
def analyze_running_correlation():
    """Analyze running correlation along pseudotime with sliding window"""

    window_size = 20

    # Sort by pseudotime
    sorted_idx = np.argsort(adata.obs['pseudotime'])

    # Your sliding window analysis here

    return results

# ==========================================
# STATEMENT 7: CYCS biphasic pattern
# ==========================================
def analyze_cycs_pattern():
    """Analyze CYCS expression pattern relative to MC1 scores"""

    if 'CYCS' in adata.var_names and 'MC1' in adata.obs.columns:
        cycs_expr = adata[:, 'CYCS'].X.flatten()
        mc1_scores = adata.obs['MC1']

        # Define groups
        low_mc1 = cycs_expr[mc1_scores < 2.5]
        high_mc1 = cycs_expr[mc1_scores >= 3.0]

        # T-test and Cohen's d
        t_stat, p_val = ttest_ind(low_mc1, high_mc1)

        # Cohen's d calculation
        pooled_std = np.sqrt(((len(low_mc1)-1)*np.std(low_mc1)**2 +
                              (len(high_mc1)-1)*np.std(high_mc1)**2) /
                             (len(low_mc1) + len(high_mc1) - 2))
        cohens_d = (np.mean(low_mc1) - np.mean(high_mc1)) / pooled_std

        print(f"Low MC1 (<2.5): {np.mean(low_mc1):.3f} ± {np.std(low_mc1):.3f}")
        print(f"High MC1 (≥3.0): {np.mean(high_mc1):.3f} ± {np.std(high_mc1):.3f}")
        print(f"t={t_stat:.3f}, p={p_val:.3e}, Cohen's d={cohens_d:.2f}")

    return results

# ==========================================
# STATEMENT 8: CYCS-VATPase correlation
# ==========================================
def analyze_cycs_vatpase():
    """Analyze correlation between CYCS and V-ATPase"""

    # Need to identify V-ATPase proteins in the dataset
    # Your analysis here

    return results

# Run all analyses
if __name__ == "__main__":
    print("="*50)
    print("Starting Finding Group Analysis")
    print("="*50)

    # Run each analysis function
    # Uncomment as you implement each one

    # analyze_ups_proteins()
    # analyze_sqstm1()
    # analyze_autophagy_proteins()
    # analyze_protein_correlations()
    # analyze_sqstm1_vdac1_global()
    # analyze_running_correlation()
    # analyze_cycs_pattern()
    # analyze_cycs_vatpase()
```
