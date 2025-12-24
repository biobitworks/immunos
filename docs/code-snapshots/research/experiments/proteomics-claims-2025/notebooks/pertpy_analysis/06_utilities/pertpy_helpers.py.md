---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/notebooks/pertpy_analysis/06_utilities/pertpy_helpers.py
relative: research/experiments/proteomics-claims-2025/notebooks/pertpy_analysis/06_utilities/pertpy_helpers.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
PertPy Helper Functions for DGE Analysis
Utility functions for differential expression analysis using PertPy
"""

import pertpy as pt
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Union, List, Optional, Tuple


def run_pydeseq2_safe(adata, design="~tau_status", contrast_column="tau_status",
                      baseline="negative", compare="positive"):
    """
    Run PyDESeq2 with fallback to traditional methods if it fails

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    design : str
        Design formula for DESeq2
    contrast_column : str
        Column to use for contrast
    baseline : str
        Baseline group
    compare : str
        Comparison group

    Returns:
    --------
    pd.DataFrame : Results dataframe with log2FC, pvalue, padj
    """
    try:
        # Use counts if available
        if 'counts' in adata.layers:
            adata.layers['log2'] = adata.X.copy()
            adata.X = adata.layers['counts'].copy()

        # Initialize PyDESeq2
        pds2 = pt.tl.PyDESeq2(adata=adata, design=design, refit_cooks=True)
        pds2.fit()

        # Test contrast
        results = pds2.test_contrasts(
            pds2.contrast(
                column=contrast_column,
                baseline=baseline,
                group_to_compare=compare
            )
        )

        print("✓ PyDESeq2 completed successfully")
        return results

    except Exception as e:
        print(f"⚠ PyDESeq2 failed: {e}")
        print("Falling back to traditional differential expression...")

        # Traditional analysis
        results_list = []
        group1 = adata.obs[contrast_column] == compare
        group2 = adata.obs[contrast_column] == baseline

        for i in range(adata.n_vars):
            expr1 = adata.X[group1, i]
            expr2 = adata.X[group2, i]

            # Calculate statistics
            log2fc = np.mean(expr1) - np.mean(expr2)
            tstat, pval_t = stats.ttest_ind(expr1, expr2)
            ustat, pval_mw = stats.mannwhitneyu(expr1, expr2, alternative='two-sided')

            # Cohen's d
            pooled_std = np.sqrt(
                ((len(expr1)-1)*np.var(expr1) + (len(expr2)-1)*np.var(expr2)) /
                (len(expr1) + len(expr2) - 2)
            )
            cohen_d = log2fc / pooled_std if pooled_std > 0 else 0

            results_list.append({
                'protein': adata.var.index[i],
                'log2FoldChange': log2fc,
                'pvalue': pval_t,
                'pvalue_mw': pval_mw,
                'stat': tstat,
                'cohen_d': cohen_d
            })

        results = pd.DataFrame(results_list)

        # Add FDR correction
        results['padj'] = multipletests(results['pvalue'], method='fdr_bh')[1]

        return results


def find_proteins_in_dataset(adata, protein_list, fuzzy_match=True):
    """
    Find proteins in dataset with fuzzy matching

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    protein_list : list
        List of protein names to find
    fuzzy_match : bool
        Whether to use fuzzy matching

    Returns:
    --------
    dict : Dictionary with 'found' and 'missing' protein lists
    """
    protein_names = adata.var['protein_name'] if 'protein_name' in adata.var else adata.var.index
    found = []
    missing = []

    for protein in protein_list:
        if protein in protein_names.tolist():
            found.append(protein)
        elif fuzzy_match:
            # Try partial matching
            matches = [p for p in protein_names if protein.upper() in p.upper() or p.upper() in protein.upper()]
            if matches:
                found.append(matches[0])
            else:
                missing.append(protein)
        else:
            missing.append(protein)

    return {'found': found, 'missing': missing}


def create_volcano_plot(results_df, title="Volcano Plot",
                        fc_threshold=0.5, pval_threshold=0.05,
                        label_top=10, save_path=None):
    """
    Create volcano plot from results

    Parameters:
    -----------
    results_df : pd.DataFrame
        Results with log2FoldChange and padj columns
    title : str
        Plot title
    fc_threshold : float
        Fold change threshold for significance
    pval_threshold : float
        P-value threshold
    label_top : int
        Number of top proteins to label
    save_path : str
        Path to save figure

    Returns:
    --------
    fig, ax : matplotlib figure and axes
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Calculate -log10(p-value)
    results_df['neg_log10_pval'] = -np.log10(results_df['padj'] + 1e-300)

    # Define colors
    colors = []
    for _, row in results_df.iterrows():
        if row['padj'] < pval_threshold and abs(row['log2FoldChange']) > fc_threshold:
            colors.append('red')  # Significant
        elif row['padj'] < pval_threshold:
            colors.append('orange')  # Significant but small effect
        elif abs(row['log2FoldChange']) > fc_threshold:
            colors.append('blue')  # Large effect only
        else:
            colors.append('gray')  # Not significant

    # Create scatter plot
    scatter = ax.scatter(results_df['log2FoldChange'],
                        results_df['neg_log10_pval'],
                        c=colors, alpha=0.6, s=50)

    # Add threshold lines
    ax.axhline(y=-np.log10(pval_threshold), color='black', linestyle='--', alpha=0.3)
    ax.axvline(x=fc_threshold, color='black', linestyle='--', alpha=0.3)
    ax.axvline(x=-fc_threshold, color='black', linestyle='--', alpha=0.3)

    # Label top proteins
    if label_top > 0:
        top_proteins = results_df.nsmallest(label_top, 'padj')
        for _, row in top_proteins.iterrows():
            if row['padj'] < pval_threshold:
                ax.annotate(row['protein'],
                           (row['log2FoldChange'], row['neg_log10_pval']),
                           fontsize=8, alpha=0.7)

    ax.set_xlabel('Log2 Fold Change', fontsize=12)
    ax.set_ylabel('-Log10(adjusted p-value)', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', alpha=0.6, label='Significant & Large Effect'),
        Patch(facecolor='orange', alpha=0.6, label='Significant'),
        Patch(facecolor='blue', alpha=0.6, label='Large Effect Only'),
        Patch(facecolor='gray', alpha=0.6, label='Not Significant')
    ]
    ax.legend(handles=legend_elements, loc='upper left')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    return fig, ax


def analyze_protein_set(adata, protein_set, set_name="Protein Set",
                        design="~tau_status", **kwargs):
    """
    Analyze a specific set of proteins

    Parameters:
    -----------
    adata : AnnData
        Full dataset
    protein_set : list
        List of proteins to analyze
    set_name : str
        Name of protein set
    design : str
        DESeq2 design formula
    **kwargs : additional arguments for run_pydeseq2_safe

    Returns:
    --------
    dict : Analysis results dictionary
    """
    # Find proteins
    protein_info = find_proteins_in_dataset(adata, protein_set)
    found = protein_info['found']

    if not found:
        print(f"No proteins from {set_name} found in dataset")
        return None

    print(f"\nAnalyzing {set_name}:")
    print(f"  Found: {len(found)}/{len(protein_set)} proteins")

    # Subset data
    protein_names = adata.var['protein_name'] if 'protein_name' in adata.var else adata.var.index
    indices = [i for i, p in enumerate(protein_names) if p in found]
    adata_subset = adata[:, indices].copy()

    # Run analysis
    results = run_pydeseq2_safe(adata_subset, design=design, **kwargs)

    # Calculate summary statistics
    n_sig = (results['padj'] < 0.05).sum()
    n_total = len(results)
    mean_log2fc = results['log2FoldChange'].mean()

    summary = {
        'set_name': set_name,
        'n_proteins_requested': len(protein_set),
        'n_proteins_found': len(found),
        'n_proteins_missing': len(protein_info['missing']),
        'n_significant': n_sig,
        'percent_significant': n_sig / n_total * 100 if n_total > 0 else 0,
        'mean_log2fc': mean_log2fc,
        'results_df': results,
        'missing_proteins': protein_info['missing']
    }

    return summary


def temporal_correlation_analysis(adata, protein_list, time_column='pseudotime'):
    """
    Analyze temporal correlations for proteins

    Parameters:
    -----------
    adata : AnnData
        Annotated data matrix
    protein_list : list
        List of proteins to analyze
    time_column : str
        Column with temporal information

    Returns:
    --------
    pd.DataFrame : Correlation results
    """
    if time_column not in adata.obs:
        print(f"Error: {time_column} not found in metadata")
        return None

    results = []
    protein_names = adata.var['protein_name'] if 'protein_name' in adata.var else adata.var.index

    for protein in protein_list:
        if protein in protein_names.tolist():
            idx = list(protein_names).index(protein)
            expr = adata.X[:, idx]

            # Remove NaN values
            valid_mask = ~(np.isnan(expr) | adata.obs[time_column].isna())

            if valid_mask.sum() > 10:
                # Spearman correlation
                corr, pval = stats.spearmanr(
                    adata.obs.loc[valid_mask, time_column],
                    expr[valid_mask]
                )

                # Linear regression for slope
                from sklearn.linear_model import LinearRegression
                X = adata.obs.loc[valid_mask, time_column].values.reshape(-1, 1)
                y = expr[valid_mask]
                lr = LinearRegression()
                lr.fit(X, y)

                results.append({
                    'protein': protein,
                    'correlation': corr,
                    'p_value': pval,
                    'slope': lr.coef_[0],
                    'intercept': lr.intercept_,
                    'direction': 'increasing' if lr.coef_[0] > 0 else 'decreasing'
                })

    if results:
        results_df = pd.DataFrame(results)
        # Add FDR correction
        results_df['p_adjusted'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
        results_df['significant'] = results_df['p_adjusted'] < 0.05
        return results_df

    return pd.DataFrame()


def evaluate_claim(observed_value, claimed_value, value_type='fold_change',
                  tolerance=0.2):
    """
    Evaluate a scientific claim against observed data

    Parameters:
    -----------
    observed_value : float
        Observed value from analysis
    claimed_value : float
        Claimed value to test
    value_type : str
        Type of value ('fold_change', 'percentage', 'correlation')
    tolerance : float
        Tolerance for considering claim supported (0.2 = 20%)

    Returns:
    --------
    tuple : (verdict, explanation)
    """
    if pd.isna(observed_value) or pd.isna(claimed_value):
        return "UNSURE", "Missing data for comparison"

    # Calculate difference
    if value_type == 'percentage':
        diff = abs(observed_value - claimed_value)
        relative_diff = diff
    else:
        diff = abs(observed_value - claimed_value)
        relative_diff = diff / abs(claimed_value) if claimed_value != 0 else float('inf')

    # Determine verdict
    if relative_diff <= tolerance:
        verdict = "SUPPORTED"
        explanation = f"Observed ({observed_value:.2f}) matches claimed ({claimed_value:.2f}) within {tolerance*100}% tolerance"
    elif relative_diff <= tolerance * 2:
        verdict = "PARTIALLY SUPPORTED"
        explanation = f"Observed ({observed_value:.2f}) differs from claimed ({claimed_value:.2f}) by {relative_diff*100:.1f}%"
    else:
        verdict = "REFUTED"
        explanation = f"Observed ({observed_value:.2f}) substantially differs from claimed ({claimed_value:.2f})"

    return verdict, explanation


# Export key functions
__all__ = [
    'run_pydeseq2_safe',
    'find_proteins_in_dataset',
    'create_volcano_plot',
    'analyze_protein_set',
    'temporal_correlation_analysis',
    'evaluate_claim'
]
```
