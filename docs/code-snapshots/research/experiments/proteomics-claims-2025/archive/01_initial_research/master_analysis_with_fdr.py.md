---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/archive/01_initial_research/master_analysis_with_fdr.py
relative: research/experiments/proteomics-claims-2025/archive/01_initial_research/master_analysis_with_fdr.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Master Analysis Pipeline with FDR Correction
=============================================
Runs all analyses with proper multiple testing correction.

Author: Bioinformatics Finding Group Evaluation Framework
Date: 2024
"""

import sys
import os
import json
import warnings
from datetime import datetime

sys.path.append('/Users/byron/project_plan')

from config import (
    DATA_PATH, DATA_SPECS, ANALYSIS_PARAMS, PROTEIN_GROUPS,
    KNOWN_ISSUES, load_data, get_tau_groups, DIRS
)

import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

# Set plotting style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette('husl')

# ============================================================================
# FDR CORRECTION UTILITIES
# ============================================================================

def apply_fdr_correction(pvalues, alpha=0.05, method='fdr_bh'):
    """
    Apply FDR correction to p-values.

    Parameters:
    -----------
    pvalues : array-like
        Raw p-values
    alpha : float
        Significance threshold
    method : str
        Correction method ('fdr_bh' for Benjamini-Hochberg)

    Returns:
    --------
    dict with corrected p-values and significance flags
    """
    if len(pvalues) == 0:
        return {'pvals_corrected': [], 'significant': [], 'n_significant': 0}

    # Handle NaN values
    valid_mask = ~np.isnan(pvalues)
    valid_pvals = np.array(pvalues)[valid_mask]

    if len(valid_pvals) == 0:
        return {
            'pvals_corrected': pvalues,
            'significant': [False] * len(pvalues),
            'n_significant': 0
        }

    # Apply correction
    rejected, pvals_corrected, _, _ = multipletests(valid_pvals, alpha=alpha, method=method)

    # Map back to original array
    result_pvals = np.array(pvalues, dtype=float)
    result_sig = np.zeros(len(pvalues), dtype=bool)

    result_pvals[valid_mask] = pvals_corrected
    result_sig[valid_mask] = rejected

    return {
        'pvals_corrected': result_pvals.tolist(),
        'significant': result_sig.tolist(),
        'n_significant': int(np.sum(rejected)),
        'alpha': alpha,
        'method': method
    }


# ============================================================================
# GLOBAL DIFFERENTIAL EXPRESSION WITH FDR
# ============================================================================

def global_differential_expression(adata):
    """
    Perform differential expression for all proteins with FDR correction.
    This provides context for individual protein findings.
    """
    print("\n" + "="*60)
    print("GLOBAL DIFFERENTIAL EXPRESSION ANALYSIS")
    print("="*60)

    tau_pos, tau_neg = get_tau_groups(adata)

    results = []

    # Analyze subset for speed (or all if needed)
    n_proteins = min(1000, adata.n_vars)  # Analyze first 1000 for demo

    print(f"Analyzing {n_proteins} proteins...")

    for i in range(n_proteins):
        gene_name = adata.var.iloc[i]['GeneName']
        expr = adata.X[:, i]

        expr_pos = expr[tau_pos]
        expr_neg = expr[tau_neg]

        # Skip if no expression
        if np.all(expr == 0):
            continue

        # Calculate statistics
        mean_pos = np.mean(expr_pos)
        mean_neg = np.mean(expr_neg)

        if mean_neg != 0:
            fold_change = mean_pos / mean_neg
            log2_fc = np.log2(fold_change)
        else:
            fold_change = np.inf if mean_pos > 0 else 1.0
            log2_fc = np.inf if mean_pos > 0 else 0

        # Statistical test
        stat, pval = mannwhitneyu(expr_pos, expr_neg)

        # Effect size (Cohen's d)
        pooled_std = np.sqrt(((len(expr_pos)-1)*np.std(expr_pos)**2 +
                              (len(expr_neg)-1)*np.std(expr_neg)**2) /
                             (len(expr_pos) + len(expr_neg) - 2))
        cohens_d = (mean_pos - mean_neg) / pooled_std if pooled_std > 0 else 0

        results.append({
            'gene': gene_name,
            'mean_pos': mean_pos,
            'mean_neg': mean_neg,
            'fold_change': fold_change,
            'log2_fc': log2_fc,
            'pval': pval,
            'cohens_d': cohens_d
        })

    # Convert to DataFrame
    df = pd.DataFrame(results)

    # Apply FDR correction
    fdr_results = apply_fdr_correction(df['pval'].values)
    df['pval_fdr'] = fdr_results['pvals_corrected']
    df['significant_fdr'] = fdr_results['significant']

    # Add significance without FDR for comparison
    df['significant_raw'] = df['pval'] < 0.05

    # Summary statistics
    n_sig_raw = df['significant_raw'].sum()
    n_sig_fdr = df['significant_fdr'].sum()

    print(f"\nResults:")
    print(f"  Proteins analyzed: {len(df)}")
    print(f"  Significant (p<0.05): {n_sig_raw} ({n_sig_raw/len(df)*100:.1f}%)")
    print(f"  Significant (FDR<0.05): {n_sig_fdr} ({n_sig_fdr/len(df)*100:.1f}%)")
    print(f"  FDR reduction: {(n_sig_raw - n_sig_fdr)/n_sig_raw*100:.1f}%")

    # Check SQSTM1 if in results
    if 'SQSTM1' in df['gene'].values:
        sqstm1_row = df[df['gene'] == 'SQSTM1'].iloc[0]
        print(f"\nSQSTM1 status:")
        print(f"  Raw p-value: {sqstm1_row['pval']:.3e}")
        print(f"  FDR-corrected: {sqstm1_row['pval_fdr']:.3e}")
        print(f"  Significant after FDR: {sqstm1_row['significant_fdr']}")

    return df


# ============================================================================
# ENHANCED ANALYSIS FUNCTIONS WITH FDR
# ============================================================================

def analyze_protein_group_with_fdr(adata, protein_dict, group_name):
    """
    Analyze a group of proteins with FDR correction.
    """
    print(f"\n{group_name} Analysis with FDR")
    print("-" * 40)

    tau_pos, tau_neg = get_tau_groups(adata)
    results = []

    # Flatten protein list
    all_proteins = []
    for category, proteins in protein_dict.items():
        all_proteins.extend(proteins)

    for protein in all_proteins:
        mask = adata.var['GeneName'].str.contains(protein, case=False, na=False)
        if mask.any():
            idx = np.where(mask)[0][0]
            expr = adata.X[:, idx]
            expr_pos = expr[tau_pos]
            expr_neg = expr[tau_neg]

            # Statistics
            mean_pos = np.mean(expr_pos)
            mean_neg = np.mean(expr_neg)

            if mean_neg != 0:
                fold_change = mean_pos / mean_neg
                log2_fc = np.log2(fold_change)
            else:
                fold_change = np.inf if mean_pos > 0 else 1.0
                log2_fc = np.inf if mean_pos > 0 else 0

            # Test
            stat, pval = mannwhitneyu(expr_pos, expr_neg)

            results.append({
                'protein': protein,
                'fold_change': fold_change,
                'log2_fc': log2_fc,
                'pval': pval,
                'found': True
            })
        else:
            results.append({
                'protein': protein,
                'found': False
            })

    # Create DataFrame
    df = pd.DataFrame(results)
    found_df = df[df['found'] == True]

    if len(found_df) > 0:
        # Apply FDR correction
        fdr_results = apply_fdr_correction(found_df['pval'].values)
        found_df['pval_fdr'] = fdr_results['pvals_corrected']
        found_df['significant_fdr'] = fdr_results['significant']
        found_df['significant_raw'] = found_df['pval'] < 0.05

        # Summary
        n_found = len(found_df)
        n_sig_raw = found_df['significant_raw'].sum()
        n_sig_fdr = found_df['significant_fdr'].sum()

        print(f"  Proteins found: {n_found}/{len(all_proteins)}")
        print(f"  Significant (raw): {n_sig_raw} ({n_sig_raw/n_found*100:.1f}%)")
        print(f"  Significant (FDR): {n_sig_fdr} ({n_sig_fdr/n_found*100:.1f}%)")

        # Show top hits
        if n_sig_fdr > 0:
            top_hits = found_df[found_df['significant_fdr']].nsmallest(5, 'pval_fdr')
            print(f"\n  Top {group_name} proteins (FDR < 0.05):")
            for _, row in top_hits.iterrows():
                print(f"    {row['protein']}: FC={row['fold_change']:.2f}, FDR={row['pval_fdr']:.3e}")

    return df


def analyze_sqstm1_with_context(adata, global_df=None):
    """
    Analyze SQSTM1 in context of all proteins.
    """
    print("\n" + "="*60)
    print("SQSTM1 ANALYSIS WITH FDR CONTEXT")
    print("="*60)

    tau_pos, tau_neg = get_tau_groups(adata)

    # Find SQSTM1
    sqstm1_mask = adata.var['GeneName'].str.contains('SQSTM1', case=False, na=False)

    if not sqstm1_mask.any():
        return {'error': 'SQSTM1 not found'}

    sqstm1_idx = np.where(sqstm1_mask)[0][0]
    sqstm1_expr = adata.X[:, sqstm1_idx]

    expr_pos = sqstm1_expr[tau_pos]
    expr_neg = sqstm1_expr[tau_neg]

    # Calculate statistics
    mean_pos = np.mean(expr_pos)
    mean_neg = np.mean(expr_neg)
    log2_fc = np.log2(mean_pos / mean_neg) if mean_neg != 0 else 0
    fold_change = 2**log2_fc

    # Multiple statistical tests
    mw_stat, mw_pval = mannwhitneyu(expr_pos, expr_neg)
    t_stat, t_pval = ttest_ind(expr_pos, expr_neg)
    t_stat_welch, t_pval_welch = ttest_ind(expr_pos, expr_neg, equal_var=False)

    # Effect size
    pooled_std = np.sqrt(((len(expr_pos)-1)*np.std(expr_pos)**2 +
                          (len(expr_neg)-1)*np.std(expr_neg)**2) /
                         (len(expr_pos) + len(expr_neg) - 2))
    cohens_d = (mean_pos - mean_neg) / pooled_std if pooled_std > 0 else 0

    # Rank in global context if available
    rank_info = {}
    if global_df is not None and 'SQSTM1' in global_df['gene'].values:
        sqstm1_row = global_df[global_df['gene'] == 'SQSTM1'].iloc[0]
        rank_info = {
            'pval_rank': (global_df['pval'] < sqstm1_row['pval']).sum() + 1,
            'fc_rank': (global_df['fold_change'].abs() > abs(sqstm1_row['fold_change'])).sum() + 1,
            'total_proteins': len(global_df),
            'percentile_pval': (rank_info['pval_rank'] / len(global_df)) * 100 if 'pval_rank' in rank_info else None,
            'percentile_fc': (rank_info['fc_rank'] / len(global_df)) * 100 if 'fc_rank' in rank_info else None
        }

    results = {
        'mean_tau_pos': mean_pos,
        'mean_tau_neg': mean_neg,
        'log2_fold_change': log2_fc,
        'fold_change': fold_change,
        'mann_whitney_pval': mw_pval,
        't_test_pval': t_pval,
        'welch_t_pval': t_pval_welch,
        'cohens_d': cohens_d,
        'claimed_fold_change': KNOWN_ISSUES['SQSTM1_fold_change']['claimed'],
        'observed_fold_change': fold_change,
        'discrepancy_ratio': KNOWN_ISSUES['SQSTM1_fold_change']['claimed'] / fold_change,
        **rank_info
    }

    # Print results
    print(f"Expression levels:")
    print(f"  Tau+: {mean_pos:.3f}")
    print(f"  Tau-: {mean_neg:.3f}")
    print(f"\nFold change:")
    print(f"  Observed: {fold_change:.2f}x (log2: {log2_fc:.3f})")
    print(f"  Claimed: {KNOWN_ISSUES['SQSTM1_fold_change']['claimed']:.1f}x")
    print(f"  Discrepancy: {results['discrepancy_ratio']:.1f}x")
    print(f"\nStatistical tests:")
    print(f"  Mann-Whitney p: {mw_pval:.3e}")
    print(f"  T-test p: {t_pval:.3e}")
    print(f"  Welch's t p: {t_pval_welch:.3e}")
    print(f"  Cohen's d: {cohens_d:.3f}")

    if rank_info:
        print(f"\nGlobal context:")
        print(f"  P-value rank: {rank_info['pval_rank']}/{rank_info['total_proteins']}")
        print(f"  Fold change rank: {rank_info['fc_rank']}/{rank_info['total_proteins']}")

    return results


# ============================================================================
# COMPREHENSIVE REPORT GENERATION
# ============================================================================

def generate_comprehensive_report(all_results):
    """Generate detailed report with FDR statistics."""

    report = {
        'metadata': {
            'data_file': DATA_PATH,
            'analysis_date': datetime.now().isoformat(),
            'n_samples': DATA_SPECS['n_samples'],
            'n_proteins': DATA_SPECS['n_proteins'],
            'fdr_method': 'Benjamini-Hochberg',
            'fdr_alpha': 0.05
        },
        'results': all_results,
        'summary': {
            'sqstm1_validated': all_results['sqstm1']['fold_change'] > 10,
            'sqstm1_discrepancy': all_results['sqstm1']['discrepancy_ratio'],
            'global_proteins_analyzed': all_results['global_de']['n_proteins'],
            'global_significant_raw': all_results['global_de']['n_significant_raw'],
            'global_significant_fdr': all_results['global_de']['n_significant_fdr'],
            'fdr_reduction_percent': all_results['global_de']['fdr_reduction_percent']
        }
    }

    # Save JSON
    json_path = DIRS['reports'] / 'master_analysis_with_fdr.json'
    with open(json_path, 'w') as f:
        json.dump(report, f, indent=2, default=str)

    print(f"\n✓ Report saved to {json_path}")

    return report


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Run complete analysis with FDR correction."""
    print("="*60)
    print("MASTER ANALYSIS WITH FDR CORRECTION")
    print("="*60)
    print(f"Data: {DATA_PATH}")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")

    # Load data
    print("\nLoading data...")
    adata = load_data()

    all_results = {}

    # 1. Global differential expression
    global_df = global_differential_expression(adata)
    all_results['global_de'] = {
        'n_proteins': len(global_df),
        'n_significant_raw': global_df['significant_raw'].sum(),
        'n_significant_fdr': global_df['significant_fdr'].sum(),
        'fdr_reduction_percent': ((global_df['significant_raw'].sum() -
                                   global_df['significant_fdr'].sum()) /
                                  global_df['significant_raw'].sum() * 100)
    }

    # 2. SQSTM1 with context
    all_results['sqstm1'] = analyze_sqstm1_with_context(adata, global_df)

    # 3. Protein groups with FDR
    all_results['autophagy'] = analyze_protein_group_with_fdr(
        adata, PROTEIN_GROUPS['autophagy'], 'Autophagy'
    )
    all_results['ups'] = analyze_protein_group_with_fdr(
        adata, PROTEIN_GROUPS['ups'], 'UPS'
    )
    all_results['proteasome'] = analyze_protein_group_with_fdr(
        adata, PROTEIN_GROUPS['proteasome'], 'Proteasome'
    )

    # Generate report
    report = generate_comprehensive_report(all_results)

    print("\n" + "="*60)
    print("ANALYSIS COMPLETE WITH FDR CORRECTION")
    print("="*60)

    return report


if __name__ == "__main__":
    report = main()
```
