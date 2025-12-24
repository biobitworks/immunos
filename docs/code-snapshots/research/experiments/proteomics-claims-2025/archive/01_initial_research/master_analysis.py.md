---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/archive/01_initial_research/master_analysis.py
relative: research/experiments/proteomics-claims-2025/archive/01_initial_research/master_analysis.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Master Analysis Pipeline for Proteomics Project
================================================
Runs all analyses using pool_processed_v2.h5ad and generates comprehensive report.

Author: Bioinformatics Finding Group Evaluation Framework
Date: 2024
"""

import sys
import os
import json
import warnings
from datetime import datetime

# Add project root to path
sys.path.append('/Users/byron/project_plan')

from config import (
    DATA_PATH, DATA_SPECS, ANALYSIS_PARAMS, PROTEIN_GROUPS,
    KNOWN_ISSUES, load_data, get_tau_groups, DIRS
)

import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

# Set plotting style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette('husl')

# ============================================================================
# MAIN ANALYSIS FUNCTIONS
# ============================================================================

def analyze_sqstm1(adata):
    """Analyze SQSTM1 upregulation - the key discrepancy."""
    print("\n" + "="*60)
    print("ANALYZING SQSTM1 UPREGULATION")
    print("="*60)

    # Get tau groups
    tau_pos, tau_neg = get_tau_groups(adata)

    # Find SQSTM1
    sqstm1_mask = adata.var['GeneName'].str.contains('SQSTM1', case=False, na=False)

    if not sqstm1_mask.any():
        return {'error': 'SQSTM1 not found'}

    sqstm1_idx = np.where(sqstm1_mask)[0][0]
    sqstm1_expr = adata.X[:, sqstm1_idx]

    # Calculate expression in each group
    expr_pos = sqstm1_expr[tau_pos]
    expr_neg = sqstm1_expr[tau_neg]

    # Statistics
    mean_pos = np.mean(expr_pos)
    mean_neg = np.mean(expr_neg)
    log2_fc = np.log2(mean_pos / mean_neg) if mean_neg != 0 else 0
    fold_change = 2**log2_fc

    # Statistical test
    stat, pval = mannwhitneyu(expr_pos, expr_neg)

    results = {
        'mean_tau_pos': mean_pos,
        'mean_tau_neg': mean_neg,
        'log2_fold_change': log2_fc,
        'fold_change': fold_change,
        'p_value': pval,
        'claimed_fold_change': KNOWN_ISSUES['SQSTM1_fold_change']['claimed'],
        'observed_fold_change': fold_change,
        'discrepancy_ratio': KNOWN_ISSUES['SQSTM1_fold_change']['claimed'] / fold_change,
        'validated': fold_change > 10
    }

    print(f"Claimed: {KNOWN_ISSUES['SQSTM1_fold_change']['claimed']:.1f}x")
    print(f"Observed: {fold_change:.1f}x")
    print(f"Log2 FC: {log2_fc:.3f}")
    print(f"P-value: {pval:.3e}")
    print(f"Status: {'✅ VALIDATED' if results['validated'] else '❌ NOT VALIDATED'}")

    return results


def analyze_autophagy_vs_ups(adata):
    """Compare autophagy and UPS protein changes."""
    print("\n" + "="*60)
    print("ANALYZING AUTOPHAGY VS UPS")
    print("="*60)

    tau_pos, tau_neg = get_tau_groups(adata)
    results = {'autophagy': {}, 'ups': {}}

    # Analyze autophagy proteins
    autophagy_proteins = []
    for category, proteins in PROTEIN_GROUPS['autophagy'].items():
        autophagy_proteins.extend(proteins)

    # Analyze UPS proteins
    ups_proteins = []
    for category, proteins in PROTEIN_GROUPS['ups'].items():
        ups_proteins.extend(proteins[:20])  # Limit to first 20 for speed

    def analyze_group(protein_list, name):
        found = 0
        significant = 0
        sig_proteins = []

        for protein in protein_list:
            mask = adata.var['GeneName'].str.contains(protein, case=False, na=False)
            if mask.any():
                found += 1
                idx = np.where(mask)[0][0]
                expr = adata.X[:, idx]
                expr_pos = expr[tau_pos]
                expr_neg = expr[tau_neg]

                _, pval = mannwhitneyu(expr_pos, expr_neg)

                if pval < ANALYSIS_PARAMS['p_value_threshold']:
                    significant += 1
                    sig_proteins.append(protein)

        return {
            'total_analyzed': len(protein_list),
            'found': found,
            'significant': significant,
            'percent_significant': (significant/found * 100) if found > 0 else 0,
            'significant_proteins': sig_proteins
        }

    results['autophagy'] = analyze_group(autophagy_proteins, 'Autophagy')
    results['ups'] = analyze_group(ups_proteins, 'UPS')

    # Determine if autophagy is specifically disrupted
    auto_pct = results['autophagy']['percent_significant']
    ups_pct = results['ups']['percent_significant']
    results['autophagy_specific'] = auto_pct > ups_pct * 1.5

    print(f"Autophagy: {results['autophagy']['significant']}/{results['autophagy']['found']} ({auto_pct:.1f}%)")
    print(f"UPS: {results['ups']['significant']}/{results['ups']['found']} ({ups_pct:.1f}%)")
    print(f"Autophagy-specific: {'✅ YES' if results['autophagy_specific'] else '❌ NO'}")

    return results


def analyze_proteasome_vatpase(adata):
    """Count and analyze proteasome and V-ATPase proteins."""
    print("\n" + "="*60)
    print("ANALYZING PROTEASOME & V-ATPASE")
    print("="*60)

    results = {'proteasome': {}, 'vatpase': {}}

    # Count proteasome proteins
    proteasome_count = 0
    for category, proteins in PROTEIN_GROUPS['proteasome'].items():
        for protein in proteins:
            mask = adata.var['GeneName'].str.contains(protein, case=False, na=False)
            if mask.any():
                proteasome_count += 1

    # Count V-ATPase proteins
    vatpase_count = 0
    for category, proteins in PROTEIN_GROUPS['vatpase'].items():
        for protein in proteins:
            mask = adata.var['GeneName'].str.contains(protein, case=False, na=False)
            if mask.any():
                vatpase_count += 1

    results['proteasome']['count'] = proteasome_count
    results['vatpase']['count'] = vatpase_count

    # Check pseudotime for sequential failure analysis
    if 'pseudotime' in adata.obs.columns:
        pseudotime = adata.obs['pseudotime'].values
        results['pseudotime_available'] = True
        results['pseudotime_range'] = [float(pseudotime.min()), float(pseudotime.max())]
    else:
        results['pseudotime_available'] = False

    print(f"Proteasome proteins: {proteasome_count}")
    print(f"V-ATPase proteins: {vatpase_count}")
    print(f"Pseudotime available: {results.get('pseudotime_available', False)}")

    return results


def analyze_mitochondrial(adata):
    """Analyze mitochondrial protein changes."""
    print("\n" + "="*60)
    print("ANALYZING MITOCHONDRIAL PROTEINS")
    print("="*60)

    tau_pos, tau_neg = get_tau_groups(adata)
    results = {}

    key_mito_proteins = ['VDAC1', 'CYCS', 'COX4I1', 'ATP5A1', 'PINK1', 'TOMM20']

    for protein in key_mito_proteins:
        mask = adata.var['GeneName'].str.contains(protein, case=False, na=False)
        if mask.any():
            idx = np.where(mask)[0][0]
            expr = adata.X[:, idx]
            expr_pos = expr[tau_pos]
            expr_neg = expr[tau_neg]

            _, pval = mannwhitneyu(expr_pos, expr_neg)

            results[protein] = {
                'found': True,
                'p_value': float(pval),
                'significant': pval < ANALYSIS_PARAMS['p_value_threshold']
            }
        else:
            results[protein] = {'found': False}

    # Count significant changes
    sig_count = sum(1 for p in results.values() if p.get('significant', False))
    total_found = sum(1 for p in results.values() if p.get('found', False))

    print(f"Mitochondrial proteins: {sig_count}/{total_found} significantly changed")
    for protein, data in results.items():
        if data.get('found'):
            status = '✅' if data.get('significant') else '❌'
            print(f"  {protein}: {status} (p={data.get('p_value', 0):.3e})")

    return results


def generate_summary_report(all_results):
    """Generate comprehensive summary report."""
    print("\n" + "="*60)
    print("GENERATING SUMMARY REPORT")
    print("="*60)

    report = {
        'metadata': {
            'data_file': DATA_PATH,
            'analysis_date': datetime.now().isoformat(),
            'n_samples': DATA_SPECS['n_samples'],
            'n_proteins': DATA_SPECS['n_proteins'],
            'n_tau_positive': DATA_SPECS['n_tau_positive'],
            'n_tau_negative': DATA_SPECS['n_tau_negative']
        },
        'results': all_results,
        'key_findings': {
            'sqstm1_validated': all_results['sqstm1']['validated'],
            'sqstm1_discrepancy': all_results['sqstm1']['discrepancy_ratio'],
            'autophagy_specific': all_results['autophagy_ups']['autophagy_specific'],
            'proteasome_count': all_results['proteasome_vatpase']['proteasome']['count'],
            'vatpase_count': all_results['proteasome_vatpase']['vatpase']['count']
        },
        'known_issues': KNOWN_ISSUES
    }

    # Save JSON report
    json_path = DIRS['reports'] / 'master_analysis_results.json'
    with open(json_path, 'w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f"✓ JSON report saved to {json_path}")

    # Create markdown report
    md_content = f"""# Master Analysis Report - pool_processed_v2.h5ad
Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}

## Dataset Overview
- **Samples**: {DATA_SPECS['n_samples']} ({DATA_SPECS['n_tau_positive']} tau+, {DATA_SPECS['n_tau_negative']} tau-)
- **Proteins**: {DATA_SPECS['n_proteins']}
- **Data file**: {DATA_PATH}

## Key Findings

### 1. SQSTM1 Upregulation {'✅' if all_results['sqstm1']['validated'] else '❌'}
- **Claimed**: {KNOWN_ISSUES['SQSTM1_fold_change']['claimed']:.1f}-fold
- **Observed**: {all_results['sqstm1']['observed_fold_change']:.1f}-fold
- **Discrepancy**: {all_results['sqstm1']['discrepancy_ratio']:.1f}x difference
- **P-value**: {all_results['sqstm1']['p_value']:.3e}

### 2. Autophagy vs UPS {'✅' if all_results['autophagy_ups']['autophagy_specific'] else '❌'}
- **Autophagy**: {all_results['autophagy_ups']['autophagy']['percent_significant']:.1f}% disrupted
- **UPS**: {all_results['autophagy_ups']['ups']['percent_significant']:.1f}% disrupted
- **Conclusion**: {'Autophagy-specific dysfunction' if all_results['autophagy_ups']['autophagy_specific'] else 'Both systems affected'}

### 3. Proteasome & V-ATPase
- **Proteasome proteins**: {all_results['proteasome_vatpase']['proteasome']['count']}
- **V-ATPase proteins**: {all_results['proteasome_vatpase']['vatpase']['count']}
- **Sequential failure testable**: {'Yes' if all_results['proteasome_vatpase'].get('pseudotime_available') else 'No'}

## Known Issues
- SQSTM1 fold change discrepancy: {KNOWN_ISSUES['SQSTM1_fold_change']['note']}
- Proteins with semicolons: {KNOWN_ISSUES['semicolon_proteins']['affected']} ({KNOWN_ISSUES['semicolon_proteins']['percentage']}%)

## Recommendations
1. Investigate SQSTM1 normalization methods
2. Apply multiple testing correction consistently
3. Validate findings with additional datasets
"""

    md_path = DIRS['reports'] / 'master_analysis_report.md'
    with open(md_path, 'w') as f:
        f.write(md_content)
    print(f"✓ Markdown report saved to {md_path}")

    return report


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Run complete analysis pipeline."""
    print("="*60)
    print("MASTER ANALYSIS PIPELINE")
    print("="*60)
    print(f"Data: {DATA_PATH}")
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")

    # Load data
    print("\nLoading data...")
    adata = load_data()

    # Run all analyses
    all_results = {}

    # 1. SQSTM1 analysis (the main discrepancy)
    all_results['sqstm1'] = analyze_sqstm1(adata)

    # 2. Autophagy vs UPS
    all_results['autophagy_ups'] = analyze_autophagy_vs_ups(adata)

    # 3. Proteasome & V-ATPase
    all_results['proteasome_vatpase'] = analyze_proteasome_vatpase(adata)

    # 4. Mitochondrial proteins
    all_results['mitochondrial'] = analyze_mitochondrial(adata)

    # Generate reports
    report = generate_summary_report(all_results)

    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print("\nKey Results:")
    print(f"  SQSTM1: {'✅ VALIDATED' if all_results['sqstm1']['validated'] else '❌ NOT VALIDATED'}")
    print(f"  Autophagy-specific: {'✅ YES' if all_results['autophagy_ups']['autophagy_specific'] else '❌ NO'}")
    print(f"  Proteins found: {all_results['proteasome_vatpase']['proteasome']['count']} proteasome, {all_results['proteasome_vatpase']['vatpase']['count']} V-ATPase")

    return report


if __name__ == "__main__":
    report = main()
```
