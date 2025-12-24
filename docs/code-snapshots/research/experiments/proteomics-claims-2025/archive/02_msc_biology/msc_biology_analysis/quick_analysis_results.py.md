---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/archive/02_msc_biology/msc_biology_analysis/quick_analysis_results.py
relative: research/experiments/proteomics-claims-2025/archive/02_msc_biology/msc_biology_analysis/quick_analysis_results.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Quick MSc Biology Analysis - Focus on key results without plots
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import json
import warnings
warnings.filterwarnings('ignore')

print("="*60)
print("MSC BIOLOGY ANALYSIS - KEY RESULTS")
print("="*60)

# Load data
print("\nLoading data...")
adata = sc.read_h5ad('/Users/byron/project_plan/data/pool_processed_v2.h5ad')
print(f"✓ Loaded {adata.n_obs} samples × {adata.n_vars} proteins")

# Get tau groups
tau_pos = adata.obs['TauStatus'] == 'positive'
tau_neg = adata.obs['TauStatus'] == 'negative'
print(f"\nTau+ neurons: {sum(tau_pos)}")
print(f"Tau- neurons: {sum(tau_neg)}")

print("\n" + "="*60)
print("KEY FINDING 1: SQSTM1 UPREGULATION")
print("="*60)

# Find and analyze SQSTM1
sqstm1_mask = adata.var['GeneName'].str.contains('SQSTM1', case=False, na=False)
if sqstm1_mask.any():
    sqstm1_idx = np.where(sqstm1_mask)[0][0]
    sqstm1_expression = adata.X[:, sqstm1_idx]

    sqstm1_tau_pos = sqstm1_expression[tau_pos]
    sqstm1_tau_neg = sqstm1_expression[tau_neg]

    mean_pos = np.mean(sqstm1_tau_pos)
    mean_neg = np.mean(sqstm1_tau_neg)

    log2_fc = np.log2(mean_pos / mean_neg) if mean_neg != 0 else 0
    fold_change = 2**log2_fc

    stat, pval = mannwhitneyu(sqstm1_tau_pos, sqstm1_tau_neg)

    print(f"\nSQSTM1 Results:")
    print(f"  Mean in Tau+: {mean_pos:.3f}")
    print(f"  Mean in Tau-: {mean_neg:.3f}")
    print(f"  Log2 Fold Change: {log2_fc:.3f}")
    print(f"  Fold Change: {fold_change:.1f}x")
    print(f"  P-value: {pval:.3e}")
    print(f"\n  {'✅ CONFIRMED' if fold_change > 10 else '❌'}: Paper claimed 10.7-fold")

print("\n" + "="*60)
print("KEY FINDING 2: AUTOPHAGY VS UPS")
print("="*60)

# Quick analysis of key proteins
autophagy_proteins = ['SQSTM1', 'NBR1', 'MAP1LC3B', 'BECN1', 'ATG5', 'ATG7']
ups_proteins = ['PSMA1', 'PSMB5', 'PSMD1', 'UBA1']

def quick_analyze(protein_list, name):
    significant = 0
    found = 0
    for protein in protein_list:
        mask = adata.var['GeneName'].str.contains(protein, case=False, na=False)
        if mask.any():
            found += 1
            idx = np.where(mask)[0][0]
            expr = adata.X[:, idx]
            _, pval = mannwhitneyu(expr[tau_pos], expr[tau_neg])
            if pval < 0.05:
                significant += 1
    return found, significant

auto_found, auto_sig = quick_analyze(autophagy_proteins, 'Autophagy')
ups_found, ups_sig = quick_analyze(ups_proteins, 'UPS')

print(f"\nAutophagy: {auto_sig}/{auto_found} proteins significantly changed")
print(f"UPS: {ups_sig}/{ups_found} proteins significantly changed")

if auto_found > 0 and ups_found > 0:
    auto_percent = auto_sig/auto_found * 100
    ups_percent = ups_sig/ups_found * 100
    print(f"\nAutophagy disruption: {auto_percent:.1f}%")
    print(f"UPS disruption: {ups_percent:.1f}%")

    if auto_percent > ups_percent * 1.5:
        print("✅ CONFIRMED: Autophagy specifically disrupted")
    else:
        print("❌ Similar disruption levels")

print("\n" + "="*60)
print("KEY FINDING 3: PROTEASOME & V-ATPASE")
print("="*60)

# Count proteasome and V-ATPase proteins
proteasome_count = 0
vatpase_count = 0

proteasome_prefixes = ['PSMA', 'PSMB', 'PSMC', 'PSMD']
for prefix in proteasome_prefixes:
    mask = adata.var['GeneName'].str.contains(f'^{prefix}\\d', regex=True, na=False)
    proteasome_count += mask.sum()

vatpase_prefixes = ['ATP6V0', 'ATP6V1']
for prefix in vatpase_prefixes:
    mask = adata.var['GeneName'].str.contains(prefix, na=False)
    vatpase_count += mask.sum()

print(f"\nProteasome subunits found: {proteasome_count}")
print(f"V-ATPase subunits found: {vatpase_count}")

# Check pseudotime if available
if 'pseudotime' in adata.obs.columns:
    pseudotime = adata.obs['pseudotime'].values
    print(f"\nPseudotime range: {pseudotime.min():.3f} to {pseudotime.max():.3f}")
    print("✓ Can test sequential failure hypothesis")

print("\n" + "="*60)
print("KEY FINDING 4: MITOCHONDRIAL MARKERS")
print("="*60)

# Check key mitochondrial proteins
mito_proteins = ['VDAC1', 'CYCS', 'COX4I1', 'ATP5A1']
mito_found = 0
mito_changed = 0

for protein in mito_proteins:
    mask = adata.var['GeneName'].str.contains(protein, case=False, na=False)
    if mask.any():
        mito_found += 1
        idx = np.where(mask)[0][0]
        expr = adata.X[:, idx]
        _, pval = mannwhitneyu(expr[tau_pos], expr[tau_neg])
        status = "changed" if pval < 0.05 else "stable"
        print(f"  {protein}: {status} (p={pval:.3e})")
        if pval < 0.05:
            mito_changed += 1

print(f"\nMitochondrial proteins: {mito_changed}/{mito_found} significantly changed")

print("\n" + "="*60)
print("FINAL SUMMARY")
print("="*60)

summary = {
    'Data_file': 'pool_processed_v2.h5ad',
    'Total_samples': int(adata.n_obs),
    'Total_proteins': int(adata.n_vars),
    'Tau_positive': int(sum(tau_pos)),
    'Tau_negative': int(sum(tau_neg)),
    'SQSTM1_fold_change': round(fold_change, 1) if 'fold_change' in locals() else 'Not found',
    'SQSTM1_validated': bool(fold_change > 10) if 'fold_change' in locals() else False,
    'Autophagy_disrupted': f"{auto_sig}/{auto_found}" if 'auto_found' in locals() else 'N/A',
    'UPS_disrupted': f"{ups_sig}/{ups_found}" if 'ups_found' in locals() else 'N/A',
    'Proteasome_proteins': proteasome_count,
    'VATPase_proteins': vatpase_count,
    'Mitochondrial_changed': f"{mito_changed}/{mito_found}"
}

# Save results
with open('/Users/byron/project_plan/msc_biology_analysis/results/quick_analysis_results.json', 'w') as f:
    json.dump(summary, f, indent=2)

print("\nResults Summary:")
for key, value in summary.items():
    print(f"  {key}: {value}")

print("\n✅ Analysis complete!")
print("Results saved to: msc_biology_analysis/results/quick_analysis_results.json")
```
