# Autophagy vs UPS Analysis - Simplified Version
## Testing: "Autophagy and UPS show differential dysfunction patterns"

**🚀 Simple approach**: Compare autophagy and UPS pathway dysfunction patterns!

---

## Setup & Data Loading

```python
# 🔧 Setup
import os
IN_COLAB = 'google.colab' in str(get_ipython())

if IN_COLAB:
    !pip install -q pertpy scanpy pandas numpy matplotlib seaborn scipy statsmodels
    from google.colab import files
    print("📁 Upload your pool_processed_v2.h5ad file:")
    uploaded = files.upload()
    data_file = list(uploaded.keys())[0]
else:
    data_file = "pool_processed_v2.h5ad"

import pertpy as pt
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

print("✅ Setup complete!")
```

## Load Data & Define Autophagy vs UPS Pathways

```python
# 📊 Load data
adata = sc.read_h5ad(data_file)
print(f"Data loaded: {adata.shape}")

# 🏷️ Setup tau status
if 'TauStatus' in adata.obs.columns:
    adata.obs['tau_status'] = adata.obs['TauStatus']
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)

print(f"Tau+: {adata.obs['tau_positive'].sum()}, Tau-: {(adata.obs['tau_positive']==0).sum()}")

# 🔄 Proteostasis Pathways
proteostasis_pathways = {
    'Autophagy_Initiation': [
        'ULK1',       # Unc-51 like autophagy activating kinase 1
        'ULK2',       # Unc-51 like autophagy activating kinase 2
        'ATG13',      # Autophagy related 13
        'RB1CC1',     # RB1 inducible coiled-coil 1 (FIP200)
        'BECN1',      # Beclin 1
        'PIK3C3',     # Phosphatidylinositol 3-kinase catalytic subunit type 3
        'PIK3R4'      # Phosphoinositide-3-kinase regulatory subunit 4
    ],
    'Autophagy_Elongation': [
        'ATG3',       # Autophagy related 3
        'ATG5',       # Autophagy related 5
        'ATG7',       # Autophagy related 7
        'ATG12',      # Autophagy related 12
        'ATG16L1',    # Autophagy related 16 like 1
        'MAP1LC3A',   # Microtubule associated protein 1 light chain 3 alpha (LC3A)
        'MAP1LC3B',   # LC3B
        'GABARAP'     # GABA type A receptor associated protein
    ],
    'Autophagy_Maturation': [
        'SQSTM1',     # Sequestosome 1 (p62)
        'NBR1',       # Neighbor of BRCA1 gene 1
        'OPTN',       # Optineurin
        'CALCOCO2',   # Calcium binding and coiled-coil domain 2 (NDP52)
        'STX17',      # Syntaxin 17
        'SNAP29',     # Synaptosome associated protein 29
        'VAMP8'       # Vesicle associated membrane protein 8
    ],
    'UPS_Core': [
        'PSMA1',      # Proteasome 20S subunit alpha 1
        'PSMA2',      # Proteasome 20S subunit alpha 2
        'PSMA7',      # Proteasome 20S subunit alpha 7
        'PSMB1',      # Proteasome 20S subunit beta 1
        'PSMB5',      # Proteasome 20S subunit beta 5
        'PSMB7',      # Proteasome 20S subunit beta 7
        'PSMC1',      # Proteasome 26S subunit, ATPase 1
        'PSMC2'       # Proteasome 26S subunit, ATPase 2
    ],
    'UPS_Regulatory': [
        'PSMD1',      # Proteasome 26S subunit, non-ATPase 1
        'PSMD2',      # Proteasome 26S subunit, non-ATPase 2
        'PSMD4',      # Proteasome 26S subunit, non-ATPase 4
        'PSMD7',      # Proteasome 26S subunit, non-ATPase 7
        'PSMD11',     # Proteasome 26S subunit, non-ATPase 11
        'PSMD12',     # Proteasome 26S subunit, non-ATPase 12
        'PSMD14'      # Proteasome 26S subunit, non-ATPase 14
    ],
    'UPS_Ubiquitination': [
        'UBA1',       # Ubiquitin like modifier activating enzyme 1
        'UBE2D1',     # Ubiquitin conjugating enzyme E2 D1
        'UBE2D3',     # Ubiquitin conjugating enzyme E2 D3
        'UBE2K',      # Ubiquitin conjugating enzyme E2 K
        'STUB1',      # STIP1 homology and U-box containing protein 1 (CHIP)
        'TRIM32',     # Tripartite motif containing 32
        'RNF2'        # Ring finger protein 2
    ]
}

total_proteins = sum(len(proteins) for proteins in proteostasis_pathways.values())
print(f"🎯 Testing {total_proteins} proteostasis pathway proteins across {len(proteostasis_pathways)} categories")

# Group by system
autophagy_categories = ['Autophagy_Initiation', 'Autophagy_Elongation', 'Autophagy_Maturation']
ups_categories = ['UPS_Core', 'UPS_Regulatory', 'UPS_Ubiquitination']
```

## Find & Analyze Proteostasis Pathways

```python
# 🔍 Find available proteins and analyze each pathway
protein_names = list(adata.var_names)
pathway_results = []

for pathway_name, proteins in proteostasis_pathways.items():
    print(f"\n📊 Analyzing {pathway_name}...")

    # Find available proteins
    found_proteins = [p for p in proteins if p in protein_names]
    print(f"  Found: {len(found_proteins)}/{len(proteins)} proteins")

    if len(found_proteins) >= 1:  # Need at least 1 protein for analysis
        # Get expression data for each protein
        protein_results = []

        for protein in found_proteins:
            protein_idx = protein_names.index(protein)
            expr = adata.X[:, protein_idx]

            # Split by tau status
            tau_pos_expr = expr[adata.obs['tau_positive'] == 1]
            tau_neg_expr = expr[adata.obs['tau_positive'] == 0]

            # Calculate statistics
            mean_pos = np.mean(tau_pos_expr)
            mean_neg = np.mean(tau_neg_expr)
            log2fc = mean_pos - mean_neg

            # T-test
            t_stat, p_val = stats.ttest_ind(tau_pos_expr, tau_neg_expr)

            protein_results.append({
                'pathway': pathway_name,
                'protein': protein,
                'log2FC': log2fc,
                'p_value': p_val,
                'mean_tau_pos': mean_pos,
                'mean_tau_neg': mean_neg,
                'system': 'Autophagy' if pathway_name in autophagy_categories else 'UPS',
                'changed': abs(log2fc) > 0.2
            })

        # Pathway-level statistics
        log2fcs = [r['log2FC'] for r in protein_results]
        mean_log2fc = np.mean(log2fcs)
        std_log2fc = np.std(log2fcs)
        changed_count = sum(r['changed'] for r in protein_results)
        changed_pct = changed_count / len(protein_results) * 100

        # Add to pathway results
        pathway_results.append({
            'pathway': pathway_name,
            'proteins_found': len(found_proteins),
            'proteins_total': len(proteins),
            'mean_log2FC': mean_log2fc,
            'std_log2FC': std_log2fc,
            'changed_count': changed_count,
            'changed_pct': changed_pct,
            'system': 'Autophagy' if pathway_name in autophagy_categories else 'UPS',
            'proteins_tested': protein_results
        })

        print(f"  Mean log2FC: {mean_log2fc:.3f} (±{std_log2fc:.3f})")
        print(f"  Changed: {changed_count}/{len(found_proteins)} ({changed_pct:.1f}%)")

# Create combined DataFrame for all proteins
all_proteins_df = []
for pathway_result in pathway_results:
    all_proteins_df.extend(pathway_result['proteins_tested'])

df = pd.DataFrame(all_proteins_df)
if len(df) > 0:
    # Add FDR correction
    from statsmodels.stats.multitest import fdrcorrection
    df['p_adjusted'] = fdrcorrection(df['p_value'])[1]
    df['significant'] = df['p_adjusted'] < 0.05

print(f"\n✅ Analysis complete: {len(df)} proteins tested across {len(pathway_results)} pathways")
```

## Compare Autophagy vs UPS Systems

```python
if len(df) > 0:
    # 📊 System-level comparison
    autophagy_df = df[df['system'] == 'Autophagy']
    ups_df = df[df['system'] == 'UPS']

    print(f"\n🔄 SYSTEM COMPARISON:")
    print(f"Autophagy proteins: {len(autophagy_df)}")
    print(f"UPS proteins: {len(ups_df)}")

    if len(autophagy_df) > 0:
        autophagy_mean_fc = autophagy_df['log2FC'].mean()
        autophagy_changed_pct = (autophagy_df['changed'].sum() / len(autophagy_df)) * 100
        autophagy_sig_pct = (autophagy_df['significant'].sum() / len(autophagy_df)) * 100
        print(f"Autophagy - Mean log2FC: {autophagy_mean_fc:.3f}, Changed: {autophagy_changed_pct:.1f}%, Significant: {autophagy_sig_pct:.1f}%")

    if len(ups_df) > 0:
        ups_mean_fc = ups_df['log2FC'].mean()
        ups_changed_pct = (ups_df['changed'].sum() / len(ups_df)) * 100
        ups_sig_pct = (ups_df['significant'].sum() / len(ups_df)) * 100
        print(f"UPS - Mean log2FC: {ups_mean_fc:.3f}, Changed: {ups_changed_pct:.1f}%, Significant: {ups_sig_pct:.1f}%")

    # Statistical comparison between systems
    if len(autophagy_df) > 0 and len(ups_df) > 0:
        # Compare log2FC distributions
        t_stat_systems, p_val_systems = stats.ttest_ind(autophagy_df['log2FC'], ups_df['log2FC'])
        print(f"System comparison p-value: {p_val_systems:.4f}")

        # Effect size (Cohen's d)
        pooled_std = np.sqrt(((len(autophagy_df)-1)*autophagy_df['log2FC'].var() +
                            (len(ups_df)-1)*ups_df['log2FC'].var()) /
                           (len(autophagy_df) + len(ups_df) - 2))
        cohens_d = (autophagy_mean_fc - ups_mean_fc) / pooled_std
        print(f"Effect size (Cohen's d): {cohens_d:.3f}")
```

## Visualize Autophagy vs UPS

```python
if len(df) > 0:
    # 📊 Create comprehensive visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

    # 1. System comparison boxplot
    systems = ['Autophagy', 'UPS']
    system_data = [df[df['system'] == system]['log2FC'] for system in systems]
    box_plot = ax1.boxplot(system_data, labels=systems, patch_artist=True)
    box_plot['boxes'][0].set_facecolor('lightblue')
    box_plot['boxes'][1].set_facecolor('lightcoral')
    ax1.set_ylabel('Log2 Fold Change')
    ax1.set_title('Autophagy vs UPS Expression Changes')
    ax1.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax1.grid(True, alpha=0.3)

    # 2. Pathway-level comparison
    pathway_names = [pr['pathway'] for pr in pathway_results]
    pathway_fcs = [pr['mean_log2FC'] for pr in pathway_results]
    colors = ['lightblue' if pr['system'] == 'Autophagy' else 'lightcoral' for pr in pathway_results]

    bars = ax2.bar(range(len(pathway_names)), pathway_fcs, color=colors, alpha=0.7)
    ax2.set_xticks(range(len(pathway_names)))
    ax2.set_xticklabels(pathway_names, rotation=45, ha='right')
    ax2.set_ylabel('Mean Log2 Fold Change')
    ax2.set_title('Pathway-Specific Changes')
    ax2.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax2.grid(True, alpha=0.3)

    # Add legend
    ax2.scatter([], [], c='lightblue', alpha=0.7, s=60, label='Autophagy')
    ax2.scatter([], [], c='lightcoral', alpha=0.7, s=60, label='UPS')
    ax2.legend()

    # 3. Volcano plot colored by system
    colors = []
    for _, row in df.iterrows():
        if row['system'] == 'Autophagy':
            colors.append('blue')
        else:
            colors.append('red')

    ax3.scatter(df['log2FC'], -np.log10(df['p_value']), c=colors, alpha=0.6, s=40)
    ax3.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    ax3.axvline(0.2, color='green', linestyle='--', alpha=0.5)
    ax3.axvline(-0.2, color='green', linestyle='--', alpha=0.5)
    ax3.axvline(0, color='black', linestyle='-', alpha=0.3)
    ax3.set_xlabel('Log2 Fold Change (Tau+ vs Tau-)')
    ax3.set_ylabel('-Log10(p-value)')
    ax3.set_title('Proteostasis Systems (Blue=Autophagy, Red=UPS)')
    ax3.grid(True, alpha=0.3)

    # 4. System dysfunction comparison
    if len(autophagy_df) > 0 and len(ups_df) > 0:
        system_stats = [
            ['Autophagy', autophagy_changed_pct, autophagy_sig_pct],
            ['UPS', ups_changed_pct, ups_sig_pct]
        ]

        x_pos = np.arange(len(system_stats))
        changed_vals = [stat[1] for stat in system_stats]
        sig_vals = [stat[2] for stat in system_stats]

        width = 0.35
        ax4.bar(x_pos - width/2, changed_vals, width, label='% Changed', color='orange', alpha=0.7)
        ax4.bar(x_pos + width/2, sig_vals, width, label='% Significant', color='purple', alpha=0.7)

        ax4.set_xlabel('Proteostasis System')
        ax4.set_ylabel('Percentage')
        ax4.set_title('System Dysfunction Comparison')
        ax4.set_xticks(x_pos)
        ax4.set_xticklabels([stat[0] for stat in system_stats])
        ax4.legend()
        ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
```

## Evaluate Claim

```python
# 🎯 CLAIM EVALUATION
print("\n" + "="*65)
print("🎯 CLAIM EVALUATION")
print("="*65)
print("Claim: 'Autophagy and UPS show differential dysfunction patterns'")
print()

if len(df) > 0 and len(autophagy_df) > 0 and len(ups_df) > 0:
    # Differential analysis
    fc_difference = abs(autophagy_mean_fc - ups_mean_fc)
    changed_difference = abs(autophagy_changed_pct - ups_changed_pct)
    sig_difference = abs(autophagy_sig_pct - ups_sig_pct)

    print(f"📊 Differential Results:")
    print(f"Log2FC difference: {fc_difference:.3f}")
    print(f"Changed % difference: {changed_difference:.1f}%")
    print(f"Significant % difference: {sig_difference:.1f}%")
    print(f"Statistical comparison p-value: {p_val_systems:.4f}")
    print(f"Effect size (Cohen's d): {cohens_d:.3f}")

    print(f"\n📈 System-specific Results:")
    print(f"Autophagy: {autophagy_changed_pct:.1f}% changed, {autophagy_sig_pct:.1f}% significant")
    print(f"UPS: {ups_changed_pct:.1f}% changed, {ups_sig_pct:.1f}% significant")

    # Overall verdict
    if p_val_systems < 0.05 and abs(cohens_d) > 0.5:
        verdict = "✅ STRONGLY SUPPORTED"
        explanation = f"Significant differential dysfunction (p={p_val_systems:.4f}, d={cohens_d:.3f})"
    elif fc_difference > 0.3 or changed_difference > 20:
        verdict = "✅ SUPPORTED"
        explanation = f"Clear differential patterns (FC diff: {fc_difference:.3f}, % diff: {changed_difference:.1f}%)"
    elif fc_difference > 0.1 or changed_difference > 10:
        verdict = "⚠️ PARTIALLY SUPPORTED"
        explanation = f"Some differential patterns observed"
    else:
        verdict = "❌ REFUTED"
        explanation = f"No clear differential dysfunction patterns"

    print(f"\n🎯 VERDICT: {verdict}")
    print(f"📝 {explanation}")

    # Biological interpretation
    print(f"\n🧬 BIOLOGICAL IMPLICATIONS:")
    if verdict.startswith("✅"):
        print("• Proteostasis systems show distinct dysfunction patterns")
        print("• Compensatory mechanisms may differ between systems")
        print("• Selective therapeutic targeting possible")

        if autophagy_mean_fc > ups_mean_fc:
            print("• Autophagy more activated/upregulated than UPS")
            print("• UPS may be more severely compromised")
            print("• Cells may rely more on autophagy for clearance")
        else:
            print("• UPS more activated/upregulated than autophagy")
            print("• Autophagy may be more severely compromised")
            print("• Cells may rely more on UPS for clearance")

    elif verdict.startswith("⚠️"):
        print("• Some differential dysfunction between systems")
        print("• Both systems may be similarly affected")
        print("• Mixed compensatory responses")
    else:
        print("• Systems appear similarly affected")
        print("• Coordinated proteostasis dysfunction")
        print("• No clear selective vulnerability")

else:
    verdict = "❌ UNSURE"
    explanation = "Insufficient proteostasis pathway proteins found for comparison"
    print(f"VERDICT: {verdict}")
    print(f"Reasoning: {explanation}")
```

## Save Results

```python
# 💾 Save results
if len(df) > 0:
    # Pathway summary
    pathway_summary = pd.DataFrame([{
        'pathway': pr['pathway'],
        'system': pr['system'],
        'proteins_found': pr['proteins_found'],
        'proteins_total': pr['proteins_total'],
        'mean_log2FC': pr['mean_log2FC'],
        'changed_pct': pr['changed_pct']
    } for pr in pathway_results])

    # System comparison summary
    if len(autophagy_df) > 0 and len(ups_df) > 0:
        system_comparison = pd.DataFrame([
            {
                'system': 'Autophagy',
                'proteins': len(autophagy_df),
                'mean_log2FC': autophagy_mean_fc,
                'changed_pct': autophagy_changed_pct,
                'significant_pct': autophagy_sig_pct
            },
            {
                'system': 'UPS',
                'proteins': len(ups_df),
                'mean_log2FC': ups_mean_fc,
                'changed_pct': ups_changed_pct,
                'significant_pct': ups_sig_pct
            }
        ])

    # Overall summary
    summary = {
        'analysis': 'Autophagy and UPS show differential dysfunction patterns',
        'verdict': verdict,
        'proteins_tested': len(df),
        'autophagy_proteins': len(autophagy_df) if 'autophagy_df' in locals() else 0,
        'ups_proteins': len(ups_df) if 'ups_df' in locals() else 0,
        'systems_comparison_pvalue': p_val_systems if 'p_val_systems' in locals() else None,
        'effect_size': cohens_d if 'cohens_d' in locals() else None
    }

    if IN_COLAB:
        df.to_csv('autophagy_ups_comparison_results.csv', index=False)
        pathway_summary.to_csv('proteostasis_pathway_summary.csv', index=False)
        if 'system_comparison' in locals():
            system_comparison.to_csv('system_comparison_summary.csv', index=False)
        pd.DataFrame([summary]).to_csv('autophagy_ups_analysis_summary.csv', index=False)
        files.download('autophagy_ups_comparison_results.csv')
        files.download('proteostasis_pathway_summary.csv')
        if 'system_comparison' in locals():
            files.download('system_comparison_summary.csv')
        print("📁 Results downloaded!")
    else:
        df.to_csv('autophagy_ups_comparison_analysis.csv', index=False)
        pathway_summary.to_csv('proteostasis_pathway_summary.csv', index=False)
        if 'system_comparison' in locals():
            system_comparison.to_csv('system_comparison_summary.csv', index=False)
        print("📁 Results saved!")

print("\n✅ Autophagy vs UPS differential analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified autophagy vs UPS analysis:
- ✅ **Tests both systems comprehensively** (autophagy initiation, elongation, maturation vs UPS core, regulatory, ubiquitination)
- ✅ **Statistical comparison** with effect size calculation
- ✅ **Clear visualizations** showing differential patterns
- ✅ **Biological interpretation** linking to selective proteostasis failure
- ✅ **Therapeutic relevance** for system-specific interventions

**Perfect for evaluating proteostasis system selectivity in neurodegeneration!** ⚖️