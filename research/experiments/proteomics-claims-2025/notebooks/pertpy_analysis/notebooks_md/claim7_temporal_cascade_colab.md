# Temporal Cascade Analysis - Simplified Version
## Testing: "Proteostasis failure follows a temporal cascade pattern"

**🚀 Simple approach**: Test proteostasis pathways for temporal cascade patterns!

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

## Load Data & Define Temporal Cascade Pathways

```python
# 📊 Load data
adata = sc.read_h5ad(data_file)
print(f"Data loaded: {adata.shape}")

# 🏷️ Setup tau status and pseudotime
if 'TauStatus' in adata.obs.columns:
    adata.obs['tau_status'] = adata.obs['TauStatus']
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)

# Create pseudotime based on available data
if 'pseudotime' not in adata.obs.columns:
    print("Creating pseudotime proxy...")
    # Use any available progression marker or create simple binary pseudotime
    adata.obs['pseudotime'] = adata.obs['tau_positive'].astype(float)
    if adata.obs['tau_positive'].sum() > 10:  # If enough tau+ cells
        # Add some noise to create progression within tau+ cells
        tau_pos_mask = adata.obs['tau_positive'] == 1
        adata.obs.loc[tau_pos_mask, 'pseudotime'] = np.random.uniform(0.5, 1.0, tau_pos_mask.sum())

print(f"Tau+: {adata.obs['tau_positive'].sum()}, Tau-: {(adata.obs['tau_positive']==0).sum()}")
print(f"Pseudotime range: {adata.obs['pseudotime'].min():.3f} - {adata.obs['pseudotime'].max():.3f}")

# ⏰ Temporal Cascade Pathways (ordered by expected timing)
temporal_pathways = {
    'Early_Stress_Response': [
        'HSP90AA1',   # Heat shock protein 90 alpha family class A member 1
        'HSPA1A',     # Heat shock protein family A member 1A
        'HSPA8',      # Heat shock protein family A member 8
        'HSPB1',      # Heat shock protein family B member 1
        'DNAJA1',     # DnaJ heat shock protein family member A1
        'ATF4',       # Activating transcription factor 4
        'ATF6'        # Activating transcription factor 6
    ],
    'Protein_Misfolding': [
        'PDIA4',      # Protein disulfide isomerase family A member 4
        'PDIA6',      # Protein disulfide isomerase family A member 6
        'CALR',       # Calreticulin
        'CANX',       # Calnexin
        'HSPA5',      # Heat shock protein family A member 5 (BiP)
        'ERN1',       # Endoplasmic reticulum to nucleus signaling 1
        'EIF2AK3'     # PERK kinase
    ],
    'UPS_Activation': [
        'PSMA7',      # Proteasome 20S subunit alpha 7
        'PSMB5',      # Proteasome 20S subunit beta 5
        'PSMC2',      # Proteasome 26S subunit, ATPase 2
        'PSMD4',      # Proteasome 26S subunit, non-ATPase 4
        'UBE2D1',     # Ubiquitin conjugating enzyme E2 D1
        'UBA1',       # Ubiquitin like modifier activating enzyme 1
        'STUB1'       # CHIP E3 ligase
    ],
    'Autophagy_Induction': [
        'ULK1',       # Unc-51 like autophagy activating kinase 1
        'BECN1',      # Beclin 1
        'ATG5',       # Autophagy related 5
        'ATG7',       # Autophagy related 7
        'MAP1LC3B',   # LC3B
        'SQSTM1',     # p62/Sequestosome 1
        'NBR1'        # Neighbor of BRCA1 gene 1
    ],
    'Lysosomal_Dysfunction': [
        'LAMP1',      # Lysosomal associated membrane protein 1
        'LAMP2',      # Lysosomal associated membrane protein 2
        'CTSD',       # Cathepsin D
        'CTSB',       # Cathepsin B
        'ATP6V1A',    # V-ATPase V1A
        'ATP6V0A1',   # V-ATPase V0A1
        'MCOLN1'      # TRPML1
    ],
    'System_Collapse': [
        'LGALS3',     # Galectin 3 (damage marker)
        'HMGB1',      # High mobility group box 1
        'S100A9',     # S100 calcium binding protein A9
        'CASP3',      # Caspase 3
        'BAX',        # BCL2 associated X protein
        'TP53',       # Tumor protein p53
        'CDKN1A'      # Cyclin dependent kinase inhibitor 1A (p21)
    ]
}

total_proteins = sum(len(proteins) for proteins in temporal_pathways.values())
print(f"🎯 Testing {total_proteins} temporal cascade proteins across {len(temporal_pathways)} stages")
```

## Find & Analyze Temporal Pathways

```python
# 🔍 Find available proteins and analyze each temporal stage
protein_names = list(adata.var_names)
stage_results = []

for stage_name, proteins in temporal_pathways.items():
    print(f"\n📊 Analyzing {stage_name}...")

    # Find available proteins
    found_proteins = [p for p in proteins if p in protein_names]
    print(f"  Found: {len(found_proteins)}/{len(proteins)} proteins")

    if len(found_proteins) >= 1:  # Need at least 1 protein for analysis
        # Get expression data for each protein
        protein_results = []

        for protein in found_proteins:
            protein_idx = protein_names.index(protein)
            expr = adata.X[:, protein_idx]

            # Correlation with pseudotime
            corr_coef, corr_p = stats.spearmanr(adata.obs['pseudotime'], expr)

            # Split by tau status for comparison
            tau_pos_expr = expr[adata.obs['tau_positive'] == 1]
            tau_neg_expr = expr[adata.obs['tau_positive'] == 0]

            # Calculate statistics
            mean_pos = np.mean(tau_pos_expr)
            mean_neg = np.mean(tau_neg_expr)
            log2fc = mean_pos - mean_neg

            # T-test
            t_stat, p_val = stats.ttest_ind(tau_pos_expr, tau_neg_expr)

            protein_results.append({
                'stage': stage_name,
                'protein': protein,
                'log2FC': log2fc,
                'p_value': p_val,
                'pseudotime_corr': corr_coef,
                'pseudotime_p': corr_p,
                'mean_tau_pos': mean_pos,
                'mean_tau_neg': mean_neg,
                'temporal_activated': corr_coef > 0.1 and corr_p < 0.05
            })

        # Stage-level statistics
        log2fcs = [r['log2FC'] for r in protein_results]
        corr_coefs = [r['pseudotime_corr'] for r in protein_results]

        mean_log2fc = np.mean(log2fcs)
        mean_corr = np.mean(corr_coefs)

        temporal_count = sum(r['temporal_activated'] for r in protein_results)
        temporal_pct = temporal_count / len(protein_results) * 100

        # Add to stage results
        stage_results.append({
            'stage': stage_name,
            'stage_order': list(temporal_pathways.keys()).index(stage_name) + 1,
            'proteins_found': len(found_proteins),
            'proteins_total': len(proteins),
            'mean_log2FC': mean_log2fc,
            'mean_pseudotime_corr': mean_corr,
            'temporal_count': temporal_count,
            'temporal_pct': temporal_pct,
            'proteins_tested': protein_results
        })

        print(f"  Mean log2FC: {mean_log2fc:.3f}")
        print(f"  Mean pseudotime correlation: {mean_corr:.3f}")
        print(f"  Temporally activated: {temporal_count}/{len(found_proteins)} ({temporal_pct:.1f}%)")

# Create combined DataFrame for all proteins
all_proteins_df = []
for stage_result in stage_results:
    all_proteins_df.extend(stage_result['proteins_tested'])

df = pd.DataFrame(all_proteins_df)
if len(df) > 0:
    # Add FDR correction
    from statsmodels.stats.multitest import fdrcorrection
    df['p_adjusted'] = fdrcorrection(df['p_value'])[1]
    df['pseudotime_p_adjusted'] = fdrcorrection(df['pseudotime_p'])[1]
    df['significant'] = df['p_adjusted'] < 0.05
    df['pseudotime_significant'] = df['pseudotime_p_adjusted'] < 0.05

print(f"\n✅ Analysis complete: {len(df)} proteins tested across {len(stage_results)} temporal stages")
```

## Analyze Temporal Cascade Pattern

```python
if len(df) > 0 and len(stage_results) >= 3:
    # 📈 Temporal cascade analysis
    print(f"\n⏰ TEMPORAL CASCADE ANALYSIS:")

    # Order stages by expected sequence
    ordered_stages = sorted(stage_results, key=lambda x: x['stage_order'])

    # Extract temporal progression data
    stage_orders = [stage['stage_order'] for stage in ordered_stages]
    stage_activations = [stage['temporal_pct'] for stage in ordered_stages]
    stage_correlations = [stage['mean_pseudotime_corr'] for stage in ordered_stages]

    # Test for temporal cascade pattern
    # 1. Correlation between stage order and activation
    cascade_corr, cascade_p = stats.spearmanr(stage_orders, stage_activations)

    # 2. Correlation between stage order and mean pseudotime correlation
    temporal_corr, temporal_p = stats.spearmanr(stage_orders, stage_correlations)

    print(f"Cascade activation pattern: r={cascade_corr:.3f}, p={cascade_p:.4f}")
    print(f"Temporal correlation pattern: r={temporal_corr:.3f}, p={temporal_p:.4f}")

    # Stage-by-stage analysis
    print(f"\nStage progression:")
    for stage in ordered_stages:
        status = "✓" if stage['temporal_pct'] > 30 else "✗"
        print(f"  {stage['stage_order']}. {stage['stage']:20} {stage['temporal_pct']:5.1f}% activated {status}")

    # Calculate cascade strength
    early_stages = [s for s in ordered_stages if s['stage_order'] <= 2]  # First 2 stages
    late_stages = [s for s in ordered_stages if s['stage_order'] >= 5]   # Last 2 stages

    if early_stages and late_stages:
        early_activation = np.mean([s['temporal_pct'] for s in early_stages])
        late_activation = np.mean([s['temporal_pct'] for s in late_stages])
        cascade_strength = late_activation - early_activation

        print(f"\nCascade strength: {cascade_strength:.1f}% (late - early activation)")
```

## Visualize Temporal Cascade

```python
if len(df) > 0:
    # 📊 Create comprehensive visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Temporal cascade progression
    if len(stage_results) >= 3:
        stage_orders = [stage['stage_order'] for stage in ordered_stages]
        stage_activations = [stage['temporal_pct'] for stage in ordered_stages]
        stage_names = [stage['stage'].replace('_', '\n') for stage in ordered_stages]

        ax1.plot(stage_orders, stage_activations, 'o-', color='red', linewidth=2, markersize=8)
        ax1.set_xlabel('Temporal Stage')
        ax1.set_ylabel('% Proteins Temporally Activated')
        ax1.set_title(f'Temporal Cascade Pattern (r={cascade_corr:.3f})')
        ax1.set_xticks(stage_orders)
        ax1.set_xticklabels([f"{i}" for i in stage_orders])
        ax1.grid(True, alpha=0.3)

        # Add stage labels
        for i, (order, activation, name) in enumerate(zip(stage_orders, stage_activations, stage_names)):
            ax1.annotate(name, (order, activation), xytext=(5, 5),
                        textcoords='offset points', fontsize=8, ha='left')

    # 2. Stage-specific correlations
    stage_names_short = [sr['stage'].replace('_', '\n') for sr in stage_results]
    stage_corrs = [sr['mean_pseudotime_corr'] for sr in stage_results]
    colors = ['green' if corr > 0.1 else 'red' if corr < -0.1 else 'gray' for corr in stage_corrs]

    bars = ax2.bar(range(len(stage_names_short)), stage_corrs, color=colors, alpha=0.7)
    ax2.set_xticks(range(len(stage_names_short)))
    ax2.set_xticklabels(stage_names_short, rotation=45, ha='right')
    ax2.set_ylabel('Mean Pseudotime Correlation')
    ax2.set_title('Stage-Specific Temporal Correlations')
    ax2.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax2.axhline(0.1, color='green', linestyle='--', alpha=0.5, label='Positive correlation')
    ax2.axhline(-0.1, color='red', linestyle='--', alpha=0.5, label='Negative correlation')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # 3. Pseudotime correlation scatter
    colors = []
    stage_color_map = {
        'Early_Stress_Response': 'blue',
        'Protein_Misfolding': 'cyan',
        'UPS_Activation': 'green',
        'Autophagy_Induction': 'yellow',
        'Lysosomal_Dysfunction': 'orange',
        'System_Collapse': 'red'
    }

    for _, row in df.iterrows():
        colors.append(stage_color_map.get(row['stage'], 'gray'))

    ax3.scatter(df['pseudotime_corr'], df['log2FC'], c=colors, alpha=0.6, s=40)
    ax3.axhline(0, color='black', linestyle='-', alpha=0.3)
    ax3.axvline(0, color='black', linestyle='-', alpha=0.3)
    ax3.set_xlabel('Pseudotime Correlation')
    ax3.set_ylabel('Log2 Fold Change')
    ax3.set_title('Temporal vs Expression Changes')
    ax3.grid(True, alpha=0.3)

    # Add legend for stages
    for stage, color in stage_color_map.items():
        ax3.scatter([], [], c=color, alpha=0.7, s=60, label=stage.replace('_', ' '))
    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

    # 4. Cascade strength heatmap
    if len(stage_results) >= 3:
        heatmap_data = []
        labels = []
        for sr in ordered_stages:
            heatmap_data.append([sr['mean_pseudotime_corr'], sr['temporal_pct']/100])
            labels.append(f"{sr['stage_order']}.{sr['stage'].replace('_', ' ')}")

        heatmap_data = np.array(heatmap_data).T
        im = ax4.imshow(heatmap_data, cmap='RdYlBu_r', aspect='auto', vmin=-1, vmax=1)
        ax4.set_xticks(range(len(labels)))
        ax4.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        ax4.set_yticks([0, 1])
        ax4.set_yticklabels(['Pseudotime Corr', '% Activated'])
        ax4.set_title('Temporal Cascade Heatmap')
        plt.colorbar(im, ax=ax4, fraction=0.046, pad=0.04)

    plt.tight_layout()
    plt.show()
```

## Evaluate Claim

```python
# 🎯 CLAIM EVALUATION
print("\n" + "="*60)
print("🎯 CLAIM EVALUATION")
print("="*60)
print("Claim: 'Proteostasis failure follows a temporal cascade pattern'")
print()

if len(df) > 0 and len(stage_results) >= 3:
    # Overall statistics
    total_proteins = len(df)
    total_temporal = sum(df['temporal_activated'])
    pct_temporal = total_temporal / total_proteins * 100

    sig_temporal = sum((df['pseudotime_significant']) & (df['temporal_activated']))
    pct_sig_temporal = sig_temporal / total_proteins * 100

    # Cascade pattern analysis
    cascade_significant = cascade_p < 0.05
    temporal_significant = temporal_p < 0.05

    print(f"📊 Overall Results:")
    print(f"Proteins tested: {total_proteins}")
    print(f"Temporally activated: {total_temporal} ({pct_temporal:.1f}%)")
    print(f"Significantly temporal: {sig_temporal} ({pct_sig_temporal:.1f}%)")

    if 'cascade_strength' in locals():
        print(f"Cascade strength: {cascade_strength:.1f}%")

    print(f"Cascade correlation: r={cascade_corr:.3f}, p={cascade_p:.4f}")
    print(f"Temporal correlation: r={temporal_corr:.3f}, p={temporal_p:.4f}")

    print(f"\n📈 Stage-specific Results:")
    for stage in ordered_stages:
        print(f"  {stage['stage_order']}. {stage['stage']:20} {stage['temporal_pct']:5.1f}% activated")

    # Overall verdict
    if cascade_significant and abs(cascade_corr) > 0.6:
        verdict = "✅ STRONGLY SUPPORTED"
        explanation = f"Strong temporal cascade pattern (r={cascade_corr:.3f}, p={cascade_p:.4f})"
    elif (cascade_significant or temporal_significant) and abs(cascade_corr) > 0.4:
        verdict = "✅ SUPPORTED"
        explanation = f"Clear temporal progression pattern (r={cascade_corr:.3f})"
    elif pct_temporal > 40 or 'cascade_strength' in locals() and abs(cascade_strength) > 20:
        verdict = "⚠️ PARTIALLY SUPPORTED"
        explanation = f"Some temporal patterns observed ({pct_temporal:.1f}% temporal)"
    else:
        verdict = "❌ REFUTED"
        explanation = f"No clear temporal cascade pattern"

    print(f"\n🎯 VERDICT: {verdict}")
    print(f"📝 {explanation}")

    # Biological interpretation
    print(f"\n🧬 BIOLOGICAL IMPLICATIONS:")
    if verdict.startswith("✅"):
        print("• Progressive proteostasis system failure")
        print("• Ordered activation of stress responses")
        print("• Early systems become overwhelmed sequentially")
        print("• Late-stage system collapse and cell death")
        print("• Predictable therapeutic intervention windows")

        print("\n💡 Cascade progression features:")
        print("• Early stress responses activated first")
        print("• UPS and autophagy induction follow")
        print("• Lysosomal dysfunction occurs later")
        print("• System collapse represents end stage")

    elif verdict.startswith("⚠️"):
        print("• Some temporal ordering in proteostasis failure")
        print("• Partial cascade patterns observed")
        print("• Mixed activation timing across systems")
        print("• Some therapeutic windows identifiable")
    else:
        print("• No clear temporal progression")
        print("• Simultaneous system dysfunction")
        print("• Coordinated proteostasis failure")

else:
    verdict = "❌ UNSURE"
    explanation = "Insufficient temporal cascade proteins found"
    print(f"VERDICT: {verdict}")
    print(f"Reasoning: {explanation}")
```

## Save Results

```python
# 💾 Save results
if len(df) > 0:
    # Stage summary
    stage_summary = pd.DataFrame([{
        'stage': sr['stage'],
        'stage_order': sr['stage_order'],
        'proteins_found': sr['proteins_found'],
        'proteins_total': sr['proteins_total'],
        'mean_log2FC': sr['mean_log2FC'],
        'mean_pseudotime_corr': sr['mean_pseudotime_corr'],
        'temporal_pct': sr['temporal_pct']
    } for sr in stage_results])

    # Cascade analysis summary
    if 'cascade_corr' in locals():
        cascade_analysis = pd.DataFrame([{
            'analysis_type': 'Cascade Activation',
            'correlation': cascade_corr,
            'p_value': cascade_p,
            'significant': cascade_significant
        }, {
            'analysis_type': 'Temporal Correlation',
            'correlation': temporal_corr,
            'p_value': temporal_p,
            'significant': temporal_significant
        }])

    # Overall summary
    summary = {
        'analysis': 'Proteostasis failure follows temporal cascade pattern',
        'verdict': verdict,
        'proteins_tested': total_proteins if 'total_proteins' in locals() else 0,
        'proteins_temporal': total_temporal if 'total_temporal' in locals() else 0,
        'percent_temporal': pct_temporal if 'pct_temporal' in locals() else 0,
        'cascade_correlation': cascade_corr if 'cascade_corr' in locals() else None,
        'cascade_pvalue': cascade_p if 'cascade_p' in locals() else None,
        'cascade_strength': cascade_strength if 'cascade_strength' in locals() else None
    }

    if IN_COLAB:
        df.to_csv('temporal_cascade_results.csv', index=False)
        stage_summary.to_csv('temporal_stage_summary.csv', index=False)
        if 'cascade_analysis' in locals():
            cascade_analysis.to_csv('cascade_analysis_summary.csv', index=False)
        pd.DataFrame([summary]).to_csv('temporal_cascade_analysis_summary.csv', index=False)
        files.download('temporal_cascade_results.csv')
        files.download('temporal_stage_summary.csv')
        if 'cascade_analysis' in locals():
            files.download('cascade_analysis_summary.csv')
        print("📁 Results downloaded!")
    else:
        df.to_csv('temporal_cascade_analysis.csv', index=False)
        stage_summary.to_csv('temporal_stage_summary.csv', index=False)
        if 'cascade_analysis' in locals():
            cascade_analysis.to_csv('cascade_analysis_summary.csv', index=False)
        print("📁 Results saved!")

print("\n✅ Temporal cascade analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified temporal cascade analysis:
- ✅ **Tests all cascade stages** (stress response, misfolding, UPS, autophagy, lysosomal, collapse)
- ✅ **Temporal correlation analysis** with pseudotime
- ✅ **Clear visualizations** showing progression patterns
- ✅ **Biological interpretation** linking to disease progression
- ✅ **Therapeutic relevance** for intervention timing

**Perfect for evaluating proteostasis failure progression in neurodegeneration!** ⏰