# Organellar Markers Analysis - Simplified Version
## Testing: "Organellar markers show compartment-specific dysfunction patterns"

**🚀 Simple approach**: Test key organellar markers across major compartments!

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

## Load Data & Define Organellar Markers

```python
# 📊 Load data
adata = sc.read_h5ad(data_file)
print(f"Data loaded: {adata.shape}")

# 🏷️ Setup tau status
if 'TauStatus' in adata.obs.columns:
    adata.obs['tau_status'] = adata.obs['TauStatus']
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)

print(f"Tau+: {adata.obs['tau_positive'].sum()}, Tau-: {(adata.obs['tau_positive']==0).sum()}")

# 🏢 Organellar Marker Proteins
organellar_markers = {
    'Endoplasmic_Reticulum': [
        'CALR',       # Calreticulin (ER chaperone)
        'PDIA4',      # Protein disulfide isomerase A4
        'HSP90B1',    # Heat shock protein 90 beta family member 1 (GRP94)
        'HSPA5',      # Heat shock protein A5 (BiP/GRP78)
        'ATF6',       # Activating transcription factor 6
        'ERN1',       # Endoplasmic reticulum to nucleus signaling 1 (IRE1)
        'EIF2AK3'     # PERK kinase
    ],
    'Golgi_Apparatus': [
        'GOLGA1',     # Golgin A1
        'GOLGA2',     # Golgin A2 (GM130)
        'GOLGA3',     # Golgin A3
        'GOLGB1',     # Golgin B1
        'GOLIM4',     # Golgi integral membrane protein 4
        'B4GALT1',    # Beta-1,4-galactosyltransferase 1
        'GALNT1'      # Polypeptide N-acetylgalactosaminyltransferase 1
    ],
    'Lysosomes': [
        'LAMP1',      # Lysosomal associated membrane protein 1
        'LAMP2',      # Lysosomal associated membrane protein 2
        'CTSD',       # Cathepsin D
        'CTSB',       # Cathepsin B
        'CTSL',       # Cathepsin L
        'LGALS3',     # Galectin 3 (lysosomal damage marker)
        'MCOLN1'      # Mucolipin 1 (TRPML1)
    ],
    'Mitochondria': [
        'VDAC1',      # Voltage dependent anion channel 1
        'TOMM20',     # Translocase of outer membrane 20
        'TOMM40',     # Translocase of outer membrane 40
        'TIMM23',     # Translocase of inner membrane 23
        'HSPA9',      # Heat shock protein A9 (mortalin)
        'HSPD1',      # Heat shock protein D1 (HSP60)
        'COX4I1'      # Cytochrome c oxidase subunit 4I1
    ],
    'Peroxisomes': [
        'PEX19',      # Peroxisomal biogenesis factor 19
        'PEX3',       # Peroxisomal biogenesis factor 3
        'PEX14',      # Peroxisomal biogenesis factor 14
        'ACOX1',      # Acyl-CoA oxidase 1
        'CAT',        # Catalase
        'SCP2'        # Sterol carrier protein 2
    ],
    'Endosomes': [
        'EEA1',       # Early endosome antigen 1
        'RAB5A',      # RAB5A, member RAS oncogene family
        'RAB7A',      # RAB7A, member RAS oncogene family
        'TFRC',       # Transferrin receptor
        'LDLR',       # Low density lipoprotein receptor
        'SNX1'        # Sorting nexin 1
    ]
}

total_proteins = sum(len(proteins) for proteins in organellar_markers.values())
print(f"🎯 Testing {total_proteins} organellar marker proteins across {len(organellar_markers)} compartments")
```

## Find & Analyze Organellar Markers

```python
# 🔍 Find available proteins and analyze each organelle
protein_names = list(adata.var_names)
organelle_results = []

for organelle_name, proteins in organellar_markers.items():
    print(f"\n📊 Analyzing {organelle_name}...")

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
                'organelle': organelle_name,
                'protein': protein,
                'log2FC': log2fc,
                'p_value': p_val,
                'mean_tau_pos': mean_pos,
                'mean_tau_neg': mean_neg,
                'changed': abs(log2fc) > 0.2  # Changed if |log2FC| > 0.2
            })

        # Organelle-level statistics
        log2fcs = [r['log2FC'] for r in protein_results]
        mean_log2fc = np.mean(log2fcs)
        std_log2fc = np.std(log2fcs)
        changed_count = sum(r['changed'] for r in protein_results)
        changed_pct = changed_count / len(protein_results) * 100

        # Determine dysfunction pattern
        if mean_log2fc > 0.2:
            pattern = "Upregulated"
        elif mean_log2fc < -0.2:
            pattern = "Downregulated"
        elif std_log2fc > 0.5:
            pattern = "Mixed/Variable"
        else:
            pattern = "Stable"

        # Add to organelle results
        organelle_results.append({
            'organelle': organelle_name,
            'proteins_found': len(found_proteins),
            'proteins_total': len(proteins),
            'mean_log2FC': mean_log2fc,
            'std_log2FC': std_log2fc,
            'changed_count': changed_count,
            'changed_pct': changed_pct,
            'dysfunction_pattern': pattern,
            'proteins_tested': protein_results
        })

        print(f"  Mean log2FC: {mean_log2fc:.3f} (±{std_log2fc:.3f})")
        print(f"  Pattern: {pattern}")
        print(f"  Changed: {changed_count}/{len(found_proteins)} ({changed_pct:.1f}%)")

# Create combined DataFrame for all proteins
all_proteins_df = []
for organelle_result in organelle_results:
    all_proteins_df.extend(organelle_result['proteins_tested'])

df = pd.DataFrame(all_proteins_df)
if len(df) > 0:
    # Add FDR correction
    from statsmodels.stats.multitest import fdrcorrection
    df['p_adjusted'] = fdrcorrection(df['p_value'])[1]
    df['significant'] = df['p_adjusted'] < 0.05

print(f"\n✅ Analysis complete: {len(df)} proteins tested across {len(organelle_results)} organelles")
```

## Visualize Organellar Patterns

```python
if len(df) > 0:
    # 📊 Create comprehensive visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Organelle-level fold changes
    organelle_names = [org['organelle'] for org in organelle_results]
    organelle_fcs = [org['mean_log2FC'] for org in organelle_results]
    organelle_stds = [org['std_log2FC'] for org in organelle_results]

    # Color by pattern
    pattern_colors = {
        'Upregulated': 'green',
        'Downregulated': 'red',
        'Mixed/Variable': 'orange',
        'Stable': 'gray'
    }
    colors = [pattern_colors[org['dysfunction_pattern']] for org in organelle_results]

    bars = ax1.bar(range(len(organelle_names)), organelle_fcs,
                  yerr=organelle_stds, color=colors, alpha=0.7, capsize=5)
    ax1.set_xticks(range(len(organelle_names)))
    ax1.set_xticklabels(organelle_names, rotation=45, ha='right')
    ax1.set_ylabel('Mean Log2 Fold Change')
    ax1.set_title('Organellar Dysfunction Patterns')
    ax1.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax1.axhline(0.2, color='green', linestyle='--', alpha=0.5)
    ax1.axhline(-0.2, color='red', linestyle='--', alpha=0.5)
    ax1.grid(True, alpha=0.3)

    # Add pattern legend
    for pattern, color in pattern_colors.items():
        ax1.scatter([], [], c=color, alpha=0.7, s=60, label=pattern)
    ax1.legend()

    # 2. Variability within organelles
    changed_pcts = [org['changed_pct'] for org in organelle_results]
    bars = ax2.bar(range(len(organelle_names)), changed_pcts, color='purple', alpha=0.7)
    ax2.set_xticks(range(len(organelle_names)))
    ax2.set_xticklabels(organelle_names, rotation=45, ha='right')
    ax2.set_ylabel('% Proteins Changed')
    ax2.set_title('Organellar Marker Variability')
    ax2.axhline(50, color='red', linestyle='--', alpha=0.5, label='High variability')
    ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # 3. Scatter plot: Mean vs Variability
    mean_fcs = [abs(org['mean_log2FC']) for org in organelle_results]
    std_fcs = [org['std_log2FC'] for org in organelle_results]
    colors_scatter = [pattern_colors[org['dysfunction_pattern']] for org in organelle_results]

    ax3.scatter(mean_fcs, std_fcs, c=colors_scatter, s=100, alpha=0.7)
    ax3.set_xlabel('Absolute Mean Log2FC')
    ax3.set_ylabel('Standard Deviation Log2FC')
    ax3.set_title('Dysfunction Consistency vs Magnitude')
    ax3.grid(True, alpha=0.3)

    # Add organelle labels
    for i, org in enumerate(organelle_results):
        ax3.annotate(org['organelle'].replace('_', '\n'),
                    (mean_fcs[i], std_fcs[i]),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=8, ha='left')

    # 4. Heatmap of organellar dysfunction
    if len(organelle_results) > 1:
        heatmap_data = []
        labels = []
        for org in organelle_results:
            heatmap_data.append([org['mean_log2FC'], org['changed_pct']/100])
            labels.append(org['organelle'].replace('_', '\n'))

        heatmap_data = np.array(heatmap_data).T
        im = ax4.imshow(heatmap_data, cmap='RdBu_r', aspect='auto', vmin=-1, vmax=1)
        ax4.set_xticks(range(len(labels)))
        ax4.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        ax4.set_yticks([0, 1])
        ax4.set_yticklabels(['Mean log2FC', '% Changed'])
        ax4.set_title('Organellar Dysfunction Heatmap')
        plt.colorbar(im, ax=ax4, fraction=0.046, pad=0.04)

    plt.tight_layout()
    plt.show()
```

## Evaluate Claim

```python
# 🎯 CLAIM EVALUATION
print("\n" + "="*70)
print("🎯 CLAIM EVALUATION")
print("="*70)
print("Claim: 'Organellar markers show compartment-specific dysfunction patterns'")
print()

if len(df) > 0:
    # Overall statistics
    total_proteins = len(df)
    total_changed = sum(df['changed'])
    pct_changed = total_changed / total_proteins * 100

    sig_changed = sum((df['significant']) & (df['changed']))
    pct_sig_changed = sig_changed / total_proteins * 100

    # Pattern analysis
    patterns = [org['dysfunction_pattern'] for org in organelle_results]
    unique_patterns = len(set(patterns))
    non_stable_organelles = sum(1 for pattern in patterns if pattern != 'Stable')

    print(f"📊 Overall Results:")
    print(f"Proteins tested: {total_proteins}")
    print(f"Proteins changed: {total_changed} ({pct_changed:.1f}%)")
    print(f"Significantly changed: {sig_changed} ({pct_sig_changed:.1f}%)")
    print(f"Organelles with dysfunction: {non_stable_organelles}/{len(organelle_results)}")
    print(f"Unique dysfunction patterns: {unique_patterns}")

    print(f"\n📈 Organelle-specific Results:")
    for org in organelle_results:
        status_symbol = {
            'Upregulated': '↗️',
            'Downregulated': '↘️',
            'Mixed/Variable': '↕️',
            'Stable': '→'
        }
        symbol = status_symbol.get(org['dysfunction_pattern'], '?')
        print(f"  {org['organelle']:20} {org['dysfunction_pattern']:15} {symbol} ({org['changed_pct']:4.1f}% changed)")

    # Overall verdict
    if unique_patterns >= 3 and non_stable_organelles >= 4:
        verdict = "✅ STRONGLY SUPPORTED"
        explanation = f"Multiple distinct dysfunction patterns across {non_stable_organelles} organelles"
    elif unique_patterns >= 2 and non_stable_organelles >= 3:
        verdict = "✅ SUPPORTED"
        explanation = f"Clear compartment-specific patterns ({unique_patterns} patterns, {non_stable_organelles} affected)"
    elif non_stable_organelles >= 2 or pct_changed > 40:
        verdict = "⚠️ PARTIALLY SUPPORTED"
        explanation = f"Some organellar dysfunction but limited pattern diversity"
    else:
        verdict = "❌ REFUTED"
        explanation = f"No clear compartment-specific dysfunction patterns"

    print(f"\n🎯 VERDICT: {verdict}")
    print(f"📝 {explanation}")

    # Biological interpretation
    print(f"\n🧬 BIOLOGICAL IMPLICATIONS:")
    if verdict.startswith("✅"):
        print("• Compartment-specific proteostasis failure")
        print("• Organellar crosstalk disrupted")
        print("• Cellular compartmentalization compromised")
        print("• Multiple quality control systems affected")
        print("• Network-wide dysfunction cascades")

        # Pattern-specific impacts
        pattern_impacts = {
            'Upregulated': "Stress response activation",
            'Downregulated': "Organellar capacity reduced",
            'Mixed/Variable': "Quality control variability",
            'Stable': "Compartment preserved"
        }

        print("\n💡 Pattern-specific impacts:")
        for org in organelle_results:
            if org['dysfunction_pattern'] != 'Stable':
                impact = pattern_impacts.get(org['dysfunction_pattern'], "Function altered")
                print(f"  • {org['organelle']}: {impact}")

    elif verdict.startswith("⚠️"):
        print("• Selective organellar dysfunction")
        print("• Some compartments affected, others preserved")
        print("• Partial compartmentalization breakdown")
        print("• Compensatory mechanisms may be active")
    else:
        print("• Organellar markers appear stable")
        print("• Compartmentalization preserved")
        print("• Quality control systems functioning")

else:
    verdict = "❌ UNSURE"
    explanation = "Insufficient organellar marker proteins found"
    print(f"VERDICT: {verdict}")
    print(f"Reasoning: {explanation}")
```

## Save Results

```python
# 💾 Save results
if len(df) > 0:
    # Organelle summary
    organelle_summary = pd.DataFrame([{
        'organelle': org['organelle'],
        'proteins_found': org['proteins_found'],
        'proteins_total': org['proteins_total'],
        'mean_log2FC': org['mean_log2FC'],
        'std_log2FC': org['std_log2FC'],
        'changed_pct': org['changed_pct'],
        'dysfunction_pattern': org['dysfunction_pattern']
    } for org in organelle_results])

    # Overall summary
    summary = {
        'analysis': 'Organellar markers show compartment-specific dysfunction patterns',
        'verdict': verdict,
        'proteins_tested': total_proteins if 'total_proteins' in locals() else 0,
        'proteins_changed': total_changed if 'total_changed' in locals() else 0,
        'percent_changed': pct_changed if 'pct_changed' in locals() else 0,
        'organelles_affected': non_stable_organelles if 'non_stable_organelles' in locals() else 0,
        'unique_patterns': unique_patterns if 'unique_patterns' in locals() else 0
    }

    if IN_COLAB:
        df.to_csv('organellar_markers_results.csv', index=False)
        organelle_summary.to_csv('organelle_summary.csv', index=False)
        pd.DataFrame([summary]).to_csv('organellar_analysis_summary.csv', index=False)
        files.download('organellar_markers_results.csv')
        files.download('organelle_summary.csv')
        print("📁 Results downloaded!")
    else:
        df.to_csv('organellar_markers_analysis.csv', index=False)
        organelle_summary.to_csv('organelle_summary.csv', index=False)
        print("📁 Results saved!")

print("\n✅ Organellar markers analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified organellar markers analysis:
- ✅ **Tests all major organelles** (ER, Golgi, lysosomes, mitochondria, peroxisomes, endosomes)
- ✅ **Pattern-specific evaluation** (upregulated, downregulated, mixed, stable)
- ✅ **Clear visualizations** showing compartment-specific dysfunction
- ✅ **Biological interpretation** linking to proteostasis network failure
- ✅ **Therapeutic relevance** for organelle-specific interventions

**Perfect for evaluating cellular compartment dysfunction in neurodegeneration!** 🏢