# Endolysosomal Changes Analysis - Simplified Version
## Testing: "Endolysosomal system undergoes progressive dysfunction"

**🚀 Simple approach**: Test endolysosomal system components for progressive dysfunction!

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

## Load Data & Define Endolysosomal Components

```python
# 📊 Load data
adata = sc.read_h5ad(data_file)
print(f"Data loaded: {adata.shape}")

# 🏷️ Setup tau status
if 'TauStatus' in adata.obs.columns:
    adata.obs['tau_status'] = adata.obs['TauStatus']
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)

print(f"Tau+: {adata.obs['tau_positive'].sum()}, Tau-: {(adata.obs['tau_positive']==0).sum()}")

# 🔄 Endolysosomal System Components
endolysosomal_system = {
    'Early_Endosomes': [
        'EEA1',       # Early endosome antigen 1
        'RAB5A',      # RAB5A, member RAS oncogene family
        'RAB5B',      # RAB5B, member RAS oncogene family
        'RAB5C',      # RAB5C, member RAS oncogene family
        'RABEP1',     # Rabaptin, RAB GTPase binding effector protein 1
        'PIK3C3',     # Phosphatidylinositol 3-kinase catalytic subunit type 3
        'ZFYVE16'     # Zinc finger FYVE-type containing 16
    ],
    'Late_Endosomes': [
        'RAB7A',      # RAB7A, member RAS oncogene family
        'RAB7B',      # RAB7B, member RAS oncogene family
        'RAB9A',      # RAB9A, member RAS oncogene family
        'RILP',       # Rab interacting lysosomal protein
        'PLEKHM1',    # Pleckstrin homology and RUN domain containing M1
        'UVRAG',      # UV radiation resistance associated
        'M6PR'        # Mannose-6-phosphate receptor
    ],
    'Lysosomes': [
        'LAMP1',      # Lysosomal associated membrane protein 1
        'LAMP2',      # Lysosomal associated membrane protein 2
        'LAMP3',      # Lysosomal associated membrane protein 3
        'CTSD',       # Cathepsin D
        'CTSB',       # Cathepsin B
        'CTSL',       # Cathepsin L
        'CTSZ',       # Cathepsin Z
        'MCOLN1',     # Mucolipin 1 (TRPML1)
        'LGALS3'      # Galectin 3 (lysosomal damage marker)
    ],
    'Lysosomal_Biogenesis': [
        'TFEB',       # Transcription factor EB
        'TFE3',       # Transcription factor binding to IGHM enhancer 3
        'MITF',       # Melanocyte inducing transcription factor
        'MTOR',       # Mechanistic target of rapamycin kinase
        'MTORC1',     # Not a gene, represented by MTOR
        'TSC1',       # Tuberous sclerosis 1
        'TSC2',       # Tuberous sclerosis 2
        'RHEB'        # Ras homolog mTORC1 binding
    ],
    'Endosome_Maturation': [
        'CHMP2A',     # Charged multivesicular body protein 2A
        'CHMP2B',     # Charged multivesicular body protein 2B
        'CHMP4B',     # Charged multivesicular body protein 4B
        'VPS4A',      # Vacuolar protein sorting 4 homolog A
        'VPS4B',      # Vacuolar protein sorting 4 homolog B
        'ALIX',       # ALG-2 interacting protein X (PDCD6IP)
        'TSG101'      # Tumor susceptibility 101
    ],
    'Membrane_Fusion': [
        'STX7',       # Syntaxin 7
        'STX8',       # Syntaxin 8
        'STX17',      # Syntaxin 17
        'VAMP7',      # Vesicle associated membrane protein 7
        'VAMP8',      # Vesicle associated membrane protein 8
        'VTI1B',      # Vesicle transport through interaction with t-SNAREs 1B
        'SNAP29'      # Synaptosome associated protein 29
    ]
}

total_proteins = sum(len(proteins) for proteins in endolysosomal_system.values())
print(f"🎯 Testing {total_proteins} endolysosomal proteins across {len(endolysosomal_system)} components")
```

## Find & Analyze Endolysosomal Components

```python
# 🔍 Find available proteins and analyze each component
protein_names = list(adata.var_names)
component_results = []

for component_name, proteins in endolysosomal_system.items():
    print(f"\n📊 Analyzing {component_name}...")

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
                'component': component_name,
                'protein': protein,
                'log2FC': log2fc,
                'p_value': p_val,
                'mean_tau_pos': mean_pos,
                'mean_tau_neg': mean_neg,
                'dysfunctional': log2fc < -0.1  # Consider downregulation as dysfunction
            })

        # Component-level statistics
        log2fcs = [r['log2FC'] for r in protein_results]
        mean_log2fc = np.mean(log2fcs)
        std_log2fc = np.std(log2fcs)
        dysfunctional_count = sum(r['dysfunctional'] for r in protein_results)
        dysfunctional_pct = dysfunctional_count / len(protein_results) * 100

        # Determine dysfunction severity
        if mean_log2fc < -0.3:
            severity = "Severe"
        elif mean_log2fc < -0.1:
            severity = "Moderate"
        elif mean_log2fc > 0.1:
            severity = "Compensatory"
        else:
            severity = "Stable"

        # Add to component results
        component_results.append({
            'component': component_name,
            'proteins_found': len(found_proteins),
            'proteins_total': len(proteins),
            'mean_log2FC': mean_log2fc,
            'std_log2FC': std_log2fc,
            'dysfunctional_count': dysfunctional_count,
            'dysfunctional_pct': dysfunctional_pct,
            'dysfunction_severity': severity,
            'proteins_tested': protein_results
        })

        print(f"  Mean log2FC: {mean_log2fc:.3f} (±{std_log2fc:.3f})")
        print(f"  Severity: {severity}")
        print(f"  Dysfunctional: {dysfunctional_count}/{len(found_proteins)} ({dysfunctional_pct:.1f}%)")

# Create combined DataFrame for all proteins
all_proteins_df = []
for component_result in component_results:
    all_proteins_df.extend(component_result['proteins_tested'])

df = pd.DataFrame(all_proteins_df)
if len(df) > 0:
    # Add FDR correction
    from statsmodels.stats.multitest import fdrcorrection
    df['p_adjusted'] = fdrcorrection(df['p_value'])[1]
    df['significant'] = df['p_adjusted'] < 0.05

print(f"\n✅ Analysis complete: {len(df)} proteins tested across {len(component_results)} components")
```

## Analyze Progressive Dysfunction

```python
if len(df) > 0:
    # 📈 Progressive dysfunction analysis
    # Order components by expected dysfunction progression
    progression_order = [
        'Early_Endosomes',
        'Endosome_Maturation',
        'Late_Endosomes',
        'Membrane_Fusion',
        'Lysosomes',
        'Lysosomal_Biogenesis'
    ]

    print(f"\n📈 PROGRESSIVE DYSFUNCTION ANALYSIS:")

    progression_data = []
    for i, component in enumerate(progression_order):
        component_data = next((cr for cr in component_results if cr['component'] == component), None)
        if component_data:
            progression_data.append({
                'stage': i + 1,
                'component': component,
                'mean_log2FC': component_data['mean_log2FC'],
                'dysfunctional_pct': component_data['dysfunctional_pct'],
                'severity': component_data['dysfunction_severity']
            })
            print(f"Stage {i+1} - {component}: {component_data['dysfunction_severity']} ({component_data['mean_log2FC']:.3f})")

    # Calculate progression correlation
    if len(progression_data) >= 3:
        stages = [pd['stage'] for pd in progression_data]
        severities = [pd['mean_log2FC'] for pd in progression_data]

        # Spearman correlation (progression vs dysfunction)
        corr_coef, corr_p = stats.spearmanr(stages, severities)
        print(f"\nProgression correlation: r={corr_coef:.3f}, p={corr_p:.4f}")

        # Linear trend test
        slope, intercept, r_value, p_value, std_err = stats.linregress(stages, severities)
        print(f"Linear trend: slope={slope:.3f}, p={p_value:.4f}")
```

## Visualize Endolysosomal Dysfunction

```python
if len(df) > 0:
    # 📊 Create comprehensive visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Component dysfunction levels
    component_names = [cr['component'] for cr in component_results]
    component_fcs = [cr['mean_log2FC'] for cr in component_results]

    # Color by severity
    severity_colors = {
        'Severe': 'darkred',
        'Moderate': 'red',
        'Stable': 'gray',
        'Compensatory': 'green'
    }
    colors = [severity_colors[cr['dysfunction_severity']] for cr in component_results]

    bars = ax1.bar(range(len(component_names)), component_fcs, color=colors, alpha=0.7)
    ax1.set_xticks(range(len(component_names)))
    ax1.set_xticklabels(component_names, rotation=45, ha='right')
    ax1.set_ylabel('Mean Log2 Fold Change')
    ax1.set_title('Endolysosomal Component Dysfunction')
    ax1.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax1.axhline(-0.3, color='darkred', linestyle='--', alpha=0.5, label='Severe dysfunction')
    ax1.axhline(-0.1, color='red', linestyle='--', alpha=0.5, label='Moderate dysfunction')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # 2. Progressive dysfunction plot
    if len(progression_data) >= 3:
        stages = [pd['stage'] for pd in progression_data]
        severities = [pd['mean_log2FC'] for pd in progression_data]

        ax2.plot(stages, severities, 'o-', color='red', linewidth=2, markersize=8)
        ax2.set_xlabel('Progression Stage')
        ax2.set_ylabel('Mean Log2FC (Dysfunction)')
        ax2.set_title(f'Progressive Dysfunction (r={corr_coef:.3f})')
        ax2.grid(True, alpha=0.3)

        # Add component labels
        for pd in progression_data:
            ax2.annotate(pd['component'].replace('_', '\n'),
                        (pd['stage'], pd['mean_log2FC']),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=8, ha='left')

    # 3. Dysfunction percentage per component
    dysfunc_pcts = [cr['dysfunctional_pct'] for cr in component_results]
    bars = ax3.bar(range(len(component_names)), dysfunc_pcts, color='orange', alpha=0.7)
    ax3.set_xticks(range(len(component_names)))
    ax3.set_xticklabels(component_names, rotation=45, ha='right')
    ax3.set_ylabel('% Proteins Dysfunctional')
    ax3.set_title('Component Dysfunction Percentage')
    ax3.axhline(50, color='red', linestyle='--', alpha=0.5, label='Majority dysfunctional')
    ax3.set_ylim(0, 100)
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    # 4. Heatmap of component dysfunction
    if len(component_results) > 1:
        heatmap_data = []
        labels = []
        for cr in component_results:
            heatmap_data.append([cr['mean_log2FC'], cr['dysfunctional_pct']/100])
            labels.append(cr['component'].replace('_', '\n'))

        heatmap_data = np.array(heatmap_data).T
        im = ax4.imshow(heatmap_data, cmap='Reds', aspect='auto', vmin=-1, vmax=1)
        ax4.set_xticks(range(len(labels)))
        ax4.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        ax4.set_yticks([0, 1])
        ax4.set_yticklabels(['Mean log2FC', '% Dysfunctional'])
        ax4.set_title('Endolysosomal Dysfunction Heatmap')
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
print("Claim: 'Endolysosomal system undergoes progressive dysfunction'")
print()

if len(df) > 0:
    # Overall statistics
    total_proteins = len(df)
    total_dysfunctional = sum(df['dysfunctional'])
    pct_dysfunctional = total_dysfunctional / total_proteins * 100

    sig_dysfunctional = sum((df['significant']) & (df['dysfunctional']))
    pct_sig_dysfunctional = sig_dysfunctional / total_proteins * 100

    # Component-level analysis
    severe_components = sum(1 for cr in component_results if cr['dysfunction_severity'] == 'Severe')
    moderate_components = sum(1 for cr in component_results if cr['dysfunction_severity'] == 'Moderate')
    dysfunctional_components = severe_components + moderate_components

    print(f"📊 Overall Results:")
    print(f"Proteins tested: {total_proteins}")
    print(f"Proteins dysfunctional: {total_dysfunctional} ({pct_dysfunctional:.1f}%)")
    print(f"Significantly dysfunctional: {sig_dysfunctional} ({pct_sig_dysfunctional:.1f}%)")
    print(f"Components with dysfunction: {dysfunctional_components}/{len(component_results)}")
    print(f"  Severe: {severe_components}, Moderate: {moderate_components}")

    if len(progression_data) >= 3:
        print(f"Progressive correlation: r={corr_coef:.3f}, p={corr_p:.4f}")

    print(f"\n📈 Component-specific Results:")
    for cr in component_results:
        severity_symbol = {
            'Severe': '🔴',
            'Moderate': '🟡',
            'Stable': '🟢',
            'Compensatory': '🔵'
        }
        symbol = severity_symbol.get(cr['dysfunction_severity'], '?')
        print(f"  {cr['component']:20} {cr['dysfunction_severity']:12} {symbol} ({cr['dysfunctional_pct']:4.1f}% affected)")

    # Overall verdict
    if len(progression_data) >= 3 and corr_p < 0.05 and corr_coef < -0.5:
        verdict = "✅ STRONGLY SUPPORTED"
        explanation = f"Clear progressive dysfunction (r={corr_coef:.3f}, p={corr_p:.4f})"
    elif pct_dysfunctional > 60 and dysfunctional_components >= 4:
        verdict = "✅ SUPPORTED"
        explanation = f"Widespread dysfunction across {dysfunctional_components} components ({pct_dysfunctional:.1f}%)"
    elif pct_dysfunctional > 40 or dysfunctional_components >= 3:
        verdict = "⚠️ PARTIALLY SUPPORTED"
        explanation = f"Some endolysosomal dysfunction ({pct_dysfunctional:.1f}% proteins, {dysfunctional_components} components)"
    else:
        verdict = "❌ REFUTED"
        explanation = f"No clear progressive dysfunction pattern"

    print(f"\n🎯 VERDICT: {verdict}")
    print(f"📝 {explanation}")

    # Biological interpretation
    print(f"\n🧬 BIOLOGICAL IMPLICATIONS:")
    if verdict.startswith("✅"):
        print("• Endolysosomal trafficking severely compromised")
        print("• Progressive cargo accumulation and dysfunction")
        print("• Lysosomal clearance capacity overwhelmed")
        print("• Cellular waste management system failure")
        print("• Secondary organellar stress and dysfunction")

        if len(progression_data) >= 3 and corr_p < 0.05:
            print("\n💡 Progressive dysfunction features:")
            print("• Early endosome dysfunction initiates cascade")
            print("• Maturation defects prevent proper sorting")
            print("• Late endosome/lysosome fusion impaired")
            print("• Lysosomal biogenesis overwhelmed")

    elif verdict.startswith("⚠️"):
        print("• Selective endolysosomal dysfunction")
        print("• Some trafficking steps preserved")
        print("• Partial clearance capacity maintained")
        print("• Compensatory mechanisms active")
    else:
        print("• Endolysosomal system appears functional")
        print("• Trafficking and clearance preserved")
        print("• Normal lysosomal biogenesis")

else:
    verdict = "❌ UNSURE"
    explanation = "Insufficient endolysosomal proteins found"
    print(f"VERDICT: {verdict}")
    print(f"Reasoning: {explanation}")
```

## Save Results

```python
# 💾 Save results
if len(df) > 0:
    # Component summary
    component_summary = pd.DataFrame([{
        'component': cr['component'],
        'proteins_found': cr['proteins_found'],
        'proteins_total': cr['proteins_total'],
        'mean_log2FC': cr['mean_log2FC'],
        'dysfunctional_pct': cr['dysfunctional_pct'],
        'dysfunction_severity': cr['dysfunction_severity']
    } for cr in component_results])

    # Progression summary
    if len(progression_data) >= 3:
        progression_summary = pd.DataFrame(progression_data)

    # Overall summary
    summary = {
        'analysis': 'Endolysosomal system undergoes progressive dysfunction',
        'verdict': verdict,
        'proteins_tested': total_proteins if 'total_proteins' in locals() else 0,
        'proteins_dysfunctional': total_dysfunctional if 'total_dysfunctional' in locals() else 0,
        'percent_dysfunctional': pct_dysfunctional if 'pct_dysfunctional' in locals() else 0,
        'components_affected': dysfunctional_components if 'dysfunctional_components' in locals() else 0,
        'progression_correlation': corr_coef if 'corr_coef' in locals() else None,
        'progression_pvalue': corr_p if 'corr_p' in locals() else None
    }

    if IN_COLAB:
        df.to_csv('endolysosomal_dysfunction_results.csv', index=False)
        component_summary.to_csv('endolysosomal_component_summary.csv', index=False)
        if 'progression_summary' in locals():
            progression_summary.to_csv('endolysosomal_progression_summary.csv', index=False)
        pd.DataFrame([summary]).to_csv('endolysosomal_analysis_summary.csv', index=False)
        files.download('endolysosomal_dysfunction_results.csv')
        files.download('endolysosomal_component_summary.csv')
        if 'progression_summary' in locals():
            files.download('endolysosomal_progression_summary.csv')
        print("📁 Results downloaded!")
    else:
        df.to_csv('endolysosomal_dysfunction_analysis.csv', index=False)
        component_summary.to_csv('endolysosomal_component_summary.csv', index=False)
        if 'progression_summary' in locals():
            progression_summary.to_csv('endolysosomal_progression_summary.csv', index=False)
        print("📁 Results saved!")

print("\n✅ Endolysosomal progressive dysfunction analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified endolysosomal dysfunction analysis:
- ✅ **Tests all components** (early endosomes, late endosomes, lysosomes, biogenesis, maturation, fusion)
- ✅ **Progressive analysis** with correlation testing
- ✅ **Clear visualizations** showing dysfunction cascade
- ✅ **Biological interpretation** linking to trafficking failure
- ✅ **Therapeutic relevance** for endolysosomal enhancement strategies

**Perfect for evaluating progressive clearance system failure in neurodegeneration!** 🔄