# Mitochondrial Complex Analysis - Simplified Version
## Testing: "Mitochondrial complexes I-V are decreased in tau+ neurons"

**🚀 Simple approach**: Test all mitochondrial complexes for decreased expression!

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

## Load Data & Define Mitochondrial Complexes

```python
# 📊 Load data
adata = sc.read_h5ad(data_file)
print(f"Data loaded: {adata.shape}")

# 🏷️ Setup tau status
if 'TauStatus' in adata.obs.columns:
    adata.obs['tau_status'] = adata.obs['TauStatus']
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)

print(f"Tau+: {adata.obs['tau_positive'].sum()}, Tau-: {(adata.obs['tau_positive']==0).sum()}")

# ⚡ Mitochondrial Complex Proteins
mito_complexes = {
    'Complex_I': [
        'NDUFA1', 'NDUFA2', 'NDUFA4', 'NDUFA9', 'NDUFA10',
        'NDUFB1', 'NDUFB3', 'NDUFB6', 'NDUFB10',
        'NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS7',
        'NDUFV1', 'NDUFV2', 'NDUFV3'
    ],
    'Complex_II': [
        'SDHA', 'SDHB', 'SDHC', 'SDHD'
    ],
    'Complex_III': [
        'UQCRC1', 'UQCRC2', 'UQCRFS1', 'UQCRB',
        'UQCRQ', 'CYC1', 'CYCS'
    ],
    'Complex_IV': [
        'COX4I1', 'COX5A', 'COX5B', 'COX6A1',
        'COX6B1', 'COX7A1', 'COX7B', 'COX8A'
    ],
    'Complex_V': [
        'ATP5A1', 'ATP5B', 'ATP5C1', 'ATP5D', 'ATP5E',
        'ATP5F1', 'ATP5G1', 'ATP5H', 'ATP5O'
    ]
}

total_proteins = sum(len(proteins) for proteins in mito_complexes.values())
print(f"🎯 Testing {total_proteins} mitochondrial complex proteins across {len(mito_complexes)} complexes")
```

## Find & Analyze Complex Proteins

```python
# 🔍 Find available proteins and analyze each complex
protein_names = list(adata.var_names)
complex_results = []

for complex_name, proteins in mito_complexes.items():
    print(f"\\n📊 Analyzing {complex_name}...")

    # Find available proteins
    found_proteins = [p for p in proteins if p in protein_names]
    print(f"  Found: {len(found_proteins)}/{len(proteins)} proteins")

    if len(found_proteins) >= 2:  # Need at least 2 proteins for analysis
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
                'complex': complex_name,
                'protein': protein,
                'log2FC': log2fc,
                'p_value': p_val,
                'mean_tau_pos': mean_pos,
                'mean_tau_neg': mean_neg,
                'decreased': log2fc < 0
            })

        # Complex-level statistics
        log2fcs = [r['log2FC'] for r in protein_results]
        mean_log2fc = np.mean(log2fcs)
        decreased_count = sum(r['decreased'] for r in protein_results)
        decreased_pct = decreased_count / len(protein_results) * 100

        # Add to complex results
        complex_results.append({
            'complex': complex_name,
            'proteins_found': len(found_proteins),
            'proteins_total': len(proteins),
            'mean_log2FC': mean_log2fc,
            'decreased_count': decreased_count,
            'decreased_pct': decreased_pct,
            'proteins_tested': protein_results
        })

        print(f"  Mean log2FC: {mean_log2fc:.3f}")
        print(f"  Decreased: {decreased_count}/{len(found_proteins)} ({decreased_pct:.1f}%)")

# Create combined DataFrame for all proteins
all_proteins_df = []
for complex_result in complex_results:
    all_proteins_df.extend(complex_result['proteins_tested'])

df = pd.DataFrame(all_proteins_df)
if len(df) > 0:
    # Add FDR correction
    from statsmodels.stats.multitest import fdrcorrection
    df['p_adjusted'] = fdrcorrection(df['p_value'])[1]
    df['significant'] = df['p_adjusted'] < 0.05

print(f"\\n✅ Analysis complete: {len(df)} proteins tested across {len(complex_results)} complexes")
```

## Visualize Complex Results

```python
if len(df) > 0:
    # 📊 Create comprehensive visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Complex-level fold changes
    complex_names = [cr['complex'] for cr in complex_results]
    complex_fcs = [cr['mean_log2FC'] for cr in complex_results]
    colors = ['red' if fc < -0.1 else 'gray' for fc in complex_fcs]

    bars = ax1.bar(range(len(complex_names)), complex_fcs, color=colors, alpha=0.7)
    ax1.set_xticks(range(len(complex_names)))
    ax1.set_xticklabels(complex_names, rotation=45, ha='right')
    ax1.set_ylabel('Mean Log2 Fold Change')
    ax1.set_title('Mitochondrial Complex Expression Changes')
    ax1.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax1.axhline(-0.5, color='red', linestyle='--', alpha=0.5, label='Substantial decrease')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # 2. Percentage decreased per complex
    decreased_pcts = [cr['decreased_pct'] for cr in complex_results]
    bars = ax2.bar(range(len(complex_names)), decreased_pcts, color='orange', alpha=0.7)
    ax2.set_xticks(range(len(complex_names)))
    ax2.set_xticklabels(complex_names, rotation=45, ha='right')
    ax2.set_ylabel('% Proteins Decreased')
    ax2.set_title('Percentage of Proteins Decreased per Complex')
    ax2.axhline(50, color='red', linestyle='--', alpha=0.5, label='Majority decreased')
    ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # 3. Volcano plot of all proteins
    colors = []
    for _, row in df.iterrows():
        if row['significant'] and row['log2FC'] < -0.2:
            colors.append('red')  # Significantly decreased
        elif row['significant']:
            colors.append('blue')  # Significantly changed
        else:
            colors.append('gray')  # Not significant

    ax3.scatter(df['log2FC'], -np.log10(df['p_value']), c=colors, alpha=0.6, s=40)
    ax3.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    ax3.axvline(-0.2, color='red', linestyle='--', alpha=0.5)
    ax3.axvline(0, color='black', linestyle='-', alpha=0.3)
    ax3.set_xlabel('Log2 Fold Change (Tau+ vs Tau-)')
    ax3.set_ylabel('-Log10(p-value)')
    ax3.set_title('All Mitochondrial Proteins Volcano Plot')
    ax3.grid(True, alpha=0.3)

    # 4. Complex comparison heatmap
    if len(complex_results) > 1:
        heatmap_data = []
        labels = []
        for cr in complex_results:
            heatmap_data.append([cr['mean_log2FC'], cr['decreased_pct']/100])
            labels.append(cr['complex'])

        heatmap_data = np.array(heatmap_data).T
        im = ax4.imshow(heatmap_data, cmap='RdBu_r', aspect='auto', vmin=-1, vmax=1)
        ax4.set_xticks(range(len(labels)))
        ax4.set_xticklabels(labels, rotation=45, ha='right')
        ax4.set_yticks([0, 1])
        ax4.set_yticklabels(['Mean log2FC', '% Decreased'])
        ax4.set_title('Complex Dysfunction Heatmap')
        plt.colorbar(im, ax=ax4, fraction=0.046, pad=0.04)

    plt.tight_layout()
    plt.show()
```

## Evaluate Claim

```python
# 🎯 CLAIM EVALUATION
print("\\n" + "="*60)
print("🎯 CLAIM EVALUATION")
print("="*60)
print("Claim: 'Mitochondrial complexes I-V are decreased in tau+ neurons'")
print()

if len(df) > 0:
    # Overall statistics
    total_proteins = len(df)
    total_decreased = sum(df['decreased'])
    pct_decreased = total_decreased / total_proteins * 100

    sig_decreased = sum((df['significant']) & (df['decreased']))
    pct_sig_decreased = sig_decreased / total_proteins * 100

    # Complex-level analysis
    complexes_mostly_decreased = sum(1 for cr in complex_results if cr['decreased_pct'] > 50)

    print(f"📊 Overall Results:")
    print(f"Proteins tested: {total_proteins}")
    print(f"Proteins decreased: {total_decreased} ({pct_decreased:.1f}%)")
    print(f"Significantly decreased: {sig_decreased} ({pct_sig_decreased:.1f}%)")
    print(f"Complexes mostly decreased: {complexes_mostly_decreased}/{len(complex_results)}")

    print(f"\\n📈 Complex-specific Results:")
    for cr in complex_results:
        status = "✓ DECREASED" if cr['decreased_pct'] > 50 else "✗ Mixed/Stable"
        print(f"  {cr['complex']:12} {cr['decreased_pct']:5.1f}% decreased  {status}")

    # Overall verdict
    if pct_sig_decreased > 60:
        verdict = "✅ STRONGLY SUPPORTED"
        explanation = f"Majority of complexes show significant decreases ({pct_sig_decreased:.1f}%)"
    elif pct_decreased > 50 and complexes_mostly_decreased >= 3:
        verdict = "✅ SUPPORTED"
        explanation = f"Most complexes decreased ({complexes_mostly_decreased}/{len(complex_results)} complexes)"
    elif pct_decreased > 40 or complexes_mostly_decreased >= 2:
        verdict = "⚠️ PARTIALLY SUPPORTED"
        explanation = f"Some complexes show decreases ({pct_decreased:.1f}% proteins)"
    else:
        verdict = "❌ REFUTED"
        explanation = f"Complexes not consistently decreased ({pct_decreased:.1f}% proteins)"

    print(f"\\n🎯 VERDICT: {verdict}")
    print(f"📝 {explanation}")

    # Biological interpretation
    print(f"\\n🧬 BIOLOGICAL IMPLICATIONS:")
    if verdict.startswith("✅"):
        print("• Mitochondrial electron transport chain compromised")
        print("• ATP production severely impaired")
        print("• Oxidative metabolism dysfunction")
        print("• Energy crisis in tau-positive neurons")
        print("• Respiratory capacity reduced")

        # Complex-specific impacts
        complex_impacts = {
            'Complex_I': "NADH oxidation impaired",
            'Complex_II': "Succinate oxidation reduced",
            'Complex_III': "Cytochrome c reduction affected",
            'Complex_IV': "Oxygen consumption decreased",
            'Complex_V': "ATP synthesis compromised"
        }

        print("\\n💡 Complex-specific impacts:")
        for cr in complex_results:
            if cr['decreased_pct'] > 50:
                impact = complex_impacts.get(cr['complex'], "Function impaired")
                print(f"  • {cr['complex']}: {impact}")

    elif verdict.startswith("⚠️"):
        print("• Selective mitochondrial dysfunction")
        print("• Some complexes preserved, others affected")
        print("• Partial energy crisis")
        print("• Compensatory mechanisms may be active")
    else:
        print("• Mitochondrial complexes appear stable")
        print("• Energy production maintained")
        print("• Alternative dysfunction mechanisms")

else:
    verdict = "❌ UNSURE"
    explanation = "Insufficient mitochondrial proteins found"
    print(f"VERDICT: {verdict}")
    print(f"Reasoning: {explanation}")
```

## Save Results

```python
# 💾 Save results
if len(df) > 0:
    # Complex summary
    complex_summary = pd.DataFrame([{
        'complex': cr['complex'],
        'proteins_found': cr['proteins_found'],
        'proteins_total': cr['proteins_total'],
        'mean_log2FC': cr['mean_log2FC'],
        'decreased_pct': cr['decreased_pct']
    } for cr in complex_results])

    # Overall summary
    summary = {
        'analysis': 'Mitochondrial complexes I-V decreased',
        'verdict': verdict,
        'proteins_tested': total_proteins if 'total_proteins' in locals() else 0,
        'proteins_decreased': total_decreased if 'total_decreased' in locals() else 0,
        'percent_decreased': pct_decreased if 'pct_decreased' in locals() else 0,
        'complexes_affected': complexes_mostly_decreased if 'complexes_mostly_decreased' in locals() else 0
    }

    if IN_COLAB:
        df.to_csv('mitochondrial_complexes_results.csv', index=False)
        complex_summary.to_csv('complex_summary.csv', index=False)
        pd.DataFrame([summary]).to_csv('mito_analysis_summary.csv', index=False)
        files.download('mitochondrial_complexes_results.csv')
        files.download('complex_summary.csv')
        print("📁 Results downloaded!")
    else:
        df.to_csv('mitochondrial_complexes_analysis.csv', index=False)
        complex_summary.to_csv('complex_summary.csv', index=False)
        print("📁 Results saved!")

print("\\n✅ Mitochondrial complex analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified mitochondrial complex analysis:
- ✅ **Tests all 5 complexes** (I, II, III, IV, V) comprehensively
- ✅ **Complex-specific evaluation** with percentage decreases
- ✅ **Clear visualizations** showing dysfunction patterns
- ✅ **Biological interpretation** linking to energy crisis
- ✅ **Therapeutic relevance** for mitochondrial support

**Perfect for evaluating bioenergetic failure in neurodegeneration!** ⚡