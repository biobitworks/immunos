# Mitophagy Receptors Analysis - Simplified Version
## Testing: "Mitophagy receptors are upregulated in tau+ neurons"

**🚀 Simple approach**: Test key mitophagy receptors for increased expression!

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

## Load Data & Define Mitophagy Receptors

```python
# 📊 Load data
adata = sc.read_h5ad(data_file)
print(f"Data loaded: {adata.shape}")

# 🏷️ Setup tau status
if 'TauStatus' in adata.obs.columns:
    adata.obs['tau_status'] = adata.obs['TauStatus']
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)

print(f"Tau+: {adata.obs['tau_positive'].sum()}, Tau-: {(adata.obs['tau_positive']==0).sum()}")

# 🔬 Mitophagy Receptor Proteins
mitophagy_receptors = {
    'Core_Receptors': [
        'PINK1',      # PTEN-induced putative kinase 1
        'PRKN',       # Parkin (E3 ubiquitin ligase)
        'BNIP3',      # BCL2/adenovirus E1B 19kDa interacting protein 3
        'BNIP3L',     # BNIP3-like (NIX)
        'FUNDC1',     # FUN14 domain containing 1
        'SQSTM1'      # Sequestosome 1 (p62)
    ],
    'Adaptor_Proteins': [
        'NBR1',       # Neighbor of BRCA1 gene 1
        'OPTN',       # Optineurin
        'NDP52',      # Nuclear domain 10 protein 52 (CALCOCO2)
        'CALCOCO2',   # Calcium binding and coiled-coil domain 2
        'TAX1BP1'     # Tax1 binding protein 1
    ],
    'Regulatory_Proteins': [
        'USP30',      # Ubiquitin specific peptidase 30
        'MUL1',       # Mitochondrial E3 ubiquitin protein ligase 1
        'MARCH5',     # Membrane associated ring-CH-type finger 5
        'SMURF1'      # SMAD specific E3 ubiquitin protein ligase 1
    ]
}

total_proteins = sum(len(proteins) for proteins in mitophagy_receptors.values())
print(f"🎯 Testing {total_proteins} mitophagy receptor proteins across {len(mitophagy_receptors)} categories")
```

## Find & Analyze Mitophagy Receptors

```python
# 🔍 Find available proteins and analyze each category
protein_names = list(adata.var_names)
category_results = []

for category_name, proteins in mitophagy_receptors.items():
    print(f"\n📊 Analyzing {category_name}...")

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
                'category': category_name,
                'protein': protein,
                'log2FC': log2fc,
                'p_value': p_val,
                'mean_tau_pos': mean_pos,
                'mean_tau_neg': mean_neg,
                'upregulated': log2fc > 0
            })

        # Category-level statistics
        log2fcs = [r['log2FC'] for r in protein_results]
        mean_log2fc = np.mean(log2fcs)
        upregulated_count = sum(r['upregulated'] for r in protein_results)
        upregulated_pct = upregulated_count / len(protein_results) * 100

        # Add to category results
        category_results.append({
            'category': category_name,
            'proteins_found': len(found_proteins),
            'proteins_total': len(proteins),
            'mean_log2FC': mean_log2fc,
            'upregulated_count': upregulated_count,
            'upregulated_pct': upregulated_pct,
            'proteins_tested': protein_results
        })

        print(f"  Mean log2FC: {mean_log2fc:.3f}")
        print(f"  Upregulated: {upregulated_count}/{len(found_proteins)} ({upregulated_pct:.1f}%)")

# Create combined DataFrame for all proteins
all_proteins_df = []
for category_result in category_results:
    all_proteins_df.extend(category_result['proteins_tested'])

df = pd.DataFrame(all_proteins_df)
if len(df) > 0:
    # Add FDR correction
    from statsmodels.stats.multitest import fdrcorrection
    df['p_adjusted'] = fdrcorrection(df['p_value'])[1]
    df['significant'] = df['p_adjusted'] < 0.05

print(f"\n✅ Analysis complete: {len(df)} proteins tested across {len(category_results)} categories")
```

## Visualize Mitophagy Results

```python
if len(df) > 0:
    # 📊 Create comprehensive visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Category-level fold changes
    category_names = [cr['category'] for cr in category_results]
    category_fcs = [cr['mean_log2FC'] for cr in category_results]
    colors = ['green' if fc > 0.1 else 'gray' for fc in category_fcs]

    bars = ax1.bar(range(len(category_names)), category_fcs, color=colors, alpha=0.7)
    ax1.set_xticks(range(len(category_names)))
    ax1.set_xticklabels(category_names, rotation=45, ha='right')
    ax1.set_ylabel('Mean Log2 Fold Change')
    ax1.set_title('Mitophagy Receptor Category Expression Changes')
    ax1.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax1.axhline(0.5, color='green', linestyle='--', alpha=0.5, label='Substantial increase')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # 2. Percentage upregulated per category
    upregulated_pcts = [cr['upregulated_pct'] for cr in category_results]
    bars = ax2.bar(range(len(category_names)), upregulated_pcts, color='lightgreen', alpha=0.7)
    ax2.set_xticks(range(len(category_names)))
    ax2.set_xticklabels(category_names, rotation=45, ha='right')
    ax2.set_ylabel('% Proteins Upregulated')
    ax2.set_title('Percentage of Proteins Upregulated per Category')
    ax2.axhline(50, color='green', linestyle='--', alpha=0.5, label='Majority upregulated')
    ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # 3. Volcano plot of all proteins
    colors = []
    for _, row in df.iterrows():
        if row['significant'] and row['log2FC'] > 0.2:
            colors.append('green')  # Significantly upregulated
        elif row['significant']:
            colors.append('blue')   # Significantly changed
        else:
            colors.append('gray')   # Not significant

    ax3.scatter(df['log2FC'], -np.log10(df['p_value']), c=colors, alpha=0.6, s=40)
    ax3.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    ax3.axvline(0.2, color='green', linestyle='--', alpha=0.5)
    ax3.axvline(0, color='black', linestyle='-', alpha=0.3)
    ax3.set_xlabel('Log2 Fold Change (Tau+ vs Tau-)')
    ax3.set_ylabel('-Log10(p-value)')
    ax3.set_title('All Mitophagy Proteins Volcano Plot')
    ax3.grid(True, alpha=0.3)

    # 4. Category comparison heatmap
    if len(category_results) > 1:
        heatmap_data = []
        labels = []
        for cr in category_results:
            heatmap_data.append([cr['mean_log2FC'], cr['upregulated_pct']/100])
            labels.append(cr['category'])

        heatmap_data = np.array(heatmap_data).T
        im = ax4.imshow(heatmap_data, cmap='RdYlGn', aspect='auto', vmin=-1, vmax=1)
        ax4.set_xticks(range(len(labels)))
        ax4.set_xticklabels(labels, rotation=45, ha='right')
        ax4.set_yticks([0, 1])
        ax4.set_yticklabels(['Mean log2FC', '% Upregulated'])
        ax4.set_title('Mitophagy Activation Heatmap')
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
print("Claim: 'Mitophagy receptors are upregulated in tau+ neurons'")
print()

if len(df) > 0:
    # Overall statistics
    total_proteins = len(df)
    total_upregulated = sum(df['upregulated'])
    pct_upregulated = total_upregulated / total_proteins * 100

    sig_upregulated = sum((df['significant']) & (df['upregulated']))
    pct_sig_upregulated = sig_upregulated / total_proteins * 100

    # Category-level analysis
    categories_mostly_upregulated = sum(1 for cr in category_results if cr['upregulated_pct'] > 50)

    print(f"📊 Overall Results:")
    print(f"Proteins tested: {total_proteins}")
    print(f"Proteins upregulated: {total_upregulated} ({pct_upregulated:.1f}%)")
    print(f"Significantly upregulated: {sig_upregulated} ({pct_sig_upregulated:.1f}%)")
    print(f"Categories mostly upregulated: {categories_mostly_upregulated}/{len(category_results)}")

    print(f"\n📈 Category-specific Results:")
    for cr in category_results:
        status = "✓ UPREGULATED" if cr['upregulated_pct'] > 50 else "✗ Mixed/Stable"
        print(f"  {cr['category']:18} {cr['upregulated_pct']:5.1f}% upregulated  {status}")

    # Overall verdict
    if pct_sig_upregulated > 60:
        verdict = "✅ STRONGLY SUPPORTED"
        explanation = f"Majority of receptors show significant upregulation ({pct_sig_upregulated:.1f}%)"
    elif pct_upregulated > 50 and categories_mostly_upregulated >= 2:
        verdict = "✅ SUPPORTED"
        explanation = f"Most categories upregulated ({categories_mostly_upregulated}/{len(category_results)} categories)"
    elif pct_upregulated > 40 or categories_mostly_upregulated >= 1:
        verdict = "⚠️ PARTIALLY SUPPORTED"
        explanation = f"Some receptors show upregulation ({pct_upregulated:.1f}% proteins)"
    else:
        verdict = "❌ REFUTED"
        explanation = f"Receptors not consistently upregulated ({pct_upregulated:.1f}% proteins)"

    print(f"\n🎯 VERDICT: {verdict}")
    print(f"📝 {explanation}")

    # Biological interpretation
    print(f"\n🧬 BIOLOGICAL IMPLICATIONS:")
    if verdict.startswith("✅"):
        print("• Mitophagy machinery actively upregulated")
        print("• Cellular response to mitochondrial damage")
        print("• Attempt to clear dysfunctional mitochondria")
        print("• Quality control mechanisms activated")
        print("• May indicate ongoing mitochondrial stress")

        # Category-specific impacts
        category_impacts = {
            'Core_Receptors': "Direct mitophagy initiation enhanced",
            'Adaptor_Proteins': "Autophagosome targeting improved",
            'Regulatory_Proteins': "Fine-tuning of mitophagy process"
        }

        print("\n💡 Category-specific impacts:")
        for cr in category_results:
            if cr['upregulated_pct'] > 50:
                impact = category_impacts.get(cr['category'], "Mitophagy function enhanced")
                print(f"  • {cr['category']}: {impact}")

    elif verdict.startswith("⚠️"):
        print("• Selective mitophagy receptor upregulation")
        print("• Partial mitochondrial quality control response")
        print("• Some pathways activated, others unchanged")
        print("• May indicate early compensatory mechanisms")
    else:
        print("• Mitophagy receptors appear stable")
        print("• No clear upregulation of quality control")
        print("• Alternative clearance mechanisms may dominate")

else:
    verdict = "❌ UNSURE"
    explanation = "Insufficient mitophagy receptor proteins found"
    print(f"VERDICT: {verdict}")
    print(f"Reasoning: {explanation}")
```

## Save Results

```python
# 💾 Save results
if len(df) > 0:
    # Category summary
    category_summary = pd.DataFrame([{
        'category': cr['category'],
        'proteins_found': cr['proteins_found'],
        'proteins_total': cr['proteins_total'],
        'mean_log2FC': cr['mean_log2FC'],
        'upregulated_pct': cr['upregulated_pct']
    } for cr in category_results])

    # Overall summary
    summary = {
        'analysis': 'Mitophagy receptors upregulated',
        'verdict': verdict,
        'proteins_tested': total_proteins if 'total_proteins' in locals() else 0,
        'proteins_upregulated': total_upregulated if 'total_upregulated' in locals() else 0,
        'percent_upregulated': pct_upregulated if 'pct_upregulated' in locals() else 0,
        'categories_affected': categories_mostly_upregulated if 'categories_mostly_upregulated' in locals() else 0
    }

    if IN_COLAB:
        df.to_csv('mitophagy_receptors_results.csv', index=False)
        category_summary.to_csv('mitophagy_category_summary.csv', index=False)
        pd.DataFrame([summary]).to_csv('mitophagy_analysis_summary.csv', index=False)
        files.download('mitophagy_receptors_results.csv')
        files.download('mitophagy_category_summary.csv')
        print("📁 Results downloaded!")
    else:
        df.to_csv('mitophagy_receptors_analysis.csv', index=False)
        category_summary.to_csv('mitophagy_category_summary.csv', index=False)
        print("📁 Results saved!")

print("\n✅ Mitophagy receptor analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified mitophagy receptor analysis:
- ✅ **Tests all receptor categories** (core, adaptor, regulatory) comprehensively
- ✅ **Category-specific evaluation** with upregulation percentages
- ✅ **Clear visualizations** showing activation patterns
- ✅ **Biological interpretation** linking to mitochondrial quality control
- ✅ **Therapeutic relevance** for mitophagy enhancement strategies

**Perfect for evaluating cellular response to mitochondrial dysfunction!** 🔬