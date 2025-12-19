# Rab GTPases Analysis - Simplified Version
## Testing: "Rab GTPases show widespread trafficking dysfunction"

**🚀 Simple approach**: Test Rab GTPase family for trafficking dysfunction patterns!

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

## Load Data & Define Rab GTPases

```python
# 📊 Load data
adata = sc.read_h5ad(data_file)
print(f"Data loaded: {adata.shape}")

# 🏷️ Setup tau status
if 'TauStatus' in adata.obs.columns:
    adata.obs['tau_status'] = adata.obs['TauStatus']
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)

print(f"Tau+: {adata.obs['tau_positive'].sum()}, Tau-: {(adata.obs['tau_positive']==0).sum()}")

# 🚛 Rab GTPase Trafficking Network
rab_gtpases = {
    'Early_Endocytic': [
        'RAB4A',      # RAB4A, member RAS oncogene family
        'RAB4B',      # RAB4B, member RAS oncogene family
        'RAB5A',      # RAB5A, member RAS oncogene family
        'RAB5B',      # RAB5B, member RAS oncogene family
        'RAB5C',      # RAB5C, member RAS oncogene family
        'RAB11A',     # RAB11A, member RAS oncogene family
        'RAB11B'      # RAB11B, member RAS oncogene family
    ],
    'Late_Endocytic': [
        'RAB7A',      # RAB7A, member RAS oncogene family
        'RAB7B',      # RAB7B, member RAS oncogene family
        'RAB9A',      # RAB9A, member RAS oncogene family
        'RAB9B',      # RAB9B, member RAS oncogene family
        'RAB22A',     # RAB22A, member RAS oncogene family
        'RAB34'       # RAB34, member RAS oncogene family
    ],
    'Golgi_ER_Trafficking': [
        'RAB1A',      # RAB1A, member RAS oncogene family
        'RAB1B',      # RAB1B, member RAS oncogene family
        'RAB2A',      # RAB2A, member RAS oncogene family
        'RAB2B',      # RAB2B, member RAS oncogene family
        'RAB6A',      # RAB6A, member RAS oncogene family
        'RAB6B',      # RAB6B, member RAS oncogene family
        'RAB18',      # RAB18, member RAS oncogene family
        'RAB33A'      # RAB33A, member RAS oncogene family
    ],
    'Secretory_Pathway': [
        'RAB3A',      # RAB3A, member RAS oncogene family
        'RAB3B',      # RAB3B, member RAS oncogene family
        'RAB3C',      # RAB3C, member RAS oncogene family
        'RAB3D',      # RAB3D, member RAS oncogene family
        'RAB8A',      # RAB8A, member RAS oncogene family
        'RAB8B',      # RAB8B, member RAS oncogene family
        'RAB10',      # RAB10, member RAS oncogene family
        'RAB27A',     # RAB27A, member RAS oncogene family
        'RAB27B'      # RAB27B, member RAS oncogene family
    ],
    'Lysosomal_Trafficking': [
        'RAB24',      # RAB24, member RAS oncogene family
        'RAB32',      # RAB32, member RAS oncogene family
        'RAB38',      # RAB38, member RAS oncogene family
        'RAB39A',     # RAB39A, member RAS oncogene family
        'RAB39B'      # RAB39B, member RAS oncogene family
    ],
    'Autophagy_Related': [
        'RAB26',      # RAB26, member RAS oncogene family
        'RAB33B',     # RAB33B, member RAS oncogene family
        'RAB12',      # RAB12, member RAS oncogene family
        'RAB13',      # RAB13, member RAS oncogene family
        'RAB35'       # RAB35, member RAS oncogene family
    ],
    'Neuronal_Specific': [
        'RAB3A',      # RAB3A (synaptic vesicles)
        'RAB5A',      # RAB5A (early endosomes in neurons)
        'RAB7A',      # RAB7A (late endosomes in neurons)
        'RAB11A',     # RAB11A (recycling in neurons)
        'RAB17',      # RAB17, member RAS oncogene family
        'RAB40B',     # RAB40B, member RAS oncogene family
        'RAB40C'      # RAB40C, member RAS oncogene family
    ]
}

total_proteins = sum(len(proteins) for proteins in rab_gtpases.values())
print(f"🎯 Testing {total_proteins} Rab GTPase proteins across {len(rab_gtpases)} trafficking categories")

# Note: Remove duplicates for total count
all_rabs = set()
for proteins in rab_gtpases.values():
    all_rabs.update(proteins)
unique_rabs = len(all_rabs)
print(f"📊 Unique Rab GTPases: {unique_rabs}")
```

## Find & Analyze Rab GTPases

```python
# 🔍 Find available proteins and analyze each trafficking category
protein_names = list(adata.var_names)
category_results = []

for category_name, proteins in rab_gtpases.items():
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
                'dysfunctional': abs(log2fc) > 0.2  # Consider any large change as dysfunction
            })

        # Category-level statistics
        log2fcs = [r['log2FC'] for r in protein_results]
        mean_log2fc = np.mean(log2fcs)
        std_log2fc = np.std(log2fcs)
        dysfunctional_count = sum(r['dysfunctional'] for r in protein_results)
        dysfunctional_pct = dysfunctional_count / len(protein_results) * 100

        # Determine trafficking status
        if abs(mean_log2fc) > 0.3:
            if mean_log2fc > 0:
                status = "Upregulated"
            else:
                status = "Downregulated"
        elif std_log2fc > 0.3:
            status = "Variable"
        else:
            status = "Stable"

        # Add to category results
        category_results.append({
            'category': category_name,
            'proteins_found': len(found_proteins),
            'proteins_total': len(proteins),
            'mean_log2FC': mean_log2fc,
            'std_log2FC': std_log2fc,
            'dysfunctional_count': dysfunctional_count,
            'dysfunctional_pct': dysfunctional_pct,
            'trafficking_status': status,
            'proteins_tested': protein_results
        })

        print(f"  Mean log2FC: {mean_log2fc:.3f} (±{std_log2fc:.3f})")
        print(f"  Status: {status}")
        print(f"  Dysfunctional: {dysfunctional_count}/{len(found_proteins)} ({dysfunctional_pct:.1f}%)")

# Create combined DataFrame for all proteins (remove duplicates)
all_proteins_df = []
seen_proteins = set()

for category_result in category_results:
    for protein_result in category_result['proteins_tested']:
        if protein_result['protein'] not in seen_proteins:
            all_proteins_df.append(protein_result)
            seen_proteins.add(protein_result['protein'])

df = pd.DataFrame(all_proteins_df)
if len(df) > 0:
    # Add FDR correction
    from statsmodels.stats.multitest import fdrcorrection
    df['p_adjusted'] = fdrcorrection(df['p_value'])[1]
    df['significant'] = df['p_adjusted'] < 0.05

print(f"\n✅ Analysis complete: {len(df)} unique Rab proteins tested across {len(category_results)} categories")
```

## Analyze Trafficking Network Dysfunction

```python
if len(df) > 0:
    # 🚛 Network-wide analysis
    print(f"\n🚛 TRAFFICKING NETWORK ANALYSIS:")

    # Overall Rab dysfunction
    total_dysfunctional = sum(df['dysfunctional'])
    pct_dysfunctional = total_dysfunctional / len(df) * 100

    sig_dysfunctional = sum((df['significant']) & (df['dysfunctional']))
    pct_sig_dysfunctional = sig_dysfunctional / len(df) * 100

    print(f"Total Rab proteins: {len(df)}")
    print(f"Dysfunctional Rabs: {total_dysfunctional} ({pct_dysfunctional:.1f}%)")
    print(f"Significantly dysfunctional: {sig_dysfunctional} ({pct_sig_dysfunctional:.1f}%)")

    # Category-specific analysis
    categories_affected = sum(1 for cr in category_results if cr['dysfunctional_pct'] > 50)
    print(f"Categories mostly affected: {categories_affected}/{len(category_results)}")

    # Distribution analysis
    upregulated = sum(1 for _, row in df.iterrows() if row['log2FC'] > 0.2)
    downregulated = sum(1 for _, row in df.iterrows() if row['log2FC'] < -0.2)
    stable = len(df) - upregulated - downregulated

    print(f"\nExpression patterns:")
    print(f"  Upregulated: {upregulated} ({upregulated/len(df)*100:.1f}%)")
    print(f"  Downregulated: {downregulated} ({downregulated/len(df)*100:.1f}%)")
    print(f"  Stable: {stable} ({stable/len(df)*100:.1f}%)")

    # Pathway-specific dysfunction scores
    pathway_scores = {}
    for cr in category_results:
        pathway_scores[cr['category']] = {
            'dysfunction_score': cr['dysfunctional_pct'],
            'expression_change': abs(cr['mean_log2FC']),
            'variability': cr['std_log2FC']
        }

    print(f"\nTop dysfunctional pathways:")
    sorted_pathways = sorted(pathway_scores.items(),
                           key=lambda x: x[1]['dysfunction_score'], reverse=True)
    for pathway, scores in sorted_pathways[:3]:
        print(f"  {pathway}: {scores['dysfunction_score']:.1f}% dysfunction")
```

## Visualize Rab GTPase Dysfunction

```python
if len(df) > 0:
    # 📊 Create comprehensive visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Category-level dysfunction
    category_names = [cr['category'] for cr in category_results]
    category_dysfunc = [cr['dysfunctional_pct'] for cr in category_results]

    # Color by trafficking status
    status_colors = {
        'Upregulated': 'green',
        'Downregulated': 'red',
        'Variable': 'orange',
        'Stable': 'gray'
    }
    colors = [status_colors[cr['trafficking_status']] for cr in category_results]

    bars = ax1.bar(range(len(category_names)), category_dysfunc, color=colors, alpha=0.7)
    ax1.set_xticks(range(len(category_names)))
    ax1.set_xticklabels(category_names, rotation=45, ha='right')
    ax1.set_ylabel('% Rabs Dysfunctional')
    ax1.set_title('Rab GTPase Trafficking Dysfunction by Category')
    ax1.axhline(50, color='red', linestyle='--', alpha=0.5, label='Majority dysfunctional')
    ax1.set_ylim(0, 100)
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Add status legend
    for status, color in status_colors.items():
        ax1.scatter([], [], c=color, alpha=0.7, s=60, label=status)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # 2. Expression change distribution
    category_fcs = [cr['mean_log2FC'] for cr in category_results]
    category_stds = [cr['std_log2FC'] for cr in category_results]

    bars = ax2.bar(range(len(category_names)), category_fcs,
                  yerr=category_stds, color=colors, alpha=0.7, capsize=5)
    ax2.set_xticks(range(len(category_names)))
    ax2.set_xticklabels(category_names, rotation=45, ha='right')
    ax2.set_ylabel('Mean Log2 Fold Change')
    ax2.set_title('Rab GTPase Expression Changes')
    ax2.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax2.axhline(0.2, color='green', linestyle='--', alpha=0.5)
    ax2.axhline(-0.2, color='red', linestyle='--', alpha=0.5)
    ax2.grid(True, alpha=0.3)

    # 3. Volcano plot of all Rabs
    colors_volcano = []
    for _, row in df.iterrows():
        if row['significant'] and abs(row['log2FC']) > 0.2:
            colors_volcano.append('red')     # Significantly dysfunctional
        elif row['significant']:
            colors_volcano.append('blue')    # Significantly changed
        else:
            colors_volcano.append('gray')    # Not significant

    ax3.scatter(df['log2FC'], -np.log10(df['p_value']), c=colors_volcano, alpha=0.6, s=40)
    ax3.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    ax3.axvline(0.2, color='green', linestyle='--', alpha=0.5)
    ax3.axvline(-0.2, color='red', linestyle='--', alpha=0.5)
    ax3.axvline(0, color='black', linestyle='-', alpha=0.3)
    ax3.set_xlabel('Log2 Fold Change (Tau+ vs Tau-)')
    ax3.set_ylabel('-Log10(p-value)')
    ax3.set_title('All Rab GTPases Volcano Plot')
    ax3.grid(True, alpha=0.3)

    # 4. Network dysfunction summary
    summary_data = [
        ['Upregulated', upregulated],
        ['Downregulated', downregulated],
        ['Stable', stable]
    ]

    sizes = [data[1] for data in summary_data]
    labels = [f"{data[0]}\n({data[1]} Rabs)" for data in summary_data]
    colors_pie = ['green', 'red', 'gray']

    ax4.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.1f%%', startangle=90)
    ax4.set_title('Rab GTPase Network Dysfunction Summary')

    plt.tight_layout()
    plt.show()
```

## Evaluate Claim

```python
# 🎯 CLAIM EVALUATION
print("\n" + "="*60)
print("🎯 CLAIM EVALUATION")
print("="*60)
print("Claim: 'Rab GTPases show widespread trafficking dysfunction'")
print()

if len(df) > 0:
    print(f"📊 Overall Results:")
    print(f"Rab proteins tested: {len(df)}")
    print(f"Dysfunctional Rabs: {total_dysfunctional} ({pct_dysfunctional:.1f}%)")
    print(f"Significantly dysfunctional: {sig_dysfunctional} ({pct_sig_dysfunctional:.1f}%)")
    print(f"Categories affected: {categories_affected}/{len(category_results)}")

    print(f"\n📈 Category-specific Results:")
    for cr in category_results:
        status_symbol = {
            'Upregulated': '↗️',
            'Downregulated': '↘️',
            'Variable': '↕️',
            'Stable': '→'
        }
        symbol = status_symbol.get(cr['trafficking_status'], '?')
        print(f"  {cr['category']:20} {cr['trafficking_status']:12} {symbol} ({cr['dysfunctional_pct']:4.1f}% dysfunctional)")

    # Overall verdict
    if pct_sig_dysfunctional > 60:
        verdict = "✅ STRONGLY SUPPORTED"
        explanation = f"Majority of Rab GTPases significantly dysfunctional ({pct_sig_dysfunctional:.1f}%)"
    elif pct_dysfunctional > 50 and categories_affected >= 4:
        verdict = "✅ SUPPORTED"
        explanation = f"Widespread dysfunction across {categories_affected} trafficking categories ({pct_dysfunctional:.1f}%)"
    elif pct_dysfunctional > 40 or categories_affected >= 3:
        verdict = "⚠️ PARTIALLY SUPPORTED"
        explanation = f"Some trafficking dysfunction ({pct_dysfunctional:.1f}% Rabs, {categories_affected} categories)"
    else:
        verdict = "❌ REFUTED"
        explanation = f"No widespread trafficking dysfunction ({pct_dysfunctional:.1f}% dysfunctional)"

    print(f"\n🎯 VERDICT: {verdict}")
    print(f"📝 {explanation}")

    # Biological interpretation
    print(f"\n🧬 BIOLOGICAL IMPLICATIONS:")
    if verdict.startswith("✅"):
        print("• Vesicular trafficking system severely compromised")
        print("• Membrane dynamics and cargo sorting disrupted")
        print("• Organellar communication networks impaired")
        print("• Synaptic function and neuronal connectivity affected")
        print("• Cellular compartmentalization breakdown")

        # Category-specific impacts
        category_impacts = {
            'Early_Endocytic': "Receptor recycling and uptake impaired",
            'Late_Endocytic': "Cargo sorting and degradation affected",
            'Golgi_ER_Trafficking': "Secretory pathway dysfunction",
            'Secretory_Pathway': "Vesicle release and secretion compromised",
            'Lysosomal_Trafficking': "Degradative pathway impaired",
            'Autophagy_Related': "Autophagosome trafficking disrupted",
            'Neuronal_Specific': "Synaptic transmission and plasticity affected"
        }

        print("\n💡 Category-specific impacts:")
        for cr in category_results:
            if cr['dysfunctional_pct'] > 50:
                impact = category_impacts.get(cr['category'], "Trafficking function impaired")
                print(f"  • {cr['category']}: {impact}")

    elif verdict.startswith("⚠️"):
        print("• Selective trafficking pathway dysfunction")
        print("• Some vesicular transport preserved")
        print("• Partial membrane dynamics impairment")
        print("• Compensatory mechanisms may be active")
    else:
        print("• Rab GTPase network appears functional")
        print("• Vesicular trafficking preserved")
        print("• Normal membrane dynamics")

else:
    verdict = "❌ UNSURE"
    explanation = "Insufficient Rab GTPase proteins found"
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
        'std_log2FC': cr['std_log2FC'],
        'dysfunctional_pct': cr['dysfunctional_pct'],
        'trafficking_status': cr['trafficking_status']
    } for cr in category_results])

    # Network summary
    network_summary = pd.DataFrame([{
        'network_metric': 'Total Rabs',
        'value': len(df)
    }, {
        'network_metric': 'Dysfunctional Rabs',
        'value': total_dysfunctional
    }, {
        'network_metric': 'Percent Dysfunctional',
        'value': pct_dysfunctional
    }, {
        'network_metric': 'Categories Affected',
        'value': categories_affected
    }, {
        'network_metric': 'Upregulated',
        'value': upregulated
    }, {
        'network_metric': 'Downregulated',
        'value': downregulated
    }])

    # Overall summary
    summary = {
        'analysis': 'Rab GTPases show widespread trafficking dysfunction',
        'verdict': verdict,
        'rabs_tested': len(df),
        'rabs_dysfunctional': total_dysfunctional,
        'percent_dysfunctional': pct_dysfunctional,
        'categories_affected': categories_affected,
        'upregulated_rabs': upregulated,
        'downregulated_rabs': downregulated
    }

    if IN_COLAB:
        df.to_csv('rab_gtpases_results.csv', index=False)
        category_summary.to_csv('rab_category_summary.csv', index=False)
        network_summary.to_csv('rab_network_summary.csv', index=False)
        pd.DataFrame([summary]).to_csv('rab_gtpases_analysis_summary.csv', index=False)
        files.download('rab_gtpases_results.csv')
        files.download('rab_category_summary.csv')
        files.download('rab_network_summary.csv')
        print("📁 Results downloaded!")
    else:
        df.to_csv('rab_gtpases_analysis.csv', index=False)
        category_summary.to_csv('rab_category_summary.csv', index=False)
        network_summary.to_csv('rab_network_summary.csv', index=False)
        print("📁 Results saved!")

print("\n✅ Rab GTPase trafficking dysfunction analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified Rab GTPase analysis:
- ✅ **Tests all trafficking categories** (early endocytic, late endocytic, Golgi-ER, secretory, lysosomal, autophagy, neuronal)
- ✅ **Network-wide dysfunction assessment** with category-specific patterns
- ✅ **Clear visualizations** showing trafficking system breakdown
- ✅ **Biological interpretation** linking to vesicular transport failure
- ✅ **Therapeutic relevance** for trafficking enhancement strategies

**Perfect for evaluating membrane trafficking system dysfunction in neurodegeneration!** 🚛