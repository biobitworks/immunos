# Cristae Organization Analysis - Simplified Version
## Testing: "Cristae organization proteins are disrupted in tau pathology"

**🚀 Simple approach**: Test proteins involved in mitochondrial cristae structure!

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

## Load Data & Define Cristae Proteins

```python
# 📊 Load data
adata = sc.read_h5ad(data_file)
print(f"Data loaded: {adata.shape}")

# 🏷️ Setup tau status
if 'TauStatus' in adata.obs.columns:
    adata.obs['tau_status'] = adata.obs['TauStatus']
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)

print(f"Tau+: {adata.obs['tau_positive'].sum()}, Tau-: {(adata.obs['tau_positive']==0).sum()}")

# 🏗️ Cristae Organization Proteins
cristae_proteins = {
    'MICOS_complex': [
        'MIC60',     # IMMT - Inner membrane protein
        'MIC19',     # CHCHD3 - MICOS subunit
        'MIC25',     # CHCHD6 - MICOS subunit
        'MIC27',     # APOOL - MICOS subunit
        'SAMM50',    # SAM complex component
        'MIC10',     # MICOS10
        'MIC13'      # QIL1
    ],
    'Cardiolipin_synthesis': [
        'CLS1',      # Cardiolipin synthase
        'TAZ',       # Tafazzin - cardiolipin remodeling
        'PGS1',      # Phosphatidylglycerophosphate synthase
        'PTPMT1',    # Protein tyrosine phosphatase
        'CRLS1'      # Cardiolipin synthase 1
    ],
    'Cristae_junctions': [
        'OPA1',      # Optic atrophy 1 - cristae maintenance
        'MFND1',     # MiD49 - fission factor
        'MFND2',     # MiD51 - fission factor
        'DNM1L',     # DRP1 - dynamin-related protein
        'FIS1',      # Fission protein 1
        'MFF'        # Mitochondrial fission factor
    ],
    'Inner_membrane_structure': [
        'IMMT',      # Inner membrane protein 60
        'CHCHD3',    # Coiled-coil-helix domain 3
        'CHCHD6',    # Coiled-coil-helix domain 6
        'CHCHD10',   # Coiled-coil-helix domain 10
        'TIMM50',    # Translocase subunit
        'TIMM23',    # Translocase subunit
        'TIMM17A'    # Translocase subunit
    ],
    'Contact_sites': [
        'VDAC1',     # Voltage-dependent anion channel
        'VDAC2',     # Voltage-dependent anion channel
        'VDAC3',     # Voltage-dependent anion channel
        'ANT1',      # Adenine nucleotide translocase
        'ANT2',      # Adenine nucleotide translocase
        'TSPO',      # Translocator protein
        'ACAD9'      # Acyl-CoA dehydrogenase
    ]
}

total_proteins = sum(len(proteins) for proteins in cristae_proteins.values())
print(f"🎯 Testing {total_proteins} cristae organization proteins across {len(cristae_proteins)} categories")
```

## Find & Analyze Cristae Proteins

```python
# 🔍 Find available proteins and analyze each category
protein_names = list(adata.var_names)
category_results = []

for category, proteins in cristae_proteins.items():
    print(f"\\n📊 Analyzing {category}...")

    # Find available proteins (try multiple name formats)
    found_proteins = []
    for protein in proteins:
        if protein in protein_names:
            found_proteins.append(protein)
        else:
            # Try alternative names
            alt_names = []
            if protein == 'MIC60':
                alt_names = ['IMMT']
            elif protein == 'MIC19':
                alt_names = ['CHCHD3']
            elif protein == 'MIC25':
                alt_names = ['CHCHD6']
            elif protein == 'ANT1':
                alt_names = ['SLC25A4']
            elif protein == 'ANT2':
                alt_names = ['SLC25A5']

            for alt_name in alt_names:
                if alt_name in protein_names:
                    found_proteins.append(alt_name)
                    break

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
                'category': category,
                'protein': protein,
                'log2FC': log2fc,
                'p_value': p_val,
                'mean_tau_pos': mean_pos,
                'mean_tau_neg': mean_neg,
                'disrupted': abs(log2fc) > 0.2  # Consider >20% change as disrupted
            })

        # Category-level statistics
        log2fcs = [r['log2FC'] for r in protein_results]
        abs_log2fcs = [abs(fc) for fc in log2fcs]
        mean_abs_log2fc = np.mean(abs_log2fcs)
        disrupted_count = sum(r['disrupted'] for r in protein_results)
        disrupted_pct = disrupted_count / len(protein_results) * 100

        # Add to category results
        category_results.append({
            'category': category,
            'proteins_found': len(found_proteins),
            'proteins_total': len(proteins),
            'mean_abs_log2FC': mean_abs_log2fc,
            'disrupted_count': disrupted_count,
            'disrupted_pct': disrupted_pct,
            'proteins_tested': protein_results
        })

        print(f"  Mean |log2FC|: {mean_abs_log2fc:.3f}")
        print(f"  Disrupted: {disrupted_count}/{len(found_proteins)} ({disrupted_pct:.1f}%)")

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

print(f"\\n✅ Analysis complete: {len(df)} proteins tested across {len(category_results)} categories")
```

## Visualize Cristae Results

```python
if len(df) > 0:
    # 📊 Create comprehensive visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Category-level disruption
    category_names = [cr['category'] for cr in category_results]
    disrupted_pcts = [cr['disrupted_pct'] for cr in category_results]
    colors = ['red' if pct > 50 else 'orange' if pct > 25 else 'gray' for pct in disrupted_pcts]

    bars = ax1.bar(range(len(category_names)), disrupted_pcts, color=colors, alpha=0.7)
    ax1.set_xticks(range(len(category_names)))
    ax1.set_xticklabels(category_names, rotation=45, ha='right')
    ax1.set_ylabel('% Proteins Disrupted')
    ax1.set_title('Cristae Organization Disruption by Category')
    ax1.axhline(50, color='red', linestyle='--', alpha=0.5, label='Majority disrupted')
    ax1.set_ylim(0, 100)
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # 2. Effect sizes per category
    effect_sizes = [cr['mean_abs_log2FC'] for cr in category_results]
    bars = ax2.bar(range(len(category_names)), effect_sizes, color='purple', alpha=0.7)
    ax2.set_xticks(range(len(category_names)))
    ax2.set_xticklabels(category_names, rotation=45, ha='right')
    ax2.set_ylabel('Mean |Log2 Fold Change|')
    ax2.set_title('Effect Sizes by Category')
    ax2.axhline(0.5, color='red', linestyle='--', alpha=0.5, label='Large effect')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # 3. Individual protein volcano plot
    colors = []
    for _, row in df.iterrows():
        if row['significant'] and abs(row['log2FC']) > 0.5:
            colors.append('red')  # Significantly disrupted
        elif row['significant']:
            colors.append('orange')  # Significantly changed
        elif abs(row['log2FC']) > 0.5:
            colors.append('blue')  # Large effect but not significant
        else:
            colors.append('gray')  # Not disrupted

    ax3.scatter(df['log2FC'], -np.log10(df['p_value']), c=colors, alpha=0.6, s=50)
    ax3.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    ax3.axvline(0.5, color='red', linestyle='--', alpha=0.5)
    ax3.axvline(-0.5, color='red', linestyle='--', alpha=0.5)
    ax3.axvline(0, color='black', linestyle='-', alpha=0.3)
    ax3.set_xlabel('Log2 Fold Change (Tau+ vs Tau-)')
    ax3.set_ylabel('-Log10(p-value)')
    ax3.set_title('Cristae Proteins Volcano Plot')
    ax3.grid(True, alpha=0.3)

    # 4. Category comparison heatmap
    if len(category_results) > 1:
        heatmap_data = []
        labels = []
        for cr in category_results:
            heatmap_data.append([cr['mean_abs_log2FC'], cr['disrupted_pct']/100])
            labels.append(cr['category'].replace('_', '\\n'))

        heatmap_data = np.array(heatmap_data).T
        im = ax4.imshow(heatmap_data, cmap='Reds', aspect='auto', vmin=0, vmax=1)
        ax4.set_xticks(range(len(labels)))
        ax4.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
        ax4.set_yticks([0, 1])
        ax4.set_yticklabels(['Effect Size', '% Disrupted'])
        ax4.set_title('Category Disruption Heatmap')
        plt.colorbar(im, ax=ax4, fraction=0.046, pad=0.04)

    plt.tight_layout()
    plt.show()

    # Print detailed results
    if df['significant'].any():
        print(f"\\n🔍 Most Disrupted Cristae Proteins:")
        sig_sorted = df[df['significant']].sort_values('p_adjusted')
        for _, row in sig_sorted.head(5).iterrows():
            direction = "↑" if row['log2FC'] > 0 else "↓"
            print(f"  {row['protein']:12} ({row['category']}) {direction} {abs(row['log2FC']):.3f} (FDR={row['p_adjusted']:.2e})")
```

## Evaluate Claim

```python
# 🎯 CLAIM EVALUATION
print("\\n" + "="*60)
print("🎯 CLAIM EVALUATION")
print("="*60)
print("Claim: 'Cristae organization proteins are disrupted in tau pathology'")
print()

if len(df) > 0:
    # Overall statistics
    total_proteins = len(df)
    total_disrupted = sum(df['disrupted'])
    pct_disrupted = total_disrupted / total_proteins * 100

    sig_disrupted = sum((df['significant']) & (df['disrupted']))
    pct_sig_disrupted = sig_disrupted / total_proteins * 100

    # Category-level analysis
    categories_mostly_disrupted = sum(1 for cr in category_results if cr['disrupted_pct'] > 40)
    mean_effect_size = df['log2FC'].abs().mean()

    print(f"📊 Overall Results:")
    print(f"Proteins tested: {total_proteins}")
    print(f"Proteins disrupted: {total_disrupted} ({pct_disrupted:.1f}%)")
    print(f"Significantly disrupted: {sig_disrupted} ({pct_sig_disrupted:.1f}%)")
    print(f"Categories affected: {categories_mostly_disrupted}/{len(category_results)}")
    print(f"Mean effect size: {mean_effect_size:.3f}")

    print(f"\\n📈 Category-specific Results:")
    for cr in category_results:
        status = "✓ DISRUPTED" if cr['disrupted_pct'] > 40 else "✗ Stable"
        print(f"  {cr['category']:20} {cr['disrupted_pct']:5.1f}% disrupted  {status}")

    # Overall verdict
    if pct_sig_disrupted > 50:
        verdict = "✅ STRONGLY SUPPORTED"
        explanation = f"Majority of cristae proteins significantly disrupted ({pct_sig_disrupted:.1f}%)"
    elif pct_disrupted > 40 and categories_mostly_disrupted >= 3:
        verdict = "✅ SUPPORTED"
        explanation = f"Cristae organization widely disrupted ({categories_mostly_disrupted}/{len(category_results)} categories)"
    elif pct_disrupted > 30 or categories_mostly_disrupted >= 2:
        verdict = "⚠️ PARTIALLY SUPPORTED"
        explanation = f"Some cristae disruption detected ({pct_disrupted:.1f}% proteins)"
    else:
        verdict = "❌ REFUTED"
        explanation = f"Cristae organization largely preserved ({pct_disrupted:.1f}% disrupted)"

    print(f"\\n🎯 VERDICT: {verdict}")
    print(f"📝 {explanation}")

    # Biological interpretation
    print(f"\\n🧬 BIOLOGICAL IMPLICATIONS:")
    if verdict.startswith("✅"):
        print("• Mitochondrial ultrastructure severely compromised")
        print("• Cristae junctions and membranes disrupted")
        print("• Respiratory complex organization affected")
        print("• Cardiolipin-dependent processes impaired")
        print("• Mitochondrial contact sites dysfunctional")

        # Category-specific impacts
        category_impacts = {
            'MICOS_complex': "Cristae junction formation impaired",
            'Cardiolipin_synthesis': "Membrane composition altered",
            'Cristae_junctions': "Inner membrane organization disrupted",
            'Inner_membrane_structure': "Protein import/assembly affected",
            'Contact_sites': "Organelle communication compromised"
        }

        print("\\n💡 Category-specific impacts:")
        for cr in category_results:
            if cr['disrupted_pct'] > 40:
                impact = category_impacts.get(cr['category'], "Function impaired")
                print(f"  • {cr['category']}: {impact}")

        print("\\n⚡ Functional consequences:")
        print("• Reduced respiratory efficiency")
        print("• Impaired calcium buffering")
        print("• Decreased membrane potential")
        print("• Enhanced ROS production")

    elif verdict.startswith("⚠️"):
        print("• Selective cristae dysfunction")
        print("• Some structural elements preserved")
        print("• Partial mitochondrial reorganization")
        print("• Compensatory mechanisms may be active")
    else:
        print("• Cristae structure appears preserved")
        print("• Mitochondrial ultrastructure maintained")
        print("• Alternative dysfunction mechanisms")

else:
    verdict = "❌ UNSURE"
    explanation = "Insufficient cristae proteins found"
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
        'mean_abs_log2FC': cr['mean_abs_log2FC'],
        'disrupted_pct': cr['disrupted_pct']
    } for cr in category_results])

    # Overall summary
    summary = {
        'analysis': 'Cristae organization proteins disrupted',
        'verdict': verdict,
        'proteins_tested': total_proteins if 'total_proteins' in locals() else 0,
        'proteins_disrupted': total_disrupted if 'total_disrupted' in locals() else 0,
        'percent_disrupted': pct_disrupted if 'pct_disrupted' in locals() else 0,
        'categories_affected': categories_mostly_disrupted if 'categories_mostly_disrupted' in locals() else 0
    }

    if IN_COLAB:
        df.to_csv('cristae_organization_results.csv', index=False)
        category_summary.to_csv('cristae_category_summary.csv', index=False)
        pd.DataFrame([summary]).to_csv('cristae_analysis_summary.csv', index=False)
        files.download('cristae_organization_results.csv')
        files.download('cristae_category_summary.csv')
        print("📁 Results downloaded!")
    else:
        df.to_csv('cristae_organization_analysis.csv', index=False)
        category_summary.to_csv('cristae_category_summary.csv', index=False)
        print("📁 Results saved!")

print("\\n✅ Cristae organization analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified cristae organization analysis:
- ✅ **Tests 5 key categories** of cristae proteins (MICOS, cardiolipin, junctions, etc.)
- ✅ **Evaluates structural disruption** with clear thresholds
- ✅ **Links to ultrastructure** and respiratory function
- ✅ **Category-specific analysis** showing selective vulnerability
- ✅ **Clinical relevance** for mitochondrial diseases

**Perfect for studying mitochondrial ultrastructural dysfunction!** 🏗️