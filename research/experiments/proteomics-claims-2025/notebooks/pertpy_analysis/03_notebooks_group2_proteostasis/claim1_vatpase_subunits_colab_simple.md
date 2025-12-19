# V-ATPase Analysis - Simplified Version
## Testing: "V-ATPase subunits show differential expression patterns"

**🚀 Simple approach**: Test lysosomal acidification machinery dysfunction!

---

## Setup & Data Loading

```python
# 🔧 Setup
import os
IN_COLAB = 'google.colab' in str(get_ipython())

if IN_COLAB:
    !pip install -q pertpy scanpy pandas numpy matplotlib seaborn scipy
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

## Load Data & Define V-ATPase Subunits

```python
# 📊 Load data
adata = sc.read_h5ad(data_file)
print(f"Data loaded: {adata.shape}")

# 🏷️ Setup tau status
if 'TauStatus' in adata.obs.columns:
    adata.obs['tau_status'] = adata.obs['TauStatus']
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)

print(f"Tau+: {adata.obs['tau_positive'].sum()}, Tau-: {(adata.obs['tau_positive']==0).sum()}")

# 🔋 V-ATPase subunits (key components)
vatpase_subunits = {
    'V0_membrane': [
        'ATP6V0A1',  # a1 - neuronal specific
        'ATP6V0A2',  # a2 - ubiquitous
        'ATP6V0B',   # b subunit
        'ATP6V0C',   # c - proton channel
        'ATP6V0D1',  # d1
        'ATP6V0E1',  # e1
        'ATP6AP1',   # Accessory protein 1
        'ATP6AP2'    # Accessory protein 2
    ],
    'V1_cytoplasmic': [
        'ATP6V1A',   # A - catalytic
        'ATP6V1B1',  # B1 - kidney
        'ATP6V1B2',  # B2 - brain
        'ATP6V1C1',  # C1 - ubiquitous
        'ATP6V1D',   # D - central rotor
        'ATP6V1E1',  # E1 - ubiquitous
        'ATP6V1F',   # F - central rotor
        'ATP6V1G1',  # G1 - ubiquitous
        'ATP6V1G2',  # G2 - brain enriched
        'ATP6V1H'    # H - regulatory
    ]
}

all_vatpase = vatpase_subunits['V0_membrane'] + vatpase_subunits['V1_cytoplasmic']
print(f"🎯 Testing {len(all_vatpase)} V-ATPase subunits")
print(f"  - V0 domain (membrane): {len(vatpase_subunits['V0_membrane'])}")
print(f"  - V1 domain (cytoplasmic): {len(vatpase_subunits['V1_cytoplasmic'])}")
```

## Find Available V-ATPase Proteins

```python
# 🔍 Search for V-ATPase proteins
protein_names = list(adata.var_names)
found_vatpase = {'V0_membrane': [], 'V1_cytoplasmic': []}
missing_vatpase = []

print("🔍 Searching for V-ATPase proteins...")

for domain, proteins in vatpase_subunits.items():
    for protein in proteins:
        # Try exact match first
        if protein in protein_names:
            found_vatpase[domain].append(protein)
        else:
            # Try partial matches (case-insensitive)
            matches = [p for p in protein_names if protein.upper() in p.upper()]
            if matches:
                found_vatpase[domain].append(matches[0])
            else:
                # Try without ATP6 prefix
                short_name = protein.replace('ATP6', '').replace('ATP', '')
                matches = [p for p in protein_names if short_name in p.upper()]
                if matches:
                    found_vatpase[domain].append(matches[0])
                else:
                    missing_vatpase.append((domain, protein))

# Summary
total_found = len(found_vatpase['V0_membrane']) + len(found_vatpase['V1_cytoplasmic'])
print(f"\\n📊 V-ATPase Search Results:")
print(f"V0 domain: {len(found_vatpase['V0_membrane'])}/{len(vatpase_subunits['V0_membrane'])} found")
print(f"V1 domain: {len(found_vatpase['V1_cytoplasmic'])}/{len(vatpase_subunits['V1_cytoplasmic'])} found")
print(f"Total: {total_found}/{len(all_vatpase)} ({total_found/len(all_vatpase)*100:.1f}%)")

if total_found > 0:
    print(f"✅ Sufficient V-ATPase proteins for analysis!")
    all_found = found_vatpase['V0_membrane'] + found_vatpase['V1_cytoplasmic']
else:
    print(f"❌ No V-ATPase proteins found!")
    all_found = []
```

## Analyze V-ATPase Expression

```python
if total_found > 0:
    # 🧮 Differential expression analysis
    results = []

    for protein in all_found:
        if protein in protein_names:
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

            # Determine domain
            if protein in found_vatpase['V0_membrane']:
                domain = 'V0'
            else:
                domain = 'V1'

            results.append({
                'protein': protein,
                'domain': domain,
                'log2FC': log2fc,
                'p_value': p_val,
                'mean_tau_pos': mean_pos,
                'mean_tau_neg': mean_neg
            })

    # Convert to DataFrame and add FDR correction
    df = pd.DataFrame(results)
    from statsmodels.stats.multitest import fdrcorrection
    df['p_adjusted'] = fdrcorrection(df['p_value'])[1]
    df['significant'] = df['p_adjusted'] < 0.05

    print(f"\\n📈 V-ATPASE ANALYSIS RESULTS")
    print(f"="*40)
    print(f"Proteins analyzed: {len(df)}")
    print(f"Significant (FDR < 0.05): {df['significant'].sum()}")

    # Domain-specific results
    for domain in ['V0', 'V1']:
        domain_df = df[df['domain'] == domain]
        if len(domain_df) > 0:
            n_sig = domain_df['significant'].sum()
            print(f"{domain} domain: {n_sig}/{len(domain_df)} significant")

else:
    print("⚠️ Cannot analyze - no V-ATPase proteins found")
    df = pd.DataFrame()
```

## Visualize Results

```python
if len(df) > 0:
    # 📊 Create comprehensive visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Volcano plot
    colors = []
    for _, row in df.iterrows():
        if row['significant']:
            colors.append('red' if row['domain'] == 'V0' else 'blue')
        else:
            colors.append('gray')

    ax1.scatter(df['log2FC'], -np.log10(df['p_value']), c=colors, alpha=0.7, s=60)
    ax1.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    ax1.axvline(0, color='black', linestyle='-', alpha=0.3)
    ax1.set_xlabel('Log2 Fold Change (Tau+ vs Tau-)')
    ax1.set_ylabel('-Log10(p-value)')
    ax1.set_title('V-ATPase Subunits Volcano Plot')
    ax1.grid(True, alpha=0.3)

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', alpha=0.7, label='V0 domain (significant)'),
        Patch(facecolor='blue', alpha=0.7, label='V1 domain (significant)'),
        Patch(facecolor='gray', alpha=0.7, label='Not significant')
    ]
    ax1.legend(handles=legend_elements, loc='upper right', fontsize=9)

    # 2. Domain comparison
    domain_stats = df.groupby('domain').agg({
        'log2FC': 'mean',
        'significant': 'sum',
        'protein': 'count'
    }).rename(columns={'protein': 'total'})

    x_pos = np.arange(len(domain_stats))
    bars = ax2.bar(x_pos, domain_stats['log2FC'],
                   color=['red', 'blue'], alpha=0.7)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels(domain_stats.index)
    ax2.set_ylabel('Mean Log2 Fold Change')
    ax2.set_title('Mean Expression Change by Domain')
    ax2.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax2.grid(True, alpha=0.3)

    # Add value labels
    for bar, val in zip(bars, domain_stats['log2FC']):
        ax2.text(bar.get_x() + bar.get_width()/2, val + 0.01,
                f'{val:.3f}', ha='center', va='bottom', fontweight='bold')

    # 3. Significance by domain
    ax3.bar(x_pos, domain_stats['significant'],
            color=['red', 'blue'], alpha=0.7)
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(domain_stats.index)
    ax3.set_ylabel('Significant Proteins')
    ax3.set_title('Differential Expression by Domain')
    ax3.grid(True, alpha=0.3)

    # 4. Individual protein effects
    if df['significant'].any():
        sig_df = df[df['significant']].sort_values('log2FC')
        y_pos = np.arange(len(sig_df))
        colors = ['red' if d == 'V0' else 'blue' for d in sig_df['domain']]

        ax4.barh(y_pos, sig_df['log2FC'], color=colors, alpha=0.7)
        ax4.set_yticks(y_pos)
        ax4.set_yticklabels(sig_df['protein'], fontsize=8)
        ax4.set_xlabel('Log2 Fold Change')
        ax4.set_title('Significant V-ATPase Subunits')
        ax4.axvline(0, color='black', linestyle='-', linewidth=0.5)
        ax4.grid(True, alpha=0.3)
    else:
        ax4.text(0.5, 0.5, 'No significant\\nproteins found',
                ha='center', va='center', transform=ax4.transAxes, fontsize=12)

    plt.tight_layout()
    plt.show()

    # 📋 Print top results
    if df['significant'].any():
        print(f"\\n🔍 Top Differentially Expressed V-ATPase Subunits:")
        sig_sorted = df[df['significant']].sort_values('p_adjusted')
        for _, row in sig_sorted.head(5).iterrows():
            direction = "↑" if row['log2FC'] > 0 else "↓"
            print(f"  {row['protein']:12} ({row['domain']}) {direction} {abs(row['log2FC']):.3f} (FDR={row['p_adjusted']:.2e})")
    else:
        print(f"\\n⚠️ No significantly altered V-ATPase subunits found")
```

## Evaluate Claim

```python
# 🎯 CLAIM EVALUATION
print("\\n" + "="*50)
print("🎯 CLAIM EVALUATION")
print("="*50)
print("Claim: 'V-ATPase subunits show differential expression patterns'")
print()

if len(df) > 0:
    n_sig = df['significant'].sum()
    n_total = len(df)
    pct_sig = n_sig / n_total * 100

    # Check for domain-specific patterns
    v0_sig = df[(df['domain'] == 'V0') & df['significant']].shape[0]
    v1_sig = df[(df['domain'] == 'V1') & df['significant']].shape[0]

    # Check for bidirectional changes
    if n_sig > 0:
        sig_df = df[df['significant']]
        n_up = (sig_df['log2FC'] > 0).sum()
        n_down = (sig_df['log2FC'] < 0).sum()
        bidirectional = n_up > 0 and n_down > 0
    else:
        bidirectional = False

    print(f"📊 Analysis Results:")
    print(f"V-ATPase subunits tested: {n_total}/{len(all_vatpase)}")
    print(f"Significantly changed: {n_sig} ({pct_sig:.1f}%)")
    print(f"V0 domain affected: {v0_sig}")
    print(f"V1 domain affected: {v1_sig}")
    print(f"Bidirectional changes: {'Yes' if bidirectional else 'No'}")

    # Verdict logic
    if n_sig >= 4 and pct_sig > 25:
        verdict = "✅ STRONGLY SUPPORTED"
        explanation = f"Multiple V-ATPase subunits show clear differential patterns ({n_sig}/{n_total})"
    elif n_sig >= 3 or pct_sig > 15:
        verdict = "⚠️ SUPPORTED"
        explanation = f"V-ATPase shows differential expression patterns ({n_sig}/{n_total})"
    elif n_sig >= 2:
        verdict = "⚠️ PARTIALLY SUPPORTED"
        explanation = f"Some V-ATPase subunits affected ({n_sig}/{n_total})"
    elif n_sig == 1:
        verdict = "⚠️ WEAKLY SUPPORTED"
        explanation = "Only one V-ATPase subunit significantly changed"
    else:
        verdict = "❌ REFUTED"
        explanation = "No significant differential expression in V-ATPase subunits"

    print(f"\\n🎯 VERDICT: {verdict}")
    print(f"📝 {explanation}")

    # Biological interpretation
    print(f"\\n🧬 BIOLOGICAL IMPLICATIONS:")
    if n_sig > 0:
        print("• V-ATPase dysfunction detected")
        print("• Lysosomal acidification likely impaired")
        print("• Consequences:")
        print("  - Reduced autophagy efficiency")
        print("  - Impaired protein degradation")
        print("  - Cathepsin inactivation")
        print("  - Autophagosome-lysosome fusion defects")

        if bidirectional:
            print("• Mixed V-ATPase response suggests complex regulation")

        if v0_sig > 0:
            print(f"• V0 domain affected: proton translocation impaired")
        if v1_sig > 0:
            print(f"• V1 domain affected: ATP hydrolysis disrupted")

        print("\\n💊 Therapeutic relevance:")
        print("• V-ATPase could be therapeutic target")
        print("• Lysosomal acidification restoration strategies")
        print("• Links to SQSTM1/p62 accumulation")
    else:
        print("• V-ATPase function appears preserved")
        print("• Lysosomal acidification likely maintained")
        print("• Other proteostasis mechanisms may be primary targets")

else:
    verdict = "❌ UNSURE"
    explanation = "V-ATPase proteins not found in dataset"
    print(f"VERDICT: {verdict}")
    print(f"Reasoning: {explanation}")
```

## Save Results

```python
# 💾 Save results
if len(df) > 0:
    summary = {
        'claim': 'V-ATPase differential expression',
        'verdict': verdict,
        'proteins_tested': len(df),
        'significant_proteins': n_sig,
        'percent_significant': f"{pct_sig:.1f}%",
        'v0_affected': v0_sig,
        'v1_affected': v1_sig,
        'bidirectional': bidirectional
    }

    if IN_COLAB:
        df.to_csv('vatpase_results.csv', index=False)
        pd.DataFrame([summary]).to_csv('vatpase_summary.csv', index=False)
        files.download('vatpase_results.csv')
        files.download('vatpase_summary.csv')
        print("📁 Results downloaded!")
    else:
        df.to_csv('vatpase_analysis.csv', index=False)
        pd.DataFrame([summary]).to_csv('vatpase_summary.csv', index=False)
        print("📁 Results saved!")

print("\\n✅ V-ATPase analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified V-ATPase analysis:
- ✅ **Tests key V-ATPase subunits** (V0 membrane + V1 cytoplasmic domains)
- ✅ **Domain-specific analysis** to understand functional impacts
- ✅ **Clear visualizations** showing expression patterns
- ✅ **Biological interpretation** linking to lysosomal dysfunction
- ✅ **Therapeutic relevance** for proteostasis restoration

**Perfect for testing lysosomal acidification dysfunction!** 🔋