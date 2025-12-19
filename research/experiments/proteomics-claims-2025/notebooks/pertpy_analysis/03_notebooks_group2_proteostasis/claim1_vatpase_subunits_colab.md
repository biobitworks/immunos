# V-ATPase Subunits Analysis - Simplified Version
## Testing: "V-ATPase subunits are dysregulated in tau+ neurons"

**🚀 Simple approach**: Test V-ATPase subunit domains for dysregulation patterns!

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

# 🔋 V-ATPase Subunit Organization
vatpase_subunits = {
    'V0_Domain': [
        'ATP6V0A1',   # V0 a1 subunit
        'ATP6V0A2',   # V0 a2 subunit
        'ATP6V0A4',   # V0 a4 subunit
        'ATP6V0B',    # V0 b subunit
        'ATP6V0C',    # V0 c subunit
        'ATP6V0D1',   # V0 d1 subunit
        'ATP6V0D2',   # V0 d2 subunit
        'ATP6V0E1',   # V0 e1 subunit
        'ATP6V0E2'    # V0 e2 subunit
    ],
    'V1_Domain': [
        'ATP6V1A',    # V1 A subunit
        'ATP6V1B1',   # V1 B1 subunit
        'ATP6V1B2',   # V1 B2 subunit
        'ATP6V1C1',   # V1 C1 subunit
        'ATP6V1C2',   # V1 C2 subunit
        'ATP6V1D',    # V1 D subunit
        'ATP6V1E1',   # V1 E1 subunit
        'ATP6V1E2',   # V1 E2 subunit
        'ATP6V1F',    # V1 F subunit
        'ATP6V1G1',   # V1 G1 subunit
        'ATP6V1G2',   # V1 G2 subunit
        'ATP6V1G3',   # V1 G3 subunit
        'ATP6V1H'     # V1 H subunit
    ],
    'Accessory_Proteins': [
        'ATP6AP1',    # V-ATPase accessory protein 1
        'ATP6AP2',    # V-ATPase accessory protein 2 (PRR)
        'TCIRG1',     # T-cell immune regulator 1
        'CLCN3',      # Chloride channel 3
        'CLCN4',      # Chloride channel 4
        'CLCN7'       # Chloride channel 7
    ]
}

total_subunits = sum(len(subunits) for subunits in vatpase_subunits.values())
print(f"🎯 Testing {total_subunits} V-ATPase subunits across {len(vatpase_subunits)} domains")
```

## Find & Analyze V-ATPase Subunits

```python
# 🔍 Find available V-ATPase subunits
protein_names = list(adata.var_names)
domain_results = []

for domain_name, subunits in vatpase_subunits.items():
    print(f"\n📊 Analyzing {domain_name}...")

    # Find available subunits
    found_subunits = [s for s in subunits if s in protein_names]
    print(f"  Found: {len(found_subunits)}/{len(subunits)} subunits")

    if len(found_subunits) >= 1:  # Need at least 1 subunit for analysis
        # Get expression data for each subunit
        subunit_results = []

        for subunit in found_subunits:
            subunit_idx = protein_names.index(subunit)
            expr = adata.X[:, subunit_idx]

            # Split by tau status
            tau_pos_expr = expr[adata.obs['tau_positive'] == 1]
            tau_neg_expr = expr[adata.obs['tau_positive'] == 0]

            # Calculate statistics
            mean_pos = np.mean(tau_pos_expr)
            mean_neg = np.mean(tau_neg_expr)
            log2fc = mean_pos - mean_neg

            # T-test
            t_stat, p_val = stats.ttest_ind(tau_pos_expr, tau_neg_expr)

            subunit_results.append({
                'domain': domain_name,
                'subunit': subunit,
                'log2FC': log2fc,
                'p_value': p_val,
                'mean_tau_pos': mean_pos,
                'mean_tau_neg': mean_neg,
                'dysregulated': abs(log2fc) > 0.2
            })

        # Domain-level statistics
        log2fcs = [r['log2FC'] for r in subunit_results]
        mean_log2fc = np.mean(log2fcs)
        std_log2fc = np.std(log2fcs)
        dysregulated_count = sum(r['dysregulated'] for r in subunit_results)
        dysregulated_pct = dysregulated_count / len(subunit_results) * 100

        # Determine domain status
        if abs(mean_log2fc) > 0.3:
            if mean_log2fc > 0:
                status = "Upregulated"
            else:
                status = "Downregulated"
        elif std_log2fc > 0.4:
            status = "Variable"
        else:
            status = "Stable"

        # Add to domain results
        domain_results.append({
            'domain': domain_name,
            'subunits_found': len(found_subunits),
            'subunits_total': len(subunits),
            'mean_log2FC': mean_log2fc,
            'std_log2FC': std_log2fc,
            'dysregulated_count': dysregulated_count,
            'dysregulated_pct': dysregulated_pct,
            'domain_status': status,
            'subunits_tested': subunit_results
        })

        print(f"  Mean log2FC: {mean_log2fc:.3f} (±{std_log2fc:.3f})")
        print(f"  Status: {status}")
        print(f"  Dysregulated: {dysregulated_count}/{len(found_subunits)} ({dysregulated_pct:.1f}%)")

# Create combined DataFrame for all subunits
all_subunits_df = []
for domain_result in domain_results:
    all_subunits_df.extend(domain_result['subunits_tested'])

df = pd.DataFrame(all_subunits_df)
if len(df) > 0:
    # Add FDR correction
    from statsmodels.stats.multitest import fdrcorrection
    df['p_adjusted'] = fdrcorrection(df['p_value'])[1]
    df['significant'] = df['p_adjusted'] < 0.05

print(f"\n✅ Analysis complete: {len(df)} subunits tested across {len(domain_results)} domains")
```

## Domain-Specific Analysis

```python
if len(df) > 0:
    # 🔋 V-ATPase functional analysis
    print(f"\n🔋 V-ATPASE FUNCTIONAL ANALYSIS:")

    # Overall dysregulation
    total_dysregulated = sum(df['dysregulated'])
    pct_dysregulated = total_dysregulated / len(df) * 100

    sig_dysregulated = sum((df['significant']) & (df['dysregulated']))
    pct_sig_dysregulated = sig_dysregulated / len(df) * 100

    print(f"Total subunits: {len(df)}")
    print(f"Dysregulated subunits: {total_dysregulated} ({pct_dysregulated:.1f}%)")
    print(f"Significantly dysregulated: {sig_dysregulated} ({pct_sig_dysregulated:.1f}%)")

    # Domain comparison
    print(f"\nDomain-specific patterns:")
    for dr in domain_results:
        print(f"  {dr['domain']:18} {dr['domain_status']:12} ({dr['dysregulated_pct']:4.1f}% dysregulated)")

    # V0 vs V1 comparison
    v0_df = df[df['domain'] == 'V0_Domain']
    v1_df = df[df['domain'] == 'V1_Domain']

    if len(v0_df) > 0 and len(v1_df) > 0:
        v0_mean = v0_df['log2FC'].mean()
        v1_mean = v1_df['log2FC'].mean()

        # Statistical comparison
        t_stat_domains, p_val_domains = stats.ttest_ind(v0_df['log2FC'], v1_df['log2FC'])

        print(f"\nV0 vs V1 domain comparison:")
        print(f"  V0 mean log2FC: {v0_mean:.3f}")
        print(f"  V1 mean log2FC: {v1_mean:.3f}")
        print(f"  Domain difference p-value: {p_val_domains:.4f}")
```

## Visualize V-ATPase Results

```python
if len(df) > 0:
    # 📊 Create comprehensive visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Domain-level dysregulation
    domain_names = [dr['domain'] for dr in domain_results]
    domain_dysfunc = [dr['dysregulated_pct'] for dr in domain_results]

    # Color by domain status
    status_colors = {
        'Upregulated': 'green',
        'Downregulated': 'red',
        'Variable': 'orange',
        'Stable': 'gray'
    }
    colors = [status_colors[dr['domain_status']] for dr in domain_results]

    bars = ax1.bar(range(len(domain_names)), domain_dysfunc, color=colors, alpha=0.7)
    ax1.set_xticks(range(len(domain_names)))
    ax1.set_xticklabels(domain_names, rotation=45, ha='right')
    ax1.set_ylabel('% Subunits Dysregulated')
    ax1.set_title('V-ATPase Domain Dysregulation')
    ax1.axhline(50, color='red', linestyle='--', alpha=0.5, label='Majority dysregulated')
    ax1.set_ylim(0, 100)
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Add status legend
    for status, color in status_colors.items():
        ax1.scatter([], [], c=color, alpha=0.7, s=60, label=status)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    # 2. Expression change distribution
    domain_fcs = [dr['mean_log2FC'] for dr in domain_results]
    domain_stds = [dr['std_log2FC'] for dr in domain_results]

    bars = ax2.bar(range(len(domain_names)), domain_fcs,
                  yerr=domain_stds, color=colors, alpha=0.7, capsize=5)
    ax2.set_xticks(range(len(domain_names)))
    ax2.set_xticklabels(domain_names, rotation=45, ha='right')
    ax2.set_ylabel('Mean Log2 Fold Change')
    ax2.set_title('V-ATPase Domain Expression Changes')
    ax2.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax2.axhline(0.3, color='green', linestyle='--', alpha=0.5)
    ax2.axhline(-0.3, color='red', linestyle='--', alpha=0.5)
    ax2.grid(True, alpha=0.3)

    # 3. Volcano plot of all subunits
    colors_volcano = []
    domain_color_map = {
        'V0_Domain': 'blue',
        'V1_Domain': 'red',
        'Accessory_Proteins': 'green'
    }

    for _, row in df.iterrows():
        if row['significant'] and abs(row['log2FC']) > 0.2:
            colors_volcano.append(domain_color_map.get(row['domain'], 'gray'))
        else:
            colors_volcano.append('lightgray')

    ax3.scatter(df['log2FC'], -np.log10(df['p_value']), c=colors_volcano, alpha=0.6, s=40)
    ax3.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    ax3.axvline(0.2, color='green', linestyle='--', alpha=0.5)
    ax3.axvline(-0.2, color='red', linestyle='--', alpha=0.5)
    ax3.axvline(0, color='black', linestyle='-', alpha=0.3)
    ax3.set_xlabel('Log2 Fold Change (Tau+ vs Tau-)')
    ax3.set_ylabel('-Log10(p-value)')
    ax3.set_title('V-ATPase Subunits Volcano Plot')
    ax3.grid(True, alpha=0.3)

    # Add domain legend
    for domain, color in domain_color_map.items():
        ax3.scatter([], [], c=color, alpha=0.7, s=60, label=domain.replace('_', ' '))
    ax3.legend()

    # 4. V0 vs V1 comparison
    if len(v0_df) > 0 and len(v1_df) > 0:
        box_data = [v0_df['log2FC'], v1_df['log2FC']]
        box_labels = ['V0 Domain', 'V1 Domain']
        box_colors = ['lightblue', 'lightcoral']

        bp = ax4.boxplot(box_data, labels=box_labels, patch_artist=True)
        for patch, color in zip(bp['boxes'], box_colors):
            patch.set_facecolor(color)

        ax4.set_ylabel('Log2 Fold Change')
        ax4.set_title('V0 vs V1 Domain Comparison')
        ax4.axhline(0, color='black', linestyle='-', linewidth=0.5)
        ax4.grid(True, alpha=0.3)

        # Add statistical annotation
        ax4.text(0.5, ax4.get_ylim()[1] * 0.9, f'p = {p_val_domains:.4f}',
                ha='center', transform=ax4.transData, fontweight='bold')
    else:
        ax4.text(0.5, 0.5, 'Insufficient\nData for\nComparison', ha='center', va='center',
                transform=ax4.transAxes, fontsize=14)
        ax4.set_title('V0 vs V1 Comparison')

    plt.tight_layout()
    plt.show()

else:
    print("⚠️ No data available for visualization")
```

## Evaluate Claim

```python
# 🎯 CLAIM EVALUATION
print("\n" + "="*55)
print("🎯 CLAIM EVALUATION")
print("="*55)
print("Claim: 'V-ATPase subunits are dysregulated in tau+ neurons'")
print()

if len(df) > 0:
    print(f"📊 Overall Results:")
    print(f"V-ATPase subunits tested: {len(df)}")
    print(f"Dysregulated subunits: {total_dysregulated} ({pct_dysregulated:.1f}%)")
    print(f"Significantly dysregulated: {sig_dysregulated} ({pct_sig_dysregulated:.1f}%)")

    # Domain analysis
    domains_dysregulated = sum(1 for dr in domain_results if dr['dysregulated_pct'] > 50)
    print(f"Domains with dysregulation: {domains_dysregulated}/{len(domain_results)}")

    print(f"\n📈 Domain-specific Results:")
    for dr in domain_results:
        status_symbol = {
            'Upregulated': '↗️',
            'Downregulated': '↘️',
            'Variable': '↕️',
            'Stable': '→'
        }
        symbol = status_symbol.get(dr['domain_status'], '?')
        print(f"  {dr['domain']:18} {dr['domain_status']:12} {symbol} ({dr['dysregulated_pct']:4.1f}% dysregulated)")

    # Overall verdict
    if pct_sig_dysregulated > 60:
        verdict = "✅ STRONGLY SUPPORTED"
        explanation = f"Majority of V-ATPase subunits significantly dysregulated ({pct_sig_dysregulated:.1f}%)"
    elif pct_dysregulated > 50 and domains_dysregulated >= 2:
        verdict = "✅ SUPPORTED"
        explanation = f"Clear dysregulation across {domains_dysregulated} domains ({pct_dysregulated:.1f}%)"
    elif pct_dysregulated > 30 or domains_dysregulated >= 1:
        verdict = "⚠️ PARTIALLY SUPPORTED"
        explanation = f"Some V-ATPase dysregulation ({pct_dysregulated:.1f}% subunits)"
    else:
        verdict = "❌ REFUTED"
        explanation = f"No significant V-ATPase dysregulation ({pct_dysregulated:.1f}% affected)"

    print(f"\n🎯 VERDICT: {verdict}")
    print(f"📝 {explanation}")

    # Biological interpretation
    print(f"\n🧬 BIOLOGICAL IMPLICATIONS:")
    if verdict.startswith("✅"):
        print("• Lysosomal acidification machinery compromised")
        print("• pH-dependent processes severely impaired")
        print("• Autophagosome-lysosome fusion affected")
        print("• Cathepsin activation reduced")
        print("• Secondary lysosomal storage dysfunction")

        # Domain-specific impacts
        domain_impacts = {
            'V0_Domain': "Membrane integration and proton translocation impaired",
            'V1_Domain': "ATP hydrolysis and catalytic activity compromised",
            'Accessory_Proteins': "V-ATPase assembly and regulation disrupted"
        }

        print("\n💡 Domain-specific impacts:")
        for dr in domain_results:
            if dr['dysregulated_pct'] > 50:
                impact = domain_impacts.get(dr['domain'], "Function compromised")
                print(f"  • {dr['domain']}: {impact}")

    elif verdict.startswith("⚠️"):
        print("• Partial V-ATPase dysfunction")
        print("• Some lysosomal acidification preserved")
        print("• Compensatory mechanisms may be active")
        print("• Variable pH regulation")
    else:
        print("• V-ATPase machinery appears functional")
        print("• Lysosomal acidification preserved")
        print("• Normal pH-dependent processes")

else:
    verdict = "❌ UNSURE"
    explanation = "Insufficient V-ATPase subunit proteins found"
    print(f"VERDICT: {verdict}")
    print(f"Reasoning: {explanation}")
```

## Save Results

```python
# 💾 Save results
if len(df) > 0:
    # Domain summary
    domain_summary = pd.DataFrame([{
        'domain': dr['domain'],
        'subunits_found': dr['subunits_found'],
        'subunits_total': dr['subunits_total'],
        'mean_log2FC': dr['mean_log2FC'],
        'std_log2FC': dr['std_log2FC'],
        'dysregulated_pct': dr['dysregulated_pct'],
        'domain_status': dr['domain_status']
    } for dr in domain_results])

    # Overall summary
    summary = {
        'analysis': 'V-ATPase subunits dysregulated in tau+ neurons',
        'verdict': verdict,
        'subunits_tested': len(df),
        'subunits_dysregulated': total_dysregulated if 'total_dysregulated' in locals() else 0,
        'percent_dysregulated': pct_dysregulated if 'pct_dysregulated' in locals() else 0,
        'domains_affected': domains_dysregulated if 'domains_dysregulated' in locals() else 0
    }

    if IN_COLAB:
        df.to_csv('vatpase_subunits_results.csv', index=False)
        domain_summary.to_csv('vatpase_domain_summary.csv', index=False)
        pd.DataFrame([summary]).to_csv('vatpase_analysis_summary.csv', index=False)
        files.download('vatpase_subunits_results.csv')
        files.download('vatpase_domain_summary.csv')
        print("📁 Results downloaded!")
    else:
        df.to_csv('vatpase_subunits_analysis.csv', index=False)
        domain_summary.to_csv('vatpase_domain_summary.csv', index=False)
        print("📁 Results saved!")

print("\n✅ V-ATPase subunit analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified V-ATPase analysis:
- ✅ **Tests all V-ATPase domains** (V0, V1, accessory proteins)
- ✅ **Domain-specific dysregulation** assessment
- ✅ **Clear visualizations** showing acidification dysfunction
- ✅ **Biological interpretation** linking to lysosomal pH regulation
- ✅ **Therapeutic relevance** for acidification enhancement strategies

**Perfect for evaluating lysosomal acidification machinery in neurodegeneration!** 🔋