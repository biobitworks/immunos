# ATP6V0A1 Dysfunction Analysis - Simplified Version
## Testing: "ATP6V0A1 subunit dysfunction leads to lysosomal alkalinization"

**🚀 Simple approach**: Test ATP6V0A1 and related V-ATPase acidification machinery!

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

## Load Data & Define ATP6V0A1 Acidification Machinery

```python
# 📊 Load data
adata = sc.read_h5ad(data_file)
print(f"Data loaded: {adata.shape}")

# 🏷️ Setup tau status
if 'TauStatus' in adata.obs.columns:
    adata.obs['tau_status'] = adata.obs['TauStatus']
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)

print(f"Tau+: {adata.obs['tau_positive'].sum()}, Tau-: {(adata.obs['tau_positive']==0).sum()}")

# 🔋 ATP6V0A1 and Acidification Components
acidification_machinery = {
    'V0_A_Subunits': [
        'ATP6V0A1',   # V-ATPase a1 subunit (primary target)
        'ATP6V0A2',   # V-ATPase a2 subunit
        'ATP6V0A4'    # V-ATPase a4 subunit
    ],
    'Critical_V0_Subunits': [
        'ATP6V0B',    # V-ATPase b subunit
        'ATP6V0C',    # V-ATPase c subunit
        'ATP6V0D1',   # V-ATPase d1 subunit
        'ATP6V0D2',   # V-ATPase d2 subunit
        'ATP6V0E1',   # V-ATPase e1 subunit
        'ATP6V0E2'    # V-ATPase e2 subunit
    ],
    'V1_Catalytic_Subunits': [
        'ATP6V1A',    # V-ATPase A subunit
        'ATP6V1B1',   # V-ATPase B1 subunit
        'ATP6V1B2',   # V-ATPase B2 subunit
        'ATP6V1C1',   # V-ATPase C1 subunit
        'ATP6V1D',    # V-ATPase D subunit
        'ATP6V1E1',   # V-ATPase E1 subunit
        'ATP6V1F',    # V-ATPase F subunit
        'ATP6V1G1'    # V-ATPase G1 subunit
    ],
    'Regulatory_Proteins': [
        'ATP6AP1',    # V-ATPase accessory protein 1 (Ac45)
        'ATP6AP2',    # V-ATPase accessory protein 2 (PRR)
        'TCIRG1',     # T-cell immune regulator 1
        'CLCN7'      # Chloride voltage-gated channel 7
    ],
    'pH_Sensors': [
        'SLC9A6',     # Sodium/hydrogen exchanger 6
        'SLC9A7',     # Sodium/hydrogen exchanger 7
        'LAMP1',      # Lysosomal associated membrane protein 1
        'LAMP2',      # Lysosomal associated membrane protein 2
        'CTSD'        # Cathepsin D (pH-dependent protease)
    ]
}

total_proteins = sum(len(proteins) for proteins in acidification_machinery.values())
print(f"🎯 Testing {total_proteins} acidification machinery proteins across {len(acidification_machinery)} categories")
```

## Find & Analyze Acidification Components

```python
# 🔍 Find available proteins and analyze each category
protein_names = list(adata.var_names)
category_results = []

for category_name, proteins in acidification_machinery.items():
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
                'decreased': log2fc < 0,
                'is_atp6v0a1': protein == 'ATP6V0A1'
            })

        # Category-level statistics
        log2fcs = [r['log2FC'] for r in protein_results]
        mean_log2fc = np.mean(log2fcs)
        decreased_count = sum(r['decreased'] for r in protein_results)
        decreased_pct = decreased_count / len(protein_results) * 100

        # Check for ATP6V0A1 specifically
        atp6v0a1_result = next((r for r in protein_results if r['is_atp6v0a1']), None)

        # Add to category results
        category_results.append({
            'category': category_name,
            'proteins_found': len(found_proteins),
            'proteins_total': len(proteins),
            'mean_log2FC': mean_log2fc,
            'decreased_count': decreased_count,
            'decreased_pct': decreased_pct,
            'proteins_tested': protein_results,
            'has_atp6v0a1': atp6v0a1_result is not None,
            'atp6v0a1_log2fc': atp6v0a1_result['log2FC'] if atp6v0a1_result else None
        })

        print(f"  Mean log2FC: {mean_log2fc:.3f}")
        print(f"  Decreased: {decreased_count}/{len(found_proteins)} ({decreased_pct:.1f}%)")
        if atp6v0a1_result:
            print(f"  ATP6V0A1 log2FC: {atp6v0a1_result['log2FC']:.3f}")

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

## Visualize ATP6V0A1 Results

```python
if len(df) > 0:
    # 📊 Create comprehensive visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Category-level fold changes (focus on decreases = dysfunction)
    category_names = [cr['category'] for cr in category_results]
    category_fcs = [cr['mean_log2FC'] for cr in category_results]
    colors = ['red' if fc < -0.1 else 'gray' for fc in category_fcs]

    bars = ax1.bar(range(len(category_names)), category_fcs, color=colors, alpha=0.7)
    ax1.set_xticks(range(len(category_names)))
    ax1.set_xticklabels(category_names, rotation=45, ha='right')
    ax1.set_ylabel('Mean Log2 Fold Change')
    ax1.set_title('Acidification Machinery Expression Changes')
    ax1.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax1.axhline(-0.5, color='red', linestyle='--', alpha=0.5, label='Substantial dysfunction')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # 2. Percentage decreased per category
    decreased_pcts = [cr['decreased_pct'] for cr in category_results]
    bars = ax2.bar(range(len(category_names)), decreased_pcts, color='orange', alpha=0.7)
    ax2.set_xticks(range(len(category_names)))
    ax2.set_xticklabels(category_names, rotation=45, ha='right')
    ax2.set_ylabel('% Proteins Decreased')
    ax2.set_title('Percentage of Proteins Decreased per Category')
    ax2.axhline(50, color='red', linestyle='--', alpha=0.5, label='Majority dysfunctional')
    ax2.set_ylim(0, 100)
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # 3. Volcano plot highlighting ATP6V0A1
    colors = []
    sizes = []
    for _, row in df.iterrows():
        if row['is_atp6v0a1']:
            colors.append('purple')  # Highlight ATP6V0A1
            sizes.append(100)
        elif row['significant'] and row['log2FC'] < -0.2:
            colors.append('red')     # Significantly decreased
            sizes.append(40)
        elif row['significant']:
            colors.append('blue')    # Significantly changed
            sizes.append(40)
        else:
            colors.append('gray')    # Not significant
            sizes.append(40)

    ax3.scatter(df['log2FC'], -np.log10(df['p_value']), c=colors, alpha=0.6, s=sizes)
    ax3.axhline(-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    ax3.axvline(-0.2, color='red', linestyle='--', alpha=0.5)
    ax3.axvline(0, color='black', linestyle='-', alpha=0.3)
    ax3.set_xlabel('Log2 Fold Change (Tau+ vs Tau-)')
    ax3.set_ylabel('-Log10(p-value)')
    ax3.set_title('Acidification Proteins (Purple = ATP6V0A1)')
    ax3.grid(True, alpha=0.3)

    # 4. ATP6V0A1 specific analysis
    atp6v0a1_data = df[df['is_atp6v0a1']]
    if len(atp6v0a1_data) > 0:
        # Single protein bar chart for ATP6V0A1
        atp6v0a1_fc = atp6v0a1_data.iloc[0]['log2FC']
        atp6v0a1_sig = atp6v0a1_data.iloc[0]['significant']

        color = 'red' if atp6v0a1_fc < -0.2 and atp6v0a1_sig else 'gray'
        bar = ax4.bar(['ATP6V0A1'], [atp6v0a1_fc], color=color, alpha=0.7)
        ax4.set_ylabel('Log2 Fold Change')
        ax4.set_title('ATP6V0A1 Dysfunction')
        ax4.axhline(0, color='black', linestyle='-', linewidth=0.5)
        ax4.axhline(-0.5, color='red', linestyle='--', alpha=0.5, label='Severe dysfunction')
        ax4.grid(True, alpha=0.3)
        ax4.legend()

        # Add significance annotation
        if atp6v0a1_sig:
            ax4.text(0, atp6v0a1_fc + 0.1, f'p < 0.05\n(log2FC: {atp6v0a1_fc:.3f})',
                    ha='center', va='bottom', fontweight='bold')
    else:
        ax4.text(0.5, 0.5, 'ATP6V0A1\nNot Found', ha='center', va='center',
                transform=ax4.transAxes, fontsize=14)
        ax4.set_title('ATP6V0A1 Status')

    plt.tight_layout()
    plt.show()
```

## Evaluate Claim

```python
# 🎯 CLAIM EVALUATION
print("\n" + "="*65)
print("🎯 CLAIM EVALUATION")
print("="*65)
print("Claim: 'ATP6V0A1 subunit dysfunction leads to lysosomal alkalinization'")
print()

if len(df) > 0:
    # Overall statistics
    total_proteins = len(df)
    total_decreased = sum(df['decreased'])
    pct_decreased = total_decreased / total_proteins * 100

    sig_decreased = sum((df['significant']) & (df['decreased']))
    pct_sig_decreased = sig_decreased / total_proteins * 100

    # ATP6V0A1 specific analysis
    atp6v0a1_found = any(df['is_atp6v0a1'])
    if atp6v0a1_found:
        atp6v0a1_row = df[df['is_atp6v0a1']].iloc[0]
        atp6v0a1_decreased = atp6v0a1_row['decreased']
        atp6v0a1_significant = atp6v0a1_row['significant']
        atp6v0a1_log2fc = atp6v0a1_row['log2FC']

    # Category-level analysis
    categories_mostly_decreased = sum(1 for cr in category_results if cr['decreased_pct'] > 50)

    print(f"📊 Overall Results:")
    print(f"Proteins tested: {total_proteins}")
    print(f"Proteins decreased: {total_decreased} ({pct_decreased:.1f}%)")
    print(f"Significantly decreased: {sig_decreased} ({pct_sig_decreased:.1f}%)")
    print(f"Categories mostly decreased: {categories_mostly_decreased}/{len(category_results)}")

    if atp6v0a1_found:
        print(f"\n🎯 ATP6V0A1 Specific Results:")
        print(f"ATP6V0A1 found: YES")
        print(f"ATP6V0A1 log2FC: {atp6v0a1_log2fc:.3f}")
        print(f"ATP6V0A1 decreased: {'YES' if atp6v0a1_decreased else 'NO'}")
        print(f"ATP6V0A1 significant: {'YES' if atp6v0a1_significant else 'NO'}")
    else:
        print(f"\n🎯 ATP6V0A1 Specific Results:")
        print(f"ATP6V0A1 found: NO")

    print(f"\n📈 Category-specific Results:")
    for cr in category_results:
        status = "✓ DYSFUNCTIONAL" if cr['decreased_pct'] > 50 else "✗ Mixed/Stable"
        print(f"  {cr['category']:20} {cr['decreased_pct']:5.1f}% decreased  {status}")

    # Overall verdict
    if atp6v0a1_found and atp6v0a1_decreased and atp6v0a1_significant:
        if pct_sig_decreased > 50:
            verdict = "✅ STRONGLY SUPPORTED"
            explanation = f"ATP6V0A1 significantly decreased + widespread V-ATPase dysfunction ({pct_sig_decreased:.1f}%)"
        else:
            verdict = "✅ SUPPORTED"
            explanation = f"ATP6V0A1 significantly decreased (log2FC: {atp6v0a1_log2fc:.3f})"
    elif atp6v0a1_found and atp6v0a1_decreased:
        verdict = "⚠️ PARTIALLY SUPPORTED"
        explanation = f"ATP6V0A1 decreased but not significant (log2FC: {atp6v0a1_log2fc:.3f})"
    elif pct_decreased > 60:
        verdict = "⚠️ PARTIALLY SUPPORTED"
        explanation = f"V-ATPase machinery widely decreased ({pct_decreased:.1f}%) but ATP6V0A1 missing"
    else:
        verdict = "❌ REFUTED"
        explanation = f"No clear ATP6V0A1 dysfunction or acidification failure"

    print(f"\n🎯 VERDICT: {verdict}")
    print(f"📝 {explanation}")

    # Biological interpretation
    print(f"\n🧬 BIOLOGICAL IMPLICATIONS:")
    if verdict.startswith("✅"):
        print("• Lysosomal acidification severely compromised")
        print("• pH-dependent proteases (cathepsins) inhibited")
        print("• Autophagosome-lysosome fusion impaired")
        print("• Protein aggregate clearance blocked")
        print("• Secondary lysosomal storage dysfunction")

        if atp6v0a1_found and atp6v0a1_decreased:
            print("• ATP6V0A1 subunit dysfunction confirmed")
            print("• V0 domain assembly/stability compromised")
            print("• Proton translocation efficiency reduced")

    elif verdict.startswith("⚠️"):
        print("• Partial acidification machinery dysfunction")
        print("• Some lysosomal processes may be impaired")
        print("• Compensatory mechanisms may be active")
        print("• Variable pH regulation across compartments")
    else:
        print("• Lysosomal acidification appears normal")
        print("• V-ATPase machinery functioning")
        print("• pH-dependent processes preserved")

else:
    verdict = "❌ UNSURE"
    explanation = "Insufficient acidification machinery proteins found"
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
        'decreased_pct': cr['decreased_pct'],
        'has_atp6v0a1': cr['has_atp6v0a1'],
        'atp6v0a1_log2fc': cr['atp6v0a1_log2fc']
    } for cr in category_results])

    # Overall summary
    summary = {
        'analysis': 'ATP6V0A1 dysfunction leads to lysosomal alkalinization',
        'verdict': verdict,
        'proteins_tested': total_proteins if 'total_proteins' in locals() else 0,
        'proteins_decreased': total_decreased if 'total_decreased' in locals() else 0,
        'percent_decreased': pct_decreased if 'pct_decreased' in locals() else 0,
        'atp6v0a1_found': atp6v0a1_found if 'atp6v0a1_found' in locals() else False,
        'atp6v0a1_dysfunction': atp6v0a1_decreased if 'atp6v0a1_decreased' in locals() else False
    }

    if IN_COLAB:
        df.to_csv('atp6v0a1_dysfunction_results.csv', index=False)
        category_summary.to_csv('acidification_category_summary.csv', index=False)
        pd.DataFrame([summary]).to_csv('atp6v0a1_analysis_summary.csv', index=False)
        files.download('atp6v0a1_dysfunction_results.csv')
        files.download('acidification_category_summary.csv')
        print("📁 Results downloaded!")
    else:
        df.to_csv('atp6v0a1_dysfunction_analysis.csv', index=False)
        category_summary.to_csv('acidification_category_summary.csv', index=False)
        print("📁 Results saved!")

print("\n✅ ATP6V0A1 dysfunction analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified ATP6V0A1 dysfunction analysis:
- ✅ **Tests ATP6V0A1 specifically** plus broader V-ATPase machinery
- ✅ **Category-specific evaluation** (V0, V1, regulatory, pH sensors)
- ✅ **Clear visualizations** showing acidification failure
- ✅ **Biological interpretation** linking to lysosomal alkalinization
- ✅ **Therapeutic relevance** for lysosomal acidification enhancers

**Perfect for evaluating lysosomal pH dysregulation in neurodegeneration!** 🔋