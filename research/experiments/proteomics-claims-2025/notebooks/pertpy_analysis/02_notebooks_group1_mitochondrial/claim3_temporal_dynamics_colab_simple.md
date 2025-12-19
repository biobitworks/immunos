# Temporal Dynamics Analysis - Simplified Version
## Testing: "Progressive mitochondrial dysfunction over time"

**🚀 Simple approach**: Analyze how key proteins change with disease progression!

---

## Setup & Data Loading

```python
# 🔧 Setup
import os
IN_COLAB = 'google.colab' in str(get_ipython())

if IN_COLAB:
    !pip install -q pertpy scanpy pandas numpy matplotlib seaborn scipy scikit-learn
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
from sklearn.linear_model import LinearRegression
import warnings
warnings.filterwarnings('ignore')

print("✅ Setup complete!")
```

## Load Data & Define Key Proteins

```python
# 📊 Load data
adata = sc.read_h5ad(data_file)
print(f"Data loaded: {adata.shape}")

# 🏷️ Setup variables
if 'TauStatus' in adata.obs.columns:
    adata.obs['tau_status'] = adata.obs['TauStatus']
if 'Pseudotime' in adata.obs.columns and 'pseudotime' not in adata.obs.columns:
    adata.obs['pseudotime'] = adata.obs['Pseudotime']

adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)

# Check pseudotime availability
has_pseudotime = 'pseudotime' in adata.obs.columns and not adata.obs['pseudotime'].isna().all()
print(f"Pseudotime available: {has_pseudotime}")

if not has_pseudotime:
    print("⚠️ Creating mock pseudotime for demonstration")
    np.random.seed(42)
    adata.obs['pseudotime'] = np.random.uniform(0, 1, adata.n_obs)

# 🧬 Key protein sets (simplified)
protein_sets = {
    'Mitochondrial Complex I': [
        'NDUFA1', 'NDUFA2', 'NDUFA4', 'NDUFA9', 'NDUFB1', 'NDUFB3',
        'NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFV1', 'NDUFV2'
    ],
    'Mitochondrial Complex V': [
        'ATP5A1', 'ATP5B', 'ATP5C1', 'ATP5D', 'ATP5F1', 'ATP5G1', 'ATP5H'
    ],
    'Autophagy': [
        'ULK1', 'BECN1', 'ATG5', 'ATG7', 'MAP1LC3A', 'MAP1LC3B',
        'SQSTM1', 'NBR1', 'OPTN'
    ],
    'Proteasome': [
        'PSMA1', 'PSMA2', 'PSMB1', 'PSMB5', 'PSMC1', 'PSMC4', 'PSMD1'
    ],
    'Heat Shock': [
        'HSPA1A', 'HSPA5', 'HSPA8', 'HSPA9', 'HSPB1', 'HSP90AA1', 'HSPD1'
    ]
}

print(f"🎯 Analyzing {sum(len(v) for v in protein_sets.values())} proteins across {len(protein_sets)} pathways")
```

## Find Proteins & Calculate Correlations

```python
# 🔍 Find available proteins and calculate temporal correlations
results = []
protein_names = list(adata.var_names)

for pathway, proteins in protein_sets.items():
    print(f"\\n📊 Analyzing {pathway}...")

    pathway_correlations = []
    found_proteins = []

    for protein in proteins:
        # Find protein (exact or partial match)
        if protein in protein_names:
            found_proteins.append(protein)
            protein_idx = protein_names.index(protein)
        else:
            # Try partial match
            matches = [p for p in protein_names if protein in p.upper()]
            if matches:
                found_proteins.append(matches[0])
                protein_idx = protein_names.index(matches[0])
            else:
                continue

        # Get expression and calculate correlation with pseudotime
        expr = adata.X[:, protein_idx]
        valid_mask = ~(np.isnan(expr) | np.isnan(adata.obs['pseudotime']))

        if valid_mask.sum() > 10:
            corr, p_val = stats.spearmanr(adata.obs['pseudotime'][valid_mask], expr[valid_mask])
            pathway_correlations.append(corr)

            # Store individual protein results
            results.append({
                'pathway': pathway,
                'protein': found_proteins[-1],
                'correlation': corr,
                'p_value': p_val,
                'direction': 'decreasing' if corr < 0 else 'increasing'
            })

    # Pathway summary
    if pathway_correlations:
        mean_corr = np.mean(pathway_correlations)
        decreasing_pct = sum(1 for c in pathway_correlations if c < 0) / len(pathway_correlations) * 100
        print(f"  Found: {len(found_proteins)}/{len(proteins)} proteins")
        print(f"  Mean correlation: {mean_corr:.3f}")
        print(f"  Decreasing over time: {decreasing_pct:.1f}%")

# Convert to DataFrame
df = pd.DataFrame(results)
if len(df) > 0:
    # Add FDR correction
    from statsmodels.stats.multitest import fdrcorrection
    df['p_adjusted'] = fdrcorrection(df['p_value'])[1]
    df['significant'] = df['p_adjusted'] < 0.05

print(f"\\n✅ Analysis complete: {len(df)} proteins analyzed")
```

## Visualize Temporal Patterns

```python
if len(df) > 0:
    # 📊 Create summary visualization
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Correlation by pathway
    pathway_stats = df.groupby('pathway').agg({
        'correlation': 'mean',
        'significant': 'sum',
        'protein': 'count'
    }).rename(columns={'protein': 'total'})

    ax1 = axes[0, 0]
    bars = ax1.bar(range(len(pathway_stats)), pathway_stats['correlation'],
                   color=['red' if x < 0 else 'blue' for x in pathway_stats['correlation']])
    ax1.set_xticks(range(len(pathway_stats)))
    ax1.set_xticklabels(pathway_stats.index, rotation=45, ha='right')
    ax1.set_ylabel('Mean Correlation with Pseudotime')
    ax1.set_title('Pathway Temporal Patterns')
    ax1.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax1.grid(True, alpha=0.3)

    # 2. Significant proteins per pathway
    ax2 = axes[0, 1]
    ax2.bar(range(len(pathway_stats)), pathway_stats['significant'],
            color='green', alpha=0.7)
    ax2.set_xticks(range(len(pathway_stats)))
    ax2.set_xticklabels(pathway_stats.index, rotation=45, ha='right')
    ax2.set_ylabel('Significant Proteins')
    ax2.set_title('Temporal Significance by Pathway')
    ax2.grid(True, alpha=0.3)

    # 3. Overall correlation distribution
    ax3 = axes[1, 0]
    ax3.hist(df['correlation'], bins=20, alpha=0.7, color='skyblue', edgecolor='black')
    ax3.axvline(0, color='red', linestyle='--', linewidth=2)
    ax3.set_xlabel('Correlation with Pseudotime')
    ax3.set_ylabel('Number of Proteins')
    ax3.set_title('Distribution of Temporal Correlations')
    ax3.grid(True, alpha=0.3)

    # 4. Example protein trajectory
    ax4 = axes[1, 1]
    # Pick a significant mitochondrial protein
    mito_proteins = df[df['pathway'].str.contains('Mitochondrial')]
    if len(mito_proteins) > 0:
        example_protein = mito_proteins.iloc[0]['protein']
        if example_protein in protein_names:
            protein_idx = protein_names.index(example_protein)
            expr = adata.X[:, protein_idx]

            # Color by tau status
            colors = ['red' if x == 1 else 'blue' for x in adata.obs['tau_positive']]
            ax4.scatter(adata.obs['pseudotime'], expr, c=colors, alpha=0.6, s=20)

            # Add trend line
            valid_mask = ~(np.isnan(expr) | np.isnan(adata.obs['pseudotime']))
            if valid_mask.sum() > 10:
                z = np.polyfit(adata.obs['pseudotime'][valid_mask], expr[valid_mask], 1)
                p = np.poly1d(z)
                x_line = np.linspace(adata.obs['pseudotime'].min(), adata.obs['pseudotime'].max(), 100)
                ax4.plot(x_line, p(x_line), 'k--', linewidth=2)

            ax4.set_xlabel('Pseudotime')
            ax4.set_ylabel('Expression')
            ax4.set_title(f'Example: {example_protein}')
            ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
```

## Analyze Mitochondrial Dysfunction

```python
# 🔬 Focus on mitochondrial proteins
print("\\n" + "="*50)
print("🔬 MITOCHONDRIAL DYSFUNCTION ANALYSIS")
print("="*50)

mito_df = df[df['pathway'].str.contains('Mitochondrial')]

if len(mito_df) > 0:
    n_mito = len(mito_df)
    n_decreasing = sum(mito_df['correlation'] < 0)
    n_significant = sum(mito_df['significant'])
    mean_corr = mito_df['correlation'].mean()

    print(f"Mitochondrial proteins analyzed: {n_mito}")
    print(f"Decreasing over time: {n_decreasing}/{n_mito} ({n_decreasing/n_mito*100:.1f}%)")
    print(f"Significantly changing: {n_significant}/{n_mito} ({n_significant/n_mito*100:.1f}%)")
    print(f"Mean temporal correlation: {mean_corr:.3f}")

    # Complex-specific breakdown
    print(f"\\n📊 By Complex:")
    for pathway in mito_df['pathway'].unique():
        complex_df = mito_df[mito_df['pathway'] == pathway]
        complex_decreasing = sum(complex_df['correlation'] < 0)
        complex_mean = complex_df['correlation'].mean()
        print(f"  {pathway}: {complex_decreasing}/{len(complex_df)} decreasing, r̄={complex_mean:.3f}")

    # Show most changed proteins
    if n_significant > 0:
        sig_mito = mito_df[mito_df['significant']].sort_values('correlation')
        print(f"\\n🔍 Most changed mitochondrial proteins:")
        for _, row in sig_mito.head(3).iterrows():
            direction = "↓" if row['correlation'] < 0 else "↑"
            print(f"  {row['protein']:10} {direction} r={row['correlation']:.3f} (p={row['p_adjusted']:.2e})")

else:
    print("⚠️ No mitochondrial proteins found for analysis")
    n_decreasing = 0
    n_mito = 0
```

## Evaluate Claim

```python
# 🎯 CLAIM EVALUATION
print("\\n" + "="*50)
print("🎯 CLAIM EVALUATION")
print("="*50)
print("Claim: 'Temporal dynamics reveal progressive mitochondrial dysfunction'")
print()

if len(mito_df) > 0:
    pct_decreasing = n_decreasing / n_mito * 100
    pct_significant = n_significant / n_mito * 100

    # Decision logic
    if pct_decreasing > 60 and pct_significant > 30:
        verdict = "✅ STRONGLY SUPPORTED"
        explanation = f"{pct_decreasing:.1f}% of mitochondrial proteins decline over time"
    elif pct_decreasing > 50 or pct_significant > 20:
        verdict = "⚠️ PARTIALLY SUPPORTED"
        explanation = f"Moderate evidence: {n_decreasing}/{n_mito} proteins declining"
    else:
        verdict = "❌ REFUTED"
        explanation = "No clear progressive mitochondrial dysfunction pattern"

    print(f"VERDICT: {verdict}")
    print(f"Reasoning: {explanation}")

    # Additional evidence
    print(f"\\n📊 Evidence Summary:")
    print(f"• Mitochondrial proteins: {n_mito} analyzed")
    print(f"• Progressive decline: {pct_decreasing:.1f}%")
    print(f"• Statistically significant: {pct_significant:.1f}%")
    print(f"• Mean correlation: {mean_corr:.3f}")

    # Biological interpretation
    print(f"\\n🧬 Biological Meaning:")
    if verdict.startswith("✅"):
        print("• Clear progressive decline in mitochondrial function")
        print("• Energy production decreases with disease progression")
        print("• Supports bioenergetic failure hypothesis")
        print("• May trigger compensatory autophagy upregulation")
    elif verdict.startswith("⚠️"):
        print("• Some mitochondrial dysfunction detected")
        print("• Selective complex vulnerability")
        print("• Mixed evidence for progressive decline")
    else:
        print("• Mitochondrial function appears stable")
        print("• No clear temporal decline pattern")
        print("• Other mechanisms may be primary drivers")

else:
    verdict = "❌ UNSURE"
    explanation = "Insufficient mitochondrial proteins for analysis"
    print(f"VERDICT: {verdict}")
    print(f"Reasoning: {explanation}")
```

## Save Results

```python
# 💾 Save results
if len(df) > 0:
    summary = {
        'claim': 'Progressive mitochondrial dysfunction',
        'verdict': verdict,
        'proteins_analyzed': len(df),
        'mito_proteins': len(mito_df) if len(mito_df) > 0 else 0,
        'mito_decreasing': n_decreasing if 'n_decreasing' in locals() else 0,
        'mean_correlation': mean_corr if 'mean_corr' in locals() else None
    }

    # Save detailed results
    if IN_COLAB:
        df.to_csv('temporal_analysis.csv', index=False)
        pd.DataFrame([summary]).to_csv('temporal_summary.csv', index=False)
        files.download('temporal_analysis.csv')
        files.download('temporal_summary.csv')
        print("📁 Results downloaded!")
    else:
        df.to_csv('temporal_dynamics_results.csv', index=False)
        pd.DataFrame([summary]).to_csv('temporal_summary.csv', index=False)
        print("📁 Results saved!")

print("\\n✅ Temporal analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified temporal analysis:
- ✅ **Tests key protein pathways** for temporal changes
- ✅ **Focuses on mitochondrial dysfunction** specifically
- ✅ **Uses robust correlation analysis** with pseudotime
- ✅ **Clear visualization** of temporal patterns
- ✅ **Objective claim evaluation** with biological context

**Perfect for studying disease progression patterns!** 📈