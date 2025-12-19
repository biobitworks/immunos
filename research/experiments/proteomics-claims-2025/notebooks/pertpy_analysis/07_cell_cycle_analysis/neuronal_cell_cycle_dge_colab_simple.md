# Cell Cycle Analysis - Simplified Version
## Testing: "Are genes differentially expressed based on neuronal cell cycle?"

**🚀 Simple approach**: Test if cell cycle re-entry affects gene expression in neurons!

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
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings('ignore')

print("✅ Setup complete!")
```

## Load Data & Define Cell Cycle Genes

```python
# 📊 Load data
adata = sc.read_h5ad(data_file)
print(f"Data loaded: {adata.shape}")

# 🏷️ Setup tau status
if 'TauStatus' in adata.obs.columns:
    adata.obs['tau_status'] = adata.obs['TauStatus']
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)

print(f"Tau+: {adata.obs['tau_positive'].sum()}, Tau-: {(adata.obs['tau_positive']==0).sum()}")

# 🧬 Key cell cycle genes (simplified sets)
cell_cycle_genes = {
    'G1_S_transition': [
        'CCND1', 'CCND2', 'CCNE1', 'CCNE2',  # Cyclins
        'CDK4', 'CDK6', 'CDK2',              # CDKs
        'E2F1', 'E2F2', 'MYC',               # Transcription factors
        'RB1', 'RBL1'                        # Tumor suppressors
    ],
    'S_phase': [
        'PCNA',                               # Replication marker
        'MCM2', 'MCM3', 'MCM4', 'MCM6',     # Replication complex
        'TYMS', 'RRM2',                      # DNA synthesis
        'POLA1', 'POLE'                      # DNA polymerases
    ],
    'G2_M_transition': [
        'CCNB1', 'CCNB2', 'CCNA2',          # Mitotic cyclins
        'CDK1',                              # Mitotic CDK
        'AURKA', 'AURKB', 'PLK1',           # Mitotic kinases
        'TOP2A'                              # DNA topology
    ],
    'Cell_cycle_exit': [
        'CDKN1A', 'CDKN1B', 'CDKN2A',       # CDK inhibitors (p21, p27, p16)
        'BTG1', 'BTG2',                      # Growth arrest
        'GAS1'                               # Growth arrest specific
    ],
    'DNA_damage': [
        'ATM', 'ATR',                        # Damage sensors
        'CHEK1', 'CHEK2',                    # Checkpoint kinases
        'TP53',                              # Tumor suppressor
        'GADD45A'                            # DNA damage response
    ],
    'Apoptosis': [
        'BAX', 'BCL2',                       # BCL2 family
        'CASP3', 'CASP9',                    # Caspases
        'TP53BP1'                            # p53 pathway
    ]
}

total_genes = sum(len(genes) for genes in cell_cycle_genes.values())
print(f"🎯 Testing {total_genes} cell cycle genes across {len(cell_cycle_genes)} categories")
```

## Find Cell Cycle Genes & Calculate Scores

```python
# 🔍 Find available cell cycle genes
protein_names = list(adata.var_names)
found_genes = {}
all_found = []

for category, genes in cell_cycle_genes.items():
    found = [g for g in genes if g in protein_names]
    # Try partial matches for missing genes
    for missing_gene in set(genes) - set(found):
        matches = [p for p in protein_names if missing_gene in p.upper()]
        if matches:
            found.append(matches[0])

    found_genes[category] = found
    all_found.extend(found)
    print(f"{category:20} {len(found):2}/{len(genes):2} genes found")

print(f"\\nTotal found: {len(set(all_found))}/{total_genes} unique genes")

# 📊 Calculate cell cycle phase scores
phase_scores = pd.DataFrame(index=adata.obs.index)

for category, genes in found_genes.items():
    if genes:
        # Get expression matrix for this gene set
        gene_indices = [protein_names.index(g) for g in genes if g in protein_names]
        if gene_indices:
            expr_matrix = adata.X[:, gene_indices]
            # Calculate mean expression as phase score
            phase_scores[category] = np.mean(expr_matrix, axis=1)
    else:
        phase_scores[category] = 0

# Add to adata
for col in phase_scores.columns:
    adata.obs[f'score_{col}'] = phase_scores[col]

print("✅ Cell cycle scores calculated!")
```

## Analyze Cell Cycle Patterns

```python
# 🧮 Compare cell cycle scores between tau groups
print("\\n📊 CELL CYCLE ANALYSIS")
print("="*50)

comparison_results = []

for category in cell_cycle_genes.keys():
    score_col = f'score_{category}'

    if score_col in adata.obs.columns:
        # Get scores for each group
        tau_pos_scores = adata.obs[adata.obs['tau_positive'] == 1][score_col]
        tau_neg_scores = adata.obs[adata.obs['tau_positive'] == 0][score_col]

        # Statistical test
        t_stat, p_val = stats.ttest_ind(tau_pos_scores, tau_neg_scores)

        # Effect size (Cohen's d)
        pooled_std = np.sqrt(((len(tau_pos_scores)-1)*tau_pos_scores.var() +
                             (len(tau_neg_scores)-1)*tau_neg_scores.var()) /
                            (len(tau_pos_scores) + len(tau_neg_scores) - 2))
        cohens_d = (tau_pos_scores.mean() - tau_neg_scores.mean()) / pooled_std if pooled_std > 0 else 0

        # Fold change
        fold_change = tau_pos_scores.mean() / (tau_neg_scores.mean() + 1e-10)

        comparison_results.append({
            'category': category,
            'tau_pos_mean': tau_pos_scores.mean(),
            'tau_neg_mean': tau_neg_scores.mean(),
            'fold_change': fold_change,
            'p_value': p_val,
            'cohens_d': cohens_d,
            'direction': 'increased' if fold_change > 1 else 'decreased'
        })

# Convert to DataFrame and add FDR correction
df = pd.DataFrame(comparison_results)
if len(df) > 0:
    from statsmodels.stats.multitest import fdrcorrection
    df['p_adjusted'] = fdrcorrection(df['p_value'])[1]
    df['significant'] = df['p_adjusted'] < 0.05

    print("Cell cycle score comparison (Tau+ vs Tau-):")
    print("-"*60)
    for _, row in df.iterrows():
        sig_marker = "*" if row['significant'] else " "
        direction = "↑" if row['fold_change'] > 1 else "↓"
        print(f"{row['category']:20} {direction} {row['fold_change']:5.2f}x  p={row['p_adjusted']:.2e} {sig_marker}")
```

## Visualize Cell Cycle Patterns

```python
if len(df) > 0:
    # 📊 Create visualization
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Fold change comparison
    colors = ['red' if sig else 'gray' for sig in df['significant']]
    bars = ax1.bar(range(len(df)), df['fold_change'], color=colors, alpha=0.7)
    ax1.set_xticks(range(len(df)))
    ax1.set_xticklabels(df['category'], rotation=45, ha='right')
    ax1.set_ylabel('Fold Change (Tau+ / Tau-)')
    ax1.set_title('Cell Cycle Score Changes in Tau+ Neurons')
    ax1.axhline(1, color='black', linestyle='-', linewidth=0.5)
    ax1.grid(True, alpha=0.3)

    # 2. Effect sizes
    ax2.bar(range(len(df)), df['cohens_d'], color=colors, alpha=0.7)
    ax2.set_xticks(range(len(df)))
    ax2.set_xticklabels(df['category'], rotation=45, ha='right')
    ax2.set_ylabel("Cohen's d")
    ax2.set_title('Effect Sizes')
    ax2.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax2.grid(True, alpha=0.3)

    # 3. PCA of cell cycle scores
    if len(phase_scores.columns) > 1:
        pca = PCA(n_components=2)
        pca_coords = pca.fit_transform(phase_scores.fillna(0))

        scatter = ax3.scatter(pca_coords[:, 0], pca_coords[:, 1],
                            c=adata.obs['tau_positive'], cmap='RdBu_r', alpha=0.6, s=20)
        ax3.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} var)')
        ax3.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} var)')
        ax3.set_title('PCA of Cell Cycle Scores')
        plt.colorbar(scatter, ax=ax3, label='Tau Status')

    # 4. Significant categories heatmap
    if df['significant'].any():
        sig_categories = df[df['significant']]['category'].tolist()
        if len(sig_categories) > 1:
            sig_scores = phase_scores[['score_' + cat.replace('score_', '') for cat in sig_categories if 'score_' + cat.replace('score_', '') in phase_scores.columns]]
            if len(sig_scores.columns) > 0:
                # Create heatmap data
                heatmap_data = []
                for tau_status in ['negative', 'positive']:
                    mask = adata.obs['tau_status'] == tau_status
                    mean_scores = sig_scores.loc[mask].mean()
                    heatmap_data.append(mean_scores.values)

                im = ax4.imshow(heatmap_data, cmap='RdBu_r', aspect='auto')
                ax4.set_xticks(range(len(sig_scores.columns)))
                ax4.set_xticklabels([col.replace('score_', '') for col in sig_scores.columns],
                                   rotation=45, ha='right')
                ax4.set_yticks([0, 1])
                ax4.set_yticklabels(['Tau-', 'Tau+'])
                ax4.set_title('Significant Categories Heatmap')
                plt.colorbar(im, ax=ax4, label='Mean Score')

    plt.tight_layout()
    plt.show()
```

## Test Individual Cell Cycle Genes

```python
# 🔬 Test individual genes for most significant category
if df['significant'].any():
    top_category = df[df['significant']].iloc[0]['category']
    print(f"\\n🔍 DETAILED ANALYSIS: {top_category}")
    print("="*50)

    category_genes = found_genes[top_category]
    if category_genes:
        gene_results = []

        for gene in category_genes[:10]:  # Limit to top 10 for brevity
            if gene in protein_names:
                gene_idx = protein_names.index(gene)
                expr = adata.X[:, gene_idx]

                tau_pos_expr = expr[adata.obs['tau_positive'] == 1]
                tau_neg_expr = expr[adata.obs['tau_positive'] == 0]

                # Statistics
                mean_pos = np.mean(tau_pos_expr)
                mean_neg = np.mean(tau_neg_expr)
                log2fc = mean_pos - mean_neg
                t_stat, p_val = stats.ttest_ind(tau_pos_expr, tau_neg_expr)

                gene_results.append({
                    'gene': gene,
                    'log2FC': log2fc,
                    'p_value': p_val
                })

        gene_df = pd.DataFrame(gene_results)
        if len(gene_df) > 0:
            gene_df['p_adjusted'] = fdrcorrection(gene_df['p_value'])[1]
            gene_df = gene_df.sort_values('p_value')

            print(f"Top {top_category} genes:")
            for _, row in gene_df.head(5).iterrows():
                direction = "↑" if row['log2FC'] > 0 else "↓"
                sig = "*" if row['p_adjusted'] < 0.05 else " "
                print(f"  {row['gene']:10} {direction} {abs(row['log2FC']):.3f}  p={row['p_value']:.2e} {sig}")
```

## Evaluate Cell Cycle Re-entry

```python
# 🎯 EVALUATION: Cell Cycle Re-entry in Neurodegeneration
print("\\n" + "="*60)
print("🎯 NEURONAL CELL CYCLE RE-ENTRY ANALYSIS")
print("="*60)

if len(df) > 0:
    # Check for aberrant cell cycle re-entry signature
    g1s_increased = df[df['category'] == 'G1_S_transition']['fold_change'].iloc[0] > 1.2 if 'G1_S_transition' in df['category'].values else False
    s_increased = df[df['category'] == 'S_phase']['fold_change'].iloc[0] > 1.2 if 'S_phase' in df['category'].values else False
    g2m_low = df[df['category'] == 'G2_M_transition']['fold_change'].iloc[0] < 1.1 if 'G2_M_transition' in df['category'].values else True
    damage_increased = df[df['category'] == 'DNA_damage']['fold_change'].iloc[0] > 1.2 if 'DNA_damage' in df['category'].values else False
    apoptosis_increased = df[df['category'] == 'Apoptosis']['fold_change'].iloc[0] > 1.2 if 'Apoptosis' in df['category'].values else False

    # Count significant alterations
    n_sig_categories = df['significant'].sum()
    pct_sig = n_sig_categories / len(df) * 100

    print("🔬 Cell Cycle Re-entry Signature:")
    print(f"✓ G1/S activation: {'YES' if g1s_increased else 'NO'}")
    print(f"✓ S-phase entry: {'YES' if s_increased else 'NO'}")
    print(f"✓ G2/M blocked: {'YES' if g2m_low else 'NO'}")
    print(f"✓ DNA damage response: {'YES' if damage_increased else 'NO'}")
    print(f"✓ Apoptosis activation: {'YES' if apoptosis_increased else 'NO'}")
    print(f"\\nSignificant categories: {n_sig_categories}/{len(df)} ({pct_sig:.1f}%)")

    # Verdict
    if g1s_increased and damage_increased and not g2m_low:
        verdict = "✅ STRONG EVIDENCE"
        explanation = "Clear abortive cell cycle re-entry detected"
    elif n_sig_categories >= 3:
        verdict = "⚠️ MODERATE EVIDENCE"
        explanation = f"Multiple cell cycle alterations detected ({n_sig_categories} categories)"
    elif n_sig_categories >= 1:
        verdict = "⚠️ WEAK EVIDENCE"
        explanation = f"Some cell cycle changes detected"
    else:
        verdict = "❌ NO EVIDENCE"
        explanation = "No significant cell cycle alterations"

    print(f"\\n🎯 VERDICT: {verdict}")
    print(f"📝 {explanation}")

    # Biological interpretation
    print(f"\\n🧬 BIOLOGICAL INTERPRETATION:")
    if verdict.startswith("✅") or verdict.startswith("⚠️"):
        print("• Neurons show aberrant cell cycle re-activation")
        print("• Post-mitotic neurons attempting division")
        print("• Leads to neuronal dysfunction and death")
        print("• Common feature of neurodegeneration")
        print("• Potential therapeutic target: CDK4/6 inhibitors")

        if damage_increased:
            print("• DNA damage response activated")
        if apoptosis_increased:
            print("• Apoptotic pathways engaged")
    else:
        print("• Cell cycle machinery appears quiescent")
        print("• Neurons maintain post-mitotic state")
        print("• No evidence of aberrant division attempts")

else:
    print("❌ Cannot evaluate - insufficient cell cycle gene data")
```

## Save Results

```python
# 💾 Save results
if len(df) > 0:
    summary = {
        'analysis': 'Neuronal cell cycle re-entry',
        'verdict': verdict if 'verdict' in locals() else 'Unknown',
        'categories_tested': len(df),
        'significant_categories': n_sig_categories if 'n_sig_categories' in locals() else 0,
        'g1s_activated': g1s_increased if 'g1s_increased' in locals() else False,
        'dna_damage_activated': damage_increased if 'damage_increased' in locals() else False,
        'genes_found': len(set(all_found))
    }

    if IN_COLAB:
        df.to_csv('cell_cycle_results.csv', index=False)
        pd.DataFrame([summary]).to_csv('cell_cycle_summary.csv', index=False)
        files.download('cell_cycle_results.csv')
        files.download('cell_cycle_summary.csv')
        print("📁 Results downloaded!")
    else:
        df.to_csv('cell_cycle_analysis.csv', index=False)
        pd.DataFrame([summary]).to_csv('cell_cycle_summary.csv', index=False)
        print("📁 Results saved!")

print("\\n✅ Cell cycle analysis complete! 🎉")
```

---

## 🎯 Summary

This simplified cell cycle analysis:
- ✅ **Tests key cell cycle pathways** for aberrant activation
- ✅ **Detects neuronal cell cycle re-entry** patterns
- ✅ **Links to neurodegeneration** mechanisms
- ✅ **Clear biological interpretation** with therapeutic implications
- ✅ **Identifies abortive cell cycle** signatures

**Perfect for studying neuronal dysfunction in disease!** 🧠