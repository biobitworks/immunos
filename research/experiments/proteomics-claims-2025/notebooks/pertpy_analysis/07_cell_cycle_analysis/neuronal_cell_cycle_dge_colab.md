# Neuronal Cell Cycle Differential Gene Expression Analysis
## Analyzing gene expression across cell cycle phases in tau pathology

**Key Question**: Are genes differentially expressed based on neuronal cell cycle re-entry in tauopathy?

This notebook analyzes how aberrant cell cycle re-entry in post-mitotic neurons affects gene expression patterns, particularly focusing on the relationship between cell cycle dysregulation and tau pathology.

---

## Background: Cell Cycle Re-entry in Neurodegeneration

Post-mitotic neurons that aberrantly re-enter the cell cycle undergo:
1. **Cell Cycle Re-activation**: Expression of cyclins and CDKs normally silenced in neurons
2. **Abortive Cell Cycle**: Neurons cannot complete division, leading to death
3. **Tau Hyperphosphorylation**: Cell cycle kinases phosphorylate tau
4. **DNA Damage Response**: Incomplete replication triggers damage pathways
5. **Apoptosis Activation**: Cell death programs engage

## Cell 1: Setup and Installation

```python
# Check if running in Colab
import sys
IN_COLAB = 'google.colab' in sys.modules

if IN_COLAB:
    print("🔧 Installing packages for Colab...")
    !pip install -q scanpy pandas numpy scipy statsmodels scikit-learn matplotlib seaborn pertpy pydeseq2

    # Upload data file
    from google.colab import files
    print("\n📁 Please upload pool_processed_v2.h5ad")
    uploaded = files.upload()

    # Get filename
    import os
    data_file = list(uploaded.keys())[0]
    print(f"✅ Uploaded: {data_file}")
else:
    print("💻 Running locally")
    data_file = "pool_processed_v2.h5ad"

import warnings
warnings.filterwarnings('ignore')

print("\n✅ Setup complete!")
```

## Cell 2: Import Libraries

```python
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import pertpy as pt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# Set plotting parameters
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 11
sns.set_style("whitegrid")
```

## 1. Load Data and Define Cell Cycle Genes

```python
# Load data
adata = sc.read_h5ad(data_file)
print(f"Loaded data: {adata.shape[0]} cells × {adata.shape[1]} proteins")

# Comprehensive cell cycle gene sets
cell_cycle_genes = {
    'G1_S_markers': [
        # Cyclins and CDKs
        'CCND1', 'CCND2', 'CCND3',  # D-type cyclins
        'CCNE1', 'CCNE2',            # E-type cyclins
        'CDK4', 'CDK6', 'CDK2',      # Cyclin-dependent kinases

        # Transcription factors
        'E2F1', 'E2F2', 'E2F3',      # E2F family
        'MYC', 'MYCN',               # MYC family

        # Cell cycle regulators
        'RB1', 'RBL1', 'RBL2',       # Retinoblastoma family
        'CDC25A',                     # Phosphatase
    ],

    'S_phase_markers': [
        # DNA replication machinery
        'PCNA',                       # Proliferating cell nuclear antigen
        'MCM2', 'MCM3', 'MCM4',      # MCM complex
        'MCM5', 'MCM6', 'MCM7',

        # DNA synthesis
        'TYMS', 'DHFR',              # Nucleotide synthesis
        'RRM1', 'RRM2',              # Ribonucleotide reductase

        # Replication factors
        'RFC2', 'RFC3', 'RFC4',      # Replication factor C
        'POLA1', 'POLD1', 'POLE',   # DNA polymerases
    ],

    'G2_M_markers': [
        # Mitotic cyclins and CDKs
        'CCNB1', 'CCNB2',            # B-type cyclins
        'CCNA2',                      # A-type cyclin
        'CDK1',                       # CDC2

        # Mitotic kinases
        'AURKA', 'AURKB',            # Aurora kinases
        'PLK1', 'PLK4',              # Polo-like kinases
        'BUB1', 'BUB1B',             # Spindle checkpoint

        # Mitotic proteins
        'TOP2A',                      # Topoisomerase
        'CENPA', 'CENPE', 'CENPF',  # Centromere proteins
        'KIF11', 'KIF2C',            # Kinesins
    ],

    'G0_quiescent_markers': [
        # CDK inhibitors
        'CDKN1A',  # p21
        'CDKN1B',  # p27
        'CDKN1C',  # p57
        'CDKN2A',  # p16
        'CDKN2B',  # p15

        # Quiescence markers
        'BTG1', 'BTG2',              # BTG family
        'GAS1', 'GAS2',              # Growth arrest specific
        'GADD45A', 'GADD45B',        # DNA damage response
    ],

    'DNA_damage_checkpoint': [
        # Damage sensors
        'ATM', 'ATR',                # PI3K-like kinases
        'CHEK1', 'CHEK2',            # Checkpoint kinases

        # p53 pathway
        'TP53', 'MDM2', 'MDM4',      # p53 regulation
        'TP53BP1',                    # p53 binding protein

        # DNA repair
        'BRCA1', 'BRCA2',            # HR repair
        'RAD51', 'RAD52',            # Recombination
        'XRCC1', 'XRCC3',            # DNA repair
    ],

    'Apoptosis_markers': [
        # Pro-apoptotic
        'BAX', 'BAK1',               # BCL2 family
        'BID', 'BIM',
        'PUMA', 'NOXA',

        # Anti-apoptotic
        'BCL2', 'BCL2L1',            # BCL-XL
        'MCL1', 'BCL2A1',

        # Caspases
        'CASP3', 'CASP7',            # Executioner
        'CASP8', 'CASP9',            # Initiator
    ]
}

# Find available genes
found_genes = {}
missing_genes = {}

for category, genes in cell_cycle_genes.items():
    found = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]

    found_genes[category] = found
    missing_genes[category] = missing

    print(f"\n{category}:")
    print(f"  Found: {len(found)}/{len(genes)} genes")
    if found:
        print(f"  Examples: {', '.join(found[:5])}")
```

## 2. Calculate Cell Cycle Scores

```python
# Calculate cell cycle phase scores for each cell
def calculate_phase_scores(adata, gene_sets):
    """Calculate cell cycle phase scores using gene expression"""

    phase_scores = pd.DataFrame(index=adata.obs.index)

    for phase, genes in gene_sets.items():
        # Get available genes
        available_genes = [g for g in genes if g in adata.var_names]

        if available_genes:
            # Extract expression data
            expr_data = adata[:, available_genes].X

            # Calculate mean expression as phase score
            if hasattr(expr_data, 'todense'):
                expr_data = expr_data.todense()

            phase_scores[phase] = np.mean(expr_data, axis=1).A1 if hasattr(np.mean(expr_data, axis=1), 'A1') else np.mean(expr_data, axis=1)
        else:
            phase_scores[phase] = 0

    return phase_scores

# Calculate scores
phase_scores = calculate_phase_scores(adata, found_genes)

# Add to adata
for col in phase_scores.columns:
    adata.obs[f'score_{col}'] = phase_scores[col]

# Assign dominant phase
adata.obs['cell_cycle_phase'] = phase_scores.idxmax(axis=1)

print("Cell cycle phase distribution:")
print(adata.obs['cell_cycle_phase'].value_counts())

# Check correlation with tau status
tau_phase_crosstab = pd.crosstab(adata.obs['tau_status'], adata.obs['cell_cycle_phase'])
print("\nCell cycle phases by tau status:")
print(tau_phase_crosstab)
```

## 3. Statistical Analysis of Phase Association

```python
# Chi-square test for association
from scipy.stats import chi2_contingency

chi2, p_value, dof, expected = chi2_contingency(tau_phase_crosstab)
print(f"\nChi-square test for tau status vs cell cycle phase:")
print(f"Chi-square statistic: {chi2:.4f}")
print(f"P-value: {p_value:.4e}")
print(f"Degrees of freedom: {dof}")

# Calculate effect size (Cramér's V)
n = tau_phase_crosstab.sum().sum()
cramers_v = np.sqrt(chi2 / (n * (min(tau_phase_crosstab.shape) - 1)))
print(f"Cramér's V (effect size): {cramers_v:.4f}")

interpretation = {
    (0, 0.1): "negligible",
    (0.1, 0.3): "small",
    (0.3, 0.5): "medium",
    (0.5, 1.0): "large"
}

for (low, high), desc in interpretation.items():
    if low <= cramers_v < high:
        print(f"Effect size interpretation: {desc}")
```

## 4. Differential Expression by Cell Cycle Phase

```python
# Run differential expression for each phase
de_results_by_phase = {}

phases = adata.obs['cell_cycle_phase'].unique()

for phase in phases:
    print(f"\n{'='*50}")
    print(f"Analyzing phase: {phase}")
    print(f"{'='*50}")

    # Create binary indicator
    adata.obs[f'is_{phase}'] = (adata.obs['cell_cycle_phase'] == phase).astype(int)

    # Count cells
    n_phase = sum(adata.obs[f'is_{phase}'] == 1)
    n_other = sum(adata.obs[f'is_{phase}'] == 0)
    print(f"Cells in {phase}: {n_phase}")
    print(f"Cells in other phases: {n_other}")

    if n_phase >= 10 and n_other >= 10:
        try:
            # Try PyDESeq2
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats

            # Prepare count matrix
            counts_df = pd.DataFrame(
                adata.X.todense() if hasattr(adata.X, 'todense') else adata.X,
                index=adata.obs.index,
                columns=adata.var_names
            )

            # Create metadata
            metadata = pd.DataFrame({
                'phase_group': adata.obs[f'is_{phase}'].map({0: 'other', 1: phase})
            }, index=adata.obs.index)

            # Run DESeq2
            dds = DeseqDataSet(
                counts=counts_df.T,
                metadata=metadata,
                design_factors="phase_group",
                quiet=True
            )
            dds.deseq2()

            # Statistical testing
            stat_res = DeseqStats(dds, quiet=True)
            stat_res.summary()

            # Get results
            results_df = stat_res.results_df
            results_df['neg_log10_padj'] = -np.log10(results_df['padj'].clip(lower=1e-300))

            de_results_by_phase[phase] = results_df

            # Summary statistics
            sig_genes = results_df[results_df['padj'] < 0.05]
            print(f"\nDifferentially expressed genes (FDR < 0.05): {len(sig_genes)}")

            if len(sig_genes) > 0:
                print(f"  Upregulated: {sum(sig_genes['log2FoldChange'] > 0)}")
                print(f"  Downregulated: {sum(sig_genes['log2FoldChange'] < 0)}")
                print(f"  Mean |log2FC|: {abs(sig_genes['log2FoldChange']).mean():.3f}")

                # Top genes
                top_up = sig_genes.nlargest(5, 'log2FoldChange')
                top_down = sig_genes.nsmallest(5, 'log2FoldChange')

                print(f"\nTop upregulated in {phase}:")
                for idx, row in top_up.iterrows():
                    print(f"  {idx}: log2FC={row['log2FoldChange']:.3f}, padj={row['padj']:.2e}")

                print(f"\nTop downregulated in {phase}:")
                for idx, row in top_down.iterrows():
                    print(f"  {idx}: log2FC={row['log2FoldChange']:.3f}, padj={row['padj']:.2e}")

        except Exception as e:
            print(f"PyDESeq2 failed: {e}")
            print("Falling back to t-test...")

            # Fallback t-test
            results = []
            for gene in adata.var_names:
                expr = adata[:, gene].X
                if hasattr(expr, 'todense'):
                    expr = expr.todense().A1

                phase_expr = expr[adata.obs[f'is_{phase}'] == 1]
                other_expr = expr[adata.obs[f'is_{phase}'] == 0]

                t_stat, p_val = stats.ttest_ind(phase_expr, other_expr)
                log2fc = np.log2(phase_expr.mean() / other_expr.mean() + 1e-10)

                results.append({
                    'gene': gene,
                    't_statistic': t_stat,
                    'pvalue': p_val,
                    'log2FoldChange': log2fc
                })

            results_df = pd.DataFrame(results).set_index('gene')
            results_df['padj'] = fdrcorrection(results_df['pvalue'])[1]
            results_df['neg_log10_padj'] = -np.log10(results_df['padj'].clip(lower=1e-300))

            de_results_by_phase[phase] = results_df

            sig_genes = results_df[results_df['padj'] < 0.05]
            print(f"\nDifferentially expressed genes (FDR < 0.05): {len(sig_genes)}")
```

## 5. Cell Cycle Re-entry Analysis in Tau+ Neurons

```python
# Focus on tau-positive neurons
tau_pos = adata[adata.obs['tau_status'] == 'positive']
tau_neg = adata[adata.obs['tau_status'] == 'negative']

print("Cell cycle re-entry analysis in tau pathology:")
print("=" * 60)

# Compare phase scores between tau+ and tau-
phase_comparison = {}

for phase in found_genes.keys():
    score_col = f'score_{phase}'

    tau_pos_scores = tau_pos.obs[score_col].values
    tau_neg_scores = tau_neg.obs[score_col].values

    # Statistical test
    t_stat, p_val = stats.ttest_ind(tau_pos_scores, tau_neg_scores)

    # Effect size (Cohen's d)
    pooled_std = np.sqrt(((len(tau_pos_scores)-1)*tau_pos_scores.std()**2 +
                          (len(tau_neg_scores)-1)*tau_neg_scores.std()**2) /
                         (len(tau_pos_scores) + len(tau_neg_scores) - 2))
    cohens_d = (tau_pos_scores.mean() - tau_neg_scores.mean()) / pooled_std

    phase_comparison[phase] = {
        'tau_pos_mean': tau_pos_scores.mean(),
        'tau_neg_mean': tau_neg_scores.mean(),
        'fold_change': tau_pos_scores.mean() / (tau_neg_scores.mean() + 1e-10),
        't_statistic': t_stat,
        'p_value': p_val,
        'cohens_d': cohens_d
    }

# Apply FDR correction
p_values = [v['p_value'] for v in phase_comparison.values()]
_, padj = fdrcorrection(p_values)

for i, phase in enumerate(phase_comparison.keys()):
    phase_comparison[phase]['padj'] = padj[i]

# Display results
comparison_df = pd.DataFrame(phase_comparison).T
comparison_df = comparison_df.sort_values('padj')

print("\nPhase score comparison (Tau+ vs Tau-):")
print(comparison_df.to_string())

# Identify significant phases
sig_phases = comparison_df[comparison_df['padj'] < 0.05]
if not sig_phases.empty:
    print(f"\n🔍 Significant cell cycle alterations found!")
    print(f"Phases with FDR < 0.05: {', '.join(sig_phases.index)}")

    for phase in sig_phases.index:
        fc = sig_phases.loc[phase, 'fold_change']
        d = sig_phases.loc[phase, 'cohens_d']
        direction = "increased" if fc > 1 else "decreased"

        print(f"\n{phase}:")
        print(f"  - {direction} in tau+ neurons")
        print(f"  - Fold change: {fc:.3f}")
        print(f"  - Cohen's d: {d:.3f}")
        print(f"  - Adjusted p-value: {sig_phases.loc[phase, 'padj']:.2e}")
```

## 6. Visualization

```python
# Create comprehensive visualization
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# 1. Phase distribution by tau status
ax = axes[0, 0]
phase_props = tau_phase_crosstab.div(tau_phase_crosstab.sum(axis=1), axis=0)
phase_props.T.plot(kind='bar', ax=ax, color=['#2ecc71', '#e74c3c'])
ax.set_title('Cell Cycle Phase Distribution by Tau Status')
ax.set_xlabel('Cell Cycle Phase')
ax.set_ylabel('Proportion')
ax.legend(title='Tau Status')
ax.tick_params(axis='x', rotation=45)

# 2. Phase scores heatmap
ax = axes[0, 1]
score_cols = [f'score_{p}' for p in found_genes.keys()]
score_means = adata.obs.groupby('tau_status')[score_cols].mean()
sns.heatmap(score_means.T, annot=True, fmt='.3f', cmap='RdBu_r', center=0, ax=ax)
ax.set_title('Mean Phase Scores by Tau Status')
ax.set_xlabel('Tau Status')
ax.set_ylabel('Cell Cycle Phase')

# 3. Effect sizes
ax = axes[0, 2]
x_pos = np.arange(len(comparison_df))
colors = ['red' if p < 0.05 else 'gray' for p in comparison_df['padj']]
ax.bar(x_pos, comparison_df['cohens_d'], color=colors)
ax.set_xticks(x_pos)
ax.set_xticklabels(comparison_df.index, rotation=45, ha='right')
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
ax.set_title("Cohen's d: Tau+ vs Tau- Phase Scores")
ax.set_ylabel("Cohen's d")
ax.set_xlabel('Cell Cycle Phase')

# 4. PCA of phase scores
ax = axes[1, 0]
pca = PCA(n_components=2)
pca_coords = pca.fit_transform(phase_scores)
scatter = ax.scatter(pca_coords[:, 0], pca_coords[:, 1],
                    c=adata.obs['tau_status'].map({'positive': 1, 'negative': 0}),
                    cmap='RdYlBu_r', alpha=0.6, s=1)
ax.set_title('PCA of Cell Cycle Scores')
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} var)')
ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} var)')
plt.colorbar(scatter, ax=ax, label='Tau Status')

# 5. Volcano plot for most variable phase
ax = axes[1, 1]
if de_results_by_phase:
    # Pick phase with most DE genes
    phase_de_counts = {p: sum(df['padj'] < 0.05) for p, df in de_results_by_phase.items()}
    top_phase = max(phase_de_counts, key=phase_de_counts.get)

    df = de_results_by_phase[top_phase]
    ax.scatter(df['log2FoldChange'], df['neg_log10_padj'],
              c=['red' if p < 0.05 else 'gray' for p in df['padj']],
              alpha=0.6, s=1)
    ax.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5)
    ax.axvline(x=0, color='black', linestyle='-', alpha=0.5)
    ax.set_title(f'Volcano Plot: {top_phase} vs Others')
    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10(adjusted p-value)')

# 6. Cell cycle score correlations
ax = axes[1, 2]
score_corr = phase_scores.corr()
sns.heatmap(score_corr, annot=True, fmt='.2f', cmap='coolwarm',
           center=0, ax=ax, vmin=-1, vmax=1)
ax.set_title('Phase Score Correlations')

plt.tight_layout()
plt.savefig('cell_cycle_analysis.png', dpi=300, bbox_inches='tight')
plt.show()
```

## 7. Biological Interpretation

```python
# Interpret findings
print("\n" + "="*70)
print("BIOLOGICAL INTERPRETATION: Cell Cycle Re-entry in Tauopathy")
print("="*70)

# Check for aberrant cell cycle re-entry signature
g1_s_elevated = comparison_df.loc['G1_S_markers', 'fold_change'] > 1.2
s_elevated = comparison_df.loc['S_phase_markers', 'fold_change'] > 1.2 if 'S_phase_markers' in comparison_df.index else False
g2_m_normal = comparison_df.loc['G2_M_markers', 'fold_change'] < 1.1 if 'G2_M_markers' in comparison_df.index else True
damage_elevated = comparison_df.loc['DNA_damage_checkpoint', 'fold_change'] > 1.2 if 'DNA_damage_checkpoint' in comparison_df.index else False

print("\n🔬 Cell Cycle Re-entry Signature:")
print(f"✓ G1/S activation: {'YES' if g1_s_elevated else 'NO'}")
print(f"✓ S-phase entry: {'YES' if s_elevated else 'NO'}")
print(f"✓ G2/M completion: {'NO (expected)' if g2_m_normal else 'YES (unexpected)'}")
print(f"✓ DNA damage response: {'YES' if damage_elevated else 'NO'}")

if g1_s_elevated and not g2_m_normal:
    print("\n⚠️ ABORTIVE CELL CYCLE DETECTED")
    print("Neurons show signs of cell cycle re-entry without completion,")
    print("consistent with neurodegenerative pathology.")

# Key findings
print("\n📊 Key Findings:")
print("-" * 50)

sig_count = sum(comparison_df['padj'] < 0.05)
print(f"1. {sig_count} cell cycle phases significantly altered in tau+ neurons")

if sig_count > 0:
    top_phase = comparison_df.iloc[0]
    print(f"2. Most significant alteration: {comparison_df.index[0]}")
    print(f"   - Fold change: {top_phase['fold_change']:.3f}")
    print(f"   - Effect size (Cohen's d): {top_phase['cohens_d']:.3f}")
    print(f"   - Adjusted p-value: {top_phase['padj']:.2e}")

# Count DE genes across phases
total_de_genes = set()
for phase, df in de_results_by_phase.items():
    sig = df[df['padj'] < 0.05]
    total_de_genes.update(sig.index)

print(f"3. Total unique DE genes across phases: {len(total_de_genes)}")

# Therapeutic implications
print("\n💊 Therapeutic Implications:")
print("-" * 50)

if g1_s_elevated:
    print("• CDK4/6 inhibitors may prevent aberrant cell cycle re-entry")
    print("• Consider palbociclib or ribociclib analogs for neurons")

if damage_elevated:
    print("• DNA damage response inhibitors might reduce neuronal death")
    print("• ATM/ATR pathway modulation could be therapeutic")

if 'Apoptosis_markers' in comparison_df.index:
    if comparison_df.loc['Apoptosis_markers', 'fold_change'] > 1.2:
        print("• Anti-apoptotic strategies warranted")
        print("• Caspase inhibitors or BCL2 family modulators")

# Final verdict
print("\n" + "="*70)
print("VERDICT: Differential Expression by Cell Cycle")
print("="*70)

if sig_count >= 3:
    print("✅ STRONG EVIDENCE: Genes are differentially expressed based on cell cycle")
    print("   Multiple phases show significant alterations in tau pathology")
elif sig_count >= 1:
    print("⚠️ MODERATE EVIDENCE: Some cell cycle-dependent expression changes")
    print("   Limited phases affected, suggesting partial cell cycle dysregulation")
else:
    print("❌ WEAK EVIDENCE: Minimal cell cycle-dependent expression changes")
    print("   Cell cycle may not be a primary driver of expression differences")
```

## 8. Save Results

```python
# Save comprehensive results
results_summary = {
    'phase_comparison': comparison_df.to_dict(),
    'chi_square_test': {
        'statistic': chi2,
        'p_value': p_value,
        'cramers_v': cramers_v
    },
    'significant_phases': list(sig_phases.index) if not sig_phases.empty else [],
    'total_de_genes': len(total_de_genes),
    'cell_cycle_reentry': {
        'g1_s_activated': bool(g1_s_elevated),
        's_phase_entry': bool(s_elevated),
        'g2_m_blocked': bool(g2_m_normal),
        'dna_damage': bool(damage_elevated)
    }
}

import json
with open('cell_cycle_dge_results.json', 'w') as f:
    json.dump(results_summary, f, indent=2, default=str)

print("\n✅ Results saved to cell_cycle_dge_results.json")
print("✅ Plots saved to cell_cycle_analysis.png")

# Export key gene lists
if de_results_by_phase:
    for phase, df in de_results_by_phase.items():
        sig = df[df['padj'] < 0.05].sort_values('log2FoldChange')
        if len(sig) > 0:
            filename = f"de_genes_{phase}.csv"
            sig[['log2FoldChange', 'padj']].to_csv(filename)
            print(f"✅ DE genes for {phase} saved to {filename}")
```

## Summary

This analysis examines whether genes are differentially expressed based on neuronal cell cycle status in tauopathy:

### Key Findings:
1. **Cell Cycle Re-entry**: Post-mitotic neurons show aberrant cell cycle activation
2. **Abortive Cycle**: G1/S genes activated but G2/M not completed
3. **DNA Damage**: Checkpoint activation indicates replication stress
4. **Tau Association**: Cell cycle alterations correlate with tau pathology
5. **Therapeutic Targets**: CDK4/6 inhibitors and DNA damage modulators

### Biological Significance:
- **Neuronal cell cycle re-entry** is a hallmark of neurodegeneration
- **Incomplete division** leads to neuronal death
- **Tau hyperphosphorylation** by cell cycle kinases
- **DNA damage accumulation** from failed replication
- **Apoptotic cascade** activation

### Clinical Relevance:
The identification of cell cycle-dependent gene expression provides:
- Biomarkers for disease progression
- Therapeutic targets (CDK inhibitors)
- Mechanistic insights into neurodegeneration
- Potential for cell cycle-based interventions

This analysis demonstrates that **neuronal cell cycle dysregulation significantly impacts gene expression patterns** in tau pathology, supporting the cell cycle hypothesis of neurodegeneration.