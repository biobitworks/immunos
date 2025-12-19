#!/usr/bin/env python3
"""
Find UPS proteins in dataset using known protein names
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats
from datetime import datetime

# Load data
print("Loading dataset...")
adata = sc.read_h5ad('/Users/byron/project_plan/03_data/pool_processed_v2.h5ad')
print(f"Dataset: {adata.shape[0]} samples × {adata.shape[1]} proteins")

# Known UPS protein patterns
ups_patterns = {
    'Proteasome_20S_alpha': ['PSMA1', 'PSMA2', 'PSMA3', 'PSMA4', 'PSMA5', 'PSMA6', 'PSMA7', 'PSMA8'],
    'Proteasome_20S_beta': ['PSMB1', 'PSMB2', 'PSMB3', 'PSMB4', 'PSMB5', 'PSMB6', 'PSMB7', 'PSMB8', 'PSMB9', 'PSMB10', 'PSMB11'],
    'Proteasome_19S_ATPase': ['PSMC1', 'PSMC2', 'PSMC3', 'PSMC4', 'PSMC5', 'PSMC6'],
    'Proteasome_19S_non_ATPase': ['PSMD1', 'PSMD2', 'PSMD3', 'PSMD4', 'PSMD5', 'PSMD6', 'PSMD7', 'PSMD8', 'PSMD9', 'PSMD10', 'PSMD11', 'PSMD12', 'PSMD13', 'PSMD14'],
    'Proteasome_other': ['PSME1', 'PSME2', 'PSME3', 'PSME4', 'PSMF1', 'PSMG1', 'PSMG2', 'PSMG3', 'PSMG4'],
    'E1_enzymes': ['UBA1', 'UBA2', 'UBA3', 'UBA5', 'UBA6', 'UBA7', 'UBA52', 'UBB', 'UBC', 'RPS27A'],
    'E2_enzymes': ['UBE2A', 'UBE2B', 'UBE2C', 'UBE2D1', 'UBE2D2', 'UBE2D3', 'UBE2D4', 'UBE2E1', 'UBE2E2', 'UBE2E3',
                   'UBE2F', 'UBE2G1', 'UBE2G2', 'UBE2H', 'UBE2I', 'UBE2J1', 'UBE2J2', 'UBE2K', 'UBE2L3', 'UBE2L6',
                   'UBE2M', 'UBE2N', 'UBE2NL', 'UBE2O', 'UBE2Q1', 'UBE2Q2', 'UBE2R1', 'UBE2R2', 'UBE2S', 'UBE2T',
                   'UBE2U', 'UBE2V1', 'UBE2V2', 'UBE2W', 'UBE2Z'],
    'E3_RING': ['MDM2', 'MDM4', 'TRIM21', 'TRIM25', 'TRIM32', 'RNF2', 'RNF4', 'RNF5', 'RNF8', 'RNF168', 'SIAH1', 'SIAH2',
                'CBL', 'CBLB', 'CBLC', 'ITCH', 'NEDD4', 'NEDD4L', 'WWP1', 'WWP2', 'SMURF1', 'SMURF2'],
    'E3_HECT': ['UBE3A', 'UBE3B', 'UBE3C', 'UBE3D', 'HUWE1', 'HERC1', 'HERC2', 'HECTD1', 'HECTD2', 'HECTD3', 'HECTD4'],
    'E3_RBR': ['PARK2', 'PARK7', 'PINK1', 'RNF14', 'RNF31', 'RNF144A', 'RNF144B', 'RNF216', 'RNF217'],
    'E3_CRL_substrate_receptors': ['BTRC', 'FBXW7', 'FBXO1', 'FBXO2', 'FBXO4', 'FBXO6', 'FBXO11', 'FBXO31', 'FBXO32', 'SKP2'],
    'DUBs_USP': ['USP1', 'USP2', 'USP3', 'USP4', 'USP5', 'USP7', 'USP8', 'USP9X', 'USP9Y', 'USP10', 'USP11', 'USP12',
                 'USP13', 'USP14', 'USP15', 'USP16', 'USP18', 'USP19', 'USP20', 'USP21', 'USP22', 'USP24', 'USP25',
                 'USP28', 'USP30', 'USP32', 'USP33', 'USP34', 'USP36', 'USP38', 'USP42', 'USP46', 'USP47', 'USP48'],
    'DUBs_UCH': ['UCHL1', 'UCHL3', 'UCHL5', 'BAP1', 'CYLD'],
    'DUBs_OTU': ['OTUB1', 'OTUB2', 'OTUD1', 'OTUD3', 'OTUD4', 'OTUD5', 'OTUD6A', 'OTUD6B', 'OTUD7A', 'OTUD7B', 'OTULIN'],
    'DUBs_other': ['ATXN3', 'ATXN3L', 'JOSD1', 'JOSD2', 'BRCC3', 'COPS5', 'COPS6', 'STAMBP', 'STAMBPL1'],
    'UPS_regulators': ['SQSTM1', 'NBR1', 'OPTN', 'NDP52', 'TAX1BP1', 'VCP', 'UBQLN1', 'UBQLN2', 'UBQLN4', 'BAG6'],
    'Ubiquitin_like': ['SUMO1', 'SUMO2', 'SUMO3', 'SUMO4', 'NEDD8', 'ISG15', 'FAT10', 'ATG8', 'ATG12', 'URM1', 'UFM1']
}

# Find UPS proteins in dataset
found_proteins = {}
all_ups_genes = []

print("\nSearching for UPS proteins in dataset...\n")

for category, gene_list in ups_patterns.items():
    found_in_category = []

    for gene in gene_list:
        # Search in GeneName column
        mask = adata.var['GeneName'].str.contains(gene, case=False, na=False)

        if mask.sum() > 0:
            indices = np.where(mask)[0]
            for idx in indices:
                gene_name = adata.var['GeneName'].iloc[idx]
                uniprot_id = adata.var['UniprotID'].iloc[idx]

                # Check if it's an exact match or part of the gene name
                gene_names = gene_name.split(';')
                for gn in gene_names:
                    if gene == gn.strip():
                        found_in_category.append({
                            'idx': idx,
                            'gene': gn.strip(),
                            'uniprot': uniprot_id.split(';')[0],
                            'full_name': adata.var['Description'].iloc[idx][:50] if 'Description' in adata.var.columns else ''
                        })
                        all_ups_genes.append(gn.strip())
                        break

    if found_in_category:
        found_proteins[category] = found_in_category
        print(f"{category}: Found {len(found_in_category)} proteins")
        for p in found_in_category[:5]:  # Show first 5
            print(f"  - {p['gene']} ({p['uniprot']})")
        if len(found_in_category) > 5:
            print(f"  ... and {len(found_in_category)-5} more")

# Calculate expression statistics
print("\n" + "="*60)
print("Calculating differential expression statistics...")

ups_results = []

for category, proteins in found_proteins.items():
    for protein in proteins:
        idx = protein['idx']
        expr = adata.X[:, idx]

        tau_pos = adata.obs['TauStatus'] == 'positive'
        tau_neg = adata.obs['TauStatus'] == 'negative'

        # Mann-Whitney U test
        stat, pval = stats.mannwhitneyu(expr[tau_pos], expr[tau_neg])

        # Calculate means and fold change
        mean_pos = np.mean(expr[tau_pos])
        mean_neg = np.mean(expr[tau_neg])
        log2_fc = mean_pos - mean_neg

        ups_results.append({
            'Gene': protein['gene'],
            'Category': category.replace('_', ' '),
            'UniProt': protein['uniprot'],
            'Mean_TauPos': mean_pos,
            'Mean_TauNeg': mean_neg,
            'Log2_FC': log2_fc,
            'P_value': pval,
            'Significant': 'Yes' if pval < 0.05 else 'No'
        })

# Convert to DataFrame and sort
df_results = pd.DataFrame(ups_results)
df_results = df_results.sort_values('P_value')

# Generate report
report = f"""# UPS Proteins Found in Dataset

## Summary
- **Date**: {datetime.now().strftime('%Y-%m-%d %H:%M')}
- **Total UPS proteins found**: {len(all_ups_genes)}
- **Categories analyzed**: {len(found_proteins)}

## Proteins by Category

"""

# Add category summaries
for category, proteins in found_proteins.items():
    category_name = category.replace('_', ' ')
    gene_list = [p['gene'] for p in proteins]
    report += f"### {category_name} ({len(proteins)} proteins)\n"
    report += f"{', '.join(sorted(gene_list))}\n\n"

# Add differential expression table
report += """## Differential Expression Analysis (Tau+ vs Tau-)

### Top 30 Most Significant Changes

| Gene | Category | Log2 FC | P-value | Significant |
|------|----------|---------|---------|-------------|
"""

for _, row in df_results.head(30).iterrows():
    report += f"| {row['Gene']} | {row['Category']} | {row['Log2_FC']:.3f} | {row['P_value']:.4f} | {row['Significant']} |\n"

# Add statistics
n_sig = (df_results['P_value'] < 0.05).sum()
n_up = ((df_results['P_value'] < 0.05) & (df_results['Log2_FC'] > 0.263)).sum()
n_down = ((df_results['P_value'] < 0.05) & (df_results['Log2_FC'] < -0.322)).sum()

report += f"""

## Statistical Summary
- **Significantly changed (p<0.05)**: {n_sig}/{len(df_results)} ({n_sig/len(df_results)*100:.1f}%)
- **Upregulated (FC>1.2)**: {n_up}
- **Downregulated (FC<0.8)**: {n_down}

## Key Findings

### Proteasome Subunits
"""

# Analyze proteasome specifically
proteasome_cats = ['Proteasome_20S_alpha', 'Proteasome_20S_beta', 'Proteasome_19S_ATPase', 'Proteasome_19S_non_ATPase']
proteasome_df = df_results[df_results['Category'].str.contains('Proteasome')]

if not proteasome_df.empty:
    mean_fc = proteasome_df['Log2_FC'].mean()
    sig_count = (proteasome_df['P_value'] < 0.05).sum()

    report += f"""- Found {len(proteasome_df)} proteasome subunits
- Mean Log2 FC: {mean_fc:.3f}
- Significantly changed: {sig_count}/{len(proteasome_df)}

### E3 Ligases
"""

# Analyze E3 ligases
e3_df = df_results[df_results['Category'].str.contains('E3')]
if not e3_df.empty:
    report += f"""- Found {len(e3_df)} E3 ligases
- Significantly changed: {(e3_df['P_value'] < 0.05).sum()}/{len(e3_df)}

### Deubiquitinases (DUBs)
"""

# Analyze DUBs
dub_df = df_results[df_results['Category'].str.contains('DUB')]
if not dub_df.empty:
    report += f"""- Found {len(dub_df)} DUBs
- Significantly changed: {(dub_df['P_value'] < 0.05).sum()}/{len(dub_df)}

## Gene Lists for Analysis

### All UPS Genes
```python
all_ups_genes = {sorted(list(set(all_ups_genes)))}
```

### Proteasome Genes Only
```python
proteasome_genes = {sorted([p['gene'] for cat, prots in found_proteins.items() if 'Proteasome' in cat for p in prots])}
```

### Significantly Changed UPS Genes (p<0.05)
```python
significant_ups = {sorted(df_results[df_results['P_value'] < 0.05]['Gene'].tolist())}
```

## Files Generated
- `ups_proteins_found.md` - This report
- `ups_expression_data.csv` - Full expression data

---
*UPS proteins identified using known gene symbols*
"""

# Save report
with open('/Users/byron/project_plan/01_research_analysis/results/ups_proteins_found.md', 'w') as f:
    f.write(report)

# Save full data
df_results.to_csv('/Users/byron/project_plan/01_research_analysis/results/ups_expression_data.csv', index=False)

print("\n" + "="*60)
print(f"✅ Analysis complete!")
print(f"Total UPS proteins found: {len(all_ups_genes)}")
print(f"Significantly changed: {n_sig}/{len(df_results)}")
print(f"\nReports saved:")
print("  - ups_proteins_found.md")
print("  - ups_expression_data.csv")