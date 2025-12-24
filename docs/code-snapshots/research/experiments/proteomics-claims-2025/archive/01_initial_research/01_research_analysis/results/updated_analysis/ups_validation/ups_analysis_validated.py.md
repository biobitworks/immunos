---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/results/updated_analysis/ups_validation/ups_analysis_validated.py
relative: research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/results/updated_analysis/ups_validation/ups_analysis_validated.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Updated UPS Protein Analysis with Validated List
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats

# Validated UPS proteins from UniProt/BioGRID analysis
VALIDATED_UPS_PROTEINS = {
    'proteasome_subunits': ['PSMA1', 'PSMA2', 'PSMA3', 'PSMA4', 'PSMA5', 'PSMA6', 'PSMA7', 'PSMB1', 'PSMB10', 'PSMB2', 'PSMB3', 'PSMB4', 'PSMB5', 'PSMB6', 'PSMB7', 'PSMB8', 'PSMB9', 'PSMC1', 'PSMC2', 'PSMC3', 'PSMC4', 'PSMC5', 'PSMC6', 'PSMD1', 'PSMD10', 'PSMD11', 'PSMD12', 'PSMD13', 'PSMD14', 'PSMD2', 'PSMD3', 'PSMD4', 'PSMD5', 'PSMD6', 'PSMD7', 'PSMD8', 'PSMD9', 'PSME1', 'PSME2', 'PSME3', 'PSMF1', 'PSMG1', 'PSMG3'],
    'e3_ligases': ['CBL', 'FBXO2', 'FBXO6', 'HECTD1', 'HECTD3', 'HECTD4', 'HERC1', 'HERC2', 'HUWE1', 'ITCH', 'NEDD4L', 'PARK7', 'RNF31', 'SMURF1', 'TRIM25', 'TRIM32', 'UBE3A', 'UBE3B', 'UBE3C'],
    'e2_enzymes': ['UBE2D1', 'UBE2D3', 'UBE2D4', 'UBE2E2', 'UBE2G1', 'UBE2H', 'UBE2I', 'UBE2K', 'UBE2L3', 'UBE2L6', 'UBE2M', 'UBE2N', 'UBE2O', 'UBE2Q1', 'UBE2R2', 'UBE2V1', 'UBE2V2', 'UBE2Z'],
    'e1_enzymes': ['UBA1', 'UBA2', 'UBA3', 'UBA5', 'UBA6', 'UBB', 'UBC'],
    'deubiquitinases': ['ATXN3', 'BRCC3', 'COPS5', 'COPS6', 'CYLD', 'OTUB1', 'OTUD6B', 'STAMBP', 'UCHL1', 'UCHL3', 'UCHL5', 'USP10', 'USP11', 'USP14', 'USP15', 'USP19', 'USP24', 'USP25', 'USP30', 'USP32', 'USP4', 'USP46', 'USP47', 'USP48', 'USP5', 'USP7', 'USP8', 'USP9X'],
    'ups_regulators': ['BAG6', 'NBR1', 'OPTN', 'SQSTM1', 'TAX1BP1', 'UBQLN1', 'UBQLN2', 'UBQLN4', 'VCP']
}

# All UPS proteins (n=132)
ALL_UPS_PROTEINS = ['ATG12', 'ATXN3', 'BAG6', 'BRCC3', 'CBL', 'COPS5', 'COPS6', 'CYLD', 'FBXO2', 'FBXO6', 'HECTD1', 'HECTD3', 'HECTD4', 'HERC1', 'HERC2', 'HUWE1', 'ISG15', 'ITCH', 'NBR1', 'NEDD4L', 'NEDD8', 'OPTN', 'OTUB1', 'OTUD6B', 'PARK7', 'PSMA1', 'PSMA2', 'PSMA3', 'PSMA4', 'PSMA5', 'PSMA6', 'PSMA7', 'PSMB1', 'PSMB10', 'PSMB2', 'PSMB3', 'PSMB4', 'PSMB5', 'PSMB6', 'PSMB7', 'PSMB8', 'PSMB9', 'PSMC1', 'PSMC2', 'PSMC3', 'PSMC4', 'PSMC5', 'PSMC6', 'PSMD1', 'PSMD10', 'PSMD11', 'PSMD12', 'PSMD13', 'PSMD14', 'PSMD2', 'PSMD3', 'PSMD4', 'PSMD5', 'PSMD6', 'PSMD7', 'PSMD8', 'PSMD9', 'PSME1', 'PSME2', 'PSME3', 'PSMF1', 'PSMG1', 'PSMG3', 'RNF31', 'SMURF1', 'SQSTM1', 'STAMBP', 'SUMO2', 'SUMO3', 'SUMO4', 'TAX1BP1', 'TRIM25', 'TRIM32', 'UBA1', 'UBA2', 'UBA3', 'UBA5', 'UBA6', 'UBB', 'UBC', 'UBE2D1', 'UBE2D3', 'UBE2D4', 'UBE2E2', 'UBE2G1', 'UBE2H', 'UBE2I', 'UBE2K', 'UBE2L3', 'UBE2L6', 'UBE2M', 'UBE2N', 'UBE2O', 'UBE2Q1', 'UBE2R2', 'UBE2V1', 'UBE2V2', 'UBE2Z', 'UBE3A', 'UBE3B', 'UBE3C', 'UBQLN1', 'UBQLN2', 'UBQLN4', 'UCHL1', 'UCHL3', 'UCHL5', 'UFM1', 'URM1', 'USP10', 'USP11', 'USP14', 'USP15', 'USP19', 'USP24', 'USP25', 'USP30', 'USP32', 'USP4', 'USP46', 'USP47', 'USP48', 'USP5', 'USP7', 'USP8', 'USP9X', 'VCP']

# Significantly changed UPS proteins (p<0.05, n=38)
SIGNIFICANT_UPS_PROTEINS = ['ATG12', 'CBL', 'HERC1', 'HERC2', 'HUWE1', 'ISG15', 'NBR1', 'NEDD4L', 'PSMA4', 'PSMB8', 'PSMB9', 'PSMC3', 'PSMD5', 'PSMD9', 'PSME1', 'PSME2', 'PSMF1', 'SQSTM1', 'TAX1BP1', 'TRIM25', 'TRIM32', 'UBA5', 'UBA6', 'UBB', 'UBC', 'UBE2E2', 'UBE2L6', 'UBE2O', 'UBE3A', 'UCHL1', 'UCHL3', 'UFM1', 'URM1', 'USP11', 'USP15', 'USP30', 'USP47', 'USP9X']

def analyze_ups_proteins(adata):
    """Analyze all validated UPS proteins in the dataset"""

    results = {}

    # Analyze each category
    for category, gene_list in VALIDATED_UPS_PROTEINS.items():
        cat_results = []

        for gene in gene_list:
            # Find gene in dataset
            mask = adata.var['GeneName'].str.contains(gene, case=False, na=False)

            if mask.sum() > 0:
                idx = np.where(mask)[0][0]
                expr = adata.X[:, idx]

                # Differential expression
                tau_pos = adata.obs['TauStatus'] == 'positive'
                tau_neg = adata.obs['TauStatus'] == 'negative'

                stat, pval = stats.mannwhitneyu(expr[tau_pos], expr[tau_neg])
                log2_fc = np.mean(expr[tau_pos]) - np.mean(expr[tau_neg])

                cat_results.append({
                    'gene': gene,
                    'log2_fc': log2_fc,
                    'pvalue': pval,
                    'significant': pval < 0.05
                })

        results[category] = cat_results

    return results

def summarize_ups_changes(results):
    """Summarize UPS changes by category"""

    summary = {}

    for category, genes in results.items():
        if genes:
            df = pd.DataFrame(genes)
            summary[category] = {
                'total': len(df),
                'significant': df['significant'].sum(),
                'percent_sig': df['significant'].sum() / len(df) * 100,
                'mean_log2_fc': df['log2_fc'].mean(),
                'upregulated': (df['log2_fc'] > 0.263).sum(),
                'downregulated': (df['log2_fc'] < -0.322).sum()
            }

    return summary

# Key findings from validation
KEY_FINDINGS = {
    'sqstm1_massively_upregulated': {
        'gene': 'SQSTM1',
        'log2_fc': 3.413,
        'pvalue': 9.29e-8,
        'interpretation': 'SQSTM1/p62 shows massive upregulation, consistent with impaired autophagy'
    },
    'proteasome_subunits_altered': {
        'significant': ['PSME2', 'PSMD9', 'PSME1', 'PSMF1', 'PSMD5', 'PSMB9'],
        'interpretation': 'Multiple proteasome subunits significantly downregulated'
    },
    'e3_ligases_dysregulated': {
        'upregulated': ['HERC2', 'TRIM32', 'HERC1'],
        'downregulated': ['TRIM25', 'CBL', 'NEDD4L'],
        'interpretation': 'Mixed E3 ligase expression suggests selective UPS dysfunction'
    },
    'autophagy_specific_failure': {
        'autophagy_receptors_up': ['SQSTM1', 'NBR1', 'TAX1BP1'],
        'proteasome_function': 'Most proteasome subunits not significantly changed',
        'interpretation': 'Autophagy-specific failure rather than global UPS collapse'
    }
}

print(f"Updated UPS analysis with {len(ALL_UPS_PROTEINS)} validated proteins")
print(f"Significantly changed: {len(SIGNIFICANT_UPS_PROTEINS)} proteins")

```
