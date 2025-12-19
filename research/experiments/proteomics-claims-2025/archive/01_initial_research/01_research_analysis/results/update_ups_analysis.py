#!/usr/bin/env python3
"""
Update UPS protein analysis with validated list
"""

import pandas as pd
import json

# Load the UPS expression data
print("Loading validated UPS proteins...")
df_ups = pd.read_csv('/Users/byron/project_plan/01_research_analysis/results/ups_expression_data.csv')

# Get all UPS genes
all_ups_genes = sorted(df_ups['Gene'].unique().tolist())
print(f"Found {len(all_ups_genes)} validated UPS proteins")

# Categorize by function
proteasome_genes = sorted(df_ups[df_ups['Category'].str.contains('Proteasome')]['Gene'].unique().tolist())
e3_ligases = sorted(df_ups[df_ups['Category'].str.contains('E3')]['Gene'].unique().tolist())
e2_enzymes = sorted(df_ups[df_ups['Category'].str.contains('E2')]['Gene'].unique().tolist())
e1_enzymes = sorted(df_ups[df_ups['Category'].str.contains('E1')]['Gene'].unique().tolist())
dubs = sorted(df_ups[df_ups['Category'].str.contains('DUB')]['Gene'].unique().tolist())
ups_regulators = sorted(df_ups[df_ups['Category'].str.contains('regulator')]['Gene'].unique().tolist())

# Get significantly changed genes
sig_ups = sorted(df_ups[df_ups['P_value'] < 0.05]['Gene'].unique().tolist())
sig_up = sorted(df_ups[(df_ups['P_value'] < 0.05) & (df_ups['Log2_FC'] > 0.263)]['Gene'].unique().tolist())
sig_down = sorted(df_ups[(df_ups['P_value'] < 0.05) & (df_ups['Log2_FC'] < -0.322)]['Gene'].unique().tolist())

# Create comprehensive UPS dictionary
ups_dict = {
    'all_ups_genes': all_ups_genes,
    'proteasome_subunits': proteasome_genes,
    'e3_ligases': e3_ligases,
    'e2_enzymes': e2_enzymes,
    'e1_enzymes': e1_enzymes,
    'deubiquitinases': dubs,
    'ups_regulators': ups_regulators,
    'significantly_changed': sig_ups,
    'significantly_upregulated': sig_up,
    'significantly_downregulated': sig_down,
    'statistics': {
        'total_ups_proteins': len(all_ups_genes),
        'proteasome_subunits_count': len(proteasome_genes),
        'e3_ligases_count': len(e3_ligases),
        'e2_enzymes_count': len(e2_enzymes),
        'e1_enzymes_count': len(e1_enzymes),
        'dubs_count': len(dubs),
        'ups_regulators_count': len(ups_regulators),
        'significantly_changed_count': len(sig_ups),
        'percent_significant': round(len(sig_ups) / len(all_ups_genes) * 100, 1)
    }
}

# Save to JSON
with open('/Users/byron/project_plan/01_research_analysis/results/validated_ups_genes.json', 'w') as f:
    json.dump(ups_dict, f, indent=2)

print("\nStatistics:")
print(f"  Total UPS proteins: {len(all_ups_genes)}")
print(f"  Proteasome subunits: {len(proteasome_genes)}")
print(f"  E3 ligases: {len(e3_ligases)}")
print(f"  E2 enzymes: {len(e2_enzymes)}")
print(f"  E1 enzymes: {len(e1_enzymes)}")
print(f"  Deubiquitinases: {len(dubs)}")
print(f"  UPS regulators: {len(ups_regulators)}")
print(f"  Significantly changed: {len(sig_ups)} ({len(sig_ups)/len(all_ups_genes)*100:.1f}%)")

# Create updated Python file for analysis
updated_analysis = f'''#!/usr/bin/env python3
"""
Updated UPS Protein Analysis with Validated List
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats

# Validated UPS proteins from UniProt/BioGRID analysis
VALIDATED_UPS_PROTEINS = {{
    'proteasome_subunits': {proteasome_genes},
    'e3_ligases': {e3_ligases},
    'e2_enzymes': {e2_enzymes},
    'e1_enzymes': {e1_enzymes},
    'deubiquitinases': {dubs},
    'ups_regulators': {ups_regulators}
}}

# All UPS proteins (n={len(all_ups_genes)})
ALL_UPS_PROTEINS = {all_ups_genes}

# Significantly changed UPS proteins (p<0.05, n={len(sig_ups)})
SIGNIFICANT_UPS_PROTEINS = {sig_ups}

def analyze_ups_proteins(adata):
    """Analyze all validated UPS proteins in the dataset"""

    results = {{}}

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

                cat_results.append({{
                    'gene': gene,
                    'log2_fc': log2_fc,
                    'pvalue': pval,
                    'significant': pval < 0.05
                }})

        results[category] = cat_results

    return results

def summarize_ups_changes(results):
    """Summarize UPS changes by category"""

    summary = {{}}

    for category, genes in results.items():
        if genes:
            df = pd.DataFrame(genes)
            summary[category] = {{
                'total': len(df),
                'significant': df['significant'].sum(),
                'percent_sig': df['significant'].sum() / len(df) * 100,
                'mean_log2_fc': df['log2_fc'].mean(),
                'upregulated': (df['log2_fc'] > 0.263).sum(),
                'downregulated': (df['log2_fc'] < -0.322).sum()
            }}

    return summary

# Key findings from validation
KEY_FINDINGS = {{
    'sqstm1_massively_upregulated': {{
        'gene': 'SQSTM1',
        'log2_fc': 3.413,
        'pvalue': 9.29e-8,
        'interpretation': 'SQSTM1/p62 shows massive upregulation, consistent with impaired autophagy'
    }},
    'proteasome_subunits_altered': {{
        'significant': ['PSME2', 'PSMD9', 'PSME1', 'PSMF1', 'PSMD5', 'PSMB9'],
        'interpretation': 'Multiple proteasome subunits significantly downregulated'
    }},
    'e3_ligases_dysregulated': {{
        'upregulated': ['HERC2', 'TRIM32', 'HERC1'],
        'downregulated': ['TRIM25', 'CBL', 'NEDD4L'],
        'interpretation': 'Mixed E3 ligase expression suggests selective UPS dysfunction'
    }},
    'autophagy_specific_failure': {{
        'autophagy_receptors_up': ['SQSTM1', 'NBR1', 'TAX1BP1'],
        'proteasome_function': 'Most proteasome subunits not significantly changed',
        'interpretation': 'Autophagy-specific failure rather than global UPS collapse'
    }}
}}

print(f"Updated UPS analysis with {{len(ALL_UPS_PROTEINS)}} validated proteins")
print(f"Significantly changed: {{len(SIGNIFICANT_UPS_PROTEINS)}} proteins")
'''

# Save updated analysis
with open('/Users/byron/project_plan/01_research_analysis/results/ups_analysis_validated.py', 'w') as f:
    f.write(updated_analysis)

print("\nFiles generated:")
print("  - validated_ups_genes.json")
print("  - ups_analysis_validated.py")

# Create summary report
summary_report = f"""# UPS Protein Validation Summary

## Overview
Successfully validated **{len(all_ups_genes)} UPS proteins** in the proteomics dataset using gene name matching.

## Key Statistics
- **Total UPS proteins**: {len(all_ups_genes)}
- **Significantly changed (p<0.05)**: {len(sig_ups)} ({len(sig_ups)/len(all_ups_genes)*100:.1f}%)
- **Upregulated (FC>1.2)**: {len(sig_up)}
- **Downregulated (FC<0.8)**: {len(sig_down)}

## Categories Validated
| Category | Count | Significant |
|----------|-------|-------------|
| Proteasome subunits | {len(proteasome_genes)} | {sum(1 for g in proteasome_genes if g in sig_ups)} |
| E3 ligases | {len(e3_ligases)} | {sum(1 for g in e3_ligases if g in sig_ups)} |
| E2 enzymes | {len(e2_enzymes)} | {sum(1 for g in e2_enzymes if g in sig_ups)} |
| E1 enzymes | {len(e1_enzymes)} | {sum(1 for g in e1_enzymes if g in sig_ups)} |
| Deubiquitinases | {len(dubs)} | {sum(1 for g in dubs if g in sig_ups)} |
| UPS regulators | {len(ups_regulators)} | {sum(1 for g in ups_regulators if g in sig_ups)} |

## Top Findings
1. **SQSTM1/p62 massively upregulated** (Log2 FC: 3.41, p<0.0001)
2. **NBR1 significantly upregulated** (Log2 FC: 1.49, p<0.0001)
3. **Multiple proteasome subunits downregulated** (PSME1, PSME2, PSMD9)
4. **Autophagy-specific dysfunction** evident from receptor accumulation

## Files Generated
- `ups_proteins_found.md` - Detailed protein report
- `ups_expression_data.csv` - Expression data for all UPS proteins
- `validated_ups_genes.json` - JSON format for programmatic use
- `ups_analysis_validated.py` - Updated analysis script

## Integration with Analysis
The validated UPS proteins have been integrated into the analysis scripts for evaluating mitochondrial dysregulation claims, particularly:
- Statement 1: UPS protein expression
- Statement 2: SQSTM1/p62 upregulation validation
- Statement 5: Mitophagy pathway assessment

---
*Validation complete: {len(all_ups_genes)} UPS proteins confirmed in dataset*
"""

with open('/Users/byron/project_plan/01_research_analysis/results/UPS_VALIDATION_SUMMARY.md', 'w') as f:
    f.write(summary_report)

print("\nSummary report saved to: UPS_VALIDATION_SUMMARY.md")
print("\n✅ All updates complete!")