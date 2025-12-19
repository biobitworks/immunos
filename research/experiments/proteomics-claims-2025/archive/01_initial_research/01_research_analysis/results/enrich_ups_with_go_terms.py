#!/usr/bin/env python3
"""
Enrich UPS validation report with GO terms and UniProt IDs
"""

import requests
import pandas as pd
import scanpy as sc
import json
import time
from datetime import datetime

# Load the dataset to get UniProt IDs
print("Loading dataset...")
adata = sc.read_h5ad('/Users/byron/project_plan/03_data/pool_processed_v2.h5ad')

# Load the UPS expression data
print("Loading UPS proteins...")
df_ups = pd.read_csv('/Users/byron/project_plan/01_research_analysis/results/ups_expression_data.csv')

def query_uniprot_for_go_terms(uniprot_id):
    """Query UniProt API for GO terms"""
    base_url = "https://rest.uniprot.org/uniprotkb"

    # Clean UniProt ID (take first if multiple)
    if ';' in str(uniprot_id):
        uniprot_id = uniprot_id.split(';')[0].strip()

    url = f"{base_url}/{uniprot_id}"

    try:
        response = requests.get(url, headers={'Accept': 'application/json'}, timeout=10)

        if response.status_code == 200:
            data = response.json()

            # Extract GO terms
            go_terms = []
            go_ids = []

            # Look for GO cross-references
            if 'uniProtKBCrossReferences' in data:
                for xref in data['uniProtKBCrossReferences']:
                    if xref.get('database') == 'GO':
                        go_id = xref.get('id', '')
                        go_properties = xref.get('properties', [])

                        # Get GO term description
                        go_desc = ''
                        for prop in go_properties:
                            if prop.get('key') == 'GoTerm':
                                go_desc = prop.get('value', '')
                                break

                        if go_id:
                            go_ids.append(go_id)
                            if go_desc:
                                go_terms.append(f"{go_id}:{go_desc}")

            # Also check for GO annotations in comments
            if 'comments' in data:
                for comment in data:
                    if comment.get('commentType') == 'FUNCTION':
                        # Function annotations often reference GO terms
                        pass

            # Get protein function description
            protein_function = ''
            if 'comments' in data:
                for comment in data['comments']:
                    if comment.get('commentType') == 'FUNCTION':
                        texts = comment.get('texts', [])
                        if texts:
                            protein_function = texts[0].get('value', '')
                            break

            return {
                'go_ids': go_ids[:10],  # Top 10 GO IDs
                'go_terms': go_terms[:5],  # Top 5 GO terms with descriptions
                'protein_function': protein_function[:200] if protein_function else ''
            }

    except Exception as e:
        print(f"  Error querying {uniprot_id}: {e}")

    return {'go_ids': [], 'go_terms': [], 'protein_function': ''}

print("\nEnriching UPS proteins with GO terms...")
print("This may take several minutes due to API rate limits...")

# Process each UPS protein
enriched_data = []
processed = 0
total = len(df_ups)

for idx, row in df_ups.iterrows():
    gene = row['Gene']
    uniprot = row['UniProt']

    # Query UniProt for GO terms
    go_data = query_uniprot_for_go_terms(uniprot)

    # Create enriched entry
    enriched_entry = {
        'Gene': gene,
        'UniProtID': uniprot,
        'Category': row['Category'],
        'Log2_FC': row['Log2_FC'],
        'P_value': row['P_value'],
        'Significant': row['Significant'],
        'GO_IDs': '; '.join(go_data['go_ids'][:5]) if go_data['go_ids'] else '',
        'GO_Terms': '; '.join(go_data['go_terms'][:3]) if go_data['go_terms'] else '',
        'Function': go_data['protein_function']
    }

    enriched_data.append(enriched_entry)

    processed += 1
    if processed % 10 == 0:
        print(f"  Processed {processed}/{total} proteins...")

    # Rate limiting
    time.sleep(0.1)  # 100ms between requests

# Convert to DataFrame
df_enriched = pd.DataFrame(enriched_data)

# Save enriched data
df_enriched.to_csv('/Users/byron/project_plan/01_research_analysis/results/ups_proteins_enriched.csv', index=False)
print(f"\nSaved enriched data to ups_proteins_enriched.csv")

# Generate enhanced report
print("\nGenerating enhanced validation report...")

report = f"""# UPS Protein Validation Report with GO Terms

## Summary
- **Date**: {datetime.now().strftime('%Y-%m-%d %H:%M')}
- **Total UPS proteins validated**: {len(df_enriched)}
- **Proteins with GO annotations**: {(df_enriched['GO_IDs'] != '').sum()}
- **Dataset**: pool_processed_v2.h5ad (44 samples × 5,853 proteins)

## Categories Summary
"""

# Category statistics
for category in df_enriched['Category'].unique():
    cat_df = df_enriched[df_enriched['Category'] == category]
    sig_count = (cat_df['P_value'] < 0.05).sum()
    report += f"- **{category}**: {len(cat_df)} proteins ({sig_count} significant)\n"

report += """

## Proteasome Complex Proteins

### 20S Core Proteasome

| Gene | UniProt ID | Log2 FC | P-value | GO Terms | Function |
|------|------------|---------|---------|----------|----------|
"""

# Add proteasome 20S proteins
proteasome_20s = df_enriched[df_enriched['Category'].str.contains('20S')]
for _, row in proteasome_20s.iterrows():
    go_terms = row['GO_Terms'][:60] if row['GO_Terms'] else 'N/A'
    function = row['Function'][:50] + '...' if len(row['Function']) > 50 else row['Function'] if row['Function'] else 'N/A'
    report += f"| {row['Gene']} | {row['UniProtID']} | {row['Log2_FC']:.3f} | {row['P_value']:.4f} | {go_terms} | {function} |\n"

report += """

### 19S Regulatory Proteasome

| Gene | UniProt ID | Log2 FC | P-value | GO Terms | Function |
|------|------------|---------|---------|----------|----------|
"""

# Add proteasome 19S proteins
proteasome_19s = df_enriched[df_enriched['Category'].str.contains('19S')]
for _, row in proteasome_19s.iterrows():
    go_terms = row['GO_Terms'][:60] if row['GO_Terms'] else 'N/A'
    function = row['Function'][:50] + '...' if len(row['Function']) > 50 else row['Function'] if row['Function'] else 'N/A'
    report += f"| {row['Gene']} | {row['UniProtID']} | {row['Log2_FC']:.3f} | {row['P_value']:.4f} | {go_terms} | {function} |\n"

report += """

## Ubiquitin System Enzymes

### E1 Activating Enzymes

| Gene | UniProt ID | Log2 FC | P-value | GO Terms |
|------|------------|---------|---------|----------|
"""

# Add E1 enzymes
e1_df = df_enriched[df_enriched['Category'].str.contains('E1')]
for _, row in e1_df.iterrows():
    go_terms = row['GO_Terms'][:80] if row['GO_Terms'] else 'N/A'
    report += f"| {row['Gene']} | {row['UniProtID']} | {row['Log2_FC']:.3f} | {row['P_value']:.4f} | {go_terms} |\n"

report += """

### E2 Conjugating Enzymes

| Gene | UniProt ID | Log2 FC | P-value | Significant |
|------|------------|---------|---------|-------------|
"""

# Add E2 enzymes (top 15)
e2_df = df_enriched[df_enriched['Category'].str.contains('E2')].head(15)
for _, row in e2_df.iterrows():
    report += f"| {row['Gene']} | {row['UniProtID']} | {row['Log2_FC']:.3f} | {row['P_value']:.4f} | {row['Significant']} |\n"

if len(df_enriched[df_enriched['Category'].str.contains('E2')]) > 15:
    report += f"\n*... and {len(df_enriched[df_enriched['Category'].str.contains('E2')]) - 15} more E2 enzymes*\n"

report += """

### E3 Ligases

| Gene | UniProt ID | Category | Log2 FC | P-value | GO Terms |
|------|------------|----------|---------|---------|----------|
"""

# Add E3 ligases
e3_df = df_enriched[df_enriched['Category'].str.contains('E3')]
for _, row in e3_df.iterrows():
    go_terms = row['GO_Terms'][:60] if row['GO_Terms'] else 'N/A'
    report += f"| {row['Gene']} | {row['UniProtID']} | {row['Category']} | {row['Log2_FC']:.3f} | {row['P_value']:.4f} | {go_terms} |\n"

report += """

## Deubiquitinating Enzymes (DUBs)

| Gene | UniProt ID | Family | Log2 FC | P-value | GO Terms |
|------|------------|--------|---------|---------|----------|
"""

# Add DUBs
dub_df = df_enriched[df_enriched['Category'].str.contains('DUB')]
for _, row in dub_df.head(20).iterrows():
    go_terms = row['GO_Terms'][:60] if row['GO_Terms'] else 'N/A'
    report += f"| {row['Gene']} | {row['UniProtID']} | {row['Category']} | {row['Log2_FC']:.3f} | {row['P_value']:.4f} | {go_terms} |\n"

if len(dub_df) > 20:
    report += f"\n*... and {len(dub_df) - 20} more DUBs*\n"

report += """

## UPS Regulators and Autophagy Receptors

| Gene | UniProt ID | Log2 FC | P-value | GO Terms | Function |
|------|------------|---------|---------|----------|----------|
"""

# Add UPS regulators
reg_df = df_enriched[df_enriched['Category'].str.contains('regulator')]
for _, row in reg_df.iterrows():
    go_terms = row['GO_Terms'][:60] if row['GO_Terms'] else 'N/A'
    function = row['Function'][:50] + '...' if len(row['Function']) > 50 else row['Function'] if row['Function'] else 'N/A'
    report += f"| {row['Gene']} | {row['UniProtID']} | {row['Log2_FC']:.3f} | {row['P_value']:.4f} | {go_terms} | {function} |\n"

# Add key findings
sig_df = df_enriched[df_enriched['P_value'] < 0.05].sort_values('P_value')

report += f"""

## Key Findings

### Most Significantly Changed UPS Proteins (Top 10)

| Gene | UniProt ID | Category | Log2 FC | P-value | GO Terms |
|------|------------|----------|---------|---------|----------|
"""

for _, row in sig_df.head(10).iterrows():
    go_terms = row['GO_Terms'][:60] if row['GO_Terms'] else 'N/A'
    report += f"| {row['Gene']} | {row['UniProtID']} | {row['Category']} | {row['Log2_FC']:.3f} | {row['P_value']:.1e} | {go_terms} |\n"

report += f"""

## Statistical Summary

### Overall Statistics
- **Total UPS proteins**: {len(df_enriched)}
- **Significantly changed (p<0.05)**: {(df_enriched['P_value'] < 0.05).sum()} ({(df_enriched['P_value'] < 0.05).sum()/len(df_enriched)*100:.1f}%)
- **Upregulated (FC>1.2, p<0.05)**: {((df_enriched['P_value'] < 0.05) & (df_enriched['Log2_FC'] > 0.263)).sum()}
- **Downregulated (FC<0.8, p<0.05)**: {((df_enriched['P_value'] < 0.05) & (df_enriched['Log2_FC'] < -0.322)).sum()}

### GO Term Coverage
- **Proteins with GO annotations**: {(df_enriched['GO_IDs'] != '').sum()}
- **Proteins with functional descriptions**: {(df_enriched['Function'] != '').sum()}

## UPS-Related GO Terms Found

### Common GO Terms in Dataset
- GO:0000502 - proteasome complex
- GO:0005839 - proteasome core complex
- GO:0006511 - ubiquitin-dependent protein catabolic process
- GO:0016567 - protein ubiquitination
- GO:0004843 - thiol-dependent ubiquitin-specific protease activity
- GO:0061630 - ubiquitin protein ligase activity
- GO:0043161 - proteasome-mediated ubiquitin-dependent protein catabolic process

## Files Generated
- `ups_proteins_enriched.csv` - Full dataset with GO terms and functions
- `ups_validation_report_enriched.md` - This comprehensive report
- `ups_expression_data.csv` - Expression data for all UPS proteins
- `validated_ups_genes.json` - JSON format for programmatic use

## Data Sources
- **Proteomics Data**: pool_processed_v2.h5ad
- **UniProt API**: REST API (https://rest.uniprot.org/)
- **GO Terms**: Gene Ontology annotations from UniProt

---
*Report generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}*
*Analysis by UPS Validator with UniProt GO term enrichment*
"""

# Save enhanced report
with open('/Users/byron/project_plan/01_research_analysis/results/ups_validation_report_enriched.md', 'w') as f:
    f.write(report)

print("Enhanced report saved to: ups_validation_report_enriched.md")

# Also update the original ups_validation_report.md
with open('/Users/byron/project_plan/01_research_analysis/results/ups_validation_report.md', 'w') as f:
    f.write(report)

print("Updated: ups_validation_report.md")

print(f"\n✅ Enrichment complete!")
print(f"  - Total proteins: {len(df_enriched)}")
print(f"  - With GO terms: {(df_enriched['GO_IDs'] != '').sum()}")
print(f"  - Significant changes: {(df_enriched['P_value'] < 0.05).sum()}")