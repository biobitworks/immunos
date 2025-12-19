#!/usr/bin/env python3
"""
Validate UPS Proteins using UniProt API
Properly handles the dataset structure with UniProtID column
"""

import requests
import pandas as pd
import numpy as np
import scanpy as sc
import json
import time
from typing import Dict, List, Set, Tuple
import logging
from datetime import datetime

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class UPSValidator:
    """Validates UPS proteins using UniProt API"""

    def __init__(self, data_path: str = '/Users/byron/project_plan/03_data/pool_processed_v2.h5ad'):
        self.data_path = data_path
        self.adata = None
        self.ups_proteins = {}
        self.base_url = "https://rest.uniprot.org/uniprotkb"

        # UPS-related GO terms
        self.ups_go_terms = {
            'GO:0043161': 'proteasome-mediated ubiquitin-dependent protein catabolic process',
            'GO:0006511': 'ubiquitin-dependent protein catabolic process',
            'GO:0000209': 'protein polyubiquitination',
            'GO:0016567': 'protein ubiquitination',
            'GO:0010498': 'proteasomal protein catabolic process',
            'GO:0019941': 'modification-dependent protein catabolic process'
        }

        # Proteasome GO terms
        self.proteasome_go_terms = {
            'GO:0000502': 'proteasome complex',
            'GO:0005839': 'proteasome core complex',
            'GO:0005838': 'proteasome regulatory particle'
        }

    def load_data(self) -> None:
        """Load the proteomics dataset"""
        logger.info(f"Loading data from {self.data_path}")
        self.adata = sc.read_h5ad(self.data_path)
        logger.info(f"Loaded: {self.adata.shape[0]} samples × {self.adata.shape[1]} proteins")

    def extract_uniprot_ids(self) -> List[Tuple[str, str, str]]:
        """Extract UniProt IDs and gene names from dataset"""
        logger.info("Extracting UniProt IDs and gene names")

        protein_list = []

        for idx in range(self.adata.n_vars):
            uniprot_ids = str(self.adata.var['UniprotID'].iloc[idx])
            gene_names = str(self.adata.var['GeneName'].iloc[idx])

            # Handle multiple IDs
            uniprot_list = [uid.strip() for uid in uniprot_ids.split(';')]
            gene_list = [g.strip() for g in gene_names.split(';')]

            # Use first UniProt ID for API query
            primary_id = uniprot_list[0]
            primary_gene = gene_list[0] if gene_list else ''

            protein_list.append((idx, primary_id, primary_gene))

        logger.info(f"Extracted {len(protein_list)} proteins for validation")
        return protein_list

    def query_uniprot_batch(self, uniprot_ids: List[str]) -> Dict:
        """Query UniProt API for multiple proteins at once"""
        # Build query for batch
        query = ' OR '.join([f'accession:{uid}' for uid in uniprot_ids])

        params = {
            'query': query,
            'format': 'json',
            'fields': 'accession,gene_names,protein_name,go,go_f,go_c,go_p',
            'size': '500'
        }

        try:
            response = requests.get(f"{self.base_url}/search", params=params, timeout=30)

            if response.status_code == 200:
                data = response.json()
                results = {}

                for entry in data.get('results', []):
                    accession = entry.get('primaryAccession', '')
                    results[accession] = entry

                return results
            else:
                logger.warning(f"API returned status {response.status_code}")
                return {}

        except Exception as e:
            logger.error(f"API request failed: {e}")
            return {}

    def check_ups_go_terms(self, protein_data: Dict) -> Tuple[bool, List[str], str]:
        """Check if protein has UPS-related GO terms"""
        if not protein_data:
            return False, [], ''

        go_terms = []

        # Extract GO annotations
        for go_ref in protein_data.get('goReferences', []):
            go_id = go_ref.get('goId', '')
            if go_id:
                go_terms.append(go_id)

        # Check for UPS terms
        all_ups_terms = set()
        all_ups_terms.update(self.ups_go_terms.keys())
        all_ups_terms.update(self.proteasome_go_terms.keys())

        matching_terms = [t for t in go_terms if t in all_ups_terms]

        # Categorize
        category = ''
        if any(t in self.proteasome_go_terms for t in matching_terms):
            category = 'Proteasome'
        elif matching_terms:
            category = 'UPS'

        return len(matching_terms) > 0, matching_terms, category

    def validate_proteins_fast(self, batch_size: int = 100) -> Dict:
        """Fast validation using batch queries"""
        protein_list = self.extract_uniprot_ids()

        logger.info(f"Validating {len(protein_list)} proteins in batches of {batch_size}")

        ups_proteins = {}

        for i in range(0, len(protein_list), batch_size):
            batch = protein_list[i:i+batch_size]
            batch_ids = [p[1] for p in batch]

            logger.info(f"Processing batch {i//batch_size + 1}: proteins {i+1}-{min(i+batch_size, len(protein_list))}")

            # Query batch
            results = self.query_uniprot_batch(batch_ids)

            # Process results
            for idx, uniprot_id, gene_name in batch:
                if uniprot_id in results:
                    protein_data = results[uniprot_id]
                    is_ups, go_terms, category = self.check_ups_go_terms(protein_data)

                    if is_ups:
                        ups_proteins[idx] = {
                            'uniprot_id': uniprot_id,
                            'gene_name': gene_name,
                            'category': category,
                            'go_terms': go_terms[:5]  # First 5 GO terms
                        }
                        logger.info(f"  ✓ {gene_name} ({uniprot_id}): {category}")

            # Rate limiting
            time.sleep(0.5)

            logger.info(f"  UPS proteins found so far: {len(ups_proteins)}")

        self.ups_proteins = ups_proteins
        return ups_proteins

    def get_expression_stats(self, idx: int) -> Dict:
        """Get expression statistics for a protein"""
        try:
            expr = self.adata.X[:, idx]
            tau_pos = self.adata.obs['TauStatus'] == 'positive'
            tau_neg = self.adata.obs['TauStatus'] == 'negative'

            from scipy import stats
            stat, pval = stats.mannwhitneyu(expr[tau_pos], expr[tau_neg])

            return {
                'mean_expr': float(np.mean(expr)),
                'tau_pos_mean': float(np.mean(expr[tau_pos])),
                'tau_neg_mean': float(np.mean(expr[tau_neg])),
                'log2_fc': float(np.mean(expr[tau_pos]) - np.mean(expr[tau_neg])),
                'pvalue': float(pval)
            }
        except:
            return {}

    def generate_report(self) -> str:
        """Generate markdown report of UPS proteins"""
        logger.info("Generating UPS validation report")

        # Categorize proteins
        proteasome = []
        ups_other = []

        for idx, info in self.ups_proteins.items():
            entry = {
                'idx': idx,
                'gene': info['gene_name'],
                'uniprot': info['uniprot_id'],
                'category': info['category'],
                'go_terms': ', '.join(info['go_terms'][:3])
            }

            # Add expression data
            expr_stats = self.get_expression_stats(idx)
            entry.update(expr_stats)

            if info['category'] == 'Proteasome':
                proteasome.append(entry)
            else:
                ups_other.append(entry)

        # Sort by gene name
        proteasome.sort(key=lambda x: x['gene'])
        ups_other.sort(key=lambda x: x['gene'])

        # Generate report
        report = f"""# UPS Protein Validation Report

## Summary
- **Date**: {datetime.now().strftime('%Y-%m-%d %H:%M')}
- **Total proteins analyzed**: {self.adata.n_vars}
- **UPS proteins identified**: {len(self.ups_proteins)}
- **Proteasome subunits**: {len(proteasome)}
- **Other UPS proteins**: {len(ups_other)}

## Proteasome Subunits ({len(proteasome)} proteins)

| Gene | UniProt | Log2 FC | P-value | GO Terms |
|------|---------|---------|---------|----------|
"""

        for p in proteasome:
            if 'log2_fc' in p:
                report += f"| {p['gene']} | {p['uniprot']} | {p['log2_fc']:.3f} | {p.get('pvalue', 1):.4f} | {p['go_terms'][:40]} |\n"
            else:
                report += f"| {p['gene']} | {p['uniprot']} | N/A | N/A | {p['go_terms'][:40]} |\n"

        report += f"""

## Other UPS Proteins ({len(ups_other)} proteins)

| Gene | UniProt | Log2 FC | P-value | GO Terms |
|------|---------|---------|---------|----------|
"""

        # Show first 30
        for p in ups_other[:30]:
            if 'log2_fc' in p:
                report += f"| {p['gene']} | {p['uniprot']} | {p['log2_fc']:.3f} | {p.get('pvalue', 1):.4f} | {p['go_terms'][:40]} |\n"
            else:
                report += f"| {p['gene']} | {p['uniprot']} | N/A | N/A | {p['go_terms'][:40]} |\n"

        if len(ups_other) > 30:
            report += f"\n*... and {len(ups_other) - 30} more UPS proteins*\n"

        # Add gene lists
        proteasome_genes = [p['gene'] for p in proteasome if p['gene']]
        ups_genes = [p['gene'] for p in ups_other if p['gene']]
        all_ups = proteasome_genes + ups_genes

        report += f"""

## Gene Lists for Analysis

### All UPS Genes ({len(all_ups)} total)
```python
all_ups_genes = {all_ups[:50]}
```

### Proteasome Genes Only ({len(proteasome_genes)} total)
```python
proteasome_genes = {proteasome_genes}
```

### Key Proteasome Subunits Found
- **20S Core (α subunits)**: {', '.join([g for g in proteasome_genes if g.startswith('PSMA')])}
- **20S Core (β subunits)**: {', '.join([g for g in proteasome_genes if g.startswith('PSMB')])}
- **19S Regulatory (ATPases)**: {', '.join([g for g in proteasome_genes if g.startswith('PSMC')])}
- **19S Regulatory (non-ATPases)**: {', '.join([g for g in proteasome_genes if g.startswith('PSMD')])}

## Statistical Summary

### Differential Expression (Tau+ vs Tau-)
- **Significantly changed (p<0.05)**: {sum(1 for p in proteasome + ups_other if p.get('pvalue', 1) < 0.05)}
- **Upregulated (FC>1.2)**: {sum(1 for p in proteasome + ups_other if p.get('log2_fc', 0) > 0.263)}
- **Downregulated (FC<0.8)**: {sum(1 for p in proteasome + ups_other if p.get('log2_fc', 0) < -0.322)}

## Files Generated
- `ups_validation_report.md` - This report
- `ups_proteins.json` - Complete protein data
- `ups_gene_list.txt` - Gene names for analysis

---
*Generated by UPS Validator using UniProt API*
"""

        return report

    def save_results(self, output_dir: str = '/Users/byron/project_plan/01_research_analysis/results/'):
        """Save validation results"""
        logger.info(f"Saving results to {output_dir}")

        # Save protein data
        with open(f"{output_dir}/ups_proteins.json", 'w') as f:
            json.dump(self.ups_proteins, f, indent=2)

        # Save gene list
        gene_list = [info['gene_name'] for info in self.ups_proteins.values() if info['gene_name']]
        with open(f"{output_dir}/ups_gene_list.txt", 'w') as f:
            f.write('\n'.join(sorted(set(gene_list))))

        # Save report
        report = self.generate_report()
        with open(f"{output_dir}/ups_validation_report.md", 'w') as f:
            f.write(report)

        logger.info("Results saved successfully")

    def run(self):
        """Run complete validation"""
        logger.info("Starting UPS validation")

        # Load data
        self.load_data()

        # Validate proteins
        self.validate_proteins_fast()

        # Save results
        self.save_results()

        logger.info(f"Complete! Found {len(self.ups_proteins)} UPS proteins")

        return self.ups_proteins


def main():
    """Main execution"""
    validator = UPSValidator()
    ups_proteins = validator.run()

    print(f"\n✅ Validation complete!")
    print(f"Total UPS proteins: {len(ups_proteins)}")
    print(f"Report saved to: ups_validation_report.md")


if __name__ == "__main__":
    main()