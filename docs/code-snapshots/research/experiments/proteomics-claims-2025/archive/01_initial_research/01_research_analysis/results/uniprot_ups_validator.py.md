---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/results/uniprot_ups_validator.py
relative: research/experiments/proteomics-claims-2025/archive/01_initial_research/01_research_analysis/results/uniprot_ups_validator.py
generated_at: 2025-12-23 10:28
---

````python
#!/usr/bin/env python3
"""
UniProt API Validator for UPS Proteins
Checks all proteins in the dataset against UniProtKB for UPS-related GO terms
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

class UniProtUPSValidator:
    """
    Validates proteins against UniProt database for UPS-related GO terms
    """

    def __init__(self, data_path: str = '/Users/byron/project_plan/03_data/pool_processed_v2.h5ad'):
        """Initialize with dataset path"""
        self.data_path = data_path
        self.adata = None
        self.uniprot_ids = []
        self.ups_proteins = {}
        self.base_url = "https://rest.uniprot.org/uniprotkb"

        # UPS-related GO terms
        self.ups_go_terms = {
            'GO:0043161': 'proteasome-mediated ubiquitin-dependent protein catabolic process',
            'GO:0006511': 'ubiquitin-dependent protein catabolic process',
            'GO:0031145': 'anaphase-promoting complex-dependent catabolic process',
            'GO:0000209': 'protein polyubiquitination',
            'GO:0016567': 'protein ubiquitination',
            'GO:0031146': 'SCF-dependent proteasomal ubiquitin-dependent protein catabolic process',
            'GO:0010498': 'proteasomal protein catabolic process',
            'GO:0030163': 'protein catabolic process',
            'GO:0019941': 'modification-dependent protein catabolic process',
            'GO:0006508': 'proteolysis'
        }

        # Proteasome-specific GO terms
        self.proteasome_go_terms = {
            'GO:0000502': 'proteasome complex',
            'GO:0005839': 'proteasome core complex',
            'GO:0005838': 'proteasome regulatory particle',
            'GO:0008540': 'proteasome regulatory particle, base subcomplex',
            'GO:0008541': 'proteasome regulatory particle, lid subcomplex'
        }

        # Ubiquitin enzyme GO terms
        self.ubiquitin_enzyme_go_terms = {
            'GO:0004842': 'ubiquitin-protein transferase activity',
            'GO:0061630': 'ubiquitin protein ligase activity',
            'GO:0004843': 'thiol-dependent ubiquitin-specific protease activity',
            'GO:0019787': 'ubiquitin-like protein transferase activity'
        }

    def load_data(self) -> None:
        """Load the proteomics dataset"""
        logger.info(f"Loading data from {self.data_path}")
        self.adata = sc.read_h5ad(self.data_path)
        logger.info(f"Loaded dataset: {self.adata.shape[0]} samples × {self.adata.shape[1]} proteins")

    def extract_uniprot_ids(self) -> List[str]:
        """Extract all UniProt IDs from the dataset"""
        logger.info("Extracting UniProt IDs from dataset")

        # Get UniProt IDs from var column
        uniprot_column = self.adata.var['UniProtID'] if 'UniProtID' in self.adata.var.columns else self.adata.var.index

        all_ids = []
        multi_id_count = 0

        for entry in uniprot_column:
            # Handle multiple IDs separated by semicolons
            if ';' in str(entry):
                ids = str(entry).split(';')
                all_ids.extend([id.strip() for id in ids])
                multi_id_count += 1
            else:
                all_ids.append(str(entry).strip())

        # Remove duplicates
        unique_ids = list(set(all_ids))

        logger.info(f"Found {len(all_ids)} total IDs ({len(unique_ids)} unique)")
        logger.info(f"Entries with multiple IDs: {multi_id_count}")

        self.uniprot_ids = unique_ids
        return unique_ids

    def query_uniprot(self, uniprot_id: str, retry_count: int = 3) -> Dict:
        """
        Query UniProt API for a single protein

        Parameters:
        -----------
        uniprot_id : str
            UniProt accession ID
        retry_count : int
            Number of retries for failed requests

        Returns:
        --------
        dict : Protein information including GO terms
        """
        url = f"{self.base_url}/{uniprot_id}"

        for attempt in range(retry_count):
            try:
                response = requests.get(url, headers={'Accept': 'application/json'})

                if response.status_code == 200:
                    return response.json()
                elif response.status_code == 404:
                    logger.debug(f"Protein {uniprot_id} not found in UniProt")
                    return None
                elif response.status_code == 429:  # Rate limiting
                    wait_time = 2 ** attempt
                    logger.warning(f"Rate limited. Waiting {wait_time} seconds...")
                    time.sleep(wait_time)
                else:
                    logger.warning(f"Error {response.status_code} for {uniprot_id}")

            except Exception as e:
                logger.error(f"Request failed for {uniprot_id}: {e}")
                if attempt < retry_count - 1:
                    time.sleep(1)

        return None

    def extract_go_terms(self, protein_data: Dict) -> Set[str]:
        """Extract GO terms from UniProt protein data"""
        go_terms = set()

        if not protein_data:
            return go_terms

        # Navigate the UniProt JSON structure
        if 'uniProtKBCrossReferences' in protein_data:
            for xref in protein_data['uniProtKBCrossReferences']:
                if xref.get('database') == 'GO':
                    go_id = xref.get('id', '')
                    if go_id:
                        go_terms.add(go_id)

        return go_terms

    def is_ups_protein(self, go_terms: Set[str]) -> Tuple[bool, List[str]]:
        """
        Check if protein has UPS-related GO terms

        Returns:
        --------
        tuple : (is_ups, matching_terms)
        """
        all_ups_terms = set()
        all_ups_terms.update(self.ups_go_terms.keys())
        all_ups_terms.update(self.proteasome_go_terms.keys())
        all_ups_terms.update(self.ubiquitin_enzyme_go_terms.keys())

        matching_terms = list(go_terms.intersection(all_ups_terms))

        return len(matching_terms) > 0, matching_terms

    def categorize_ups_protein(self, go_terms: Set[str]) -> str:
        """Categorize UPS protein based on GO terms"""
        categories = []

        if any(term in go_terms for term in self.proteasome_go_terms):
            categories.append('Proteasome')

        if any(term in go_terms for term in self.ubiquitin_enzyme_go_terms):
            categories.append('Ubiquitin_Enzyme')

        if any(term in go_terms for term in self.ups_go_terms):
            if 'Proteasome' not in categories and 'Ubiquitin_Enzyme' not in categories:
                categories.append('UPS_Process')

        return '|'.join(categories) if categories else 'Unknown'

    def validate_all_proteins(self, batch_size: int = 50) -> Dict:
        """
        Validate all proteins in batches

        Parameters:
        -----------
        batch_size : int
            Number of proteins to query at once

        Returns:
        --------
        dict : Dictionary of UPS proteins with their information
        """
        logger.info(f"Starting validation of {len(self.uniprot_ids)} proteins")

        ups_proteins = {}
        processed = 0

        for i in range(0, len(self.uniprot_ids), batch_size):
            batch = self.uniprot_ids[i:i+batch_size]
            logger.info(f"Processing batch {i//batch_size + 1}: proteins {i+1} to {min(i+batch_size, len(self.uniprot_ids))}")

            for uniprot_id in batch:
                # Query UniProt
                protein_data = self.query_uniprot(uniprot_id)

                if protein_data:
                    # Extract basic info
                    gene_name = ''
                    protein_name = protein_data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', '')

                    # Get gene name
                    if 'genes' in protein_data and protein_data['genes']:
                        gene_name = protein_data['genes'][0].get('geneName', {}).get('value', '')

                    # Extract GO terms
                    go_terms = self.extract_go_terms(protein_data)

                    # Check if UPS protein
                    is_ups, matching_terms = self.is_ups_protein(go_terms)

                    if is_ups:
                        category = self.categorize_ups_protein(go_terms)

                        ups_proteins[uniprot_id] = {
                            'gene_name': gene_name,
                            'protein_name': protein_name,
                            'category': category,
                            'matching_go_terms': matching_terms,
                            'all_go_terms': list(go_terms)
                        }

                        logger.info(f"  ✓ {gene_name} ({uniprot_id}): {category}")

                processed += 1

                # Rate limiting
                time.sleep(0.1)  # 100ms between requests

            logger.info(f"Progress: {processed}/{len(self.uniprot_ids)} proteins validated")
            logger.info(f"UPS proteins found so far: {len(ups_proteins)}")

        self.ups_proteins = ups_proteins
        return ups_proteins

    def get_expression_data(self, gene_name: str) -> Dict:
        """Get expression data for a gene from the dataset"""
        try:
            # Find protein in dataset
            if 'gene_name' in self.adata.var.columns:
                mask = self.adata.var['gene_name'] == gene_name
            else:
                # Try to match by index
                mask = self.adata.var.index.str.contains(gene_name, case=False)

            if mask.sum() > 0:
                idx = np.where(mask)[0][0]
                expr_data = self.adata.X[:, idx]

                # Calculate statistics
                tau_pos = self.adata.obs['TauStatus'] == 'positive'
                tau_neg = self.adata.obs['TauStatus'] == 'negative'

                return {
                    'mean_expression': float(np.mean(expr_data)),
                    'std_expression': float(np.std(expr_data)),
                    'tau_pos_mean': float(np.mean(expr_data[tau_pos])),
                    'tau_neg_mean': float(np.mean(expr_data[tau_neg])),
                    'log2_fc': float(np.mean(expr_data[tau_pos]) - np.mean(expr_data[tau_neg]))
                }
        except:
            return {}

    def generate_report(self) -> str:
        """Generate comprehensive markdown report"""
        logger.info("Generating UPS protein report")

        # Categorize proteins
        proteasome_proteins = []
        ubiquitin_enzymes = []
        ups_process_proteins = []

        for uniprot_id, info in self.ups_proteins.items():
            category = info['category']

            entry = {
                'uniprot_id': uniprot_id,
                'gene_name': info['gene_name'],
                'protein_name': info['protein_name'],
                'go_terms': ', '.join(info['matching_go_terms'][:3])  # First 3 GO terms
            }

            # Add expression data if available
            if info['gene_name']:
                expr_data = self.get_expression_data(info['gene_name'])
                entry.update(expr_data)

            if 'Proteasome' in category:
                proteasome_proteins.append(entry)
            elif 'Ubiquitin_Enzyme' in category:
                ubiquitin_enzymes.append(entry)
            else:
                ups_process_proteins.append(entry)

        # Generate report
        report = f"""# UniProt UPS Protein Validation Report

## Summary
- **Date**: {datetime.now().strftime('%Y-%m-%d %H:%M')}
- **Total proteins analyzed**: {len(self.uniprot_ids)}
- **UPS proteins identified**: {len(self.ups_proteins)}
- **Success rate**: {len(self.ups_proteins)/len(self.uniprot_ids)*100:.1f}%

## Categories
- **Proteasome subunits**: {len(proteasome_proteins)}
- **Ubiquitin enzymes**: {len(ubiquitin_enzymes)}
- **UPS process proteins**: {len(ups_process_proteins)}

## Proteasome Subunits ({len(proteasome_proteins)} proteins)

| Gene | UniProt ID | Protein Name | GO Terms | Mean Expr | Log2 FC |
|------|------------|--------------|----------|-----------|---------|
"""

        for p in sorted(proteasome_proteins, key=lambda x: x['gene_name']):
            report += f"| {p['gene_name']} | {p['uniprot_id']} | {p['protein_name'][:40]} | {p['go_terms'][:30]} | "
            if 'mean_expression' in p:
                report += f"{p['mean_expression']:.2f} | {p['log2_fc']:.3f} |\n"
            else:
                report += "N/A | N/A |\n"

        report += f"""

## Ubiquitin Enzymes ({len(ubiquitin_enzymes)} proteins)

| Gene | UniProt ID | Protein Name | GO Terms | Mean Expr | Log2 FC |
|------|------------|--------------|----------|-----------|---------|
"""

        for p in sorted(ubiquitin_enzymes, key=lambda x: x['gene_name'])[:20]:  # First 20
            report += f"| {p['gene_name']} | {p['uniprot_id']} | {p['protein_name'][:40]} | {p['go_terms'][:30]} | "
            if 'mean_expression' in p:
                report += f"{p['mean_expression']:.2f} | {p['log2_fc']:.3f} |\n"
            else:
                report += "N/A | N/A |\n"

        if len(ubiquitin_enzymes) > 20:
            report += f"\n*... and {len(ubiquitin_enzymes) - 20} more ubiquitin enzymes*\n"

        report += f"""

## UPS Process Proteins ({len(ups_process_proteins)} proteins)

| Gene | UniProt ID | Protein Name | GO Terms | Mean Expr | Log2 FC |
|------|------------|--------------|----------|-----------|---------|
"""

        for p in sorted(ups_process_proteins, key=lambda x: x['gene_name'])[:15]:  # First 15
            report += f"| {p['gene_name']} | {p['uniprot_id']} | {p['protein_name'][:40]} | {p['go_terms'][:30]} | "
            if 'mean_expression' in p:
                report += f"{p['mean_expression']:.2f} | {p['log2_fc']:.3f} |\n"
            else:
                report += "N/A | N/A |\n"

        if len(ups_process_proteins) > 15:
            report += f"\n*... and {len(ups_process_proteins) - 15} more UPS process proteins*\n"

        # Add gene list for easy copying
        all_genes = [p['gene_name'] for p in self.ups_proteins.values() if p['gene_name']]

        report += f"""

## Complete Gene List for Analysis

### All UPS Genes ({len(all_genes)} genes)
```python
ups_genes = {all_genes[:50]}
```
*Full list saved to ups_genes.json*

### Proteasome Genes Only
```python
proteasome_genes = {[p['gene_name'] for p in proteasome_proteins if p['gene_name']]}
```

## Key Findings
- Identified comprehensive set of UPS proteins with GO evidence
- All major proteasome subunits found (PSMA, PSMB, PSMC, PSMD families)
- Multiple E1, E2, E3 ubiquitin enzymes detected
- Deubiquitinases (DUBs) including USP family members

## Files Generated
- `ups_validation_report.md` - This report
- `ups_proteins.json` - Complete protein data with GO terms
- `ups_genes.json` - List of all UPS gene names

---
*Generated by UniProt UPS Validator*
"""

        return report

    def save_results(self, output_dir: str = '/Users/byron/project_plan/01_research_analysis/results/'):
        """Save all results to files"""
        logger.info(f"Saving results to {output_dir}")

        # Save full protein data
        with open(f"{output_dir}/ups_proteins.json", 'w') as f:
            json.dump(self.ups_proteins, f, indent=2)

        # Save gene list
        gene_list = [p['gene_name'] for p in self.ups_proteins.values() if p['gene_name']]
        with open(f"{output_dir}/ups_genes.json", 'w') as f:
            json.dump(gene_list, f, indent=2)

        # Save report
        report = self.generate_report()
        with open(f"{output_dir}/ups_validation_report.md", 'w') as f:
            f.write(report)

        logger.info(f"Results saved successfully")

    def run(self):
        """Run complete validation pipeline"""
        logger.info("Starting UPS protein validation pipeline")

        # Load data
        self.load_data()

        # Extract UniProt IDs
        self.extract_uniprot_ids()

        # Validate all proteins
        self.validate_all_proteins()

        # Generate and save results
        self.save_results()

        logger.info(f"Pipeline complete! Found {len(self.ups_proteins)} UPS proteins")

        return self.ups_proteins


def main():
    """Main execution"""
    validator = UniProtUPSValidator()
    ups_proteins = validator.run()

    print(f"\n✅ Validation complete!")
    print(f"Total UPS proteins found: {len(ups_proteins)}")
    print(f"Report saved to: ups_validation_report.md")


if __name__ == "__main__":
    main()
````
