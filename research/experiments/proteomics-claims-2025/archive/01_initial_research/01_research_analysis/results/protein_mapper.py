#!/usr/bin/env python3
"""
Protein Mapper Utility
Maps gene names to protein indices in the dataset
"""

import scanpy as sc
import pandas as pd

class ProteinMapper:
    """Maps gene names to dataset indices"""

    def __init__(self, adata):
        self.adata = adata
        self._create_mapping()

    def _create_mapping(self):
        """Create gene name to index mapping"""
        self.gene_to_idx = {}

        for idx in self.adata.var_names:
            gene_name = self.adata.var.loc[idx, 'GeneName']
            if pd.notna(gene_name):
                # Handle multiple gene names separated by semicolon
                genes = gene_name.split(';')
                for gene in genes:
                    gene = gene.strip()
                    if gene not in self.gene_to_idx:
                        self.gene_to_idx[gene] = []
                    self.gene_to_idx[gene].append(idx)

    def get_protein_index(self, gene_name):
        """Get protein index(es) for a gene name"""
        if gene_name in self.gene_to_idx:
            return self.gene_to_idx[gene_name][0]  # Return first match

        # Try partial match
        for gene, indices in self.gene_to_idx.items():
            if gene_name in gene or gene in gene_name:
                return indices[0]

        return None

    def get_all_indices(self, gene_name):
        """Get all protein indices for a gene name"""
        if gene_name in self.gene_to_idx:
            return self.gene_to_idx[gene_name]
        return []

    def find_proteins(self, protein_list):
        """Find multiple proteins and return found/missing"""
        found = {}
        missing = []

        for protein in protein_list:
            idx = self.get_protein_index(protein)
            if idx:
                found[protein] = idx
            else:
                missing.append(protein)

        return found, missing

    def get_protein_info(self, gene_name):
        """Get full protein information"""
        idx = self.get_protein_index(gene_name)
        if idx:
            info = {
                'index': idx,
                'uniprot_id': self.adata.var.loc[idx, 'UniprotID'],
                'uniprot_name': self.adata.var.loc[idx, 'UniprotName'],
                'gene_name': self.adata.var.loc[idx, 'GeneName'],
                'description': self.adata.var.loc[idx, 'Description']
            }
            return info
        return None


def load_data_with_mapper(data_path):
    """Load data and create protein mapper"""
    adata = sc.read_h5ad(data_path)
    mapper = ProteinMapper(adata)
    return adata, mapper


# Commonly searched proteins in the analyses
PROTEOSTASIS_PROTEINS = {
    'v_atpase': ['ATP6V1A', 'ATP6V1B2', 'ATP6V0A1', 'ATP6V0D1', 'ATP6V1E1',
                 'ATP6V1H', 'ATP6V1C1', 'ATP6V1D', 'ATP6V1F', 'ATP6V1G1'],
    'retromer': ['VPS35', 'VPS29', 'VPS26A', 'VPS26B'],
    'lysosomes': ['LAMP1', 'LAMP2', 'CTSD', 'CTSB', 'CTSL'],
    'stress': ['HSPA5', 'HSP90AA1', 'HSPA8', 'HSPB1', 'DNAJB1'],
    'transport': ['RAB5A', 'RAB7A', 'RAB11A', 'SNX1', 'SNX2']
}

MITOCHONDRIAL_PROTEINS = {
    'ups': ['UBE2D3', 'UBE2N', 'UBE2K', 'PSMA1', 'PSMA2', 'PSMA3',
            'PSMB1', 'PSMB5', 'PSMD11', 'PSMD14'],
    'autophagy': ['SQSTM1', 'BECN1', 'ATG5', 'ATG7', 'LC3B', 'GABARAPL2'],
    'mitochondria': ['VDAC1', 'VDAC2', 'CYCS', 'COX4I1', 'TOMM20',
                     'PINK1', 'PRKN', 'MFN1', 'MFN2', 'OPA1'],
    'oxidative': ['SOD1', 'SOD2', 'CAT', 'GPX1', 'PRDX1', 'PRDX3']
}

if __name__ == "__main__":
    # Test the mapper
    data_path = '/Users/byron/project_plan/03_data/pool_processed_v2.h5ad'
    adata, mapper = load_data_with_mapper(data_path)

    print("Testing Protein Mapper")
    print("=" * 50)

    # Test key proteins
    test_proteins = ['SQSTM1', 'ATP6V0A1', 'BECN1', 'VPS35', 'UBE2D3']
    found, missing = mapper.find_proteins(test_proteins)

    print(f"\nSearched for {len(test_proteins)} proteins:")
    print(f"Found: {len(found)}")
    print(f"Missing: {len(missing)}")

    for protein, idx in found.items():
        info = mapper.get_protein_info(protein)
        print(f"\n{protein}:")
        print(f"  Index: {idx}")
        print(f"  UniProt: {info['uniprot_id']}")
        print(f"  Description: {info['description'][:60]}...")