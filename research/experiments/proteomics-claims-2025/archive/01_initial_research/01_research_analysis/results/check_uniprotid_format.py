#!/usr/bin/env python3
"""
Check the format of UniProt IDs in the dataset to understand semicolon usage
"""

import scanpy as sc
import pandas as pd

# Load the data
print("Loading data...")
adata = sc.read_h5ad('/Users/byron/project_plan/data/pool_processed_v2.h5ad')

print(f"\nDataset shape: {adata.shape}")
print(f"Samples: {adata.n_obs}, Proteins: {adata.n_vars}")

# Check the protein (variable) information
print("\n" + "="*60)
print("CHECKING PROTEIN ANNOTATIONS")
print("="*60)

# Show column names in the var dataframe
print("\nColumns in adata.var:")
print(adata.var.columns.tolist())

# Check first few entries
print("\nFirst 10 protein entries:")
print(adata.var.head(10))

# Check if UniprotID column exists and examine it
if 'UniprotID' in adata.var.columns:
    print("\n" + "="*60)
    print("UNIPROT ID FORMAT ANALYSIS")
    print("="*60)

    # Get sample of UniProt IDs
    uniprot_ids = adata.var['UniprotID'].head(20)
    print("\nFirst 20 UniProt IDs:")
    for i, uid in enumerate(uniprot_ids):
        print(f"{i:3}: {uid}")

    # Check for semicolons
    ids_with_semicolon = adata.var['UniprotID'].str.contains(';', na=False)
    count_with_semicolon = ids_with_semicolon.sum()

    print(f"\nUniProt IDs containing semicolons: {count_with_semicolon}/{len(adata.var)} ({count_with_semicolon/len(adata.var)*100:.1f}%)")

    if count_with_semicolon > 0:
        print("\nExamples of UniProt IDs with semicolons:")
        examples = adata.var[ids_with_semicolon]['UniprotID'].head(10)
        for uid in examples:
            print(f"  - {uid}")
            # Split and show what the semicolon separates
            parts = uid.split(';')
            if len(parts) > 1:
                print(f"    Parts: {parts}")

        print("\n" + "="*60)
        print("EXPLANATION OF SEMICOLONS IN UNIPROT IDs")
        print("="*60)
        print("""
Semicolons in UniProt IDs typically indicate:

1. **Multiple isoforms**: Same protein, different splice variants
   Example: P04637;P04637-2 (p53 isoforms)

2. **Multiple UniProt entries**: Protein maps to multiple database entries
   Example: When a protein has been independently sequenced multiple times

3. **Gene duplications**: Paralogous proteins with high similarity
   Example: Multiple histone variants

4. **Database merging artifacts**: Historical reasons from database consolidation

This is NORMAL in proteomics data and indicates that the mass spec peptides
could match multiple protein forms or entries.
        """)

# Check UniprotName column if it exists
if 'UniprotName' in adata.var.columns:
    print("\n" + "="*60)
    print("UNIPROT NAME FORMAT")
    print("="*60)

    print("\nFirst 10 UniProt Names:")
    for i, name in enumerate(adata.var['UniprotName'].head(10)):
        print(f"{i:3}: {name}")

    # Check for semicolons in names too
    names_with_semicolon = adata.var['UniprotName'].str.contains(';', na=False)
    if names_with_semicolon.sum() > 0:
        print(f"\nUniProt Names with semicolons: {names_with_semicolon.sum()}")

# Check GeneName column
if 'GeneName' in adata.var.columns:
    print("\n" + "="*60)
    print("GENE NAME FORMAT")
    print("="*60)

    print("\nFirst 10 Gene Names:")
    for i, name in enumerate(adata.var['GeneName'].head(10)):
        print(f"{i:3}: {name}")

    # Check for semicolons in gene names
    genes_with_semicolon = adata.var['GeneName'].str.contains(';', na=False)
    if genes_with_semicolon.sum() > 0:
        print(f"\nGene Names with semicolons: {genes_with_semicolon.sum()}")
        print("Examples:")
        for gene in adata.var[genes_with_semicolon]['GeneName'].head(5):
            print(f"  - {gene}")

# Show how to handle these in analysis
print("\n" + "="*60)
print("HOW TO HANDLE IN ANALYSIS")
print("="*60)
print("""
When analyzing proteins with semicolon-separated IDs:

1. For most analyses: Use the first ID (primary isoform)
   uniprot_primary = adata.var['UniprotID'].str.split(';').str[0]

2. For isoform-specific analysis: Keep all IDs
   all_ids = adata.var['UniprotID'].str.split(';').explode()

3. For gene-level analysis: Use GeneName instead
   genes = adata.var['GeneName']

4. For pathway analysis: Map all isoforms to same gene
   This ensures you don't miss proteins due to isoform differences
""")