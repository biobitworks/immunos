# Data Preparation for PertPy DGE Analysis (Colab Version)
## Self-Contained Notebook for Google Colab

This notebook prepares the pool_processed_v2.h5ad dataset for PertPy-based differential expression analysis.
**Colab-optimized**: All protein sets are embedded directly in the notebook - no external files needed!

---

## Cell 1: Google Colab Setup

```python
# Google Colab Setup (skip if running locally)
import os
IN_COLAB = 'COLAB_GPU' in os.environ

if IN_COLAB:
    print("Running in Google Colab")
    # Mount Google Drive if needed
    from google.colab import drive
    drive.mount('/content/drive')

    # Install required packages
    !pip install -q scanpy pertpy pydeseq2

    # Upload file prompt
    from google.colab import files
    print("\nPlease upload pool_processed_v2.h5ad when prompted:")
    uploaded = files.upload()
    data_path = 'pool_processed_v2.h5ad'
else:
    print("Running locally")
    # Adjust path for local environment
    data_path = '../../data/pool_processed_v2.h5ad'
```

## Cell 2: Import Required Packages

```python
# Import required packages
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# Set plotting parameters
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white')
sns.set_style('whitegrid')

print("Packages imported successfully")
```

## Cell 3: Define All Protein Sets (Embedded)

All protein sets are defined here - no external files needed!

```python
# COMPREHENSIVE PROTEIN SET DEFINITIONS
# All 132 UPS proteins and other key protein sets embedded directly

protein_sets = {
    # Complete 132 UPS proteins list
    'ups_proteins': [
        # Proteasome Core Subunits (20S)
        'PSMA1', 'PSMA2', 'PSMA3', 'PSMA4', 'PSMA5', 'PSMA6', 'PSMA7',  # Alpha subunits
        'PSMB1', 'PSMB2', 'PSMB3', 'PSMB4', 'PSMB5', 'PSMB6', 'PSMB7', 'PSMB8', 'PSMB9', 'PSMB10',  # Beta subunits

        # Proteasome Regulatory Subunits (19S)
        'PSMC1', 'PSMC2', 'PSMC3', 'PSMC4', 'PSMC5', 'PSMC6',  # ATPases
        'PSMD1', 'PSMD2', 'PSMD3', 'PSMD4', 'PSMD5', 'PSMD6', 'PSMD7', 'PSMD8',
        'PSMD9', 'PSMD10', 'PSMD11', 'PSMD12', 'PSMD13', 'PSMD14',  # Non-ATPases

        # Alternative Proteasome Caps
        'PSME1', 'PSME2', 'PSME3', 'PSMF1', 'PSMG1', 'PSMG3',

        # E1 Ubiquitin-Activating Enzymes
        'UBA1', 'UBA2', 'UBA3', 'UBA5', 'UBA6', 'UBB', 'UBC',

        # E2 Ubiquitin-Conjugating Enzymes
        'UBE2A', 'UBE2B', 'UBE2C', 'UBE2D1', 'UBE2D2', 'UBE2D3', 'UBE2D4',
        'UBE2E1', 'UBE2E2', 'UBE2E3', 'UBE2G1', 'UBE2G2', 'UBE2H', 'UBE2I',
        'UBE2J1', 'UBE2J2', 'UBE2K', 'UBE2L3', 'UBE2L6', 'UBE2M', 'UBE2N',
        'UBE2O', 'UBE2Q1', 'UBE2R2', 'UBE2S', 'UBE2V1', 'UBE2V2', 'UBE2Z',

        # E3 Ubiquitin Ligases
        'CBL', 'FBXO2', 'FBXO6', 'FBXW7', 'HECTD1', 'HECTD3', 'HECTD4',
        'HERC1', 'HERC2', 'HUWE1', 'ITCH', 'MDM2', 'NEDD4', 'NEDD4L',
        'PARK2', 'PARK7', 'RNF31', 'SMURF1', 'SMURF2', 'TRIM25', 'TRIM32',
        'UBE3A', 'UBE3B', 'UBE3C', 'UBR4', 'VHL',

        # Deubiquitinating Enzymes (DUBs)
        'ATXN3', 'BRCC3', 'COPS5', 'COPS6', 'CYLD', 'OTUB1', 'OTUD6B',
        'STAMBP', 'UCHL1', 'UCHL3', 'UCHL5', 'USP4', 'USP5', 'USP7',
        'USP8', 'USP9X', 'USP10', 'USP11', 'USP14', 'USP15', 'USP19',
        'USP24', 'USP25', 'USP30', 'USP32', 'USP46', 'USP47', 'USP48',

        # UPS Regulators and Adaptors
        'BAG6', 'NBR1', 'OPTN', 'SQSTM1', 'TAX1BP1',
        'UBQLN1', 'UBQLN2', 'UBQLN4', 'VCP',

        # Alternative Modifiers
        'ATG12', 'ISG15', 'NEDD8', 'SUMO1', 'SUMO2', 'SUMO3', 'SUMO4', 'UFM1', 'URM1'
    ],

    # Mitochondrial proteins
    'mitochondrial_complex_I': [
        'NDUFA1', 'NDUFA2', 'NDUFA3', 'NDUFA4', 'NDUFA5', 'NDUFA6', 'NDUFA7',
        'NDUFA8', 'NDUFA9', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13',
        'NDUFB1', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7',
        'NDUFB8', 'NDUFB9', 'NDUFB10', 'NDUFB11', 'NDUFS1', 'NDUFS2', 'NDUFS3',
        'NDUFS4', 'NDUFS5', 'NDUFS6', 'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2'
    ],

    'mitochondrial_complex_II': [
        'SDHA', 'SDHB', 'SDHC', 'SDHD'
    ],

    'mitochondrial_complex_III': [
        'UQCRC1', 'UQCRC2', 'UQCRFS1', 'UQCRB', 'UQCRQ', 'UQCRH'
    ],

    'mitochondrial_complex_IV': [
        'COX4I1', 'COX4I2', 'COX5A', 'COX5B', 'COX6A1', 'COX6A2',
        'COX6B1', 'COX6C', 'COX7A1', 'COX7A2', 'COX7B', 'COX7C'
    ],

    'mitochondrial_complex_V': [
        'ATP5A1', 'ATP5B', 'ATP5C1', 'ATP5D', 'ATP5E', 'ATP5F1',
        'ATP5G1', 'ATP5G2', 'ATP5G3', 'ATP5H', 'ATP5I', 'ATP5J',
        'ATP5J2', 'ATP5L', 'ATP5O'
    ],

    # V-ATPase subunits (24 total)
    'vatpase_V0_domain': [
        'ATP6V0A1', 'ATP6V0A2', 'ATP6V0A4',
        'ATP6V0B', 'ATP6V0C', 'ATP6V0D1', 'ATP6V0D2',
        'ATP6V0E1', 'ATP6V0E2', 'ATP6AP1', 'ATP6AP2'
    ],

    'vatpase_V1_domain': [
        'ATP6V1A', 'ATP6V1B1', 'ATP6V1B2',
        'ATP6V1C1', 'ATP6V1C2', 'ATP6V1D',
        'ATP6V1E1', 'ATP6V1E2', 'ATP6V1F',
        'ATP6V1G1', 'ATP6V1G2', 'ATP6V1G3', 'ATP6V1H'
    ],

    # Autophagy proteins
    'autophagy_core': [
        'ATG3', 'ATG4A', 'ATG4B', 'ATG5', 'ATG7', 'ATG9A', 'ATG12',
        'ATG13', 'ATG14', 'ATG16L1', 'BECN1', 'ULK1', 'ULK2'
    ],

    'autophagy_receptors': [
        'SQSTM1', 'NBR1', 'OPTN', 'TAX1BP1', 'CALCOCO2', 'TOLLIP'
    ],

    'autophagy_LC3_family': [
        'MAP1LC3A', 'MAP1LC3B', 'MAP1LC3B2', 'MAP1LC3C',
        'GABARAP', 'GABARAPL1', 'GABARAPL2'
    ],

    # Lysosomal proteins
    'lysosomal_markers': [
        'LAMP1', 'LAMP2', 'LAMP3', 'CTSD', 'CTSL', 'CTSB', 'CTSZ',
        'LGMN', 'HEXA', 'HEXB', 'GBA', 'GLA', 'GAA', 'ARSA'
    ],

    # Endosomal proteins
    'endosomal_markers': [
        'EEA1', 'RAB5A', 'RAB5B', 'RAB5C', 'RAB7A', 'RAB7B',
        'RAB9A', 'RAB11A', 'RAB11B', 'VPS35', 'VPS26A', 'VPS26B', 'VPS29'
    ],

    # Retromer complex
    'retromer_complex': [
        'VPS35', 'VPS26A', 'VPS26B', 'VPS29',
        'SNX1', 'SNX2', 'SNX3', 'SNX5', 'SNX6', 'SNX27'
    ],

    # Heat shock proteins
    'heat_shock_proteins': [
        'HSPA1A', 'HSPA1B', 'HSPA2', 'HSPA4', 'HSPA5', 'HSPA6', 'HSPA8', 'HSPA9',
        'HSPB1', 'HSPB2', 'HSPB3', 'HSPB6', 'HSPB7', 'HSPB8',
        'HSP90AA1', 'HSP90AB1', 'HSP90B1', 'HSPD1', 'HSPE1'
    ],

    # Temporal analysis sets
    'temporal_early_response': [
        'HSPA1A', 'HSPA5', 'HSPB1', 'HSP90AA1', 'HSP90AB1'
    ],

    'temporal_autophagy_early': [
        'BECN1', 'ATG5', 'ATG7', 'ATG12'
    ],

    'temporal_autophagy_late': [
        'SQSTM1', 'NBR1', 'OPTN', 'TAX1BP1'
    ]
}

# Print summary of embedded protein sets
print("Embedded Protein Sets:")
print("="*50)
for set_name, proteins in protein_sets.items():
    print(f"{set_name:30} {len(proteins):3} proteins")

total_unique = len(set(p for proteins in protein_sets.values() for p in proteins))
print("="*50)
print(f"Total unique proteins defined: {total_unique}")
```

## Cell 4: Load Proteomics Data

```python
# Load the proteomics dataset
adata = sc.read_h5ad(data_path)

print(f"Dataset shape: {adata.shape}")
print(f"Number of cells/samples: {adata.n_obs}")
print(f"Number of proteins: {adata.n_vars}")
print("\nMetadata columns:")
print(adata.obs.columns.tolist())
print("\nProtein annotation columns:")
print(adata.var.columns.tolist())
```

## Cell 5: Explore Key Variables for DGE

```python
# Check tau status distribution
print("Tau Status Distribution:")
tau_column = 'TauStatus' if 'TauStatus' in adata.obs.columns else 'tau_status'
print(adata.obs[tau_column].value_counts())
print("\n" + "="*50)

# Check MC1 score distribution
mc1_column = 'MC1' if 'MC1' in adata.obs.columns else 'mc1_score'
if mc1_column in adata.obs.columns:
    print(f"\nMC1 Score Statistics:")
    print(adata.obs[mc1_column].describe())

# Check pseudotime distribution
pseudotime_column = 'Pseudotime' if 'Pseudotime' in adata.obs.columns else 'pseudotime'
if pseudotime_column in adata.obs.columns:
    print(f"\nPseudotime Statistics:")
    print(adata.obs[pseudotime_column].describe())
```

## Cell 6: Prepare Design Matrix for PertPy

```python
# Standardize column names for consistency
column_mapping = {
    'TauStatus': 'tau_status',
    'MC1': 'mc1_score',
    'Pseudotime': 'pseudotime',
    'Age': 'age_at_death'
}

for old_name, new_name in column_mapping.items():
    if old_name in adata.obs.columns and new_name not in adata.obs.columns:
        adata.obs[new_name] = adata.obs[old_name]

# Ensure tau_status is categorical
if 'tau_status' in adata.obs.columns:
    adata.obs['tau_status'] = pd.Categorical(adata.obs['tau_status'])
    print(f"Tau status categories: {adata.obs['tau_status'].cat.categories.tolist()}")

# Create binary tau variable for cleaner analysis
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)
print(f"\nTau positive samples: {adata.obs['tau_positive'].sum()}")
print(f"Tau negative samples: {(adata.obs['tau_positive'] == 0).sum()}")
```

## Cell 7: Data Quality Checks

```python
# Check for missing values
print("Missing values in expression matrix:")
missing_expr = np.isnan(adata.X).sum()
print(f"Total missing: {missing_expr} ({missing_expr / adata.X.size * 100:.2f}%)")

# Check for zero inflation
zero_count = (adata.X == 0).sum()
print(f"\nZero values: {zero_count} ({zero_count / adata.X.size * 100:.2f}%)")

# Check expression value range (should be log2 transformed)
print(f"\nExpression value range:")
print(f"Min: {np.nanmin(adata.X):.3f}")
print(f"Max: {np.nanmax(adata.X):.3f}")
print(f"Mean: {np.nanmean(adata.X):.3f}")
print(f"Median: {np.nanmedian(adata.X):.3f}")

# Check if data appears to be log-transformed
if np.nanmax(adata.X) < 50:  # Typical for log-transformed data
    print("\n✓ Data appears to be log-transformed")
else:
    print("\n⚠ Data may not be log-transformed")
```

## Cell 8: Prepare Protein Annotations

```python
# Check protein name format
print("Sample protein names:")
print(adata.var.index[:5].tolist())

# If there's a gene name column, use it
if 'GeneName' in adata.var.columns:
    print("\nUsing GeneName for protein identification")
    adata.var['protein_name'] = adata.var['GeneName']
elif 'gene_name' in adata.var.columns:
    adata.var['protein_name'] = adata.var['gene_name']
else:
    # Use index as protein name
    adata.var['protein_name'] = adata.var.index

# Create a clean protein list
protein_list = adata.var['protein_name'].tolist()
print(f"\nTotal proteins available: {len(protein_list)}")
print(f"Sample proteins: {protein_list[:5]}")
```

## Cell 9: Check Protein Set Availability

```python
# Check availability of protein sets
protein_availability = {}

print("Protein Set Availability in Dataset:")
print("="*70)
print(f"{'Protein Set':30} {'Found':>10} {'Total':>10} {'Coverage':>10}")
print("-"*70)

for set_name, proteins in protein_sets.items():
    available = [p for p in proteins if p in protein_list]
    protein_availability[set_name] = available
    coverage = len(available) / len(proteins) * 100
    print(f"{set_name:30} {len(available):>10} {len(proteins):>10} {coverage:>9.1f}%")

    if coverage < 50 and len(proteins) > 5:
        missing = set(proteins) - set(available)
        print(f"  ⚠ Low coverage - Missing: {list(missing)[:3]}...")

print("="*70)
```

## Cell 10: Create PertPy-Compatible Data Structure

```python
# Ensure data is in dense format for PertPy
if hasattr(adata.X, 'toarray'):
    print("Converting sparse matrix to dense...")
    adata.X = adata.X.toarray()

# Create a copy for PertPy analysis
adata_pertpy = adata.copy()

# Add raw counts if not present (PertPy/DESeq2 expects counts)
if 'raw_counts' not in adata_pertpy.layers:
    # If data is log2 transformed, reverse it for DESeq2
    print("Creating pseudo-raw counts from log2 data...")
    # Convert log2 to linear scale and multiply by scaling factor
    adata_pertpy.layers['counts'] = np.power(2, adata_pertpy.X) * 1000
    # Round to integers for DESeq2
    adata_pertpy.layers['counts'] = np.round(adata_pertpy.layers['counts']).astype(int)
else:
    adata_pertpy.layers['counts'] = adata_pertpy.layers['raw_counts']

print(f"PertPy-ready data shape: {adata_pertpy.shape}")
print(f"Counts layer added: {adata_pertpy.layers['counts'].shape}")
```

## Cell 11: Helper Functions for Protein Set Analysis

```python
def get_protein_subset(adata, protein_list, set_name="protein_set"):
    """
    Create an AnnData subset containing only specified proteins.

    Parameters:
    -----------
    adata : AnnData
        Full dataset
    protein_list : list
        List of protein names to include
    set_name : str
        Name for this protein set (for reporting)

    Returns:
    --------
    AnnData subset with only specified proteins
    """
    # Find available proteins
    available = [p for p in protein_list if p in adata.var_names or p in adata.var.get('protein_name', [])]

    if len(available) == 0:
        print(f"⚠ No proteins from {set_name} found in dataset")
        return None

    print(f"Creating subset for {set_name}: {len(available)}/{len(protein_list)} proteins found")

    # Get indices
    if 'protein_name' in adata.var.columns:
        indices = [i for i, p in enumerate(adata.var['protein_name']) if p in available]
    else:
        indices = [i for i, p in enumerate(adata.var_names) if p in available]

    # Create subset
    adata_subset = adata[:, indices].copy()

    return adata_subset

# Test the function
test_subset = get_protein_subset(adata_pertpy, protein_sets['autophagy_receptors'], "Autophagy Receptors")
if test_subset:
    print(f"Test subset created: {test_subset.shape}")
```

## Cell 12: Save Prepared Data (Optional for Local Use)

```python
# Only save if not in Colab (Colab users will work with in-memory data)
if not IN_COLAB:
    output_file = 'prepared_for_pertpy.h5ad'
    adata_pertpy.write_h5ad(output_file)
    print(f"Prepared data saved to: {output_file}")
else:
    print("Running in Colab - data kept in memory")

# Create summary statistics
summary = {
    'n_samples': adata_pertpy.n_obs,
    'n_proteins': adata_pertpy.n_vars,
    'tau_positive': int(adata_pertpy.obs['tau_positive'].sum()),
    'tau_negative': int((adata_pertpy.obs['tau_positive'] == 0).sum()),
    'has_pseudotime': 'pseudotime' in adata_pertpy.obs.columns,
    'has_mc1': 'mc1_score' in adata_pertpy.obs.columns,
    'data_format': 'log2_transformed',
    'counts_layer': 'counts' in adata_pertpy.layers
}

print("\n" + "="*50)
print("Data Preparation Summary:")
print("="*50)
for key, value in summary.items():
    print(f"{key}: {value}")
```

## Cell 13: Visualize Data Distribution

```python
# Create visualization of data distribution
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Expression distribution
axes[0, 0].hist(adata_pertpy.X.flatten(), bins=50, alpha=0.7, color='blue')
axes[0, 0].set_xlabel('Log2 Expression')
axes[0, 0].set_ylabel('Frequency')
axes[0, 0].set_title('Global Expression Distribution')

# Plot 2: Sample-wise mean expression
sample_means = np.nanmean(adata_pertpy.X, axis=1)
tau_colors = ['red' if x == 1 else 'blue' for x in adata_pertpy.obs['tau_positive']]
axes[0, 1].scatter(range(len(sample_means)), sample_means, c=tau_colors, alpha=0.6)
axes[0, 1].set_xlabel('Sample Index')
axes[0, 1].set_ylabel('Mean Log2 Expression')
axes[0, 1].set_title('Sample-wise Mean Expression (Red=Tau+, Blue=Tau-)')

# Plot 3: Protein coverage
protein_coverage = np.sum(~np.isnan(adata_pertpy.X), axis=0) / adata_pertpy.n_obs
axes[1, 0].hist(protein_coverage, bins=50, alpha=0.7, color='green')
axes[1, 0].set_xlabel('Coverage (fraction of samples)')
axes[1, 0].set_ylabel('Number of Proteins')
axes[1, 0].set_title('Protein Coverage Distribution')

# Plot 4: Tau status vs pseudotime (if available)
if 'pseudotime' in adata_pertpy.obs.columns:
    for tau_val, color, label in [(0, 'blue', 'Tau-'), (1, 'red', 'Tau+')]:
        mask = adata_pertpy.obs['tau_positive'] == tau_val
        axes[1, 1].scatter(adata_pertpy.obs.loc[mask, 'pseudotime'],
                          adata_pertpy.obs.loc[mask, 'mc1_score'] if 'mc1_score' in adata_pertpy.obs.columns else np.zeros(mask.sum()),
                          c=color, label=label, alpha=0.6)
    axes[1, 1].set_xlabel('Pseudotime')
    axes[1, 1].set_ylabel('MC1 Score')
    axes[1, 1].set_title('Disease Progression Markers')
    axes[1, 1].legend()
else:
    axes[1, 1].text(0.5, 0.5, 'Pseudotime not available',
                    ha='center', va='center', transform=axes[1, 1].transAxes)

plt.tight_layout()
plt.show()

print("\n✓ Data preparation complete!")
print("Ready for PertPy differential expression analysis")
```

## Summary

This **Colab-optimized** notebook has prepared the proteomics data for PertPy analysis:

1. ✓ **All protein sets embedded directly** - no external files needed
2. ✓ Loaded pool_processed_v2.h5ad dataset
3. ✓ Standardized metadata columns
4. ✓ Created design matrix variables (tau_status, pseudotime, mc1_score)
5. ✓ Generated count matrix for DESeq2 analysis
6. ✓ Defined comprehensive protein sets (132 UPS proteins + others)
7. ✓ Created helper functions for protein subsetting

The data is now ready for differential expression analysis using PertPy's PyDESeq2 implementation.

### Next Steps
Use the embedded protein sets and helper functions in downstream analysis notebooks:
```python
# Example usage in analysis notebooks:
adata_ups = get_protein_subset(adata_pertpy, protein_sets['ups_proteins'], "UPS Proteins")
```

---

*This notebook is fully self-contained and optimized for Google Colab - just upload pool_processed_v2.h5ad and run!*