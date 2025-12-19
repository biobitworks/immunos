# Data Documentation

## Overview
Centralized data management for the proteomics analysis framework. Contains the primary dataset and supplementary resources for neurodegeneration research.

## Primary Dataset

### pool_processed_v2.h5ad
**Main proteomics dataset for all analyses**

#### Specifications
- **Format**: AnnData (h5ad) - Hierarchical Data Format 5
- **Dimensions**: 44 samples × 5,853 proteins
- **Source**: Alzheimer's disease brain tissue proteomics
- **Processing**: Log2-transformed expression values
- **Size**: ~15 MB

#### Data Structure
```python
# AnnData object structure
adata.X          # Expression matrix (44 × 5,853)
adata.obs        # Sample metadata (44 rows)
adata.var        # Protein information (5,853 rows)
adata.uns        # Unstructured metadata
```

#### Sample Metadata (adata.obs)
| Column | Type | Description |
|--------|------|-------------|
| PatientID | str | Unique patient identifier |
| age | float | Patient age at death |
| PMI | float | Post-mortem interval (hours) |
| tau_positive | bool | Tau pathology status |
| MC1 | float | MC1 antibody reactivity score |
| pseudotime | float | Disease progression metric |
| sample_id | str | Unique sample identifier |

#### Protein Information (adata.var)
| Column | Type | Description |
|--------|------|-------------|
| protein_name | str | Gene symbol/protein name |
| uniprot_id | str | UniProt accession number |
| description | str | Protein description |
| mean_expression | float | Mean expression across samples |
| variance | float | Expression variance |
| n_cells | int | Number of samples with detection |

## Supplementary Data

### supplementary_data/
Additional datasets and reference files:

#### protein_annotations.csv
- Detailed protein annotations from UniProt
- GO terms, pathways, disease associations
- Domain information and modifications

#### clinical_metadata.csv
- Extended clinical information
- Cognitive scores, neuropathology staging
- Medication history

#### reference_proteome.fasta
- Human reference proteome sequences
- Used for protein identification validation
- UniProt release information

## Data Access

### Loading the Main Dataset
```python
import scanpy as sc

# Load the AnnData object
adata = sc.read_h5ad('pool_processed_v2.h5ad')

# Quick data inspection
print(f"Samples: {adata.n_obs}")
print(f"Proteins: {adata.n_vars}")
print(f"Tau-positive samples: {adata.obs['tau_positive'].sum()}")
```

### Accessing Expression Data
```python
# Get expression for specific protein
protein = 'SQSTM1'
if protein in adata.var_names:
    expression = adata[:, protein].X.flatten()

# Get expression matrix as DataFrame
import pandas as pd
expression_df = pd.DataFrame(
    adata.X,
    index=adata.obs_names,
    columns=adata.var_names
)
```

### Filtering and Subsetting
```python
# Filter by tau status
tau_positive = adata[adata.obs['tau_positive'] == True]
tau_negative = adata[adata.obs['tau_positive'] == False]

# Filter proteins by expression
highly_expressed = adata[:, adata.var['mean_expression'] > 5]

# Get specific protein groups
ups_proteins = ['UBE2D3', 'UBE2N', 'UBE2K', 'PSMA1', 'PSMB5']
ups_subset = adata[:, [p for p in ups_proteins if p in adata.var_names]]
```

## Data Quality

### Quality Metrics
- **Missing Values**: <0.1% (handled during preprocessing)
- **Normalization**: Log2-transformed, batch-corrected
- **Quality Control**: Samples passing QC thresholds
- **Protein Coverage**: >95% coverage across samples

### Validation Performed
- ✅ Sample identity verification
- ✅ Batch effect assessment
- ✅ Outlier detection and removal
- ✅ Protein annotation validation
- ✅ Expression distribution checks

## Data Usage Guidelines

### Best Practices
1. **Always check protein existence** before accessing:
   ```python
   if protein_name in adata.var_names:
       # Safe to access
   ```

2. **Use consistent identifiers** (gene symbols vs UniProt IDs)

3. **Account for log-transformation** in statistical analyses

4. **Consider covariates** (age, PMI) in differential expression

5. **Apply multiple testing correction** for proteome-wide analyses

### Common Issues & Solutions

#### Missing Proteins
```python
# Find alternative names
def find_protein_alternatives(adata, protein):
    alternatives = [p for p in adata.var_names
                   if protein.lower() in p.lower()]
    return alternatives
```

#### Memory Management
```python
# Process in chunks for large operations
chunk_size = 1000
for i in range(0, adata.n_vars, chunk_size):
    chunk = adata[:, i:min(i+chunk_size, adata.n_vars)]
    # Process chunk
```

## Data Updates

### Version History
- **v2.0** (2024-09): Current version with complete annotations
- **v1.0** (2024-08): Initial processed dataset

### Update Process
1. Raw data validation
2. Quality control filtering
3. Normalization and batch correction
4. Annotation update from UniProt
5. Format conversion to h5ad

## References

### Data Sources
- **Proteomics Data**: [Original Publication Reference]
- **UniProt Annotations**: Release 2024_03
- **Gene Ontology**: GO database release 2024-09
- **Clinical Data**: [Clinical Study Reference]

## Navigation

- [[../INDEX|Back to Index]]
- [[../01_research_analysis/README|Research Analysis]]
- [[../04_documentation/README|Documentation]]