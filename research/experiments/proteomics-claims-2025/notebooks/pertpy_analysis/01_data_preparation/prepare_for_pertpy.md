# Data Preparation for PertPy DGE Analysis

## Overview
This notebook prepares the pool_processed_v2.h5ad proteomics dataset for differential expression analysis using PertPy and PyDESeq2.

## Workflow

### 1. Import Required Packages
```python
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns
```

### 2. Load Proteomics Data
```python
adata = sc.read_h5ad('../../data/pool_processed_v2.h5ad')
```

**Dataset Information:**
- **Samples**: 44 neurons (mini-pools of 10 neurons each)
- **Proteins**: ~5000 measured proteins
- **Metadata**: tau_status, pseudotime, MC1 score, age at death
- **Format**: Log2 transformed expression values

### 3. Explore Key Variables

#### Tau Status Distribution
- Tau-positive neurons: Disease state
- Tau-negative neurons: Control state
- Binary classification for primary comparisons

#### MC1 Score
- Quantification of misfolded tau
- Continuous variable for covariate adjustment

#### Pseudotime
- Disease progression metric
- Used for temporal analyses

### 4. Prepare Design Matrix

Standardized column names for consistency:
```python
column_mapping = {
    'TauStatus': 'tau_status',
    'MC1': 'mc1_score',
    'Pseudotime': 'pseudotime',
    'Age': 'age_at_death'
}
```

Create binary tau variable:
```python
adata.obs['tau_positive'] = (adata.obs['tau_status'] == 'positive').astype(int)
```

### 5. Data Quality Checks

#### Missing Values
- Check for NaN in expression matrix
- Typical: <1% missing values

#### Zero Inflation
- Assess percentage of zero values
- Expected in single-cell/low-input proteomics

#### Expression Range
- Verify log2 transformation (values typically -5 to 15)
- Check for outliers

### 6. Prepare Protein Annotations

Map protein identifiers:
```python
if 'GeneName' in adata.var.columns:
    adata.var['protein_name'] = adata.var['GeneName']
```

Create clean protein list for searching:
- Standardized gene names
- Handle multiple isoforms
- Remove contaminants

### 7. Create PertPy-Compatible Structure

#### Convert to Dense Matrix
```python
if hasattr(adata.X, 'toarray'):
    adata.X = adata.X.toarray()
```

#### Add Count Layer for DESeq2
PyDESeq2 expects count data:
```python
# If data is log2 transformed, create pseudo-counts
adata.layers['counts'] = np.power(2, adata.X) * 1000
adata.layers['counts'] = np.round(adata.layers['counts']).astype(int)
```

### 8. Define Protein Sets for Analysis

Key protein categories:
```python
protein_sets = {
    'ups_proteins': [132 proteins],      # Ubiquitin-proteasome system
    'mitochondrial': [50+ proteins],     # Energy metabolism
    'autophagy': [20+ proteins],         # Protein degradation
    'vatpase': [24 subunits]            # Lysosomal acidification
}
```

### 9. Save Prepared Data

Output files:
- `prepared_for_pertpy.h5ad` - Main analysis file
- `protein_sets.json` - Protein group definitions
- `data_summary.csv` - Dataset statistics

### 10. Visualization

Generated plots:
1. **Expression Distribution**: Global histogram
2. **Sample Means**: By tau status
3. **Protein Coverage**: Detection rates
4. **Disease Markers**: Pseudotime vs MC1

---

## Key Parameters

### Data Specifications
- **Expression format**: Log2 transformed
- **Sample size**: 44 neurons
- **Tau groups**: 22 positive, 22 negative
- **Protein coverage**: >80% detected per sample

### Quality Metrics
- **Missing data**: <1%
- **Zero values**: <20% typical
- **Dynamic range**: ~20 orders of magnitude

### Design Considerations
- **Primary contrast**: Tau+ vs Tau-
- **Covariates**: Pseudotime, MC1 score
- **Batch effects**: Not detected

---

## Output Summary

### Files Created
1. `prepared_for_pertpy.h5ad`
   - Ready for PyDESeq2 analysis
   - Contains counts and log2 layers
   - Standardized metadata

2. `protein_sets.json`
   - Curated protein groups
   - Used across all analyses

3. `data_distribution_overview.png`
   - Quality control visualization
   - Sample distributions

---

## Usage in Downstream Analysis

### Loading Prepared Data
```python
adata = sc.read_h5ad('prepared_for_pertpy.h5ad')

# Access layers
log2_expr = adata.layers['log2']  # Original values
counts = adata.layers['counts']    # For DESeq2
```

### Using Protein Sets
```python
import json
with open('protein_sets.json', 'r') as f:
    protein_sets = json.load(f)

ups_proteins = protein_sets['ups_proteins']
```

---

## Troubleshooting

### Common Issues

1. **Memory errors with large datasets**
   - Solution: Process in chunks
   - Use sparse matrices where possible

2. **Protein name mismatches**
   - Check gene name variations
   - Use fuzzy matching functions

3. **PyDESeq2 convergence issues**
   - Verify count matrix is integer
   - Check for extreme outliers
   - Ensure sufficient replicates

---

## Next Steps

After data preparation:
1. Run individual claim analyses
2. Apply PyDESeq2 for differential expression
3. Generate volcano plots
4. Evaluate scientific claims

---

## Dependencies

- scanpy >= 1.9.0
- anndata >= 0.9.0
- pandas >= 1.5.0
- numpy >= 1.23.0
- matplotlib >= 3.5.0
- seaborn >= 0.12.0

---

*Data preparation is critical for reliable differential expression analysis. This notebook ensures compatibility with PertPy/PyDESeq2 while maintaining data integrity.*