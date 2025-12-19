# PertPy-based Differential Gene Expression Analysis

## Overview
This directory contains comprehensive differential expression analysis using PertPy's PyDESeq2 implementation for all claims from the proteomics dataset.

## Directory Structure
```
pertpy_dge_analysis/
├── 01_data_preparation/     # Data loading and preparation
├── 02_group1_mitochondrial/ # Group 1 claims (8 notebooks)
├── 03_group2_proteostasis/  # Group 2 claims (8 notebooks)
├── 04_comprehensive_analysis/ # Summary analyses
├── 05_statistical_reports/  # Results CSV files
└── 06_utilities/           # Helper functions
```

## Installation
```bash
pip install -r requirements.txt
```

## Key Features
- **Modern DGE Methods**: Uses PyDESeq2 (Python implementation of DESeq2)
- **Comprehensive Coverage**: All 16 claims analyzed systematically
- **Professional Visualizations**: Volcano plots, heatmaps, fold change plots
- **Robust Statistics**: FDR correction, effect sizes, covariate adjustment
- **Reproducible**: All analyses in Jupyter notebooks

## Claims Analyzed

### Group 1 - Mitochondrial Dysregulation
1. **UPS Proteins**: No significant alterations (132 proteins tested)
2. **SQSTM1**: Massive upregulation analysis
3. **Temporal Dynamics**: Disease progression patterns
4. **Mitochondrial Proteins**: Complex I-V analysis
5. **Cristae Organization**: OPA1, MIC60, MIC19
6. **Sliding Window**: Temporal expression waves
7. **Mitophagy Receptors**: PINK1, Parkin, receptors
8. **Parkin-Independent**: Alternative pathways

### Group 2 - Proteostasis Failure
1. **V-ATPase Subunits**: Differential expression
2. **ATP6V0A1**: Specific validation
3. **Organellar Markers**: Compartment analysis
4. **Retromer Complex**: VPS proteins
5. **Autophagy vs UPS**: Comparative disruption
6. **Endolysosomal**: Pathway-specific changes
7. **Temporal Cascade**: Sequential failure
8. **Rab GTPases**: Trafficking analysis

## Usage

### 1. Prepare Data
```python
# Run data preparation notebook
jupyter notebook 01_data_preparation/prepare_for_pertpy.ipynb
```

### 2. Run Individual Claims
```python
# Example: Claim 1 - UPS Proteins
jupyter notebook 02_group1_mitochondrial/claim1_ups_proteins.ipynb
```

### 3. Generate Summary
```python
# Run comprehensive analysis
jupyter notebook 04_comprehensive_analysis/all_claims_summary.ipynb
```

## PertPy Advantages
- **DESeq2 Power**: Gold standard differential expression
- **Handles Covariates**: Pseudotime, MC1 score adjustment
- **Interaction Terms**: Complex experimental designs
- **Built-in Visualization**: Professional plots included
- **AnnData Compatible**: Works with h5ad format

## Output Files
- Individual claim results: `05_statistical_reports/claim*_results.csv`
- Significant proteins: `05_statistical_reports/significant_proteins_*.csv`
- Volcano plots: Saved in each notebook directory
- Summary statistics: `05_statistical_reports/all_claims_summary.csv`

## Key Findings Preview
- **Claim 1 (UPS)**: Comprehensive 132-protein analysis
- **Claim 2 (SQSTM1)**: Validated with exact fold change
- **All Claims**: Systematic evaluation with consistent methodology

## Requirements
- Python 3.8+
- PertPy >= 0.5.0
- PyDESeq2 >= 0.4.0
- Scanpy >= 1.9.0
- See requirements.txt for complete list

## Contact
For questions about the analysis or PertPy implementation, refer to the individual notebooks which contain detailed documentation and interpretation.

---
*Analysis Date: December 2024*
*Dataset: pool_processed_v2.h5ad*
*Method: PertPy PyDESeq2*