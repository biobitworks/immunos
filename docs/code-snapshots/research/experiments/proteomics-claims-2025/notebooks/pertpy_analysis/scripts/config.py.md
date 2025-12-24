---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/notebooks/pertpy_analysis/scripts/config.py
relative: research/experiments/proteomics-claims-2025/notebooks/pertpy_analysis/scripts/config.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Central Configuration for Proteomics Analysis Project
======================================================
Single source of truth for all data paths and analysis parameters.

Author: Bioinformatics Finding Group Evaluation Framework
Date: 2024
Data: pool_processed_v2.h5ad - Alzheimer's disease proteomics
"""

import os
from pathlib import Path

# ============================================================================
# DATA CONFIGURATION
# ============================================================================

# Primary data file - all analyses must use this
DATA_PATH = '/Users/byron/project_plan/data/pool_processed_v2.h5ad'

# Verify data exists
if not os.path.exists(DATA_PATH):
    raise FileNotFoundError(f"Data file not found: {DATA_PATH}")

# Data specifications (from actual data)
DATA_SPECS = {
    'n_samples': 44,
    'n_proteins': 5853,
    'n_tau_positive': 22,
    'n_tau_negative': 22,
    'tau_column': 'TauStatus',  # Correct column name (not 'tau_status')
    'tau_positive_value': 'positive',  # Correct value (not 'tau+')
    'tau_negative_value': 'negative',  # Correct value (not 'tau-')
    'mc1_column': 'MC1',
    'pseudotime_column': 'pseudotime',
    'patient_id_column': 'PatientID',
    'age_column': 'Age at death',
    'pmi_column': 'PMI hours'
}

# ============================================================================
# PROJECT STRUCTURE
# ============================================================================

# Project root
PROJECT_ROOT = Path('/Users/byron/project_plan')

# Key directories
DIRS = {
    'data': PROJECT_ROOT / 'data',
    'results': PROJECT_ROOT / 'results',
    'notebooks': PROJECT_ROOT / 'notebooks',
    'figures': PROJECT_ROOT / 'results' / 'figures',
    'reports': PROJECT_ROOT / 'results' / 'reports',
    'group1': PROJECT_ROOT / '01_research_analysis' / 'group1_mitochondrial',
    'group2': PROJECT_ROOT / '01_research_analysis' / 'group2_proteostasis',
    'ai_automation': PROJECT_ROOT / '01_research_analysis' / 'ai_automation'
}

# Create directories if they don't exist
for dir_path in DIRS.values():
    dir_path.mkdir(parents=True, exist_ok=True)

# ============================================================================
# ANALYSIS PARAMETERS
# ============================================================================

ANALYSIS_PARAMS = {
    # Statistical thresholds
    'p_value_threshold': 0.05,
    'fdr_method': 'benjamini_hochberg',
    'min_fold_change': 1.5,

    # Effect size thresholds
    'cohen_d_small': 0.2,
    'cohen_d_medium': 0.5,
    'cohen_d_large': 0.8,

    # Sample size requirements
    'min_samples_per_group': 10,

    # Visualization
    'figure_dpi': 300,
    'figure_format': 'png',

    # Multiple testing correction
    'apply_fdr': True,
    'fdr_alpha': 0.05
}

# ============================================================================
# PROTEIN GROUPS OF INTEREST
# ============================================================================

PROTEIN_GROUPS = {
    'proteasome': {
        '20S_alpha': ['PSMA1', 'PSMA2', 'PSMA3', 'PSMA4', 'PSMA5', 'PSMA6', 'PSMA7'],
        '20S_beta': ['PSMB1', 'PSMB2', 'PSMB3', 'PSMB4', 'PSMB5', 'PSMB6', 'PSMB7'],
        '19S_regulatory': ['PSMC1', 'PSMC2', 'PSMC3', 'PSMC4', 'PSMC5', 'PSMC6',
                          'PSMD1', 'PSMD2', 'PSMD3', 'PSMD4']
    },

    'autophagy': {
        'receptors': ['SQSTM1', 'NBR1', 'OPTN'],
        'core_machinery': ['MAP1LC3B', 'GABARAP', 'GABARAPL1', 'GABARAPL2'],
        'initiation': ['BECN1', 'ATG5', 'ATG7', 'ATG12', 'ATG16L1'],
        'regulators': ['MTOR', 'ULK1', 'TFEB']
    },

    'vatpase': {
        'V0_domain': ['ATP6V0A1', 'ATP6V0A2', 'ATP6V0B', 'ATP6V0C', 'ATP6V0D1', 'ATP6V0E1'],
        'V1_domain': ['ATP6V1A', 'ATP6V1B1', 'ATP6V1B2', 'ATP6V1C1', 'ATP6V1D', 'ATP6V1E1']
    },

    'mitochondrial': {
        'import': ['TOMM20', 'TOMM40', 'TIMM23'],
        'dynamics': ['MFN1', 'MFN2', 'DRP1', 'OPA1', 'FIS1'],
        'respiration': ['CYCS', 'COX4I1', 'ATP5A1', 'NDUFS1', 'SDHA'],
        'quality_control': ['PINK1', 'PRKN', 'VDAC1', 'VDAC2']
    },

    'ups': {
        'ubiquitin': ['UBB', 'UBC', 'UBA52', 'RPS27A'],
        'E1_enzymes': ['UBA1', 'UBA2', 'UBA3', 'UBA6', 'UBA7'],
        'E2_enzymes': ['UBE2A', 'UBE2B', 'UBE2C', 'UBE2D1', 'UBE2D2', 'UBE2D3'],
        'E3_ligases': ['MDM2', 'TRIM21', 'RNF2', 'PARK2', 'FBXW7'],
        'DUBs': ['USP7', 'USP14', 'UCHL1', 'UCHL3', 'UCHL5']
    }
}

# ============================================================================
# KNOWN ISSUES & DISCREPANCIES
# ============================================================================

KNOWN_ISSUES = {
    'SQSTM1_fold_change': {
        'claimed': 10.7,
        'observed': 1.3,
        'note': 'Significant discrepancy - needs investigation',
        'possible_causes': [
            'Different normalization method',
            'Subset analysis in paper',
            'Log transformation differences',
            'Statistical model differences'
        ]
    },

    'semicolon_proteins': {
        'affected': 56,
        'percentage': 1.0,
        'note': 'Proteins with multiple UniProt IDs (isoforms)',
        'examples': ['UBB;UBC', 'MAP1LC3B;MAP1LC3B2']
    }
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def load_data():
    """Load the proteomics data with proper error handling."""
    import scanpy as sc

    try:
        adata = sc.read_h5ad(DATA_PATH)
        print(f"✓ Data loaded: {adata.n_obs} samples × {adata.n_vars} proteins")
        return adata
    except Exception as e:
        raise RuntimeError(f"Failed to load data: {e}")

def get_tau_groups(adata):
    """Get tau positive and negative sample indices."""
    tau_col = DATA_SPECS['tau_column']
    tau_pos_val = DATA_SPECS['tau_positive_value']
    tau_neg_val = DATA_SPECS['tau_negative_value']

    tau_pos = adata.obs[tau_col] == tau_pos_val
    tau_neg = adata.obs[tau_col] == tau_neg_val

    return tau_pos, tau_neg

def validate_protein(adata, protein_name):
    """Check if protein exists in dataset."""
    return protein_name in adata.var['GeneName'].values

# ============================================================================
# REPORTING CONFIGURATION
# ============================================================================

REPORT_CONFIG = {
    'title': 'Proteomics Analysis Report - pool_processed_v2.h5ad',
    'author': 'Bioinformatics Finding Group Evaluation Framework',
    'date_format': '%Y-%m-%d',
    'include_methods': True,
    'include_statistics': True,
    'include_visualizations': True,
    'output_formats': ['html', 'pdf', 'md']
}

# ============================================================================
# EXPORT CONFIGURATION
# ============================================================================

__all__ = [
    'DATA_PATH',
    'DATA_SPECS',
    'PROJECT_ROOT',
    'DIRS',
    'ANALYSIS_PARAMS',
    'PROTEIN_GROUPS',
    'KNOWN_ISSUES',
    'load_data',
    'get_tau_groups',
    'validate_protein',
    'REPORT_CONFIG'
]

if __name__ == "__main__":
    print("=" * 60)
    print("PROTEOMICS ANALYSIS CONFIGURATION")
    print("=" * 60)
    print(f"Data path: {DATA_PATH}")
    print(f"Data exists: {os.path.exists(DATA_PATH)}")
    print(f"Project root: {PROJECT_ROOT}")
    print(f"\nData specifications:")
    for key, value in DATA_SPECS.items():
        print(f"  {key}: {value}")
    print(f"\nKnown issues:")
    for issue, details in KNOWN_ISSUES.items():
        print(f"  {issue}: {details['note']}")
    print("=" * 60)
```
