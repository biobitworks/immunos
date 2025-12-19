# Research Analysis Notebooks

## Overview
Core research deliverables implementing comprehensive statistical analyses for evaluating biological claims about proteostasis failure and mitochondrial dysfunction in neurodegeneration.

## Available Formats

Each notebook is available in two formats:
- **`.ipynb`** - Jupyter notebook format for running analyses
- **`.md`** - Markdown format for viewing in Obsidian

## Notebooks

### 1. Sequential Failure of Proteostasis Mechanisms
- **Jupyter**: [[Sequential_Failure_of_Proteostasis_Mechanisms_Notebook.ipynb]]
- **Markdown**: [[Sequential_Failure_of_Proteostasis_Mechanisms_Notebook.md]]
- **Purpose**: Evaluates 8 claims about proteostasis system failure
- **Key Methods**: Covariate-controlled DE, segmented regression, temporal analysis

### 2. Late-Stage Mitochondrial Dysregulation and Mitophagy Failure
- **Jupyter**: [[Late_Stage_Mitochondrial_Dysregulation_and_Mitophagy_Failure_Notebook.ipynb]]
- **Markdown**: [[Late_Stage_Mitochondrial_Dysregulation_and_Mitophagy_Failure_Notebook.md]]
- **Purpose**: Evaluates 8 claims about mitochondrial dysfunction
- **Key Methods**: UPS analysis, sliding window, biphasic detection

### 3. Master Analysis Framework
- **Jupyter**: [[Master_Analysis_Notebook.ipynb]]
- **Markdown**: [[Master_Analysis_Notebook.md]]
- **Purpose**: Comprehensive pipeline for all 16 biological claims
- **Coverage**: Both finding groups with unified analysis

## Quick Access Links

### For Obsidian Users
- [[Sequential_Failure_of_Proteostasis_Mechanisms_Notebook|📓 Proteostasis Analysis (Markdown)]]
- [[Late_Stage_Mitochondrial_Dysregulation_and_Mitophagy_Failure_Notebook|📓 Mitochondrial Analysis (Markdown)]]
- [[Master_Analysis_Notebook|📓 Master Analysis (Markdown)]]

### For Jupyter Users
```bash
# Navigate to notebooks directory
cd /Users/byron/project_plan/01_research_analysis/notebooks/

# Launch Jupyter
jupyter notebook

# Select and run desired notebook
```

## Notebook Contents

### Biological Claims Evaluated

#### Proteostasis Failure (8 claims)
1. V-ATPase and proton pump disruption
2. ATP6V0A1 upregulation patterns
3. Loss of organellar identity markers
4. Retromer complex dysfunction
5. SOS response activation
6. Segmented progression breakpoints
7. Temporal failure ordering
8. Proteostasis network collapse

#### Mitochondrial Dysfunction (8 claims)
1. UPS protein differential expression
2. SQSTM1/p62 upregulation validation
3. BECN1-SQSTM1 correlations
4. BECN1 reduction patterns
5. Mitophagy pathway impairment
6. Sliding window temporal patterns
7. Biphasic expression profiles
8. Tau-mitochondrial correlations

## Statistical Methods

All notebooks implement:
- **Multiple Testing Correction**: Benjamini-Hochberg FDR, Bonferroni
- **Effect Size Calculation**: Cohen's d, fold changes
- **Covariate Control**: Age, PMI, PatientID adjustment
- **Temporal Analysis**: Sliding windows, segmented regression
- **Pattern Detection**: Biphasic patterns, breakpoints
- **Visualization**: Publication-quality plots

## Data Requirements

- **Input**: `pool_processed_v2.h5ad` (44 samples × 5,853 proteins)
- **Format**: AnnData with log2-transformed expression values
- **Location**: `/Users/byron/project_plan/03_data/`

## Navigation

- [[../../INDEX|Back to Main Index]]
- [[../README|Research Analysis Section]]
- [[../../03_data/README|Data Documentation]]

---

*Last Updated: 2024-09-28 | Markdown versions available for Obsidian viewing*