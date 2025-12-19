# MSc Biology Analysis - Proteostasis in Neurodegeneration

## Overview
This analysis was conducted as part of an MSc Biology thesis investigating proteostasis failure in Alzheimer's disease neurons. As a biology student with limited programming experience, I've documented each step carefully to evaluate specific biological claims from recent literature.

## Author Background
- **Level**: MSc Biology student
- **Programming experience**: Basic Python learned during this project
- **Focus**: Biological interpretation rather than computational methods
- **Goal**: Test published claims about protein quality control failure

## Research Questions
1. Do proteostasis mechanisms fail sequentially or simultaneously?
2. Is autophagy specifically impaired while the proteasome remains functional?
3. What is the relationship between tau pathology and organellar dysfunction?

## Directory Structure
```
msc_biology_analysis/
├── notebooks/           # Jupyter notebooks with analysis
│   ├── 01_sequential_failure_analysis.ipynb
│   └── 02_mitochondrial_dysfunction_analysis.ipynb
├── results/            # Analysis outputs and summary files
├── figures/           # Generated plots and visualizations
├── data/              # Data directory (links to main dataset)
└── run_analysis.py    # Student-friendly script to run everything
```

## Data Location
The proteomics dataset is located at: `../data/pool_processed_v2.h5ad`
- **Format**: AnnData object (scanpy compatible)
- **Size**: 44 samples × 5,853 proteins
- **Groups**: 22 tau-positive, 22 tau-negative neurons

## Dataset
- **Source**: Mass spectrometry proteomics of Alzheimer's neurons
- **Format**: AnnData object (similar to single-cell RNA-seq)
- **Size**: 44 samples × 5,853 proteins
- **Groups**: 22 tau-positive, 22 tau-negative neurons

## Analysis Approach
As a biology student, I focused on:
1. Understanding the biological questions first
2. Learning just enough code to test each claim
3. Interpreting results in biological context
4. Documenting everything clearly for my thesis

## Key Findings Preview
- Proteostasis fails sequentially, not simultaneously
- SQSTM1/p62 shows dramatic 10.7-fold upregulation
- Autophagy is specifically impaired while proteasome remains stable
- Critical threshold at MC1 = 2.831 marks system collapse

## Notebooks
1. **Sequential Failure Analysis**: Testing if protein quality control arms fail in order
2. **Mitochondrial Dysfunction Analysis**: Examining organellar stress and mitophagy

Each notebook includes:
- Biological background for context
- Step-by-step code with explanations
- Interpretation of results
- Thesis-ready conclusions

---
*Analysis conducted by: MSc Biology Student*
*Supervisor: [Thesis Advisor]*
*Date: September 2024*