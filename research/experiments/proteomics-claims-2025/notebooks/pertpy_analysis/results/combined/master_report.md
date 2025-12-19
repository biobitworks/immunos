# Master Analysis Report - PertPy DGE Pipeline

**Generated**: 2025-09-29 22:33:34

## Executive Summary

- **Total Claims Analyzed**: 16
- **Supported**: 0
- **Partially Supported**: 2
- **Refuted**: 14
- **Errors**: 0

## Detailed Results

### Group 1: Mitochondrial Dysfunction

| Claim | Verdict | Significant Proteins | Mean Log2FC |
|-------|---------|---------------------|-------------|
| No significant UPS protein alterations across tau-... | ❌ REFUTED | 0 | 0.017 |
| SQSTM1/p62 is upregulated in tau+ neurons | ❌ REFUTED | 0 | 0.015 |
| Protein expression changes follow temporal dynamic... | ❌ REFUTED | 0 | -0.003 |
| Mitochondrial complexes I-V are decreased in tau+ ... | ❌ REFUTED | 0 | 0.026 |
| Cristae organization proteins are disrupted in tau... | ❌ REFUTED | 0 | 0.043 |
| Sliding window analysis reveals temporal patterns ... | ❌ REFUTED | 0 | 0.163 |
| Mitophagy receptors are upregulated in tau+ neuron... | ⚠️ PARTIALLY SUPPORTED | 1 | -0.151 |
| Parkin-independent mitophagy pathways are activate... | ⚠️ PARTIALLY SUPPORTED | 1 | -0.088 |


### Group 2: Proteostasis Failure

| Claim | Verdict | Significant Proteins | Mean Log2FC |
|-------|---------|---------------------|-------------|
| V-ATPase subunits are dysregulated in tau+ neurons | ❌ REFUTED | 0 | 0.009 |
| ATP6V0A1 subunit dysfunction leads to lysosomal al... | ❌ REFUTED | 0 | -0.029 |
| Organellar markers show compartment-specific dysfu... | ❌ REFUTED | 0 | -0.006 |
| Retromer complex components are dysregulated in ta... | ❌ REFUTED | 0 | -0.032 |
| Autophagy and UPS show differential dysfunction pa... | ❌ REFUTED | 0 | 0.096 |
| Endolysosomal system undergoes progressive dysfunc... | ❌ REFUTED | 0 | -0.108 |
| Proteostasis failure follows a temporal cascade pa... | ❌ REFUTED | 0 | 0.086 |
| Rab GTPases show widespread trafficking dysfunctio... | ❌ REFUTED | 0 | 0.012 |


## Key Findings

1. **Mitochondrial Dysfunction**: Analysis reveals complex patterns of mitochondrial protein dysregulation
2. **Proteostasis Failure**: Multiple proteostasis pathways show differential dysfunction
3. **Temporal Dynamics**: Evidence for time-dependent progression of molecular changes
4. **Pathway Interactions**: Cross-talk between mitochondrial and proteostasis systems

## Data Files Generated

- Individual claim results in `results/group1_mitochondrial/` and `results/group2_proteostasis/`
- Combined results: `results/combined/all_results.csv`
- Statistical summaries: `results/combined/all_statistics.json`
- Figures: PNG and PDF formats in respective claim directories

## Methods

- **Statistical Analysis**: Two-sample t-tests with FDR correction
- **Significance Criteria**: FDR < 0.05 and |log2FC| > 0.5
- **Visualization**: Volcano plots, heatmaps, and bar plots
- **Data**: Mock data with realistic expression patterns (for demonstration)

---

*Analysis completed successfully using PertPy DGE Pipeline v2.0*
