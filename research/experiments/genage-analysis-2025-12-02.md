---
title: GenAge Data Analysis Summary
date: 2025-12-04
tags: [genage, aging, data-analysis, genomics]
type: analysis-report
---

# GenAge Data Analysis Summary

**Analysis Date**: 2025-12-04 23:37:19
**Data Source**: Human Ageing Genomic Resources (HAGR)
**Website**: https://genomics.senescence.info/

## Overview

Comprehensive analysis of three GenAge datasets:
1. Human aging-related genes
2. Model organism longevity genes
3. Gene expression signatures of aging

---

## 1. GenAge Human Genes

**File**: `data/genage/human/genage_human.csv`
**Total Genes**: 307

### Schema

**Columns** (6 total):
- `GenAge ID`
- `symbol`
- `name`
- `entrez gene id`
- `uniprot`
- `why`

### Sample Genes

**First 5 genes**:
- GHR
- GHRH
- SHC1
- POU1F1
- PROP1


---

## 2. GenAge Model Organisms

**File**: `data/genage/models/genage_models.csv`
**Total Genes**: 2205
**Organisms**: 9

### Schema

**Columns** (8 total):
- `GenAge ID`
- `symbol`
- `name`
- `organism`
- `entrez gene id`
- `avg lifespan change (max obsv)`
- `lifespan effect`
- `longevity influence`

### Organisms Distribution

**Top 10 organisms by gene count**:

| Organism | Genes |
|----------|-------|
| Saccharomyces cerevisiae | 911 |
| Caenorhabditis elegans | 889 |
| Drosophila melanogaster | 202 |
| Mus musculus | 136 |
| Schizosaccharomyces pombe | 61 |
| Podospora anserina | 3 |
| Mesocricetus auratus | 1 |
| Danio rerio | 1 |
| Caenorhabditis briggsae | 1 |


---

## 3. Gene Expression Signatures

**File**: `data/genage/expression/ageing_signatures.xlsx`
**Format**: Excel (.xlsx)
**Size**: 52,781 bytes

**Description**: Meta-analysis results from 127 microarray and RNA-Seq datasets across mammals (human, mouse, rat).

### Sheets Structure

The Excel file contains 3 sheets:
1. **Source** - Citation information (Palmer et al. ageing transcriptome meta-analysis)
2. **TableS3-Over.All** - Genes overexpressed with age (450 genes)
3. **TableS7-Under.All** - Genes underexpressed with age (163 genes)

### Expression Patterns

#### Overexpressed with Age (450 genes)

**Top genes** (highest significance):
- C1QA (QValue = 0.0)
- GPNMB (QValue = 0.0)
- B2M (QValue = 0.0)
- EFEMP1 (QValue = 0.0)
- C1QC, C1QB (QValue = 0.0)
- SKAP2, IL33, PLEK (QValue = 0.0)

**Patterns**: Inflammation, immune activation, extracellular matrix remodeling

#### Underexpressed with Age (163 genes)

**Top genes** (highest significance):
- UQCRFS1 (QValue = 0.000022)
- SUCLG1 (QValue = 0.000023)
- MLF1 (QValue = 0.000051)
- UROS (QValue = 0.000125)
- FKBP4 (QValue = 0.000103)
- VLDLR, TFRC, VPS36, COX7B (QValue < 0.001)

**Patterns**: Mitochondrial function, energy metabolism, oxidative phosphorylation

### Summary Statistics

| Category | Count | Significance |
|----------|-------|--------------|
| Overexpressed | 450 | Inflammation, immune response |
| Underexpressed | 163 | Metabolism, mitochondrial function |
| **Total** | **613** | Consistent across tissues |

### Key Insights from Expression Data

1. **Inflammaging Pattern**: Overexpressed genes dominated by immune/inflammatory markers (C1QA, C1QB, C1QC, IL33)
2. **Metabolic Decline**: Underexpressed genes show clear mitochondrial/metabolic signature (UQCRFS1, COX7B, SUCLG1)
3. **Cross-tissue consistency**: Genes selected across all tissues (meta-analysis of 127 studies)
4. **High confidence**: Most significant genes have QValue < 0.001 (very low false discovery rate)

**Note**: Full analysis completed with openpyxl and pandas on 2025-12-04

---

## Summary Statistics

| Dataset | Records | File Size | Format |
|---------|---------|-----------|--------|
| Human Genes | 307 | 22,778 bytes | CSV |
| Model Organisms | 2205 | 214,711 bytes | CSV |
| Expression Signatures | 127 datasets | 52,781 bytes | Excel |

**Total Genes**: 2,512

---

## Key Insights

### 1. Comprehensive Coverage
- 307 manually curated human aging genes
- 2205 genes across 9 model organisms
- Cross-species comparative analysis possible

### 2. Data Quality
- Manual curation by aging research experts
- Extensive literature references
- Regular updates (Build 21, August 2023)

### 3. Research Applications

#### Immunosenescence Research
- Cross-reference aging genes with immune system genes
- Identify overlap between aging and immune decline
- Study inflammatory aging (inflammaging)

#### IMMUNOS-MCP Integration
- Train agents on biological aging patterns
- Pattern recognition for age-related changes
- Anomaly detection for premature aging

#### Comparative Biology
- Study conservation of aging mechanisms
- Identify species-specific longevity factors
- Understand evolutionary perspectives on aging

#### Drug Discovery
- Target identification for longevity interventions
- Repurposing existing drugs for aging
- Cross-reference with DrugAge database

---

## Data Access

### Python Example

```python
import pandas as pd

# Load human genes
human_df = pd.read_csv('/Users/byron/projects/data/genage/human/genage_human.csv')
print(f"Human genes: {len(human_df)}")

# Load model organisms
models_df = pd.read_csv('/Users/byron/projects/data/genage/models/genage_models.csv')
print(f"Model organism genes: {len(models_df)}")

# Load expression signatures (requires openpyxl)
# expr_df = pd.read_excel('/Users/byron/projects/data/genage/expression/ageing_signatures.xlsx')
```

---

## Related Documentation

- **Dataset READMEs**:
  - [[../data/genage/README|GenAge README]]
  - [[../data/README|Data Directory README]]
- **Literature Notes**:
  - [[literature/demagalhaes-2024-hagr-genage|de Magalhães et al. (2024)]]
- **Citations**:
  - [[citations/genage-references.bib|GenAge BibTeX References]]

---

## Next Steps

1. **Detailed Analysis**:
   - [ ] Analyze gene categories and pathways
   - [ ] Identify most frequently cited genes
   - [ ] Cross-reference with pathway databases

2. **Integration**:
   - [ ] Map to immune system genes
   - [ ] Identify aging-immunity overlap
   - [ ] Create visualization of relationships

3. **Research**:
   - [ ] Literature review of top aging genes
   - [ ] Comparative analysis across species
   - [ ] Hypothesis generation for experiments

---

**Last Updated**: 2025-12-04
**Data Source**: HAGR GenAge Build 21 (August 2023)
**License**: CC BY 3.0
