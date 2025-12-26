# GenAge Database - Human Ageing Genomic Resources

## Summary
Contains 3 subdirectories and 0 files. Key subfolders: expression/, human/, models/.


**Source**: Human Ageing Genomic Resources (HAGR)
**Website**: https://genomics.senescence.info/
**Download Page**: https://genomics.senescence.info/download.html
**Downloaded**: 2025-12-02
**Build**: 21 (released August 28, 2023)

## Overview

GenAge is a curated database of genes associated with aging and longevity in humans and model organisms. Part of the Human Ageing Genomic Resources (HAGR) collection maintained by the Integrative Genomics of Ageing Group.

## License

**Creative Commons Attribution 3.0 Unported License**
- Free for commercial, educational, and research use
- Attribution required in publications
- Full data export and reuse permitted

## Citation

**Primary (2024 Update)**:
```
de Magalhães, J. P., Abidi, Z., Dos Santos, G. A., Avelar, R. A., Barardo, D.,
Chatsirisupachai, K., Clark, P., De-Souza, E. A., Johnson, E. J., Lopes, I.,
Novoa, G., Senez, L., Talay, A., Thornton, D., & To, P. K. P. (2024).
Human Ageing Genomic Resources: updates on key databases in ageing research.
Nucleic Acids Research, 52(D1), D900–D908. https://doi.org/10.1093/nar/gkad927
```

**BibTeX**: See `/research/citations/genage-references.bib`

## Datasets

### 1. GenAge Human (`human/`)

**Description**: Manually curated database of human genes associated with aging
**Records**: 307 genes (Build 21)
**Format**: Tab-delimited ASCII / CSV
**Files**:
- `genage_human.zip` - Original compressed download
- `genage_human.csv` - Extracted data
- `metadata.json` - Schema and download information

**Typical Columns** (exact schema in metadata.json):
- Gene Symbol
- Gene Name
- Entrez Gene ID
- Function/Description
- Evidence for aging association
- References

### 2. GenAge Model Organisms (`models/`)

**Description**: Genes associated with longevity and aging in model species
**Records**: 2,205 genes (Build 21)
**Organisms**: Yeast, C. elegans, D. melanogaster, mice, rats, and others
**Format**: Tab-delimited ASCII / CSV
**Files**:
- `genage_models.zip` - Original compressed download
- `genage_models.csv` - Extracted data
- `metadata.json` - Schema and download information

**Typical Columns**:
- Gene Symbol
- Gene Name
- Organism
- Entrez Gene ID
- Observations
- Longevity influence
- References

### 3. Gene Expression Signatures (`expression/`)

**Description**: Meta-analysis of gene expression changes with aging across mammals
**Source**: Analysis of 127 microarray and RNA-Seq datasets
**Species**: Humans, mice, rats
**Format**: Zipped supplementary data
**Files**:
- `signatures_supplement.zip` - Original download
- `signatures_data.csv` - Extracted expression data
- `metadata.json` - Study information

**Content**: Genes consistently over- or under-expressed with age across tissues

## Usage Examples

### Python - Load Human Genes

```python
import pandas as pd

# Load GenAge human genes
df = pd.read_csv('/Users/byron/projects/data/genage/human/genage_human.csv',
                 sep='\t')

print(f"Total genes: {len(df)}")
print(df.head())
```

### Python - Load Model Organisms

```python
import pandas as pd

# Load model organism genes
df = pd.read_csv('/Users/byron/projects/data/genage/models/genage_models.csv',
                 sep='\t')

# Filter by organism
worm_genes = df[df['organism'] == 'Caenorhabditis elegans']
print(f"C. elegans longevity genes: {len(worm_genes)}")
```

## Related Resources

- **HAGR Homepage**: https://genomics.senescence.info/
- **Help Documentation**: https://genomics.senescence.info/help.html
- **Release Notes**: https://genomics.senescence.info/genes/release.html
- **Statistics**: https://genomics.senescence.info/genes/stats.php

## Obsidian Documentation

- **Literature Note**: `/research/literature/demagalhaes-2024-hagr-genage.md`
- **Dataset Docs**:
  - `/research/datasets/genage-human.md`
  - `/research/datasets/genage-models.md`
  - `/research/datasets/genage-expression.md`

## Integration Examples

- **Immunosenescence Research**: Cross-reference aging genes with immune system genes
- **IMMUNOS-MCP**: Train agents on biological aging patterns
- **Comparative Analysis**: Study aging across species
- **Drug Discovery**: Identify potential longevity interventions

---

**Last Updated**: 2025-12-02
**Build**: 21 (August 2023)
**Total Genes**: 2,512 (307 human + 2,205 model organisms)

## Directory Map
```
genage/
├── expression/
├── human/
└── models/
```
