---
title: Dataset Documentation
tags: [datasets, documentation, data-catalog]
type: index
---

# Dataset Documentation

This directory contains detailed documentation for all datasets used in research projects. Each dataset has a dedicated note with schema, usage examples, and integration information.

## Purpose

- **Schema Documentation**: Column definitions, data types, constraints
- **Usage Examples**: Code snippets for loading and analyzing data
- **Research Context**: How datasets relate to experiments and literature
- **Quality Metadata**: Source, curation process, license information
- **Integration Guides**: How to use with IMMUNOS-MCP and other projects

---

## Available Datasets

### GenAge - Human Ageing Genomic Resources

**Source**: https://genomics.senescence.info/
**Downloaded**: 2025-12-02
**Build**: 21 (August 2023)

1. **[[genage-human|GenAge Human Genes]]**
   - 307 manually curated human aging genes
   - Tab-delimited CSV format
   - 22 KB size

2. **[[genage-models|GenAge Model Organisms]]**
   - 2,205 longevity/aging genes across 9 species
   - Tab-delimited CSV format
   - 215 KB size

3. **[[genage-expression|Gene Expression Signatures]]**
   - Meta-analysis from 127 datasets
   - Excel (.xlsx) format
   - 53 KB size

**Citation**: de Magalhães et al. (2024) Nucleic Acids Research
**License**: CC BY 3.0
**Documentation**: [[../literature/demagalhaes-2024-hagr-genage|Literature Note]]

---

## Dataset Categories

### Aging & Longevity
- [[genage-human|GenAge Human]]
- [[genage-models|GenAge Model Organisms]]
- [[genage-expression|Gene Expression Signatures]]

### Medical Datasets
*(To be added)*

### Benchmark Datasets
- [[mahnob-hci|MAHNOB-HCI]] (multimodal emotion recognition; access request required)

---

## Quick Reference

### Loading Data (Python)

```python
import pandas as pd

# GenAge Human Genes
human_df = pd.read_csv('/Users/byron/projects/data/genage/human/genage_human.csv',
                       sep='\t')

# GenAge Model Organisms
models_df = pd.read_csv('/Users/byron/projects/data/genage/models/genage_models.csv',
                        sep='\t')

# Gene Expression Signatures (requires openpyxl)
# expr_df = pd.read_excel('/Users/byron/projects/data/genage/expression/ageing_signatures.xlsx')
```

### Common Operations

```python
# Get gene count
print(f"Total genes: {len(human_df)}")

# List available columns
print(f"Columns: {list(human_df.columns)}")

# Search for specific gene
gene = human_df[human_df['symbol'] == 'TP53']

# Cross-reference datasets
human_genes = set(human_df['symbol'])
models_genes = set(models_df['symbol'])
overlap = human_genes.intersection(models_genes)
```

---

## Documentation Template

Each dataset note should include:

1. **Overview**: Source, URL, download date, description
2. **File Information**: Location, format, size, records
3. **Schema**: Table of columns with types and descriptions
4. **Sample Data**: First few records or representative examples
5. **Data Quality**: Curation process, inclusion criteria
6. **Usage Examples**: Code snippets in Python/R
7. **Research Applications**: How to use for specific research questions
8. **Integration**: How to use with IMMUNOS-MCP or other tools
9. **Related Datasets**: Links to complementary data
10. **Citation**: How to cite the dataset
11. **License**: Terms of use
12. **Metadata**: Build version, update frequency, status

---

## Related Documentation

- **Data Storage**: [[../../data/README|Data Directory]]
- **Literature Notes**: [[../literature/README|Literature]]
- **Experiments**: [[../experiments/README|Experiments]]
- **Citations**: [[../citations/README|Citations]]

---

## Adding New Datasets

### Process

1. **Download & Store**: Place in `/data/[category]/[dataset-name]/`
2. **Create Documentation**: Use template above
3. **Update This Index**: Add entry to appropriate category
4. **Create Literature Note**: If from publication
5. **Add to Citations**: Create BibTeX entry
6. **Link to Experiments**: Connect to analysis notes

### Naming Convention

- **File**: `dataset-name.md` (lowercase, hyphens)
- **Title**: Human-readable dataset name
- **Tags**: Include dataset type, domain, source

---

**Last Updated**: 2025-12-02
**Total Datasets**: 3 (GenAge collection)
**Total Records**: 2,512 genes
