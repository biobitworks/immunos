# Project Datasets

Centralized storage for datasets used across all projects in this workspace.

## Directory Structure

```
data/
├── genage/              # Human Ageing Genomic Resources (HAGR)
│   ├── human/          # Human aging-related genes (307 genes)
│   ├── models/         # Model organism longevity genes (2,205 genes)
│   └── expression/     # Gene expression signatures of aging
├── medical/            # Medical datasets (future)
├── benchmarks/         # ML benchmark datasets (future)
└── experimental/       # Temporary/test datasets (future)
```

## Current Datasets

### GenAge Database

**Source**: Human Ageing Genomic Resources (HAGR)
**Website**: https://genomics.senescence.info/
**License**: Creative Commons Attribution 3.0 Unported License
**Downloaded**: 2025-12-02

**Databases Included**:
1. **GenAge Human** - 307 manually curated aging-related genes in humans
2. **GenAge Model Organisms** - 2,205 longevity/aging genes across model species
3. **Gene Expression Signatures** - Meta-analysis from 127 datasets across mammals

**Citation**:
```
de Magalhães, J. P., et al. (2024). Human Ageing Genomic Resources:
updates on key databases in ageing research. Nucleic Acids Research,
52(D1), D900–D908. https://doi.org/10.1093/nar/gkad927
```

**Documentation**: See `/research/datasets/genage-*.md` for detailed schemas and usage

## Usage Guidelines

1. **Version Control**: Small datasets (<25MB) are committed to git
2. **Large Files**: Files >25MB excluded via .gitignore
3. **Metadata**: Each dataset includes metadata.json with download date, source, schema
4. **Documentation**: Cross-referenced in Obsidian vault at `/research/datasets/`

## Adding New Datasets

1. Create subdirectory under appropriate category
2. Add README.md with source, license, citation
3. Create metadata.json with schema information
4. Document in `/research/datasets/dataset-name.md`
5. Update this README with dataset description

## Related Documentation

- **Citations**: `/research/citations/` - BibTeX files for academic references
- **Literature Notes**: `/research/literature/` - Paper summaries and analyses
- **Datasets Docs**: `/research/datasets/` - Detailed dataset documentation
- **Experiments**: `/research/experiments/` - Analyses using these datasets

---

**Last Updated**: 2025-12-02
**Total Datasets**: 1 source (3 datasets)
**Total Size**: ~10-25 MB (estimated)
