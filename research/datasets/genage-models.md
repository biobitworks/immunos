---
type: dataset
tags: [data, aging, genomics, genage, model-organisms, longevity]
source: genomics.senescence.info
downloaded: 2025-12-02
build: 21
format: CSV
records: 2205
organisms: 9
---

# GenAge Model Organisms Dataset

## Overview

**Source**: Human Ageing Genomic Resources (HAGR) - GenAge Database
**URL**: https://genomics.senescence.info/genes/models.html
**Download URL**: https://genomics.senescence.info/genes/models_genes.zip
**Build**: 21 (released August 28, 2023)
**Downloaded**: 2025-12-02

**Description**: Database of 2,205 genes associated with longevity and aging across 9 model organisms. Each gene has experimental evidence for affecting lifespan, aging rate, or age-related phenotypes.

---

## File Information

**Location**: `[[../../data/genage/models/genage_models.csv|/data/genage/models/genage_models.csv]]`
**Format**: Tab-delimited CSV
**Size**: 214,711 bytes (~215 KB)
**Records**: 2,205 genes
**Organisms**: 9 species
**Compression**: Originally distributed as `.zip` file

---

## Schema

| Column | Type | Description |
|--------|------|-------------|
| `GenAge ID` | Integer | Unique identifier in GenAge database |
| `symbol` | String | Gene symbol in organism-specific nomenclature |
| `name` | String | Full gene name/description |
| `organism` | String | Scientific name of organism |
| `entrez gene id` | Integer | NCBI Entrez Gene ID (when available) |
| `avg lifespan change (max obsv)` | Float | Maximum observed lifespan change (%) |
| `lifespan effect` | String | Increase/Decrease/Unclear/Mixed |
| `longevity influence` | String | Pro-longevity/Anti-longevity/Mixed/Unclear |

---

## Organisms Coverage

| Organism | Common Name | Genes | Percentage |
|----------|-------------|-------|------------|
| *Saccharomyces cerevisiae* | Budding yeast | 911 | 41.3% |
| *Caenorhabditis elegans* | Nematode worm | 889 | 40.3% |
| *Drosophila melanogaster* | Fruit fly | 202 | 9.2% |
| *Mus musculus* | House mouse | 136 | 6.2% |
| *Schizosaccharomyces pombe* | Fission yeast | 61 | 2.8% |
| *Podospora anserina* | Filamentous fungus | 3 | 0.1% |
| *Mesocricetus auratus* | Syrian hamster | 1 | <0.1% |
| *Danio rerio* | Zebrafish | 1 | <0.1% |
| *Caenorhabditis briggsae* | Nematode (relative of *C. elegans*) | 1 | <0.1% |
| **Total** | | **2,205** | **100%** |

**Key Insight**: Yeast (*S. cerevisiae*) and worms (*C. elegans*) account for 81.6% of entries, reflecting extensive research in these organisms.

---

## Sample Data

**Representative Genes**:

### *C. elegans* - daf-2 (Insulin/IGF-1 receptor)
- **Effect**: Decrease in daf-2 increases lifespan up to 100%
- **Influence**: Pro-longevity (when reduced)
- **Conserved**: Human ortholog is INSR/IGF1R

### *S. cerevisiae* - SIR2 (Silent information regulator 2)
- **Effect**: Overexpression increases replicative lifespan
- **Influence**: Pro-longevity
- **Conserved**: Human ortholog is SIRT1 (sirtuin family)

### *D. melanogaster* - InR (Insulin-like receptor)
- **Effect**: Reduced signaling increases lifespan
- **Influence**: Pro-longevity (when reduced)
- **Conserved**: Similar to *C. elegans* daf-2

### *M. musculus* - Gh (Growth hormone)
- **Effect**: Deficiency increases lifespan up to 40%
- **Influence**: Pro-longevity (when reduced)
- **Conserved**: Highly conserved in mammals

---

## Usage Examples

### Python - Load and Explore

```python
import pandas as pd

# Load dataset
df = pd.read_csv('/Users/byron/projects/data/genage/models/genage_models.csv',
                 sep='\t')

print(f"Total longevity genes: {len(df)}")
print(f"Organisms: {df['organism'].nunique()}")

# Count genes per organism
organism_counts = df['organism'].value_counts()
print("\nGenes per organism:")
print(organism_counts)
```

### Python - Filter by Organism

```python
# Get all C. elegans longevity genes
worm_genes = df[df['organism'] == 'Caenorhabditis elegans']
print(f"C. elegans genes: {len(worm_genes)}")

# Get all mouse genes
mouse_genes = df[df['organism'] == 'Mus musculus']
print(f"Mouse genes: {len(mouse_genes)}")

# Get all genes that increase lifespan
pro_longevity = df[df['lifespan effect'] == 'Increase']
print(f"Genes that increase lifespan: {len(pro_longevity)}")
```

### Python - Cross-Species Analysis

```python
# Find genes conserved across species
gene_symbols = df.groupby('symbol')['organism'].apply(list)
multi_species = gene_symbols[gene_symbols.apply(len) > 1]

print(f"Genes found in multiple organisms: {len(multi_species)}")
for gene, organisms in multi_species.items():
    print(f"  {gene}: {len(organisms)} species")
```

---

## Research Applications

### 1. Comparative Longevity Research

**Goal**: Identify conserved aging mechanisms across species

```python
# Compare lifespan effects across organisms
for org in df['organism'].unique():
    org_df = df[df['organism'] == org]
    increase = len(org_df[org_df['lifespan effect'] == 'Increase'])
    decrease = len(org_df[org_df['lifespan effect'] == 'Decrease'])
    print(f"{org}:")
    print(f"  Pro-longevity: {increase}")
    print(f"  Anti-longevity: {decrease}")
```

### 2. Pathway Analysis

**Goal**: Identify biological pathways enriched in longevity genes

- Group by organism
- Extract Entrez Gene IDs
- Submit to pathway analysis (KEGG, Reactome, GO)
- Compare pathways across species

### 3. Translational Research

**Goal**: Identify mouse genes for mammalian studies

```python
# Get mouse longevity genes for translation to human
mouse_df = df[df['organism'] == 'Mus musculus']

# Filter for pro-longevity genes
pro_long_mouse = mouse_df[mouse_df['longevity influence'].str.contains('Pro', na=False)]

# Export for cross-reference with human genes
pro_long_mouse[['symbol', 'name', 'lifespan effect']].to_csv(
    'mouse_pro_longevity_genes.csv', index=False
)
```

### 4. Drug Target Identification

**Goal**: Find druggable targets validated in model organisms

- Cross-reference with DrugAge database
- Identify genes with known modulators
- Prioritize by effect size and conservation

---

## Key Insights

### 1. Model Organism Distribution

**Heavily Studied**:
- *S. cerevisiae* (yeast): Simple, fast-growing, genetic tools
- *C. elegans* (worm): Short lifespan, transparent, easy genetics

**Moderately Studied**:
- *D. melanogaster* (fly): Complex multicellular, aging phenotypes
- *M. musculus* (mouse): Mammalian model, translational relevance

**Understudied**:
- Other organisms: Limited entries, specialized research

### 2. Conservation of Aging Pathways

**Conserved Mechanisms**:
- Insulin/IGF-1 signaling (daf-2, InR, Gh)
- Sir2/Sirtuin pathway (SIR2, SIRT1)
- mTOR signaling
- Mitochondrial function

**Implication**: Fundamental aging mechanisms are evolutionarily ancient

### 3. Lifespan Effects

Many genes show:
- **Lifespan increase** when reduced/deleted
- **Examples**: Insulin receptor, growth hormone, TOR
- **Interpretation**: Slowing growth/metabolism extends lifespan

---

## Integration with IMMUNOS-MCP

### Comparative Pattern Recognition

**B Cell Training**:
```python
# Train on species-specific aging patterns
from src.agents.bcell_agent import BCellAgent
from src.core.antigen import Antigen, DataType

# Group by organism
for organism, group in df.groupby('organism'):
    antigens = []
    for _, row in group.iterrows():
        antigen = Antigen(
            data=f"{row['symbol']}: {row['longevity influence']}",
            data_type=DataType.TEXT,
            class_label=row['lifespan effect'],
            identifier=f"{organism}_{row['symbol']}"
        )
        antigens.append(antigen)

    # Train organism-specific B Cell
    bcell = BCellAgent(f'bcell_{organism}')
    bcell.train(antigens)
```

### Cross-Species Analysis

**Research Questions**:
1. Can AI detect conserved aging patterns across species?
2. Do different organisms share aging gene signatures?
3. Can worm/fly data predict mouse/human aging genes?

---

## Related Datasets

- **Human Genes**: [[genage-human|GenAge Human]] (307 genes)
- **Expression**: [[genage-expression|Gene Expression Signatures]]
- **Literature**: [[../literature/demagalhaes-2024-hagr-genage|de Magalhães et al. (2024)]]

---

## Citation

**Data Citation**:
```
Human Ageing Genomic Resources (HAGR). (2024). GenAge Model Organisms Dataset.
Build 21. Retrieved from https://genomics.senescence.info/genes/models.html
```

**Primary Publication**:
```
de Magalhães, J. P., et al. (2024). Human Ageing Genomic Resources: updates on
key databases in ageing research. Nucleic Acids Research, 52(D1), D900–D908.
```

**BibTeX**: [[../citations/genage-references.bib|genage-references.bib]]

---

## Experiments Using This Dataset

- [[../experiments/genage-analysis-2025-12-02|GenAge Analysis (2025-12-02)]]

---

## Metadata Summary

- **Records**: 2,205 genes
- **Organisms**: 9 species
- **Build**: 21 (August 2023)
- **Downloaded**: 2025-12-02
- **Size**: 215 KB
- **Format**: CSV (tab-delimited)
- **License**: CC BY 3.0

---

**Last Updated**: 2025-12-02
**Status**: Downloaded and analyzed
