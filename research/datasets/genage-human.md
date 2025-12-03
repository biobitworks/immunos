---
type: dataset
tags: [data, aging, genomics, genage, human-genes]
source: genomics.senescence.info
downloaded: 2025-12-02
build: 21
format: CSV
records: 307
---

# GenAge Human Genes Dataset

## Overview

**Source**: Human Ageing Genomic Resources (HAGR) - GenAge Database
**URL**: https://genomics.senescence.info/genes/human.html
**Download URL**: https://genomics.senescence.info/genes/human_genes.zip
**Build**: 21 (released August 28, 2023)
**Downloaded**: 2025-12-02

**Description**: Manually curated database of 307 genes associated with human aging. Each gene has been identified through literature review and experimental evidence as playing a role in the aging process.

---

## File Information

**Location**: `[[../../data/genage/human/genage_human.csv|/data/genage/human/genage_human.csv]]`
**Format**: Tab-delimited CSV
**Size**: 22,778 bytes
**Records**: 307 genes
**Compression**: Originally distributed as `.zip` file

---

## Schema

| Column | Type | Description |
|--------|------|-------------|
| `GenAge ID` | Integer | Unique identifier in GenAge database |
| `symbol` | String | Standard gene symbol (e.g., "GHR", "GHRH") |
| `name` | String | Full gene name |
| `entrez gene id` | Integer | NCBI Entrez Gene identifier for cross-referencing |
| `uniprot` | String | UniProt protein database identifier |
| `why` | Text | Rationale for inclusion in database with citations |

---

## Sample Data

**First 5 genes**:
1. **GHR** - Growth hormone receptor
2. **GHRH** - Growth hormone releasing hormone
3. **SHC1** - SHC adaptor protein 1
4. **POU1F1** - POU class 1 homeobox 1
5. **PROP1** - PROP paired-like homeobox 1

---

## Data Quality

### Curation Process
- **Manual Review**: Each gene curated by aging biology experts
- **Evidence-Based**: Requires experimental evidence from literature
- **Referenced**: Each entry includes citation justification in "why" field
- **Peer-Reviewed**: Database updates undergo expert review

### Inclusion Criteria
- Experimental evidence in humans or model organisms
- Association with aging phenotypes or lifespan
- Published in peer-reviewed literature
- Validated by multiple studies (preferred)

---

## Usage Examples

### Python - Load and Explore

```python
import pandas as pd

# Load dataset
df = pd.read_csv('/Users/byron/projects/data/genage/human/genage_human.csv',
                 sep='\t')

print(f"Total aging-related genes: {len(df)}")
print(f"Columns: {list(df.columns)}")

# Sample genes
print("\nFirst 10 genes:")
print(df[['symbol', 'name']].head(10))

# Search for specific gene
gene = "TP53"
if gene in df['symbol'].values:
    info = df[df['symbol'] == gene].iloc[0]
    print(f"\n{gene}: {info['name']}")
    print(f"Rationale: {info['why'][:200]}...")
```

### Python - Cross-Reference with Other Data

```python
# Get list of aging gene symbols
aging_genes = set(df['symbol'].tolist())

# Cross-reference with your gene list
your_genes = ['TP53', 'APOE', 'IL6', 'TNF', 'INS']
overlap = aging_genes.intersection(your_genes)

print(f"Genes in both lists: {overlap}")
```

### Python - Export to Different Format

```python
# Export as JSON
df.to_json('/path/to/genage_human.json', orient='records', indent=2)

# Export subset (e.g., genes with specific keywords)
immune_related = df[df['name'].str.contains('immune|inflammation|cytokine',
                                            case=False, na=False)]
immune_related.to_csv('/path/to/aging_immune_genes.csv', index=False)
```

---

## Research Applications

### 1. Immunosenescence Research

**Goal**: Identify genes involved in both aging and immune function

```python
# Find aging genes related to immunity
immune_keywords = ['immune', 'inflammation', 'cytokine', 'interleukin',
                   'lymphocyte', 'T cell', 'B cell']

immune_aging_genes = df[df['name'].str.contains('|'.join(immune_keywords),
                                                 case=False, na=False)]

print(f"Aging genes with immune function: {len(immune_aging_genes)}")
```

### 2. Pathway Enrichment Analysis

**Goal**: Identify biological pathways enriched in aging genes

- Extract Entrez Gene IDs
- Submit to pathway analysis tools (DAVID, GSEA, etc.)
- Identify overrepresented pathways

### 3. Comparative Genomics

**Goal**: Compare human aging genes to model organisms

- Cross-reference with [[genage-models|Model Organisms dataset]]
- Identify conserved aging genes
- Study species-specific aging mechanisms

### 4. Drug Target Identification

**Goal**: Find druggable targets in aging pathways

- Cross-reference with DrugBank
- Identify known modulators
- Predict novel interventions

---

## Integration with IMMUNOS-MCP

### Training Data

**Pattern Recognition**:
- Safe aging patterns vs pathological aging
- Normal gene expression vs age-related changes
- Healthy aging trajectories

**B Cell Agent**:
```python
# Train on aging gene patterns
from src.agents.bcell_agent import BCellAgent
from src.core.antigen import Antigen, DataType

# Load aging genes as patterns
aging_df = pd.read_csv('data/genage/human/genage_human.csv', sep='\t')

# Create antigens for normal aging genes
antigens = []
for _, row in aging_df.iterrows():
    antigen = Antigen(
        data=f"{row['symbol']}: {row['name']}",
        data_type=DataType.TEXT,
        class_label='normal_aging',
        identifier=row['symbol']
    )
    antigens.append(antigen)

# Train B Cell
bcell = BCellAgent('aging_bcell')
bcell.train(antigens)
```

### Anomaly Detection

**NK Cell Agent**:
- Detect genes abnormally expressed for age
- Identify premature aging signatures
- Flag potential age-related disease risk

### Research Questions

1. Can immune system algorithms detect aging patterns?
2. Do aging genes cluster in biological networks?
3. Can we predict biological age from gene expression?
4. What genes are involved in both aging and immune decline?

---

## Related Datasets

- **Model Organisms**: [[genage-models|GenAge Model Organisms]] (2,205 genes)
- **Expression Signatures**: [[genage-expression|Gene Expression Signatures]] (127 datasets)
- **DrugAge**: Compounds affecting lifespan (separate database)
- **LongevityMap**: Human longevity variants (separate database)

---

## Citation

**Data Citation**:
```
Human Ageing Genomic Resources (HAGR). (2024). GenAge Human Genes Dataset.
Build 21. Retrieved from https://genomics.senescence.info/genes/human.html
```

**Primary Publication**:
```
de Magalhães, J. P., et al. (2024). Human Ageing Genomic Resources: updates on
key databases in ageing research. Nucleic Acids Research, 52(D1), D900–D908.
https://doi.org/10.1093/nar/gkad927
```

**BibTeX**: [[../citations/genage-references.bib|genage-references.bib]]

---

## Experiments Using This Dataset

- [[../experiments/genage-analysis-2025-12-02|GenAge Analysis (2025-12-02)]] - Initial exploratory analysis

---

## License

**Creative Commons Attribution 3.0 Unported License**
- Free for commercial, educational, and research use
- Attribution required in publications
- Full export and reuse permitted

---

## Metadata Summary

- **Records**: 307 human genes
- **Build**: 21 (August 2023)
- **Downloaded**: 2025-12-02
- **Size**: 22 KB
- **Format**: CSV (tab-delimited)
- **Quality**: Manually curated, peer-reviewed
- **Update Frequency**: Irregular (Build 21 released August 2023)

---

**Last Updated**: 2025-12-02
**Status**: Downloaded and analyzed
**Next Steps**: Cross-reference with immune pathways, train IMMUNOS-MCP agents
