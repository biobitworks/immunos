---
type: dataset
tags: [data, aging, genomics, genage, gene-expression, transcriptomics, meta-analysis]
source: genomics.senescence.info
downloaded: 2025-12-02
format: Excel
datasets: 127
---

# Gene Expression Signatures of Ageing

## Overview

**Source**: Human Ageing Genomic Resources (HAGR) - GenAge Database
**URL**: https://genomics.senescence.info/gene_expression/signatures.html
**Download URL**: https://genomics.senescence.info/gene_expression/signatures_supplement.zip
**Downloaded**: 2025-12-02

**Description**: Meta-analysis results from 127 microarray and RNA-Seq datasets identifying genes with consistent expression changes during aging across mammals (human, mouse, rat).

---

## File Information

**Location**: `[[../../data/genage/expression/ageing_signatures.xlsx|/data/genage/expression/ageing_signatures.xlsx]]`
**Format**: Excel (.xlsx)
**Size**: 52,781 bytes (~53 KB)
**Datasets Analyzed**: 127 expression studies
**Species**: Human, mouse, rat
**Technologies**: Microarray, RNA-Seq

---

## Background

### Meta-Analysis Approach

**Original Publication**: de Magalhães, Curado, & Church (2009) Bioinformatics

**Methodology**:
1. **Data Collection**: 127 independent aging studies
2. **Species Coverage**: Human, mouse, rat tissues
3. **Technology**: Both microarray and RNA-Seq platforms
4. **Analysis**: Identify genes consistently changing with age
5. **Validation**: Cross-species and cross-tissue comparisons

**Key Finding**: Despite tissue- and species-specific changes, common aging signatures exist

---

## Data Content

### What's Included

**Gene Lists**:
- Genes **up-regulated** with aging (increased expression)
- Genes **down-regulated** with aging (decreased expression)
- Statistical significance for each gene
- Cross-species comparisons

**Metadata**:
- Dataset sources and references
- Tissue types analyzed
- Species information
- Technology platforms

### Expected Patterns

**Up-regulated with age**:
- Inflammatory genes (inflammaging)
- Stress response genes
- DNA damage response
- Lysosomal genes

**Down-regulated with age**:
- Mitochondrial genes
- Energy metabolism
- Protein synthesis
- DNA repair genes

---

## Usage Examples

### Python - Load Excel File

**Requirements**: `pip install openpyxl pandas`

```python
import pandas as pd

# Load Excel file
xlsx_path = '/Users/byron/projects/data/genage/expression/ageing_signatures.xlsx'

# Read all sheets
excel_file = pd.ExcelFile(xlsx_path)
print(f"Sheets in file: {excel_file.sheet_names}")

# Load specific sheet (adjust sheet name as needed)
df = pd.read_excel(xlsx_path, sheet_name=0)

print(f"Rows: {len(df)}")
print(f"Columns: {list(df.columns)}")
print(df.head())
```

### Python - Identify Aging Signatures

```python
# Assuming columns exist for gene symbol and direction
# (Actual column names depend on file structure)

# Get up-regulated genes
if 'direction' in df.columns:
    up_regulated = df[df['direction'] == 'up']
    down_regulated = df[df['direction'] == 'down']

    print(f"Up-regulated genes: {len(up_regulated)}")
    print(f"Down-regulated genes: {len(down_regulated)}")

# Get genes with high consistency across datasets
if 'consistency_score' in df.columns:
    highly_consistent = df[df['consistency_score'] > 0.8]
    print(f"Highly consistent genes: {len(highly_consistent)}")
```

### Python - Cross-Reference with GenAge

```python
# Load GenAge human genes
genage_df = pd.read_csv('/Users/byron/projects/data/genage/human/genage_human.csv',
                        sep='\t')

# Find overlap between aging genes and expression signatures
genage_genes = set(genage_df['symbol'])
expression_genes = set(df['symbol']) if 'symbol' in df.columns else set()

overlap = genage_genes.intersection(expression_genes)
print(f"Genes in both GenAge and expression signatures: {len(overlap)}")
print(f"Examples: {list(overlap)[:10]}")
```

---

## Research Applications

### 1. Aging Biomarker Discovery

**Goal**: Identify genes as potential aging biomarkers

- Genes consistently changing with age
- Cross-tissue and cross-species validation
- Potential blood-based biomarkers

**Approach**:
- Filter for high consistency across datasets
- Prioritize genes measurable in accessible tissues
- Validate in independent cohorts

### 2. Pathway Enrichment

**Goal**: Identify biological processes affected by aging

```python
# Extract gene lists for pathway analysis
up_genes = df[df['direction'] == 'up']['symbol'].tolist()
down_genes = df[df['direction'] == 'down']['symbol'].tolist()

# Submit to:
# - DAVID: https://david.ncifcrf.gov/
# - Enrichr: https://maayanlab.cloud/Enrichr/
# - GSEA: https://www.gsea-msigdb.org/
```

### 3. Drug Repurposing

**Goal**: Find compounds that reverse aging gene expression

- Identify aging expression signature
- Query Connectivity Map (CMap) or LINCS
- Find compounds reversing the signature
- Cross-reference with DrugAge database

### 4. Immunosenescence Research

**Goal**: Study immune system aging at transcriptional level

```python
# Find immune-related genes in aging signatures
immune_keywords = ['immune', 'inflammation', 'cytokine', 'interleukin',
                   'interferon', 'lymphocyte', 'leukocyte']

# Filter for immune genes
if 'gene_name' in df.columns:
    immune_aging = df[df['gene_name'].str.contains('|'.join(immune_keywords),
                                                    case=False, na=False)]
    print(f"Immune genes changing with age: {len(immune_aging)}")
```

---

## Integration with IMMUNOS-MCP

### Pattern Recognition

**B Cell Training**:
```python
# Train on normal vs abnormal aging expression patterns
from src.agents.bcell_agent import BCellAgent
from src.core.antigen import Antigen, DataType

# Create antigens from expression signatures
antigens = []

for _, row in df.iterrows():
    # Encode gene expression pattern
    pattern = f"{row['symbol']}: {row['direction']} with age"

    antigen = Antigen(
        data=pattern,
        data_type=DataType.TEXT,
        class_label='normal_aging',
        identifier=row['symbol']
    )
    antigens.append(antigen)

# Train B Cell on normal aging patterns
bcell = BCellAgent('aging_expression_bcell')
bcell.train(antigens)
```

### Anomaly Detection

**NK Cell Application**:
- Detect abnormal expression patterns
- Identify accelerated aging signatures
- Flag potential pathological aging

**Use Case**:
```python
# Detect if patient expression profile shows accelerated aging
patient_profile = load_patient_expression_data()

# Compare to normal aging signatures
from src.agents.nk_cell_agent import NKCellAgent

nk_cell = NKCellAgent('expression_nk')
nk_cell.train(normal_aging_patterns)

# Detect anomalies
result = nk_cell.detect_anomaly(patient_profile)
if result.is_anomalous:
    print("Warning: Accelerated aging signature detected")
```

---

## Biological Insights

### Hallmarks of Aging (Transcriptional Level)

Based on meta-analysis findings:

1. **Mitochondrial Dysfunction**
   - Down-regulation of OXPHOS genes
   - Reduced energy metabolism

2. **Genomic Instability**
   - Up-regulation of DNA damage response
   - Changes in DNA repair genes

3. **Proteostasis**
   - Changes in chaperone expression
   - Up-regulation of protein degradation

4. **Inflammaging**
   - Up-regulation of inflammatory genes
   - Increased cytokine expression

5. **Cellular Senescence**
   - Expression of senescence markers
   - SASP (senescence-associated secretory phenotype) genes

---

## Limitations & Considerations

### Technical Limitations

1. **Platform Differences**: Microarray vs RNA-Seq have different sensitivities
2. **Batch Effects**: Different studies have different biases
3. **Tissue Specificity**: Some changes are tissue-specific
4. **Species Differences**: Not all findings translate across species

### Biological Considerations

1. **Causation**: Expression changes may be consequence, not cause, of aging
2. **Heterogeneity**: Aging affects different individuals differently
3. **Reversibility**: Some changes may be reversible with interventions
4. **Complexity**: Aging is multifactorial, not single-pathway

---

## Related Datasets

- **Human Genes**: [[genage-human|GenAge Human]] (307 aging genes)
- **Model Organisms**: [[genage-models|Model Organisms]] (2,205 longevity genes)
- **Literature**: [[../literature/demagalhaes-2024-hagr-genage|de Magalhães et al. (2024)]]
- **Original Paper**: de Magalhães, Curado, & Church (2009) Bioinformatics

---

## Citation

**Data Citation**:
```
Human Ageing Genomic Resources (HAGR). (2024). Gene Expression Signatures of Ageing.
Meta-analysis from 127 datasets. Retrieved from
https://genomics.senescence.info/gene_expression/signatures.html
```

**Original Meta-Analysis**:
```
de Magalhães, J. P., Curado, J., & Church, G. M. (2009). Meta-analysis of
age-related gene expression profiles identifies common signatures of aging.
Bioinformatics, 25(7), 875-881.
```

**BibTeX**: [[../citations/genage-references.bib|genage-references.bib]]

---

## Experiments Using This Dataset

- [[../experiments/genage-analysis-2025-12-02|GenAge Analysis (2025-12-02)]]

---

## Next Steps

1. **Load and Explore**:
   - [ ] Install openpyxl: `pip install openpyxl`
   - [ ] Load Excel file and inspect sheets
   - [ ] Document actual schema

2. **Analysis**:
   - [ ] Identify top aging signatures
   - [ ] Cross-reference with GenAge genes
   - [ ] Pathway enrichment analysis

3. **Integration**:
   - [ ] Train IMMUNOS-MCP agents on patterns
   - [ ] Cross-species comparison
   - [ ] Immunosenescence focus

---

## Metadata Summary

- **Datasets**: 127 expression studies
- **Species**: Human, mouse, rat
- **Technologies**: Microarray, RNA-Seq
- **Downloaded**: 2025-12-02
- **Size**: 53 KB
- **Format**: Excel (.xlsx)
- **License**: CC BY 3.0

---

**Last Updated**: 2025-12-02
**Status**: Downloaded (detailed analysis pending)
**Note**: Requires openpyxl for full analysis
