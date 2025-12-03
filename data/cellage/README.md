# CellAge Database - Cell Senescence Genes

**Source**: Human Ageing Genomic Resources (HAGR)
**Website**: https://genomics.senescence.info/cells/
**Downloaded**: 2025-12-02
**Build**: 3 (released April 22, 2023)

## Overview

CellAge is a database of genes associated with cellular senescence - the state of permanent cell cycle arrest that accumulates with age and contributes to aging phenotypes.

## License

**Creative Commons Attribution 3.0 Unported License**
- Free for commercial, educational, and research use
- Attribution required in publications
- Full data export and reuse permitted

## Datasets

### 1. CellAge Database (`cellage3.tsv`)

**Description**: Curated database of 950 genes associated with cellular senescence
**Records**: 950 genes (Build 3)
**Format**: Tab-delimited (TSV)
**File Size**: ~72 KB

**Columns**:
- **Entrez ID**: NCBI Gene database identifier
- **Gene symbol**: Standard gene name
- **Gene name**: Full descriptive name
- **Cancer Cell**: Whether senescence observed in cancer cells (Yes/No)
- **Type of senescence**: Category (Stress-induced, Oncogene-induced, Replicative, etc.)
- **Senescence Effect**: Whether gene induces or inhibits senescence
- **Reference**: PubMed ID for literature citation

**Sample Data**:
```
Entrez ID  Gene symbol  Gene name                          Cancer  Type              Effect    PubMed
22848      AAK1         AP2 associated kinase 1            No      Unclear           Induces   26583757
5243       ABCB1        ATP binding cassette subfamily B1  Yes     Stress-induced    Induces   10825123
```

### 2. Cell Senescence Gene Expression Signatures (`signatures1.csv`)

**Description**: Meta-analysis of gene expression changes in cellular senescence
**Records**: 1,259 genes
**Format**: CSV (semicolon-delimited)
**File Size**: ~96 KB

**Columns**:
- **gene_symbol**: Gene name
- **gene_name**: Full descriptive name
- **entrez_id**: NCBI Gene ID
- **total**: Total number of datasets showing differential expression
- **ovevrexp**: Number of datasets with overexpression
- **underexp**: Number of datasets with underexpression
- **p_value**: Statistical significance

**Sample Data**:
```
gene_symbol  gene_name                               entrez_id  total  overexp  underexp  p_value
ABCA3        ATP binding cassette subfamily A3       21         17     1        0         0.004279
ACTL6A       actin like 6A                           86         19     0        1         0.00000024
```

## Usage Examples

### Python - Load CellAge Database

```python
import csv

# Load CellAge genes
with open('/Users/byron/projects/data/cellage/cellage3.tsv', 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    cellage = list(reader)

print(f"Total senescence genes: {len(cellage)}")

# Filter by senescence effect
inducers = [g for g in cellage if g['Senescence Effect'] == 'Induces']
inhibitors = [g for g in cellage if g['Senescence Effect'] == 'Inhibits']

print(f"Genes that induce senescence: {len(inducers)}")
print(f"Genes that inhibit senescence: {len(inhibitors)}")
```

### Python - Load Expression Signatures

```python
import csv

# Load expression signatures (semicolon-delimited)
with open('/Users/byron/projects/data/cellage/signatures1.csv', 'r') as f:
    reader = csv.DictReader(f, delimiter=';')
    signatures = list(reader)

print(f"Total genes with expression signatures: {len(signatures)}")

# Find consistently overexpressed genes
overexpressed = [g for g in signatures if int(g['ovevrexp']) > 10]
print(f"Genes overexpressed in >10 datasets: {len(overexpressed)}")
```

## Research Applications

### 1. Aging-Senescence Connection
- Cross-reference with GenAge aging genes
- Identify overlap between aging and senescence pathways

### 2. Immunosenescence Research
- Study immune cell senescence with age
- Identify senescence markers in immune cells

### 3. Therapeutic Targets
- Senolytic drug development (drugs that kill senescent cells)
- Senostatic approaches (inhibit senescence progression)

## Citation

**Primary Publication**:
```
Avelar, R. A., Ortega, J. G., Tacutu, R., Tyler, E. J., Bennett, D., Binetti, P.,
Budovsky, A., Chatsirisupachai, K., Johnson, E., Murray, A., et al. (2020).
A multidimensional systems biology analysis of cellular senescence in aging and disease.
Genome Biology, 21(1), 91. https://doi.org/10.1186/s13059-020-01990-9
```

**Data Citation**:
```
Human Ageing Genomic Resources (HAGR). (2023). CellAge Database. Build 3.
Retrieved from https://genomics.senescence.info/cells/
```

## Related Datasets

- **GenAge Human**: [[../genage/human/|Human aging genes]] (307 genes)
- **GenAge Models**: [[../genage/models/|Model organism longevity genes]] (2,205 genes)
- **Gene Expression**: [[../genage/expression/|Aging expression signatures]]

## Integration Opportunities

### Cross-Reference Analysis
```python
# Find genes in both CellAge and GenAge
import csv

# Load both datasets
with open('data/genage/human/genage_human.csv', 'r') as f:
    genage_genes = {row['symbol'] for row in csv.DictReader(f)}

with open('data/cellage/cellage3.tsv', 'r') as f:
    cellage_genes = {row['Gene symbol'] for row in csv.DictReader(f, delimiter='\t')}

overlap = genage_genes.intersection(cellage_genes)
print(f"Genes in both GenAge and CellAge: {len(overlap)}")
```

---

**Last Updated**: 2025-12-02
**Build**: 3 (April 2023)
**Total Genes**: 950 (CellAge) + 1,259 (Signatures)
