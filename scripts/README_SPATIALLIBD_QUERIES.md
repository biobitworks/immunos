# spatialLIBD Query System for Aging Biology

Query brain spatial transcriptomics data for aging-related genes **without downloading 2GB dataset**.

## Quick Start

### Python

```bash
cd /Users/byron/projects
python scripts/query_spatiallibd.py
```

Or in Python:

```python
from scripts.query_spatiallibd import SpatialLIBDQuery

query = SpatialLIBDQuery()
query.list_questions()
query.ask_question("q1")  # Layer-specific aging genes
```

### R

```bash
cd /Users/byron/projects
Rscript scripts/query_spatiallibd.R
```

Or in R:

```r
source("scripts/query_spatiallibd.R")
list_aging_questions()
run_aging_question("q1")  # Layer-specific aging genes
```

## Predefined Questions

| ID | Question | Genes |
|----|----------|-------|
| **q1** | Are aging genes layer-specific in brain? | TP53, FOXO3, IGF1R, SIRT1 |
| **q2** | Where are senescence genes expressed in DLPFC? | CDKN2A, TP53, RB1, CDKN1A, BCL6 |
| **q3** | Do longevity genes show distinct patterns? | FOXO3, APOE, SIRT1, TERT |
| **q4** | Are immune genes enriched in specific layers? | IL6, TNF, NFKB1, CD68 |
| **q5** | How is BCL6 (senescence inhibitor) distributed? | BCL6 |

## Custom Queries

### Python
```python
# Query custom genes
query.query_genes(["TP53", "CDKN2A", "BCL6"])

# Query top 20 GenAge genes
query.query_genage_top(n=20)

# Query CellAge senescence genes
query.query_cellage_sample(n=15)
```

### R
```r
# Query custom genes
query_custom_genes(c("TP53", "CDKN2A", "BCL6"))

# Load GenAge and query
genage <- read.csv("data/genage/human/genage_human.csv")
top_genes <- genage$symbol[1:20]
query_custom_genes(top_genes)
```

## How It Works

1. **Loads your aging/senescence genes** from GenAge and CellAge datasets
2. **Generates web URLs** for each gene in spatialLIBD web interface
3. **Opens in browser** (optional) to visualize spatial expression
4. **No download required** - uses web interface, avoids 2GB data download

## What You Can See

For each gene, the web interface shows:
- **Spatial expression map**: Where gene is expressed in tissue
- **Layer enrichment**: Which cortical layers (L1-L6, WM) show highest expression
- **Multi-sample view**: Expression across 3 subjects
- **Statistics**: Mean expression, detection rate, layer-specific values
- **Export**: High-resolution figures for publications

## Benefits

✅ **No download** - Avoids 2GB spatialLIBD dataset download
✅ **Fast** - Query multiple genes in seconds
✅ **Integrated** - Automatically loads GenAge/CellAge
✅ **Extensible** - Easy to add new questions
✅ **Interactive** - Opens results in web browser

## Example Output

```
======================================================================
QUESTION q1
======================================================================

Are aging genes layer-specific in the brain?

Description: Check if GenAge aging-related genes show preferential
             expression in specific cortical layers
Genes: TP53, FOXO3, IGF1R, SIRT1
----------------------------------------------------------------------

spatialLIBD Aging Biology Query
======================================================================

Querying 4 genes in brain spatial transcriptomics data

Method: Web interface (no download required) ✓

  • TP53       → http://spatial.libd.org/spatialLIBD/?gene=TP53
  • FOXO3      → http://spatial.libd.org/spatialLIBD/?gene=FOXO3
  • IGF1R      → http://spatial.libd.org/spatialLIBD/?gene=IGF1R
  • SIRT1      → http://spatial.libd.org/spatialLIBD/?gene=SIRT1

----------------------------------------------------------------------
What you can explore in web interface:
  - Spatial expression patterns across tissue
  - Layer-specific enrichment (L1-L6, WM)
  - Expression in different subjects
  - Export high-resolution figures
----------------------------------------------------------------------

Opening first gene in browser...
```

## Files

- **Python**: `scripts/query_spatiallibd.py` (340 lines)
- **R**: `scripts/query_spatiallibd.R` (370 lines)
- **Documentation**: `research/databases/spatiallibd.md`
- **This README**: `scripts/README_SPATIALLIBD_QUERIES.md`

## Related Documentation

- [[../research/databases/spatiallibd|spatialLIBD Database Documentation]]
- [[../research/datasets/genage|GenAge Dataset]]
- [[../research/datasets/cellage|CellAge Dataset]]
- [[../research/literature/shvarts-2002-bcl6-senescence|BCL6 Senescence Paper]]

---

**Created**: 2025-12-03
**Purpose**: Query brain spatial transcriptomics without storage cost
**Storage saved**: 2 GB (vs. downloading full dataset)
