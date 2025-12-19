---
name: "RegulomeDB"
type: "database"
category: "regulatory-variants"
maintained_by: "Stanford University"
website: "https://www.regulomedb.org/"
license: "Open Access"
status: "active"
last_release: "2024"
data_size: "~10 GB (web-based)"
added_date: "2025-12-05"
tags:
  - database
  - regulatory-variants
  - gene-regulation
  - snp
  - gwas
---

# RegulomeDB

## Quick Summary

Database for annotating regulatory variants in the human genome. Integrates ENCODE data, ChIP-seq, DNase-seq, and eQTL data to score variants based on regulatory potential.

## Overview

**Purpose**: Predict regulatory consequences of DNA variants (SNPs, indels)

**Key Features**:
- Variant scoring (1-7 scale, 1 = highest regulatory evidence)
- ENCODE integration
- eQTL mapping
- ChIP-seq/DNase-seq data
- GWAS variant annotation

**Data**:
- Regulatory elements from ENCODE
- Expression QTLs (eQTLs)
- Chromatin states
- Transcription factor binding sites

## Access

**Web Interface**: https://www.regulomedb.org/regulome-search
**Search**: By rsID, genomic position, or gene

**Example Query**:
```
# Search SNP
rsID: rs7412  # APOE variant

# Result: Regulatory score + evidence
```

## Relevance to Aging Biology

**Use Cases**:
- Annotate aging GWAS variants
- Map GenAge gene regulatory regions
- Identify functional aging SNPs
- eQTL analysis of longevity genes

**Integration**:
- Cross-reference GenAge genes with regulatory variants
- Identify causal variants in aging pathways
- Prioritize SNPs for functional studies

## Citation

```bibtex
@article{boyle2012regulomedb,
  title={Annotation of functional variation in personal genomes using RegulomeDB},
  author={Boyle, Alan P and others},
  journal={Genome Research},
  volume={22},
  pages={1790--1797},
  year={2012}
}
```

## Priority

**Moderate** - Useful for variant annotation if pursuing GWAS/genomics aging research

---

**Created**: 2025-12-05
**Status**: Documented - Web portal access only
