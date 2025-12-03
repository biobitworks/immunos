#!/usr/bin/env python3
"""
Analyze GenAge Datasets

Analyzes human genes, model organism genes, and gene expression signatures
to generate summary statistics and insights.

Usage:
    python3 scripts/analyze-genage.py
"""

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from datetime import datetime


def load_csv(file_path: Path, delimiter=',') -> list[dict]:
    """Load CSV file and return list of dictionaries."""
    print(f"Loading: {file_path.name}")
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            data = list(reader)
        print(f"✓ Loaded {len(data)} records")
        return data
    except Exception as e:
        print(f"✗ Error: {e}")
        return []


def analyze_human_genes(data_path: Path) -> dict:
    """Analyze GenAge human genes dataset."""
    print(f"\n{'='*70}")
    print("Analyzing GenAge Human Genes")
    print(f"{'='*70}\n")

    csv_path = data_path / "genage_human.csv"
    data = load_csv(csv_path)

    if not data:
        return {}

    # Basic statistics
    stats = {
        "total_genes": len(data),
        "columns": list(data[0].keys()) if data else [],
        "sample_genes": [row.get('symbol', row.get('gene symbol', '')) for row in data[:5]]
    }

    # Print summary
    print(f"\nDataset Statistics:")
    print(f"  Total genes: {stats['total_genes']}")
    print(f"  Columns: {len(stats['columns'])}")
    print(f"\nColumn names:")
    for col in stats['columns']:
        print(f"  - {col}")

    print(f"\nSample genes (first 5):")
    for gene in stats['sample_genes']:
        print(f"  - {gene}")

    return stats


def analyze_model_organisms(data_path: Path) -> dict:
    """Analyze GenAge model organisms dataset."""
    print(f"\n{'='*70}")
    print("Analyzing GenAge Model Organisms")
    print(f"{'='*70}\n")

    csv_path = data_path / "genage_models.csv"
    data = load_csv(csv_path)

    if not data:
        return {}

    # Find organism column (might be 'organism', 'species', or similar)
    organism_col = None
    for col in data[0].keys():
        if 'organism' in col.lower() or 'species' in col.lower():
            organism_col = col
            break

    # Count organisms
    organism_counts = Counter()
    if organism_col:
        for row in data:
            org = row.get(organism_col, 'Unknown')
            if org:
                organism_counts[org] += 1

    # Basic statistics
    stats = {
        "total_genes": len(data),
        "columns": list(data[0].keys()) if data else [],
        "organism_column": organism_col,
        "num_organisms": len(organism_counts),
        "organisms": dict(organism_counts.most_common())
    }

    # Print summary
    print(f"\nDataset Statistics:")
    print(f"  Total genes: {stats['total_genes']}")
    print(f"  Columns: {len(stats['columns'])}")
    print(f"  Organisms: {stats['num_organisms']}")

    print(f"\nColumn names:")
    for col in stats['columns']:
        print(f"  - {col}")

    print(f"\nGenes per organism:")
    for org, count in organism_counts.most_common(10):
        print(f"  - {org}: {count} genes")
    if len(organism_counts) > 10:
        print(f"  ... and {len(organism_counts) - 10} more organisms")

    return stats


def analyze_expression_signatures(data_path: Path) -> dict:
    """Analyze gene expression signatures dataset."""
    print(f"\n{'='*70}")
    print("Analyzing Gene Expression Signatures")
    print(f"{'='*70}\n")

    xlsx_path = data_path / "ageing_signatures.xlsx"

    # For Excel files, we need openpyxl or pandas, which may not be installed
    # For now, just report file information
    print(f"File: {xlsx_path.name}")
    print(f"Size: {xlsx_path.stat().st_size:,} bytes")
    print(f"\nNote: This is an Excel file (.xlsx)")
    print(f"To analyze in detail, install: pip install openpyxl pandas")
    print(f"Or open in: Excel, LibreOffice, or Google Sheets")

    stats = {
        "file_name": xlsx_path.name,
        "file_size_bytes": xlsx_path.stat().st_size,
        "format": "Excel (.xlsx)",
        "note": "Contains meta-analysis results from 127 datasets"
    }

    return stats


def create_summary_report(human_stats: dict, models_stats: dict, expr_stats: dict) -> str:
    """Create markdown summary report."""
    report = f"""---
title: GenAge Data Analysis Summary
date: {datetime.now().strftime('%Y-%m-%d')}
tags: [genage, aging, data-analysis, genomics]
type: analysis-report
---

# GenAge Data Analysis Summary

**Analysis Date**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
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
**Total Genes**: {human_stats.get('total_genes', 'N/A')}

### Schema

**Columns** ({len(human_stats.get('columns', []))} total):
"""

    for col in human_stats.get('columns', []):
        report += f"- `{col}`\n"

    report += f"""
### Sample Genes

**First 5 genes**:
"""

    for gene in human_stats.get('sample_genes', []):
        report += f"- {gene}\n"

    report += f"""

---

## 2. GenAge Model Organisms

**File**: `data/genage/models/genage_models.csv`
**Total Genes**: {models_stats.get('total_genes', 'N/A')}
**Organisms**: {models_stats.get('num_organisms', 'N/A')}

### Schema

**Columns** ({len(models_stats.get('columns', []))} total):
"""

    for col in models_stats.get('columns', []):
        report += f"- `{col}`\n"

    report += f"""
### Organisms Distribution

**Top 10 organisms by gene count**:

| Organism | Genes |
|----------|-------|
"""

    organisms = models_stats.get('organisms', {})
    for i, (org, count) in enumerate(sorted(organisms.items(), key=lambda x: x[1], reverse=True)[:10]):
        report += f"| {org} | {count} |\n"

    if len(organisms) > 10:
        report += f"\n*... and {len(organisms) - 10} more organisms*\n"

    report += f"""

---

## 3. Gene Expression Signatures

**File**: `data/genage/expression/ageing_signatures.xlsx`
**Format**: {expr_stats.get('format', 'N/A')}
**Size**: {expr_stats.get('file_size_bytes', 0):,} bytes

**Description**: Meta-analysis results from 127 microarray and RNA-Seq datasets across mammals (human, mouse, rat).

**Note**: {expr_stats.get('note', 'N/A')}

---

## Summary Statistics

| Dataset | Records | File Size | Format |
|---------|---------|-----------|--------|
| Human Genes | {human_stats.get('total_genes', 'N/A')} | {Path('/Users/byron/projects/data/genage/human/genage_human.csv').stat().st_size:,} bytes | CSV |
| Model Organisms | {models_stats.get('total_genes', 'N/A')} | {Path('/Users/byron/projects/data/genage/models/genage_models.csv').stat().st_size:,} bytes | CSV |
| Expression Signatures | 127 datasets | {expr_stats.get('file_size_bytes', 0):,} bytes | Excel |

**Total Genes**: {human_stats.get('total_genes', 0) + models_stats.get('total_genes', 0):,}

---

## Key Insights

### 1. Comprehensive Coverage
- {human_stats.get('total_genes', 0)} manually curated human aging genes
- {models_stats.get('total_genes', 0)} genes across {models_stats.get('num_organisms', 0)} model organisms
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
print(f"Human genes: {{len(human_df)}}")

# Load model organisms
models_df = pd.read_csv('/Users/byron/projects/data/genage/models/genage_models.csv')
print(f"Model organism genes: {{len(models_df)}}")

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

**Last Updated**: {datetime.now().strftime('%Y-%m-%d')}
**Data Source**: HAGR GenAge Build 21 (August 2023)
**License**: CC BY 3.0
"""

    return report


def main():
    """Main analysis function."""
    print("=" * 70)
    print("GenAge Data Analysis")
    print("=" * 70)

    data_dir = Path("/Users/byron/projects/data/genage")

    # Analyze each dataset
    human_stats = analyze_human_genes(data_dir / "human")
    models_stats = analyze_model_organisms(data_dir / "models")
    expr_stats = analyze_expression_signatures(data_dir / "expression")

    # Create summary report
    print(f"\n{'='*70}")
    print("Generating Summary Report")
    print(f"{'='*70}\n")

    report = create_summary_report(human_stats, models_stats, expr_stats)
    report_path = Path("/Users/byron/projects/research/experiments/genage-analysis-2025-12-02.md")
    report_path.parent.mkdir(parents=True, exist_ok=True)

    with open(report_path, 'w') as f:
        f.write(report)

    print(f"✓ Report saved: {report_path}")
    print(f"\nView report in Obsidian:")
    print(f"  {report_path.relative_to(Path.cwd())}")
    print()


if __name__ == "__main__":
    main()
