# Research Tools and Software

This directory documents software tools, platforms, and computational methods used in aging biology research.

## Purpose

Track and document:
- Bioinformatics tools and pipelines
- Statistical analysis software
- Machine learning frameworks
- Visualization platforms
- Lab automation software
- Research management tools

## Tool Categories

### Bioinformatics
- **Sequence analysis** - BLAST, Clustal, MAFFT
- **Variant calling** - GATK, bcftools, VCFtools
- **RNA-seq analysis** - DESeq2, edgeR, STAR
- **Pathway analysis** - GSEA, Reactome, KEGG

### Machine Learning
- **Frameworks** - TensorFlow, PyTorch, scikit-learn
- **Specialized** - AlphaFold, ESM (protein language models)
- **AI benchmarks** - LAB-Bench (biology Q&A evaluation)

### Data Analysis
- **Statistical** - R, Python (pandas, numpy, scipy)
- **Notebooks** - Jupyter, Observable, Quarto
- **Visualization** - matplotlib, seaborn, Cytoscape

### Research Management
- **Literature** - Zotero, Mendeley, Papers
- **Notes** - Obsidian, Notion, Roam Research
- **Project management** - GitHub Projects, Asana

### Lab Software
- **Protocols** - protocols.io, Benchling
- **Data management** - LIMS systems
- **Collaboration** - Slack, Microsoft Teams

## File Naming Convention

Use lowercase with hyphens: `tool-name.md`

Examples:
- `lab-bench-benchmark.md`
- `alphafold-protein-structure.md`
- `deseq2-rnaseq.md`
- `zotero-citations.md`

## Profile Template

Use the standard template at `/templates/database-tool-profile.md` which includes:
- YAML frontmatter with metadata
- Tool overview and purpose
- Installation and setup
- Usage examples and tutorials
- Integration with our workflow
- Performance benchmarks (if applicable)
- Licensing and cost
- Community and support

## Integration Examples

### LAB-Bench (AI Biology Benchmark)
```python
from datasets import load_dataset

# Load LAB-Bench dataset from HuggingFace
dataset = load_dataset("futurehouse/lab-bench", "FigQA")
print(f"Total questions: {len(dataset['train'])}")
```

### AlphaFold Integration
```python
# Predict protein structure
from alphafold import run_alphafold
structure = run_alphafold(sequence="MKTV...", output_dir="predictions/")
```

### DESeq2 for RNA-seq
```R
# Differential expression analysis
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData, colData, design = ~ condition)
dds <- DESeq(dds)
results <- results(dds, contrast=c("condition", "old", "young"))
```

## Example Queries

### Dataview: List all ML tools
```dataview
TABLE category, language, last_used
FROM "research/tools"
WHERE category = "machine-learning"
```

### Dataview: Find tools by programming language
```dataview
LIST purpose, url
FROM "research/tools"
WHERE language = "Python"
```

## Current Tools

- **LAB-Bench** - AI evaluation benchmark for biology (HuggingFace dataset)
- **Obsidian** - Knowledge management (this system!)
- **Zotero + Better BibTeX** - Citation management
- *(more to be added)*

---

**Directory**: `/research/tools/`
**Template**: `/templates/database-tool-profile.md`
**Script**: *N/A (manual creation)*
**Last Updated**: 2025-12-03
