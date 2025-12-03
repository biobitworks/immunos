---
name: "spatialLIBD"
type: "database"
category: "spatial-transcriptomics"
maintained_by: "Lieber Institute for Brain Development (LIBD)"
website: "https://research.libd.org/spatialLIBD/"
repository: "https://github.com/LieberInstitute/spatialLIBD"
license: "Artistic-2.0"
programming_language: "R"
installation: "bioconductor"
status: "active"
last_release: "2024"
documentation_quality: "excellent"
data_size: "~2 GB (spot-level data)"
added_date: "2025-12-03"
last_updated: "2025-12-03"
tags:
  - database
  - spatial-transcriptomics
  - brain
  - neuroscience
  - bioconductor
  - 10x-visium
---

# spatialLIBD

## Quick Summary

R/Bioconductor package and web application for exploring spatial transcriptomics data from human dorsolateral prefrontal cortex (DLPFC). Provides ~2 GB dataset of 33,538 genes across 47,681 spatial spots from brain tissue using 10x Genomics Visium platform.

## Overview

### Purpose and Scope
spatialLIBD enables visualization and analysis of spatially-resolved transcriptomics data, specifically focusing on layer-specific gene expression patterns in the human brain. The package provides both interactive web tools and programmatic R access to the data.

### Key Features
- **Interactive web application**: Browse spatial gene expression at http://spatial.libd.org/spatialLIBD/
- **R/Bioconductor integration**: Programmatic data access and analysis
- **Layer annotation**: Manual annotation of cortical layers and white matter
- **Gene visualization**: Spatial expression patterns across tissue sections
- **Quality control**: Built-in QC metrics and filtering
- **Multi-sample**: Data from 3 subjects with replicates

### Version Information
- **Current Version**: Available through Bioconductor (updated regularly)
- **Data Source**: Maynard, Collado-Torres et al., Nature Neuroscience, 2021
- **Platform**: 10x Genomics Visium Spatial Gene Expression
- **Status**: Actively maintained

## Technical Details

### Data Statistics

**⚠️ DATA SIZE TRACKING - DO NOT DOWNLOAD WITHOUT PROJECT JUSTIFICATION**

| Component | Size | Records | Notes |
|-----------|------|---------|-------|
| **Spot-level data** | **~2.04 GB** | 47,681 spots | Main dataset - LARGE |
| Genes | - | 33,538 genes | Full transcriptome |
| Subjects | - | 3 individuals | 2 males, 1 female, ages 30-46 |
| Images | - | 12 images | 4 per subject (spatial replicates) |

**Storage recommendation**: Use web interface or stream specific genes rather than downloading full dataset.

### Data Schema

#### Tissue Coverage
**Brain Region**: Dorsolateral Prefrontal Cortex (DLPFC)
**Layers Covered**:
- Layer 1 (molecular layer)
- Layer 2 (external granular)
- Layer 3 (external pyramidal)
- Layer 4 (internal granular)
- Layer 5 (internal pyramidal)
- Layer 6 (multiform layer)
- White Matter (WM)

#### Experimental Design
- **3 control subjects**: All neurotypical
- **4 tissue sections per subject**:
  - 2 spatially adjacent replicates at position 0
  - 2 spatially adjacent replicates at +300 μm
- **Total**: 12 tissue sections

#### Data Structure (SpatialExperiment object)
```r
SpatialExperiment
├── Assays: counts, logcounts
├── Row Data: gene annotations (33,538 genes)
├── Column Data: spot metadata (47,681 spots)
│   ├── spatial coordinates (x, y)
│   ├── layer annotations
│   ├── sample info
│   └── QC metrics
└── Spatial Data: images and coordinates
```

### File Formats
- **R object**: SpatialExperiment (Bioconductor class)
- **Images**: High-resolution H&E tissue images
- **Coordinates**: Spatial position of each spot
- **Expression**: Count matrices (raw and normalized)

## Access Methods

### Web Interface (Recommended for Exploration)
**URL**: http://spatial.libd.org/spatialLIBD/

**Features**:
- Browse gene expression spatially
- View layer annotations
- Compare across samples
- Export visualizations
- NO download required ✅

**Use cases**:
- Quick gene lookups
- Visual exploration
- Generating figures
- Validating hypotheses

### R/Bioconductor Package

**⚠️ WARNING: Full download is ~2 GB. Use fetch_data() only if needed for specific project.**

**Installation**:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("spatialLIBD")
```

**Load pre-computed data** (downloads ~2 GB):
```r
library("spatialLIBD")

# Download full dataset - ONLY if project requires!
spe <- fetch_data(type = "spe")

# Alternative: Access specific genes via web API (smaller)
# [API methods to be documented if needed]
```

**Visualize specific gene** (after loading):
```r
# Spatial expression plot
vis_gene(
    spe = spe,
    geneid = "ENSG00000000003",  # Gene ID
    spatial = TRUE
)

# Layer enrichment
layer_stat_cor_plot(spe, geneid = "ENSG00000000003")
```

### Alternative: Query-based Access

**For space efficiency**: Use web interface or develop targeted queries rather than downloading full dataset.

```r
# Hypothetical targeted query (to be developed)
# Query specific genes without full download
get_spatial_expression(
    genes = c("SNAP25", "MOBP", "PCP4"),
    samples = c("151507", "151508"),
    return_format = "data.frame"
)
```

## Integration with Our Research

### Relevance to Aging Biology

#### Brain Aging Context
**Potential relevance**: DLPFC is affected in aging and neurodegenerative diseases.

**Aging-related questions**:
- Do GenAge genes show layer-specific expression in brain?
- Are senescence markers (CellAge) enriched in specific layers?
- Spatial patterns of aging-related genes?

#### Cross-Reference Opportunities

**GenAge genes in brain**:
```r
# Load GenAge human genes
genage <- read.csv('/Users/byron/projects/data/genage/human/genage_human.csv')

# IF we download spatialLIBD (project-dependent):
# Check which GenAge genes are expressed in DLPFC
# aging_genes_spatial <- subset(rowData(spe), gene_name %in% genage$symbol)
```

**CellAge senescence genes**:
- Are cellular senescence markers spatially localized?
- Do senescent cells accumulate in specific cortical layers?
- Connection to brain aging phenotypes?

### Project Considerations

**⚠️ Download Decision Matrix**:

| Project Type | Download? | Rationale |
|--------------|-----------|-----------|
| **Quick gene lookup** | ❌ NO | Use web interface |
| **Exploratory analysis** | ❌ NO | Use web interface |
| **Aging gene enrichment** | ⚠️ MAYBE | Small gene list query via web |
| **Layer-specific aging** | ✅ YES | Requires full data analysis |
| **Spatial aging patterns** | ✅ YES | Requires spatial coordinates |
| **Integration with IMMUNOS** | ❌ NO | Not directly relevant |

### Current Priority: **LOW**
**Rationale**:
- Not core aging database (GenAge/CellAge)
- 2 GB size vs. limited aging-specific utility
- Web interface sufficient for initial exploration
- Brain-specific (vs. general aging biology)

**Future potential**:
- If developing brain aging project → revisit
- If studying tissue-specific senescence → useful
- If investigating spatial aging patterns → valuable

## Data Provenance

### Original Study
**Paper**: Maynard, K.R., Collado-Torres, L., et al. (2021)
**Title**: "Transcriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex"
**Journal**: Nature Neuroscience
**DOI**: 10.1038/s41593-020-00787-0
**PMID**: 33558695

### Key Findings from Paper
- Identified layer-specific marker genes
- Discovered layer-enriched gene expression patterns
- Validated with Immunofluorescence
- Created reference dataset for spatial transcriptomics

### Subjects
- **N = 3** neurotypical control subjects
- **Demographics**: 2 males, 1 female
- **Age range**: 30-46 years
- **Tissue**: Postmortem DLPFC

### Technology
- **Platform**: 10x Genomics Visium Spatial Gene Expression
- **Resolution**: ~55 μm spot diameter
- **Coverage**: ~10,000 genes detected per spot
- **Depth**: ~5,000-10,000 reads per spot

## Citation and Attribution

### Citing the Dataset

**BibTeX**:
```bibtex
@article{maynard2021spatial,
  author = {Maynard, Kristen R. and Collado-Torres, Leonardo and others},
  title = {Transcriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex},
  journal = {Nature Neuroscience},
  volume = {24},
  pages = {425--436},
  year = {2021},
  doi = {10.1038/s41593-020-00787-0}
}
```

### Citing the Software

```bibtex
@article{spatiallibd2022,
  author = {Pardo, Brenda and others},
  title = {spatialLIBD: an R/Bioconductor package to visualize spatially-resolved transcriptomics data},
  journal = {BMC Genomics},
  volume = {23},
  pages = {434},
  year = {2022},
  doi = {10.1186/s12864-022-08601-w}
}
```

## Performance and Limitations

### Advantages
- High-quality manual layer annotations
- Interactive web interface (no download needed)
- Well-documented R package
- Established reference dataset
- Active maintenance

### Limitations
- **Large dataset size** (~2 GB) - storage considerations
- **Brain-specific** - not generalizable to other tissues
- **Limited subjects** (N=3) - sample size constraints
- **Neurotypical only** - no disease comparisons in base dataset
- **Specific region** (DLPFC) - not whole brain
- **Adult ages** (30-46) - limited aging time course

### Storage Impact Assessment

**If downloaded to local system**:
- **Size**: 2.04 GB
- **Impact**: Significant (20% of current data/ directory)
- **Benefit**: Brain-specific spatial analysis
- **Alternative**: Web interface + targeted queries

**Current data/ directory**:
- GenAge: ~250 KB
- CellAge: ~170 KB
- **Total current**: ~420 KB
- **With spatialLIBD**: ~2.04 GB (4,857x increase)

**Recommendation**: DO NOT download unless specific brain aging project justifies the storage cost.

## Documentation and Support

### Official Documentation
- **Website**: https://research.libd.org/spatialLIBD/
- **Vignettes**: http://research.libd.org/spatialLIBD/articles/
- **Bioconductor page**: https://www.bioconductor.org/packages/release/data/experiment/html/spatialLIBD.html
- **Paper**: https://doi.org/10.1186/s12864-022-08601-w

### Code Repository
- **GitHub**: https://github.com/LieberInstitute/spatialLIBD
- **Issues**: https://github.com/LieberInstitute/spatialLIBD/issues
- **Contributions**: Active development

### Support Channels
- **Bioconductor Support**: https://support.bioconductor.org/
- **GitHub Issues**: For bug reports and feature requests
- **LIBD Website**: Contact forms available

## Query System (NEW!)

### Automated Query Scripts

**⭐ We created custom query scripts to ask aging biology questions without downloading 2GB dataset!**

**Location**:
- Python: `/Users/byron/projects/scripts/query_spatiallibd.py`
- R: `/Users/byron/projects/scripts/query_spatiallibd.R`

### Quick Start: Python

```python
from scripts.query_spatiallibd import SpatialLIBDQuery

# Initialize query system
query = SpatialLIBDQuery()

# List available questions
query.list_questions()

# Ask predefined question
query.ask_question("q1")  # Are aging genes layer-specific?

# Custom gene query
query.query_genes(["TP53", "FOXO3", "BCL6"])

# Query top GenAge genes
query.query_genage_top(n=10)
```

### Quick Start: R

```r
# Load query functions
source("scripts/query_spatiallibd.R")

# List available questions
list_aging_questions()

# Ask predefined question
run_aging_question("q1")  # Are aging genes layer-specific?

# Custom gene query
query_custom_genes(c("TP53", "FOXO3", "BCL6"))
```

### Predefined Aging Biology Questions

**q1**: Are aging genes layer-specific in the brain?
- Genes: TP53, FOXO3, IGF1R, SIRT1

**q2**: Where are cellular senescence genes expressed in DLPFC?
- Genes: CDKN2A, TP53, RB1, CDKN1A, BCL6

**q3**: Do longevity genes show distinct spatial patterns?
- Genes: FOXO3, APOE, SIRT1, TERT

**q4**: Are immune/inflammation genes enriched in specific layers?
- Genes: IL6, TNF, NFKB1, CD68

**q5**: How is BCL6 (senescence inhibitor) distributed in brain?
- Genes: BCL6
- Follows up on [[../literature/shvarts-2002-bcl6-senescence|Shvarts 2002 paper]]

### Benefits of Query Scripts

✅ **No download** - Uses web interface (avoids 2GB download)
✅ **Automated** - Generate URLs for multiple genes
✅ **Integrated** - Loads GenAge and CellAge genes automatically
✅ **Documented** - Clear questions and descriptions
✅ **Extensible** - Easy to add new questions

## Usage Examples

### Example 1: Web Interface Exploration (Recommended)

**Goal**: Check if specific aging-related gene is layer-enriched

**Steps**:
1. Visit http://spatial.libd.org/spatialLIBD/
2. Enter gene name (e.g., "TP53", "FOXO3")
3. View spatial expression pattern
4. Check layer enrichment statistics
5. Export figure if needed

**No download required** ✅
**Time**: < 5 minutes

### Example 1b: Using Query Script (Even Faster!)

**Goal**: Check multiple aging genes at once

```bash
# Python
python scripts/query_spatiallibd.py

# Then in Python:
query.ask_question("q1")  # Opens URLs for aging genes
```

**Time**: < 1 minute ⚡

### Example 2: Programmatic Access (If Downloaded)

**Goal**: Analyze multiple GenAge genes for spatial patterns

```r
library(spatialLIBD)

# WARNING: Downloads ~2 GB
spe <- fetch_data(type = "spe")

# Load aging genes
genage_genes <- c("TP53", "FOXO3", "SIRT1", "IGF1R")

# Visualize each gene
for (gene in genage_genes) {
    vis_gene(spe, gene, spatial = TRUE)
}
```

### Example 3: Targeted Query (Space-Efficient)

**Goal**: Get expression data for specific genes without full download

```r
# Use web API or subset functionality
# (Specific implementation TBD based on project needs)

# Pseudo-code for efficient access:
# results <- query_spatiallibd_api(
#     genes = genage_top10,
#     samples = "all",
#     return_summary = TRUE
# )
```

## Related Resources

### Similar Spatial Transcriptomics Databases
- **Human Cell Atlas** - Multi-tissue spatial data
- **Allen Brain Atlas** - Extensive brain gene expression
- **10x Genomics datasets** - Public Visium data repository

### Our Related Databases
- [[genage|GenAge]] - Check if aging genes show spatial patterns
- [[cellage|CellAge]] - Senescence markers in brain tissue
- [[monarch-kg|Monarch Initiative]] - Brain phenotypes and disease

## Notes

### Strategic Importance
**Moderate** for aging biology research:
- Brain-specific, not general aging
- Valuable IF pursuing brain aging specifically
- Excellent data quality and tools
- But large storage cost (2 GB vs 420 KB current)

### Storage Management
**Key principle**: Track data size relative to project utility

**Current strategy**:
- ✅ Lightweight databases: GenAge (307 genes, 250 KB), CellAge (950 genes, 170 KB)
- ⚠️ Medium databases: Consider on project basis
- ❌ Large databases (>1 GB): Require clear justification

**spatialLIBD classification**: Large database
**Action**: Use web interface until specific brain aging project emerges

### Potential Future Projects
**IF we develop these, revisit download decision**:
1. Brain-specific aging patterns
2. Layer-specific senescence accumulation
3. Spatial distribution of aging genes
4. Neurodegenerative disease connections
5. Cross-tissue aging comparison (brain vs other organs)

### Web Interface Advantage
- No storage cost
- Always up-to-date
- Interactive exploration
- Sufficient for gene lookups
- Can generate publication figures

**Decision**: Defer download, use web interface for exploration

### Last Updated
2025-12-03 - Initial documentation with storage management strategy

---

**Created**: 2025-12-03
**Last Modified**: 2025-12-03
**Status**: Documented - Web interface only (no local download)
**Storage Decision**: DO NOT DOWNLOAD (2 GB) unless brain aging project justifies
**Priority**: Low (brain-specific, not core aging biology)
