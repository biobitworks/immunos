---
name: "BICCN"
type: "database"
category: "neuroscience-cell-census"
maintained_by: "Brain Initiative Cell Census Network"
website: "https://biccn.org/"
repository: "https://github.com/brain-initiative"
license: "CC-BY-4.0"
programming_language: "Multiple"
installation: "web-portal"
status: "active"
last_release: "2024"
documentation_quality: "excellent"
data_size: "~100+ GB (multi-modal datasets)"
added_date: "2025-12-05"
last_updated: "2025-12-05"
tags:
  - database
  - neuroscience
  - brain
  - single-cell
  - spatial-transcriptomics
  - cell-atlas
---

# BICCN - Brain Initiative Cell Census Network

## Quick Summary

Comprehensive census of brain cell types across species using multi-modal approaches including single-cell genomics, spatial transcriptomics, epigenomics, morphology, and connectivity. Part of the NIH BRAIN Initiative aiming to classify all cell types in the mammalian brain.

## Overview

### Purpose and Scope
BICCN integrates multiple data modalities to create a comprehensive atlas of brain cell types, their molecular signatures, spatial distributions, morphologies, and connectivity patterns. The network combines data from human, mouse, and non-human primate brains.

### Key Features
- **Multi-modal data**: Single-cell RNA-seq, ATAC-seq, spatial transcriptomics, Patch-seq, MERFISH
- **Cross-species**: Human, mouse, marmoset, macaque
- **Cell type taxonomy**: Hierarchical classification of brain cell types
- **Spatial mapping**: Location and distribution of cell types
- **Morphology**: 3D reconstructions of neuron morphology
- **Connectivity**: Circuit mapping and synaptic connections
- **Public data portal**: Open access to all datasets
- **Interactive tools**: Browser-based data exploration

### Version Information
- **Current Status**: Ongoing (Phase II)
- **Data Source**: NIH BRAIN Initiative funded consortium
- **Update Frequency**: Continuous data releases
- **Status**: Actively maintained and expanding

## Technical Details

### Data Statistics

**⚠️ DATA SIZE TRACKING - LARGE MULTI-MODAL DATASETS**

| Component | Size | Records | Notes |
|-----------|------|---------|-------|
| **Total datasets** | **100+ GB** | Varies by modality | Very large |
| Single-cell RNA-seq | ~50 GB | Millions of cells | Across species |
| Spatial transcriptomics | ~20 GB | Multiple platforms | MERFISH, Visium, etc. |
| Epigenomics | ~15 GB | ATAC-seq, methylation | Chromatin accessibility |
| Morphology | ~10 GB | Thousands of neurons | 3D reconstructions |
| Connectivity | ~5 GB | Projection mapping | Anterograde/retrograde tracing |

**Storage recommendation**: Use web portal and targeted downloads. Full dataset download not recommended for general users.

### Data Modalities

#### 1. Single-Cell Transcriptomics
- **scRNA-seq**: Gene expression profiles of individual cells
- **snRNA-seq**: Single-nucleus RNA-seq for frozen tissue
- **Coverage**: Whole brain regions, multiple species
- **Cell count**: Millions of profiled cells

#### 2. Spatial Transcriptomics
- **MERFISH**: Multiplexed error-robust fluorescence in situ hybridization
- **Visium**: 10x Genomics spatial platform
- **smFISH**: Single-molecule FISH
- **Resolution**: Subcellular to tissue-scale

#### 3. Epigenomics
- **scATAC-seq**: Single-cell chromatin accessibility
- **DNA methylation**: CpG methylation patterns
- **Histone marks**: ChIP-seq for regulatory elements

#### 4. Morphology
- **Patch-seq**: Simultaneous electrophysiology, morphology, and transcriptomics
- **Reconstructions**: 3D neuron morphology
- **Classification**: Morphological cell types

#### 5. Connectivity
- **Tracing**: Anterograde and retrograde viral tracing
- **Projection mapping**: Circuit connectivity
- **Synaptic partners**: Connection specificity

## Access Methods

### Web Portal (Recommended)
**URL**: https://biccn.org/

**NeMO Archive** (Main data portal):
- URL: https://nemoarchive.org/
- Browse datasets by species, region, modality
- Interactive visualization tools
- Download individual datasets

**Allen Brain Map Integration**:
- URL: https://portal.brain-map.org/
- BICCN data integrated with Allen Institute resources
- Brain region atlases
- Gene expression browsers

### Data Download

**NeMO Portal**:
```bash
# Browse and download from https://nemoarchive.org/
# Requires registration (free)
# Select datasets by:
# - Species (human, mouse, NHP)
# - Brain region (cortex, hippocampus, etc.)
# - Data modality (RNA-seq, spatial, etc.)
```

**Programmatic Access**:
```python
# Example using AWS CLI (datasets hosted on S3)
aws s3 ls s3://nemo-public/biccn/ --no-sign-request

# Download specific dataset
aws s3 cp s3://nemo-public/biccn/[dataset-path] . --no-sign-request
```

## Integration with Our Research

### Relevance to Aging Biology

#### Brain Aging Context
**Potential relevance**: Understanding brain cell types can inform age-related changes in brain cell populations.

**Aging-related questions**:
- Do specific brain cell types show age-related transcriptomic changes?
- Are GenAge genes expressed in particular brain cell types?
- Do senescence markers (CellAge) appear in aging brain cells?
- Cell type-specific vulnerability to aging?

#### Cross-Reference Opportunities

**GenAge genes in brain cell types**:
```python
# Hypothetical analysis
import scanpy as sc

# Load BICCN single-cell data (if downloaded)
adata = sc.read_h5ad('biccn_brain_scRNAseq.h5ad')

# Load GenAge genes
genage_genes = ['TP53', 'FOXO3', 'IGF1R', 'SIRT1', ...]

# Check which cell types express aging genes
sc.pl.dotplot(adata, genage_genes, groupby='cell_type')
```

**CellAge senescence in brain**:
- Are cellular senescence markers expressed in brain?
- Which cell types become senescent (neurons, glia)?
- Spatial distribution of senescent cells in aging brain?

### Current Priority: **LOW**
**Rationale**:
- Not core aging database (GenAge/CellAge)
- Very large size (100+ GB) vs. limited aging-specific utility
- Brain-specific (not general aging biology)
- Web portal sufficient for exploratory queries
- More relevant to neuroscience than longevity research

**Future potential**:
- If developing brain aging project → high value
- If studying neurodegeneration → relevant
- If investigating cell type-specific aging → useful
- Cross-reference with [[spatiallibd|spatialLIBD]] for spatial context

## Data Provenance

### BICCN Consortium

**Funding**: NIH BRAIN Initiative

**Participating Institutions** (50+ centers):
- Allen Institute for Brain Science
- Broad Institute
- Salk Institute
- Cold Spring Harbor Laboratory
- Many others

**Key Publications**:
- BICCN (2021). "A multimodal cell census and atlas of the mammalian primary motor cortex"
- *Nature* 598: 86-102
- DOI: 10.1038/s41586-021-03950-0

### Technology Platforms
- **Single-cell sequencing**: 10x Genomics, Smart-seq, sci-RNA-seq
- **Spatial**: MERFISH, Visium, STARmap, seqFISH
- **Imaging**: Two-photon microscopy, electron microscopy
- **Connectivity**: Viral tracing, MAPseq

## Citation and Attribution

### Citing BICCN Data

**Primary Citation**:
```bibtex
@article{biccn2021census,
  author = {BICCN and others},
  title = {A multimodal cell census and atlas of the mammalian primary motor cortex},
  journal = {Nature},
  volume = {598},
  pages = {86--102},
  year = {2021},
  doi = {10.1038/s41586-021-03950-0}
}
```

**Data Citation**:
```
Data from the Brain Initiative Cell Census Network (BICCN),
accessed via NeMO Archive (https://nemoarchive.org/)
```

### License
**CC-BY-4.0** - Free for academic and commercial use with attribution

## Performance and Limitations

### Advantages
- **Comprehensive**: Multiple data modalities integrated
- **High quality**: Rigorous quality control
- **Public access**: All data freely available
- **Well-documented**: Extensive metadata and protocols
- **Cross-species**: Comparative analysis enabled
- **Active development**: Ongoing data releases

### Limitations
- **Very large datasets** (100+ GB) - storage/computational requirements
- **Brain-specific** - not generalizable to other tissues
- **Complexity** - requires significant bioinformatics expertise
- **Young field** - cell type classifications still evolving
- **Species focus** - mostly mouse, some human and NHP
- **Regional coverage** - not all brain regions equally covered

### Storage Impact Assessment

**If downloaded**:
- **Full dataset**: 100+ GB (not recommended)
- **Single modality**: 10-50 GB
- **Single dataset**: 1-10 GB
- **Impact**: Very significant storage requirement

**Current data/ directory**:
- GenAge/CellAge: ~420 KB
- **With BICCN subset**: 10-100 GB (23,800x - 238,000x increase)

**Recommendation**: DO NOT download full dataset. Use web portal for queries. Download only specific datasets if brain aging project emerges.

## Documentation and Support

### Official Documentation
- **Main site**: https://biccn.org/
- **NeMO Archive**: https://nemoarchive.org/
- **Data standards**: https://biccn.org/data-standards/
- **Publications**: https://biccn.org/publications/

### Interactive Tools
- **Allen Brain Map**: https://portal.brain-map.org/
- **Cell type browsers**: Interactive exploration of cell types
- **Spatial viewers**: Visualize spatial transcriptomics data

### Support Channels
- **Email**: nemo@som.umaryland.edu (NeMO Archive)
- **Help desk**: https://nemoarchive.org/resources/help.html
- **Tutorials**: https://nemoarchive.org/resources/tutorials.html

## Related Resources

### Similar Cell Atlas Projects
- **[[spatiallibd|spatialLIBD]]** - Brain spatial transcriptomics (more focused, 2 GB)
- **Human Cell Atlas** - Multi-tissue cell census
- **Allen Brain Atlas** - Gene expression atlases
- **Tabula Muris** - Mouse cell atlas

### Our Databases
- **[[../datasets/genage|GenAge]]** - Check aging gene expression by cell type
- **[[../datasets/cellage|CellAge]]** - Senescence markers in brain cells

### Complementary Resources
- **[[biogrid|BioGRID]]** - Protein interactions (cell type context)
- **[[monarch-kg|Monarch Initiative]]** - Brain phenotypes and diseases

## Usage Examples

### Example 1: Web Portal Exploration (Recommended)

**Goal**: Check if aging genes are expressed in specific brain cell types

**Steps**:
1. Visit https://nemoarchive.org/
2. Browse datasets by "mouse" or "human"
3. Select "cortex" or region of interest
4. Choose "single-cell RNA-seq" dataset
5. Use Allen Brain Map viewer to search genes
6. View expression across cell types

**No download required** ✅
**Time**: 10-15 minutes

### Example 2: Targeted Dataset Download

**Goal**: Analyze GenAge gene expression in specific cell types

```bash
# Register at https://nemoarchive.org/ (free)
# Browse datasets and get S3 path

# Download specific dataset
aws s3 cp s3://nemo-public/biccn/grant/[specific-dataset].h5ad . --no-sign-request

# Analyze in Python
import scanpy as sc
adata = sc.read_h5ad('dataset.h5ad')

# Filter for aging genes
aging_genes = ['TP53', 'FOXO3', 'IGF1R', 'SIRT1']
adata_aging = adata[:, adata.var_names.isin(aging_genes)]

# Plot by cell type
sc.pl.dotplot(adata_aging, aging_genes, groupby='cell_type')
```

## Notes

### Strategic Importance
**Low-Moderate** for current aging biology research:
- Brain-specific, not general aging
- Valuable IF pursuing brain aging research
- Excellent resource but massive storage cost
- Web portal sufficient for exploration

### When to Revisit
**Consider downloading if**:
1. Developing brain aging project
2. Studying age-related neurodegeneration
3. Investigating cell type-specific aging in brain
4. Cross-referencing with spatialLIBD data
5. Senescence in brain cell populations

### Comparison with spatialLIBD
| Feature | BICCN | spatialLIBD |
|---------|-------|-------------|
| **Scope** | Whole brain census | DLPFC focused |
| **Size** | 100+ GB | 2 GB |
| **Modalities** | Multi-modal | Spatial transcriptomics |
| **Species** | Multiple | Human |
| **Subjects** | Hundreds | 3 |
| **Best for** | Cell type discovery | Spatial patterns |

**Recommendation for aging research**: Start with spatialLIBD (smaller, focused). Use BICCN if need cell type classification or multi-modal data.

### Last Updated
2025-12-05 - Initial documentation

---

**Created**: 2025-12-05
**Last Modified**: 2025-12-05
**Status**: Documented - Web portal recommended (avoid full download)
**Storage Decision**: DO NOT DOWNLOAD full dataset (100+ GB)
**Priority**: Low (brain-specific, massive storage requirement)
**Use Case**: Brain aging research only
