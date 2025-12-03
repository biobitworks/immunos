---
name: "longevity-db (HuggingFace Organization)"
type: "data-repository"
location: "Online (HuggingFace Hub)"
founded: "Unknown"
website: "https://huggingface.co/longevity-db"
focus_areas:
  - "Aging and longevity datasets"
  - "Open science data sharing"
  - "Machine learning for aging research"
  - "Longevity biomarkers"
  - "Aging biology databases"
size: "unknown"
funding: "unknown"
status: "active"
collaboration_status: "potential"
added_date: "2025-12-03"
last_updated: "2025-12-03"
tags:
  - organization
  - data-repository
  - aging-biology
  - huggingface
  - longevity
  - open-science
---

# longevity-db (HuggingFace Organization)

## Quick Summary

HuggingFace organization dedicated to hosting and sharing aging and longevity-related datasets for machine learning and data science applications. Provides centralized access to aging biology data in ML-friendly formats.

## Overview

### Mission
Appears to be focused on making aging and longevity datasets accessible to the machine learning and data science communities through HuggingFace Hub infrastructure.

### Purpose
- Host aging biology datasets in standardized formats
- Enable ML/AI applications for aging research
- Promote open science and data sharing
- Facilitate reproducible research

### Platform
**HuggingFace Hub**: Industry-standard platform for:
- Dataset hosting and versioning
- Easy Python/R access via `datasets` library
- Built-in dataset viewers
- Documentation and metadata
- Community engagement

## Available Datasets

### Current Dataset Status
*To be investigated*

**Potential datasets to explore**:
- Aging-related gene databases
- Longevity biomarker data
- Clinical aging studies
- Model organism lifespan data
- Transcriptomic aging signatures

**Search needed**: Browse https://huggingface.co/longevity-db to identify available datasets

### Dataset Formats
HuggingFace datasets typically provided in:
- Parquet format (efficient columnar storage)
- CSV (for smaller datasets)
- JSON (for structured data)
- Arrow format (for fast access)

**Access pattern**:
```python
from datasets import load_dataset

# Load aging dataset from longevity-db
dataset = load_dataset("longevity-db/dataset-name")
```

## Connection to Our Research

### Potential for Our Datasets

#### Upload GenAge and CellAge
**High priority**: We should upload our processed [[../datasets/genage|GenAge]] and [[../datasets/cellage|CellAge]] datasets to make them accessible to ML community.

**Benefits**:
- Wider visibility and usage
- Standardized access methods
- Version control
- Community engagement
- Citations and attribution

**Datasets to upload**:
1. **GenAge Human** (307 genes)
2. **GenAge Model Organisms** (2,205 genes)
3. **GenAge Expression Signatures**
4. **CellAge Senescence Genes** (950 genes)
5. **CellAge Expression Signatures** (1,259 genes)

#### Organization Membership
**Consider**: Join as contributor to longevity-db organization or create parallel organization

**Options**:
1. Upload to `longevity-db/genage-human`, `longevity-db/cellage` (if accepted as contributor)
2. Create our own org: `aging-research/genage`, `aging-research/cellage`
3. Upload to personal account then request transfer

### Data Integration Opportunities

1. **Discover related datasets**: Find other aging datasets in longevity-db organization
2. **Benchmark comparisons**: Compare GenAge/CellAge with other aging gene sets
3. **Meta-analyses**: Integrate multiple aging datasets for broader analyses
4. **ML applications**: Enable easy access for AI/ML researchers working on aging

### Collaboration Potential

1. **Data sharing**: Contribute our processed datasets
2. **Standardization**: Adopt common data schemas
3. **Documentation**: Share data provenance and processing pipelines
4. **Community building**: Engage with other aging biology data scientists

## HuggingFace Integration Plan

### Phase 1: Account and Authentication
```bash
# Install HuggingFace CLI
pip install huggingface_hub

# Authenticate
huggingface-cli login
# Enter token from https://huggingface.co/settings/tokens
```

### Phase 2: Dataset Preparation
```python
from datasets import Dataset
import pandas as pd

# Load our data
genage = pd.read_csv('data/genage/human/genage_human.csv')

# Convert to HuggingFace Dataset
dataset = Dataset.from_pandas(genage)

# Add metadata
dataset_dict = {
    'name': 'GenAge Human',
    'description': 'Human Ageing Genomic Resources - 307 aging-related genes',
    'version': 'Build 21',
    'license': 'CC-BY-3.0',
    'citation': 'de Magalhães et al. (2024) Nucleic Acids Research'
}
```

### Phase 3: Upload to HuggingFace
```python
from huggingface_hub import HfApi

api = HfApi()

# Create repository
api.create_repo(
    repo_id="longevity-db/genage-human",
    repo_type="dataset",
    private=False
)

# Upload dataset
dataset.push_to_hub("longevity-db/genage-human")
```

### Phase 4: Documentation
Create README.md with:
- Dataset description
- Data schema
- Usage examples
- Citation information
- License details
- Links to source (genomics.senescence.info)

## Technical Integration

### Access Methods

#### Python
```python
from datasets import load_dataset

# Load GenAge from longevity-db (once uploaded)
genage = load_dataset("longevity-db/genage-human")

# Access as pandas
df = genage['train'].to_pandas()

# Filter by criteria
aging_genes = df[df['longevity_influence'] == 'Pro-Longevity']
```

#### R (via reticulate)
```r
library(reticulate)
datasets <- import("datasets")

# Load dataset
genage <- datasets$load_dataset("longevity-db/genage-human")

# Convert to R dataframe
genage_df <- genage$`train`$to_pandas()
```

### Integration with Our Workflow

**Current workflow**:
```
genomics.senescence.info → local files → our analysis
```

**Enhanced workflow with HuggingFace**:
```
genomics.senescence.info → local files → HuggingFace → broader community
                                      ↓
                                  our analysis
```

## Community and Collaboration

### HuggingFace Community Features
- **Dataset discussions**: Community Q&A
- **Pull requests**: Community contributions
- **Stars**: Popularity metric
- **Downloads**: Usage tracking
- **Model cards**: Linked ML models using datasets

### Engagement Opportunities
- Respond to questions about datasets
- Accept community improvements
- Link to papers and analyses
- Collaborate with ML researchers

## Next Steps

### Immediate Actions
- [ ] Browse existing longevity-db datasets
- [ ] Contact longevity-db organization (if exists)
- [ ] Create HuggingFace account (if not existing)
- [ ] Generate HuggingFace API token
- [ ] Plan GenAge/CellAge upload strategy

### Data Upload Preparation
- [ ] Clean and validate GenAge datasets
- [ ] Clean and validate CellAge datasets
- [ ] Create comprehensive README files
- [ ] Prepare dataset cards with metadata
- [ ] Test upload process with small dataset

### Documentation
- [ ] Write usage examples
- [ ] Document data provenance
- [ ] Create citation guidelines
- [ ] Link to source databases

## External Links

- **HuggingFace Profile**: https://huggingface.co/longevity-db
- **HuggingFace Datasets**: https://huggingface.co/datasets
- **Datasets Library Docs**: https://huggingface.co/docs/datasets/
- **HuggingFace Hub Python**: https://huggingface.co/docs/huggingface_hub/

## Notes

### Strategic Importance
**Very High**: HuggingFace is the standard platform for ML datasets. Having our aging biology data there will:
- Increase visibility and usage
- Enable ML applications
- Establish our data as community resource
- Facilitate citations and attribution
- Support reproducible research

### Benefits of HuggingFace Hub

**For data providers (us)**:
- Free hosting
- Version control
- Download statistics
- Community engagement
- Professional presentation

**For data users**:
- Easy access via Python/R
- Standardized formats
- Built-in dataset viewers
- No need to download manually
- Always latest version

### Collaboration Model
**Open collaboration**: Anyone can use datasets, contribute improvements, create derivatives (with attribution)

### Last Updated
2025-12-03 - Initial planning for HuggingFace integration

---

**Created**: 2025-12-03
**Last Modified**: 2025-12-03
**Status**: Planning phase - High priority for dataset uploads
