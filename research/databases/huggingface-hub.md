---
name: "HuggingFace Hub"
type: "platform"
category: "data-repository"
maintained_by: "HuggingFace Inc."
website: "https://huggingface.co"
repository: "https://github.com/huggingface/datasets"
license: "Apache-2.0"
programming_language: "Python"
installation: "pip"
status: "active"
last_release: "Continuous"
documentation_quality: "excellent"
added_date: "2025-12-03"
last_updated: "2025-12-03"
tags:
  - database
  - data-repository
  - machine-learning
  - open-science
  - aging-biology
---

# HuggingFace Hub

## Quick Summary

Leading open platform for hosting, discovering, and collaborating on machine learning models and datasets. Industry-standard infrastructure for sharing ML-ready data with built-in versioning, documentation, and community features.

## Overview

### Purpose and Scope
HuggingFace Hub serves as the central repository for the machine learning community to share models, datasets, and demos (Spaces). For our purposes, we focus on the **Datasets** component for hosting aging biology data.

### Key Features
- **Free hosting**: Unlimited public repositories
- **Version control**: Git/Git LFS integration
- **Easy access**: Python `datasets` library for seamless loading
- **Built-in viewers**: Automatic dataset visualization in browser
- **Documentation**: Dataset cards with rich markdown support
- **Community**: Discussions, pull requests, collaborative improvement
- **API access**: Programmatic upload/download
- **Format agnostic**: Supports CSV, JSON, Parquet, Arrow, and more

### Version Information
- **Platform**: Continuously updated
- **datasets library**: Latest stable version via pip
- **huggingface_hub library**: For upload/management
- **Stable**: Yes

## Technical Details

### Data Storage Architecture

#### Repository Structure
Each dataset is a Git repository:
```
dataset-repo/
├── README.md              # Dataset card with metadata
├── data/
│   ├── train.csv         # Or other formats
│   ├── test.csv
│   └── validation.csv
├── .gitattributes        # LFS configuration
└── dataset_infos.json    # Auto-generated metadata
```

#### Storage Formats
**Recommended**: Parquet (columnar, compressed, efficient)
**Supported**: CSV, TSV, JSON, JSONL, Arrow, SQL databases

#### Size Limits
- **Regular files**: Up to 50MB per file
- **LFS files**: Unlimited (automatically handled for files >10MB)
- **Total size**: No hard limit on public datasets

### Data Access Methods

#### Python: datasets Library

**Installation**:
```bash
pip install datasets
```

**Loading dataset**:
```python
from datasets import load_dataset

# Load our GenAge dataset (once uploaded)
genage = load_dataset("longevity-db/genage-human")

# Access data
df = genage['train'].to_pandas()  # Convert to pandas
genes = list(genage['train']['symbol'])  # Access specific column
```

**Streaming large datasets**:
```python
# For large datasets, stream without downloading all
genage = load_dataset("longevity-db/genage-human", streaming=True)
for record in genage['train']:
    process(record)
```

#### Python: Direct Download

```python
from huggingface_hub import hf_hub_download

# Download specific file
file_path = hf_hub_download(
    repo_id="longevity-db/genage-human",
    filename="data/genage_human.csv",
    repo_type="dataset"
)
```

#### Command Line

```bash
# Clone entire dataset repository
git clone https://huggingface.co/datasets/longevity-db/genage-human

# Use Git LFS for large files
git lfs install
git clone https://huggingface.co/datasets/longevity-db/genage-human
```

#### R: Via reticulate

```r
library(reticulate)

# Use Python datasets library
datasets <- import("datasets")
genage <- datasets$load_dataset("longevity-db/genage-human")

# Convert to R dataframe
genage_df <- genage$`train`$to_pandas()
```

## Uploading Our Datasets

### Phase 1: Setup and Authentication

#### Create Account
1. Visit https://huggingface.co/join
2. Verify email
3. Complete profile

#### Generate API Token
1. Go to https://huggingface.co/settings/tokens
2. Create new token with "write" access
3. Copy token (treat as password)

#### Install Libraries
```bash
pip install datasets huggingface_hub
```

#### Authenticate
```bash
# Method 1: Interactive login
huggingface-cli login
# Paste token when prompted

# Method 2: Environment variable
export HUGGING_FACE_HUB_TOKEN="your_token_here"
```

**Verify authentication**:
```python
from huggingface_hub import whoami
print(whoami())
```

### Phase 2: Prepare GenAge Datasets

#### Script: `scripts/prepare_genage_for_hf.py`

```python
from datasets import Dataset, DatasetDict
import pandas as pd

# Load GenAge human data
genage_human = pd.read_csv('/Users/byron/projects/data/genage/human/genage_human.csv')

# Convert to HuggingFace Dataset
dataset = Dataset.from_pandas(genage_human)

# Create DatasetDict (standard structure)
dataset_dict = DatasetDict({
    'train': dataset  # Use 'train' split by convention
})

# Add metadata
dataset_dict['train'].info.description = """
GenAge Human: Human Ageing Genomic Resources

Manually curated database of 307 genes associated with human aging and longevity.
Build 21, downloaded 2025-12-02 from genomics.senescence.info.

Citation: de Magalhães et al. (2024) Nucleic Acids Research 52(D1):D900-D908.
License: CC-BY-3.0
"""

# Save locally first for review
dataset_dict.save_to_disk('genage_human_hf')
```

#### Create README.md (Dataset Card)

**Location**: `data/genage/human/README_HF.md`

```markdown
---
license: cc-by-3.0
task_categories:
- feature-extraction
tags:
- biology
- aging
- longevity
- genomics
- genes
size_categories:
- n<1K
---

# GenAge Human: Human Ageing Genomic Resources

## Dataset Description

GenAge Human is a manually curated database of genes associated with human aging and longevity.

### Dataset Summary

- **Curated by**: Integrative Genomics of Ageing Group, University of Liverpool
- **Language(s)**: English (gene annotations)
- **License**: Creative Commons Attribution 3.0 Unported License
- **Genes**: 307 human aging-related genes
- **Build**: 21 (2024)

### Supported Tasks

- Gene function analysis
- Aging biomarker discovery
- Longevity pathway analysis
- Cross-species aging comparisons

## Dataset Structure

### Data Fields

- `genage_id`: GenAge database identifier
- `symbol`: Gene symbol (HGNC)
- `name`: Full gene name
- `entrez_gene_id`: NCBI Gene database ID
- `uniprot`: UniProt accession
- `why`: Reason for inclusion in GenAge
- `band`: Cytogenetic location
- `functional_category`: Biological process
- `longevity_influence`: Pro/anti-longevity effect

### Data Splits

Single split containing all 307 genes (conventionally named "train").

## Dataset Creation

### Source Data

Downloaded from: https://genomics.senescence.info/genes/human.html
Date: 2025-12-02
Build: 21

### Annotations

All annotations are from the original GenAge database, manually curated from literature.

## Citation

**BibTeX**:
```bibtex
@article{demagalhaes2024hagr,
  author = {de Magalh{\~a}es, Jo{\~a}o Pedro and Tacutu, Robi and Bunu{\c{s}}, Alexandru N. and Apopei, Gennady A.},
  title = {Human Ageing Genomic Resources: updates on key databases in ageing research},
  journal = {Nucleic Acids Research},
  volume = {52},
  number = {D1},
  pages = {D900--D908},
  year = {2024},
  doi = {10.1093/nar/gkad927}
}
```

## Additional Information

### Licensing
Creative Commons Attribution 3.0 Unported License - Free for commercial and research use with attribution.

### Contact
- **Maintainer**: [Your name/organization]
- **Source**: Human Ageing Genomic Resources (https://genomics.senescence.info)
```

### Phase 3: Upload to HuggingFace

#### Option 1: Upload via Python

```python
from huggingface_hub import HfApi, create_repo

api = HfApi()

# Create repository (one time)
create_repo(
    repo_id="longevity-db/genage-human",
    repo_type="dataset",
    private=False
)

# Upload dataset
dataset_dict.push_to_hub(
    "longevity-db/genage-human",
    commit_message="Initial upload of GenAge Human Build 21"
)

# Upload README separately
api.upload_file(
    path_or_fileobj="data/genage/human/README_HF.md",
    path_in_repo="README.md",
    repo_id="longevity-db/genage-human",
    repo_type="dataset"
)
```

#### Option 2: Upload via CLI

```bash
# Create repo
huggingface-cli repo create longevity-db/genage-human --type dataset

# Upload files
cd data/genage/human
huggingface-cli upload longevity-db/genage-human . .
```

### Phase 4: Datasets to Upload

#### Priority 1: Core Aging Databases
- [x] Plan GenAge Human (307 genes)
- [ ] GenAge Model Organisms (2,205 genes)
- [ ] CellAge (950 genes)
- [ ] CellAge Signatures (1,259 genes)

#### Priority 2: Expression Data
- [ ] GenAge Expression Signatures
- [ ] CellAge Expression Profiles

#### Priority 3: Integrated Datasets
- [ ] GenAge-CellAge Overlap Analysis
- [ ] Cross-species Aging Gene Comparison

## Integration with Our Workflow

### Current Data Flow
```
genomics.senescence.info
    ↓ manual download
/Users/byron/projects/data/
    ↓ analysis scripts
Results
```

### Enhanced Data Flow with HuggingFace
```
genomics.senescence.info
    ↓
/Users/byron/projects/data/ (local cache)
    ↓
HuggingFace Hub (public repository)
    ↓
ML Community + Our Analysis + Collaborators
```

### Dual Storage Strategy
**Local**: `/Users/byron/projects/data/` - Master copies, full data
**HuggingFace**: Public ML-ready versions for community

**Benefits**:
- Local: Full control, all formats, immediate access
- HuggingFace: Public sharing, version control, community engagement

### Usage in Our Scripts

```python
# In our analysis scripts, load from HuggingFace
from datasets import load_dataset

def load_aging_data():
    """Load aging datasets from HuggingFace Hub."""
    genage_human = load_dataset("longevity-db/genage-human")['train'].to_pandas()
    cellage = load_dataset("longevity-db/cellage")['train'].to_pandas()

    return genage_human, cellage

# Use in analysis
genage, cellage = load_aging_data()
overlap = set(genage['symbol']) & set(cellage['Gene symbol'])
print(f"Genes in both databases: {len(overlap)}")
```

## GitHub Actions Integration

### Auto-Sync Workflow

**File**: `.github/workflows/sync-to-huggingface.yml`

```yaml
name: Sync Datasets to HuggingFace

on:
  push:
    branches: [main]
    paths:
      - 'data/**'
  workflow_dispatch:  # Manual trigger

jobs:
  sync:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
        with:
          lfs: true

      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.10'

      - name: Install dependencies
        run: |
          pip install datasets huggingface_hub pandas

      - name: Sync to HuggingFace
        env:
          HF_TOKEN: ${{ secrets.HF_TOKEN }}
        run: |
          python scripts/sync_to_huggingface.py
```

**Script**: `scripts/sync_to_huggingface.py`

```python
import os
from datasets import Dataset
import pandas as pd

token = os.environ['HF_TOKEN']

# Load and upload GenAge
genage = pd.read_csv('data/genage/human/genage_human.csv')
Dataset.from_pandas(genage).push_to_hub(
    "longevity-db/genage-human",
    token=token,
    commit_message="Auto-sync from GitHub"
)

# Repeat for other datasets...
```

## Performance and Limitations

### Performance Characteristics
- **Download speed**: Fast (CDN-backed)
- **Load time**: Instant for cached datasets
- **Streaming**: Efficient for large datasets
- **API rate limits**: Generous, rarely hit

### Known Limitations
- **Authentication**: Required for uploading, not for downloading public datasets
- **Git LFS**: Large files require LFS, which adds complexity
- **Search**: Dataset discovery can be improved
- **Versioning**: Git-based, can be complex for non-technical users

### Advantages Over Alternatives

| Feature | HuggingFace | Zenodo | Figshare | GitHub | Kaggle |
|---------|-------------|--------|----------|--------|--------|
| **ML Integration** | ⭐⭐⭐⭐⭐ | ⭐ | ⭐ | ⭐⭐ | ⭐⭐⭐⭐ |
| **Version Control** | ⭐⭐⭐⭐⭐ | ⭐⭐ | ⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐ |
| **DOI Assignment** | ⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | ⭐ | ⭐ |
| **Ease of Use** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐ | ⭐⭐⭐⭐ | ⭐⭐⭐ |
| **Cost** | Free | Free | Free | Free | Free |

**Recommendation**: Use **HuggingFace** for ML/data science community + **Zenodo** for DOI/long-term archival

## Citation and Attribution

### Citing HuggingFace Hub

```
Lhoest, Q., et al. (2021). Datasets: A Community Library for Natural Language Processing.
Proceedings of the 2021 Conference on Empirical Methods in Natural Language Processing:
System Demonstrations, 175-184.
```

### Citing Our Datasets on HuggingFace

```
[Your Name/Organization]. (2025). GenAge Human: Human Ageing Genomic Resources.
HuggingFace Datasets. https://huggingface.co/datasets/longevity-db/genage-human

Original source: de Magalhães et al. (2024) Nucleic Acids Research 52(D1):D900-D908.
```

## Documentation and Support

### Official Documentation
- **Main docs**: https://huggingface.co/docs/datasets/
- **Upload guide**: https://huggingface.co/docs/datasets/upload_dataset
- **Dataset cards**: https://huggingface.co/docs/hub/datasets-cards
- **API reference**: https://huggingface.co/docs/huggingface_hub/

### Community Support
- **Forum**: https://discuss.huggingface.co/
- **Discord**: https://hf.co/join/discord
- **GitHub Issues**: https://github.com/huggingface/datasets/issues

### Learning Resources
- **Quickstart**: https://huggingface.co/docs/datasets/quickstart
- **Tutorial**: https://huggingface.co/docs/datasets/tutorial
- **Examples**: https://github.com/huggingface/datasets/tree/main/examples

## Next Steps

### Immediate Actions
- [x] Research HuggingFace Hub capabilities
- [ ] Create HuggingFace account
- [ ] Generate API token
- [ ] Test authentication
- [ ] Prepare GenAge Human dataset
- [ ] Create comprehensive dataset card
- [ ] Upload GenAge Human

### Short Term (This Week)
- [ ] Upload CellAge dataset
- [ ] Upload GenAge Model Organisms
- [ ] Test loading from HuggingFace
- [ ] Update analysis scripts to use HuggingFace
- [ ] Set up GitHub Actions sync

### Medium Term (This Month)
- [ ] Upload expression signatures
- [ ] Create integrated datasets
- [ ] Write blog post about datasets
- [ ] Engage with HuggingFace community
- [ ] Track dataset usage and downloads

## Related Resources

### Similar Platforms
- **Zenodo** - For DOI assignment and long-term archival
- **Figshare** - Academic data sharing
- **Kaggle Datasets** - ML competition datasets
- **GitHub** - Code and small datasets

### Our Organizations to Target
- [[../organizations/longevity-db-huggingface|longevity-db]] - Aging biology datasets organization
- Create our own org if longevity-db collaboration doesn't work out

## Notes

### Why HuggingFace for Aging Biology?

**Traditional approach**: Download from source → local storage → analysis
**Problem**: Data scattered, not ML-ready, hard to discover

**HuggingFace approach**: Centralized → standardized → easy access
**Benefits**: Discoverable, versioned, community-driven, ML-optimized

### Strategic Value
**Very High**: HuggingFace is becoming the standard for dataset sharing in ML/data science. Having aging biology data there will:
- Increase visibility 10-100x
- Enable ML applications
- Facilitate reproducibility
- Build community
- Establish us as data providers

### Best Practices
1. **Comprehensive README**: Rich dataset cards with examples
2. **Clear licensing**: CC-BY-3.0 for openness
3. **Proper attribution**: Cite original sources
4. **Version control**: Use semantic versioning
5. **Community engagement**: Respond to issues/discussions

### Last Updated
2025-12-03 - Initial documentation and planning

---

**Created**: 2025-12-03
**Last Modified**: 2025-12-03
**Status**: Planning phase - Ready for implementation
**Priority**: Very High
