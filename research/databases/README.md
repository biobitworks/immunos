# Database and Data Repository Profiles

## Summary
Contains 0 subdirectories and 5 files.


This directory documents external databases, data repositories, and platforms that provide aging biology datasets and resources.

## Purpose

Track and document:
- Public aging biology databases (GenAge, CellAge, etc.)
- Data hosting platforms (HuggingFace, GitHub, Zenodo)
- Protein interaction databases (STRING-DB, BioGRID)
- Gene expression repositories (GEO, ArrayExpress)
- Clinical databases (dbGaP, UK Biobank)
- Integration methods and APIs

## Database Categories

### Aging-Specific Databases
- **HAGR** (Human Ageing Genomic Resources) - GenAge, CellAge, etc.
- **DrugAge** - Compounds affecting longevity
- **LongevityMap** - Human longevity genetic variants

### General Biological Databases
- **NCBI** - Gene, PubMed, GEO
- **UniProt** - Protein sequences and annotations
- **STRING-DB** - Protein-protein interactions
- **Reactome** - Pathway databases

### Data Hosting Platforms
- **HuggingFace Hub** - ML datasets and models
- **GitHub** - Version control and code/data sharing
- **Zenodo** - Research data archiving (DOI assignment)
- **figshare** - Research outputs and datasets

### Clinical and Population Data
- **UK Biobank** - Large-scale population health data
- **dbGaP** - Genotype-phenotype associations
- **ClinVar** - Clinical variant interpretations

## File Naming Convention

Use lowercase with hyphens: `database-name.md`

Examples:
- `huggingface-hub.md`
- `string-db.md`
- `hagr-databases.md`
- `uk-biobank.md`

## Profile Template

Use the standard template at `/templates/database-tool-profile.md` which includes:
- YAML frontmatter with metadata
- Database overview and scope
- Data types and schemas
- Access methods (API, download, web interface)
- Authentication requirements
- Usage examples and code snippets
- Integration with our projects
- Citation information

## Integration Strategies

### HuggingFace Hub
```python
# Upload dataset to HuggingFace
from huggingface_hub import HfApi
api = HfApi()
api.upload_folder(
    folder_path="data/genage",
    repo_id="username/genage-database",
    repo_type="dataset"
)
```

### STRING-DB API
```python
# Query protein interactions
import requests
url = "https://string-db.org/api/json/network"
params = {"identifiers": "GHR", "species": 9606}
response = requests.get(url, params=params)
```

### NCBI E-utilities
```bash
# Fetch gene information
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=2690&retmode=xml"
```

## Example Queries

### Dataview: List all databases with APIs
```dataview
TABLE type, access_method, last_accessed
FROM "research/databases"
WHERE contains(access_method, "API")
```

### Dataview: Find aging-specific databases
```dataview
LIST focus_area, url
FROM "research/databases"
WHERE contains(focus_area, "aging")
```

## Current Databases

- **HuggingFace Hub** - ML dataset hosting platform
- **HAGR** - Human Ageing Genomic Resources (GenAge, CellAge)
- **STRING-DB** - Protein interaction network
- *(more to be added)*

---

**Directory**: `/research/databases/`
**Template**: `/templates/database-tool-profile.md`
**Script**: *N/A (manual creation)*
**Last Updated**: 2025-12-03

## Directory Map
```
databases/
├── biccn.md
├── biogrid.md
├── huggingface-hub.md
├── regulomedb.md
└── spatiallibd.md
```
