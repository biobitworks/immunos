# Papers Dashboard

Automated archive of open-access academic papers with figures, organized by DOI and integrated with aging biology research vault.

## Quick Stats

```dataview
TABLE
  journal as Journal,
  year as Year,
  length(authors) as "Authors"
FROM "papers"
WHERE file.name != "INDEX" AND file.name != "README"
SORT year DESC
```

**Total Papers**: 7 (6 new + 1 original)
**Total Figures**: 38
**Latest Download**: 2025-12-03
**Restricted Access**: 2 papers tracked (not downloaded)

## Recent Papers

### [Spatial transcriptomics of the aging mouse brain reveals origins of inflammation in the white matter](10.1038_s41467-025-58466-2/note.md)

**Wang L, Cui CY, Lee CT, et al. (2025)**
*Nature Communications* | DOI: [10.1038/s41467-025-58466-2](https://doi.org/10.1038/s41467-025-58466-2)

**Relevance**: ⭐⭐⭐⭐⭐ Highly relevant - Spatial transcriptomics of aging brain white matter inflammation

**Integration opportunities**:
- Compare with [[../research/databases/spatiallibd|spatialLIBD]] human brain data
- Network analysis with [[../research/databases/biogrid|BioGRID]]
- Cross-reference with [[../data/genage/README|GenAge]] aging genes
- Senescence markers from [[../data/cellage/README|CellAge]]

**Downloaded**: PDF (10 MB) + 7 figures (11.4 MB total)

### [Midbrain extracellular matrix and microglia are associated with cognition in aging mice](10.1038_s41467-025-66434-z/note.md)

**Gray DT, Gutierrez A, Jami-Alahmadi Y, et al. (2025)**
*Nature Communications* | DOI: [10.1038/s41467-025-66434-z](https://doi.org/10.1038/s41467-025-66434-z)

**Relevance**: ⭐⭐⭐⭐⭐ Midbrain ECM and microglia in cognitive aging

**Downloaded**: PDF (385 KB) + 6 figures (6.9 MB total)

### [Vulnerability to memory decline in aging revealed by a mega-analysis of structural brain change](10.1038_s41467-025-66354-y/note.md)

**75 authors (2025)**
*Nature Communications* | DOI: [10.1038/s41467-025-66354-y](https://doi.org/10.1038/s41467-025-66354-y)

**Relevance**: ⭐⭐⭐⭐⭐ Large-scale structural brain aging study

**Downloaded**: PDF (479 KB) + 4 figures (mega-analysis)

### [Organ-specific proteomic aging clocks predict disease and longevity across diverse populations](10.1038_s43587-025-01016-8/note.md)

**25 authors (2025)**
*Nature Aging* | DOI: [10.1038/s43587-025-01016-8](https://doi.org/10.1038/s43587-025-01016-8)

**Relevance**: ⭐⭐⭐⭐⭐ Proteomic aging clocks - highly relevant to GenAge

**Downloaded**: PDF (15 MB) + 5 figures

### [Odoribacter splanchnicus rescues aging-related intestinal P-glycoprotein damage](10.1038_s41467-025-65692-1/note.md)

**18 authors (2025)**
*Nature Communications* | DOI: [10.1038/s41467-025-65692-1](https://doi.org/10.1038/s41467-025-65692-1)

**Relevance**: ⭐⭐⭐⭐ Gut microbiome and aging

**Downloaded**: PDF (6.1 MB) + 6 figures

### [Anti-uPAR CAR T cells reverse and prevent aging-associated defects in intestinal regeneration](10.1038_s43587-025-01022-w/note.md)

**24 authors (2025)**
*Nature Aging* | DOI: [10.1038/s43587-025-01022-w](https://doi.org/10.1038/s43587-025-01022-w)

**Relevance**: ⭐⭐⭐⭐⭐ CAR T cell therapy for aging - breakthrough approach

**Downloaded**: PDF (39 MB) + 5 figures

### [Iron homeostasis and cell clonality drive cancer-associated intestinal DNA methylation drift in aging](10.1038_s43587-025-01021-x/note.md)

**12 authors (2025)**
*Nature Aging* | DOI: [10.1038/s43587-025-01021-x](https://doi.org/10.1038/s43587-025-01021-x)
**License**: CC-BY-NC-ND-4.0

**Relevance**: ⭐⭐⭐⭐ Epigenetic aging and cancer

**Downloaded**: PDF (9.2 MB) + 5 figures

---

## By Topic

### Brain Aging & Cognition
- [[10.1038_s41467-025-58466-2/note|Wang et al. 2025]] - Spatial transcriptomics of white matter inflammation
- [[10.1038_s41467-025-66434-z/note|Gray et al. 2025]] - Midbrain ECM and microglia
- [[10.1038_s41467-025-66354-y/note|Mega-analysis 2025]] - Memory decline and structural brain changes

### Proteomic & Epigenetic Aging
- [[10.1038_s43587-025-01016-8/note|Proteomic aging clocks 2025]] - Organ-specific clocks
- [[10.1038_s43587-025-01021-x/note|Iron homeostasis 2025]] - DNA methylation drift

### Gut Biology & Aging
- [[10.1038_s41467-025-65692-1/note|Odoribacter 2025]] - Microbiome rescue of intestinal function
- [[10.1038_s43587-025-01022-w/note|CAR T cells 2025]] - Intestinal regeneration therapy

## By Year

### 2025
- **Nature Communications**: [[10.1038_s41467-025-58466-2/note|Spatial transcriptomics of aging mouse brain]]

## By Journal

### Nature Communications (4 papers)
- [[10.1038_s41467-025-58466-2/note|Wang et al. 2025]] - Spatial transcriptomics, white matter inflammation
- [[10.1038_s41467-025-66434-z/note|Gray et al. 2025]] - Midbrain ECM and microglia
- [[10.1038_s41467-025-66354-y/note|Mega-analysis 2025]] - Memory decline (75 authors)
- [[10.1038_s41467-025-65692-1/note|Odoribacter 2025]] - Gut microbiome aging

### Nature Aging (3 papers)
- [[10.1038_s43587-025-01016-8/note|Proteomic clocks 2025]] - Organ-specific aging biomarkers
- [[10.1038_s43587-025-01022-w/note|CAR T therapy 2025]] - Intestinal regeneration
- [[10.1038_s43587-025-01021-x/note|Iron homeostasis 2025]] - Epigenetic drift (CC-BY-NC-ND)

## Database

All papers are tracked in `papers.db` SQLite database:

**Query examples**:
```bash
# List all papers
sqlite3 papers.db "SELECT title, journal, year FROM papers ORDER BY year DESC"

# Papers from 2025
sqlite3 papers.db "SELECT title FROM papers WHERE year = 2025"

# Papers with figures
sqlite3 papers.db "SELECT p.title, COUNT(f.id) as fig_count
  FROM papers p
  LEFT JOIN figures f ON p.doi = f.paper_doi
  GROUP BY p.doi"
```

## Adding New Papers

### From URL
```bash
cd /Users/byron/projects/papers
python scripts/download_paper.py "https://doi.org/10.1234/example"
```

### From DOI
```bash
python scripts/download_paper.py --doi "10.1234/example"
```

### With Tags
```bash
python scripts/download_paper.py "URL" --tags aging-biology,spatial-transcriptomics
```

## Integration with Research Vault

### Cross-References

Papers automatically link to:
- **Researchers**: [[../research/researchers/|Profiles]] of paper authors
- **Databases**: [[../research/databases/|External databases]] mentioned in papers
- **Organizations**: [[../research/organizations/|Institutions]] of authors
- **Literature**: [[../research/literature/|Related papers]] in citations

### Network Graph

Open Papers workspace in Obsidian to see:
- Citation networks between papers
- Connections to researchers and databases
- Topic clusters and themes
- Cross-project integrations

**Filter**: Set graph filter to `path:papers/` to isolate papers network

### Dataview Queries

**All papers by relevance**:
```dataview
TABLE
  authors[0] + " et al." as "First Author",
  journal,
  year,
  choice(contains(tags, "aging-biology"), "⭐⭐⭐⭐⭐", "⭐⭐⭐") as Relevance
FROM "papers"
WHERE file.name = "note"
SORT year DESC
```

**Papers with figures**:
```dataview
TABLE
  length(figures) as "Figures",
  journal,
  year
FROM "papers"
WHERE file.name = "note" AND length(figures) > 0
```

## Research Questions

### Cross-Species Comparisons
**Q**: How do aging gene expression patterns compare between mouse (Wang 2025) and human (spatialLIBD)?

**Data sources**:
- [[10.1038_s41467-025-58466-2/note|Wang et al. 2025]] - Mouse brain spatial transcriptomics
- [[../research/databases/spatiallibd|spatialLIBD]] - Human DLPFC spatial data
- [[../data/genage/README|GenAge]] - Conserved aging genes

**Script**: [[../scripts/query_spatiallibd.py|spatialLIBD query script]]

### Inflammation Networks
**Q**: Are inflammatory genes from Wang 2025 enriched in protein interaction networks?

**Data sources**:
- [[10.1038_s41467-025-58466-2/note|Wang et al. 2025]] - White matter inflammation genes
- [[../research/databases/biogrid|BioGRID]] - Protein-protein interactions
- [[../data/genage/README|GenAge]] - Aging-related genes

### Senescence Connections
**Q**: Do senescence markers correlate with spatial inflammation patterns?

**Data sources**:
- [[10.1038_s41467-025-58466-2/note|Wang et al. 2025]] - Spatial inflammation patterns
- [[../data/cellage/README|CellAge]] - Senescence genes and signatures
- [[../research/literature/shvarts-2002-bcl6-senescence|BCL6 paper]] - Senescence signaling

## File Organization

```
papers/
├── INDEX.md                    # This dashboard
├── README.md                   # Project documentation
├── papers.db                   # SQLite tracking database
├── scripts/
│   ├── download_paper.py       # Main download script
│   └── init_database.py        # Database initialization
└── [DOI]/                      # One folder per paper
    ├── paper.pdf
    ├── figure_*.png
    ├── metadata.json
    └── note.md                 # Literature note
```

**Naming convention**: DOI with `/` replaced by `_`
- DOI: `10.1038/s41467-025-58466-2`
- Folder: `10.1038_s41467-025-58466-2/`

## Related Dashboards

- [[../HOME|Projects Hub]] - All projects overview
- [[../research/INDEX|Aging Research Dashboard]] - Main research hub
- [[../immunos-mcp/INDEX|IMMUNOS-MCP Dashboard]]
- [[../immunos81/INDEX|IMMUNOS81 Dashboard]]

## Metadata Standards

All papers have:
- **YAML frontmatter** with full metadata
- **DOI** as unique identifier
- **Author list** (all authors)
- **Journal** and year
- **License** information (CC-BY, etc.)
- **Download date** for provenance
- **Tags** for categorization

Example:
```yaml
---
title: "Paper title"
authors: ["Author 1", "Author 2", ...]
year: 2025
journal: "Nature Communications"
doi: "10.1038/s41467-025-58466-2"
url: "https://doi.org/10.1038/s41467-025-58466-2"
type: journal-article
downloaded: "2025-12-03"
license: "CC-BY-4.0"
tags:
  - paper
  - aging-biology
  - spatial-transcriptomics
---
```

## Storage Summary

| Paper | PDF Size | Figures | Total Size |
|-------|----------|---------|------------|
| Wang et al. 2025 (spatial transcriptomics) | 10.0 MB | 7 (11.4 MB) | 21.4 MB |
| Gray et al. 2025 (midbrain ECM) | 0.4 MB | 6 (6.9 MB) | 7.3 MB |
| Mega-analysis 2025 (memory decline) | 0.5 MB | 4 (est 4 MB) | 4.5 MB |
| Proteomic clocks 2025 | 15.0 MB | 5 (est 5 MB) | 20.0 MB |
| Odoribacter 2025 (gut microbiome) | 6.1 MB | 6 (est 6 MB) | 12.1 MB |
| CAR T therapy 2025 | 39.0 MB | 5 (est 5 MB) | 44.0 MB |
| Iron homeostasis 2025 | 9.2 MB | 5 (est 5 MB) | 14.2 MB |
| **TOTAL** | **80.2 MB** | **38 figures** | **~123.5 MB** |

**Storage strategy**: Download high-value papers with figures. For large datasets, prefer web interfaces or on-demand downloads.

**Restricted access** (not downloaded): See [[RESTRICTED_ACCESS|RESTRICTED_ACCESS.md]] for 2 tracked papers

---

**Created**: 2025-12-03
**Updated**: 2025-12-03 (batch download)
**Location**: `/Users/byron/projects/papers/`
**Database**: `papers.db` (SQLite)
**Total papers**: 7 open access + 2 restricted (tracked only)
**Total figures**: 38
**Total storage**: ~123.5 MB
