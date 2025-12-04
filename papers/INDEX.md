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

**Total Papers**: 1
**Total Figures**: 7
**Latest Download**: 2025-12-03

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

---

## By Topic

### Aging Biology
- [[10.1038_s41467-025-58466-2/note|Wang et al. 2025 - Brain aging spatial transcriptomics]]

### Spatial Transcriptomics
- [[10.1038_s41467-025-58466-2/note|Wang et al. 2025 - Mouse brain white matter]]

### Inflammation
- [[10.1038_s41467-025-58466-2/note|Wang et al. 2025 - White matter inflammation origins]]

## By Year

### 2025
- **Nature Communications**: [[10.1038_s41467-025-58466-2/note|Spatial transcriptomics of aging mouse brain]]

## By Journal

### Nature Communications
- [[10.1038_s41467-025-58466-2/note|Wang et al. 2025]] - Spatial transcriptomics, aging, white matter

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
| Wang et al. 2025 | 10 MB | 7 (11.4 MB) | 21.4 MB |
| **TOTAL** | **10 MB** | **7 figures** | **21.4 MB** |

**Storage strategy**: Download high-value papers with figures. For large datasets, prefer web interfaces or on-demand downloads.

---

**Created**: 2025-12-03
**Location**: `/Users/byron/projects/papers/`
**Database**: `papers.db` (SQLite)
**Total papers**: 1
**Total figures**: 7
