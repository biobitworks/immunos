# Academic Papers Archive

Automated system for downloading open-access papers with figures, organizing by DOI, and integrating with Obsidian research vault.

## Quick Start

```bash
cd /Users/byron/projects/papers

# Download single paper from DOI
python scripts/download_paper.py 10.1038/s41467-025-58466-2

# Download from URL
python scripts/download_paper.py https://doi.org/10.1038/s41467-025-58466-2

# Batch download multiple papers
python scripts/download_paper.py \
  10.1038/s41467-025-66434-z \
  10.1038/s41467-025-66354-y \
  10.1038/s43587-025-01016-8

# Mixed URLs and DOIs
python scripts/download_paper.py \
  https://doi.org/10.1038/s41467-025-58466-2 \
  10.1038/s43587-025-01016-8
```

## Project Structure

```
papers/
├── INDEX.md                    # Obsidian dashboard
├── README.md                   # This file
├── papers.db                   # SQLite tracking database
├── scripts/
│   ├── download_paper.py       # Main download script
│   └── init_database.py        # Initialize database
└── [DOI]/                      # One folder per paper
    ├── paper.pdf
    ├── figure_*.{png,jpg,svg}
    ├── metadata.json
    └── note.md                 # Literature note
```

## Folder Naming

Papers are organized by DOI with `/` replaced by `_`:
- DOI: `10.1093/nar/gkad927` → Folder: `10.1093_nar_gkad927/`
- DOI: `10.1186/s13059-020-01990-9` → Folder: `10.1186_s13059-020-01990-9/`

## Features

### Automated Download
- ✅ **Batch downloads** - Download multiple papers at once
- ✅ **License checking** - Automatically checks open access status
- ✅ PDF download from multiple publishers
- ✅ Figure extraction (all formats)
- ✅ Metadata generation via Crossref API
- ✅ Literature note creation
- ✅ SQLite database tracking
- ✅ Restricted paper tracking (non-open access)

### Supported Publishers
- ✅ Nature family (Nature, Nature Communications, Nature Aging, etc.)
- 🚧 PubMed Central (PMC) - planned
- 🚧 bioRxiv / medRxiv - planned
- 🚧 arXiv - planned
- 🚧 Direct PDF URLs - planned

### License Support
**Automatically downloads**:
- ✅ CC-BY (Creative Commons Attribution)
- ✅ CC-BY-NC (Non-commercial)
- ✅ CC-BY-SA (Share-alike)
- ✅ CC-BY-NC-ND (Non-commercial, no derivatives)
- ✅ CC0 (Public domain)

**Tracks but does not download**:
- ❌ Paywalled papers
- ❌ TDM-only licenses
- ❌ Subscription-only access

Restricted papers are tracked in `RESTRICTED_ACCESS.md` for reference.

### Metadata Tracking

Each paper has `metadata.json`:
```json
{
  "doi": "10.1093/nar/gkad927",
  "title": "Paper title",
  "authors": ["Author 1", "Author 2"],
  "journal": "Journal name",
  "year": 2024,
  "url": "https://...",
  "downloaded": "2025-12-03",
  "figures": ["figure_1.png", "figure_2.png"],
  "note_path": "10.1093_nar_gkad927/note.md"
}
```

### SQLite Database

Track all papers in `papers.db`:
- Papers table (DOI, title, authors, etc.)
- Figures table (linked to papers)
- Query capabilities for searching

### Literature Notes

Auto-generated markdown notes:
- YAML frontmatter for Obsidian
- Abstract section
- Figure references
- Links to PDF and images
- Wiki-links to related research

## Database Schema

```sql
CREATE TABLE papers (
    doi TEXT PRIMARY KEY,
    title TEXT NOT NULL,
    authors TEXT,
    journal TEXT,
    year INTEGER,
    url TEXT,
    pdf_path TEXT,
    download_date TEXT,
    note_path TEXT
);

CREATE TABLE figures (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    paper_doi TEXT,
    figure_num INTEGER,
    filename TEXT,
    file_path TEXT,
    FOREIGN KEY (paper_doi) REFERENCES papers(doi)
);
```

## Obsidian Integration

### Dashboard

Open [[INDEX|INDEX.md]] to see:
- All downloaded papers
- Search and filter capabilities
- Links to literature notes
- Statistics

### Literature Notes

Notes are compatible with research vault:
- Cross-reference with researchers
- Tag with relevant topics
- Link to databases and datasets
- Appear in network graph

## Usage Examples

### Example 1: Single Paper Download

```bash
cd /Users/byron/projects/papers

python scripts/download_paper.py 10.1038/s41467-025-58466-2
```

**Output**:
```
======================================================================
Checking licenses for 1 papers...
======================================================================

✅ 10.1038/s41467-025-58466-2
   Spatial transcriptomics of the aging mouse brain reveals origins...
   License: https://creativecommons.org/licenses/by/4.0

======================================================================
License Check Summary:
  ✅ Open Access: 1 papers
  ❌ Restricted: 0 papers
======================================================================

Downloading 1 open access papers...
```

**Creates**:
- `10.1038_s41467-025-58466-2/paper.pdf` (10 MB)
- `10.1038_s41467-025-58466-2/figure_*.png` (7 figures)
- `10.1038_s41467-025-58466-2/metadata.json`
- `10.1038_s41467-025-58466-2/note.md`
- Entry in `papers.db` database

### Example 2: Batch Download with License Checking

```bash
python scripts/download_paper.py \
  10.1038/s41467-025-66434-z \
  10.1038/s41467-025-66354-y \
  10.1038/s43587-025-01016-8 \
  10.1038/s43587-024-00616-0 \
  10.1038/s43587-025-01043-5
```

**Output**:
```
======================================================================
Checking licenses for 5 papers...
======================================================================

✅ 10.1038/s41467-025-66434-z
   Midbrain extracellular matrix and microglia are associated with...
   License: https://creativecommons.org/licenses/by/4.0

✅ 10.1038/s41467-025-66354-y
   Vulnerability to memory decline in aging revealed by a mega-ana...
   License: https://creativecommons.org/licenses/by/4.0

✅ 10.1038/s43587-025-01016-8
   Organ-specific proteomic aging clocks predict disease and longe...
   License: https://creativecommons.org/licenses/by/4.0

❌ 10.1038/s43587-024-00616-0 - RESTRICTED
   Nature of epigenetic aging from a single-cell perspective
   License: https://www.springernature.com/gp/researchers/text-and-data-mining

❌ 10.1038/s43587-025-01043-5 - RESTRICTED
   Protein restriction reshapes aging across organs
   License: https://www.springernature.com/gp/researchers/text-and-data-mining

======================================================================
License Check Summary:
  ✅ Open Access: 3 papers
  ❌ Restricted: 2 papers
======================================================================

Downloading 3 open access papers...
[Downloads only the 3 open access papers]

⚠️  2 papers were RESTRICTED and not downloaded
✓ Updated RESTRICTED_ACCESS.md
```

### Example 3: Mixed URLs and DOIs

```bash
python scripts/download_paper.py \
  https://doi.org/10.1038/s41467-025-58466-2 \
  10.1038/s43587-025-01016-8 \
  https://www.nature.com/articles/s41467-025-66354-y
```

The script automatically extracts DOIs from URLs and processes them.

### Example 4: Query Database

```bash
# List all papers
sqlite3 papers.db "SELECT title, journal, year FROM papers ORDER BY year DESC"

# Count papers by journal
sqlite3 papers.db "SELECT journal, COUNT(*) FROM papers GROUP BY journal"

# Papers with most figures
sqlite3 papers.db "
  SELECT p.title, COUNT(f.id) as fig_count
  FROM papers p
  LEFT JOIN figures f ON p.doi = f.paper_doi
  GROUP BY p.doi
  ORDER BY fig_count DESC
"
```

### Example 5: Python Database Query

```python
import sqlite3
import json

conn = sqlite3.connect('papers.db')
cursor = conn.cursor()

# Find all papers from 2025
cursor.execute("SELECT title, doi, authors FROM papers WHERE year = 2025")
for title, doi, authors_json in cursor.fetchall():
    authors = json.loads(authors_json)
    print(f"{title}")
    print(f"  DOI: {doi}")
    print(f"  Authors: {len(authors)} authors")
    print()
```

## Requirements

**Python packages**:
```bash
pip install requests beautifulsoup4 pypdf pillow
```

**Optional**:
- `pdfplumber` - Better PDF text extraction
- `scholarly` - Google Scholar API

## File Organization

### Good Practices

✅ **Do**:
- Use official DOI when available
- Include all figures
- Add descriptive tags
- Create literature notes

❌ **Don't**:
- Download paywalled content
- Rename files manually
- Edit metadata.json by hand
- Duplicate papers

## Troubleshooting

### Paper Won't Download

**Issue**: URL not supported
**Solution**: Try direct PDF URL or DOI-based download

### Figures Missing

**Issue**: Publisher-specific format
**Solution**: Manual download, then add to folder

### Database Locked

**Issue**: SQLite file in use
**Solution**: Close other connections to papers.db

## Current Limitations & Known Issues

⚠️ **Important**: As of December 2025, SpringerNature has implemented aggressive bot detection that prevents automated downloads.

### Status
- **6 papers successfully downloaded** (19%) - Downloaded before protections engaged
- **26 papers failed** (81%) - Blocked by bot detection and authentication requirements
- See [DOWNLOAD_STATUS.md](DOWNLOAD_STATUS.md) and [INVESTIGATION_SUMMARY.md](INVESTIGATION_SUMMARY.md) for details

### Specific Issues

1. **PDF Access - ALL PAYWALLED**
   - All Nature journal PDFs redirect to authentication portal
   - Requires institutional login despite "open access" licenses
   - **Workaround**: Use institution VPN + manual download

2. **Figure Downloads - Bot Detection**
   - SpringerNature CDN blocks automated requests after ~10-20 downloads
   - IP/session gets blocked for 24-48 hours
   - **Workaround**: Manual downloads or official API key

3. **Rate Limiting**
   - Even with 2s delays between requests, blocks occur
   - Both Python (requests/subprocess) and bash curl affected
   - **Workaround**: Wait 24-48hrs between batch attempts

### Recommended Approaches

1. **Institutional Access** (Best)
   - Use university VPN
   - Download through authenticated browser
   - Full access to PDFs and figures

2. **SpringerNature API**
   - Apply for API key: https://dev.springernature.com/
   - Legitimate bulk download access
   - Proper rate limiting and authentication

3. **Alternative Sources**
   - PubMed Central (PMC) mirrors
   - Author preprints (arXiv/bioRxiv)
   - Contact authors directly

### Scripts Available

- `scripts/download_paper.py` - Original script (currently blocked)
- `scripts/batch_download_figures.sh` - Batch retry script with rate limiting
- `scripts/download_figures_only.sh` - Single paper figure downloader

**Note**: These scripts work but will trigger rate limits after ~10-20 requests.

## Integration with Research Vault

### Link to Existing Notes

```markdown
See [[research/literature/shvarts-2002-bcl6-senescence|BCL6 paper]]
```

### Tag Organization

Use consistent tags:
- #aging-biology
- #cellular-senescence
- #database
- #methodology
- #review

### Cross-Reference

Link papers to:
- [[research/researchers/|Researchers]] (if author tracked)
- [[research/databases/|Databases]] (if paper describes dataset)
- [[research/organizations/|Organizations]] (author institutions)

## Future Enhancements

- [ ] Batch download from list
- [ ] Citation network visualization
- [ ] PDF annotation import
- [ ] Automatic citation generation
- [ ] Integration with Zotero
- [ ] Full-text search

## Related Documentation

- [[INDEX|Papers Dashboard]]
- [[../research/INDEX|Aging Research Dashboard]]
- [[../HOME|Projects Hub]]

---

**Created**: 2025-12-03
**Location**: `/Users/byron/projects/papers/`
**Database**: `papers.db` (SQLite)
**License**: For personal research use
