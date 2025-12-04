# Academic Papers Archive

Automated system for downloading open-access papers with figures, organizing by DOI, and integrating with Obsidian research vault.

## Quick Start

```bash
cd /Users/byron/projects/papers

# Download paper from URL
python scripts/download_paper.py "https://pmc.ncbi.nlm.nih.gov/articles/PMC7737760/"

# Download by DOI
python scripts/download_paper.py --doi "10.1093/nar/gkad927"

# With custom tags
python scripts/download_paper.py "URL" --tags aging-biology,database
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
- ✅ PDF download from multiple publishers
- ✅ Figure extraction (all formats)
- ✅ Metadata generation
- ✅ Literature note creation

### Supported Publishers
- PubMed Central (PMC)
- bioRxiv / medRxiv
- arXiv
- Direct PDF URLs
- *(More can be added)*

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

### Example 1: Download GenAge Paper

```bash
python scripts/download_paper.py \
  "https://academic.oup.com/nar/article/52/D1/D900/7449493" \
  --tags aging-biology,database,genage
```

**Creates**:
- `10.1093_nar_gkad927/paper.pdf`
- `10.1093_nar_gkad927/figure_*.png`
- `10.1093_nar_gkad927/metadata.json`
- `10.1093_nar_gkad927/note.md`

### Example 2: Download from PubMed Central

```bash
python scripts/download_paper.py \
  "https://pmc.ncbi.nlm.nih.gov/articles/PMC7737760/"
```

### Example 3: Query Database

```python
import sqlite3

conn = sqlite3.connect('papers.db')
cursor = conn.cursor()

# Find all papers from 2024
cursor.execute("SELECT title, doi FROM papers WHERE year = 2024")
for title, doi in cursor.fetchall():
    print(f"{title} ({doi})")
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
