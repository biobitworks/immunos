# IMMUNOS Inbox

## Summary
Contains 0 subdirectories and 1 files.


Automated intake system for processing downloaded files and routing them to appropriate destinations.

## Quick Start

```bash
# Preview what will happen
python3 ../scripts/immunos_inbox.py --dry-run --verbose

# Process inbox files
python3 ../scripts/immunos_inbox.py --verbose

# Process without journal logging
python3 ../scripts/immunos_inbox.py --no-journal
```

## File Routing

The inbox processor automatically routes files based on type and filename patterns:

### PDFs
- **Destination**: `/papers/{doi-or-name}/`
- **DOI Extraction**: 
  - From filename: `paper_DOI_10.1234-5678.pdf`
  - From content: Searches first 50KB for DOI patterns
  - Fallback: Uses filename without extension
- **Logged to**: Publications log and daily journal

### Images
- **Patterns**:
  - `*figure*`, `*fig*` → `/assets/figures/`
  - `*diagram*`, `*chart*` → `/assets/diagrams/`
  - `*screenshot*`, `*screen*` → `/assets/screenshots/`
  - `*logo*`, `*icon*` → `/assets/logos/`
  - Default → `/assets/images/`

### Datasets
- **Patterns**:
  - `*research*`, `*experiment*` → `/research/datasets/`
  - `*raw*`, `*processed*` → `/data/`
  - Default → `/data/`
- **Formats**: `.csv`, `.json`, `.jsonl`, `.parquet`, `.tsv`, `.xlsx`

### Citations
- **Destination**: `/research/citations/`
- **Formats**: `.bib`, `.ris`, `.enw`, `.nbib`

## Ignored Files

The following files are automatically ignored:
- `.DS_Store`
- `README.md`
- `.gitkeep`
- `Thumbs.db`

## Workflow

1. Download files to this inbox directory
2. Run `python3 ../scripts/immunos_inbox.py --verbose`
3. Files are automatically:
   - Routed to appropriate destinations
   - Logged to daily journal
   - Added to publications log (if PDFs with DOIs)
   - Directories created as needed

## Logging

### Daily Journal
Entries added to `/daily/YYYY-MM-DD.md`:
```markdown
## Inbox Intake (HH:MM)

- PDF: `paper.pdf` (DOI: 10.1234/example) → `/papers/10.1234_example/paper.pdf`
- Image: `figure_1.png` → `/assets/figures/figure_1.png`
- Dataset: `research_data.csv` → `/research/datasets/research_data.csv`
```

### Publications Log
PDFs with DOIs added to `/docs/reference/publications-log.md`:
```markdown
## YYYY-MM-DD
- paper.pdf
  DOI: https://doi.org/10.1234/example
  Path: /papers/10.1234_example/paper.pdf
  Intake: 2025-12-24T10:30:00
```

## Command Line Options

```
--inbox INBOX          Inbox directory (default: /Users/byron/projects/inbox)
--projects PROJECTS    Projects root (default: /Users/byron/projects)
--dry-run              Preview without moving files
--verbose              Show detailed processing info
--no-journal           Skip daily journal logging
```

## Examples

### Daily Workflow
```bash
# Morning: Check what's in inbox
ls -lh

# Preview processing
python3 ../scripts/immunos_inbox.py --dry-run --verbose

# Process files
python3 ../scripts/immunos_inbox.py --verbose
```

### Research Paper Workflow
```bash
# Download paper to inbox with DOI in filename
# e.g., smith2025_DOI_10.1234-example.2025.pdf

# Process (automatically extracts DOI and creates /papers/10.1234_example.2025/)
python3 ../scripts/immunos_inbox.py --verbose

# Check daily journal for confirmation
cat ../daily/$(date +%Y-%m-%d).md
```

### Batch Import
```bash
# Move multiple files to inbox
cp ~/Downloads/*.pdf .

# Dry run to verify routing
python3 ../scripts/immunos_inbox.py --dry-run --verbose

# Process all
python3 ../scripts/immunos_inbox.py --verbose
```

## Dependencies

Uses only Python standard library:
- `pathlib` - Path operations
- `re` - DOI pattern matching
- `json` - Data serialization
- `datetime` - Timestamps
- `argparse` - Command line interface

No external dependencies required.

## Notes

- **DOI Extraction**: Uses basic text search on first 50KB of PDF. For better extraction, consider installing `PyPDF2` or `pdfplumber`
- **Dry Run**: Always test with `--dry-run` first when processing important files
- **Backups**: Original files are moved (not copied). Ensure you have backups of important files
- **Directory Creation**: Destination directories are created automatically as needed
- **File Conflicts**: If destination file exists, operation will fail (file won't be overwritten)

## Troubleshooting

### DOI not detected
- Ensure DOI is in filename: `paper_DOI_10.1234-example.pdf`
- Check PDF contains DOI in first page
- Verify DOI format: `10.####/...`

### File not moved
- Check file extension is recognized
- Verify no file exists at destination
- Run with `--verbose` to see error messages

### Journal not updated
- Check `/daily/YYYY-MM-DD.md` exists
- Verify you're not using `--no-journal`
- Run without `--dry-run`

---

**Part of**: IMMUNOS Research Integrity System  
**Script**: `/scripts/immunos_inbox.py`  
**Last Updated**: 2025-12-24

## Directory Map
```
inbox/
└── cursor_projects_intro.md
```
