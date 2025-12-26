# Citation Management

## Summary
Contains 0 subdirectories and 3 files.


This directory contains bibliographic data in formats compatible with academic reference managers (Zotero, Mendeley, EndNote) and writing tools (LaTeX, Pandoc, Word).

## File Formats

### BibTeX (.bib)
- **Standard format** for academic citations
- Compatible with: LaTeX, Zotero, Mendeley, EndNote, Pandoc
- Human-readable and version-control friendly
- Recommended for maximum compatibility

### CSL JSON (.json)
- **Citation Style Language** format
- Modern standard for web tools and Pandoc
- Used internally by most reference managers
- Good for programmatic access

## Files

### `library.bib`
Master BibTeX file containing all references across all projects.

**Recommended Workflow**:
1. Manage references in Zotero with Better BibTeX plugin
2. Export library to this file with "Keep updated" enabled
3. Zotero automatically syncs changes to this file

### Project-Specific BibTeX Files

- `genage-references.bib` - GenAge database and aging research
- `cellage-references.bib` - CellAge database and cellular senescence research
- `immunos-references.bib` - IMMUNOS-MCP and artificial immune systems
- *(add more as needed)*

### `bibliography.json`
Optional CSL JSON version of the master library for tools that prefer JSON.

### `attachments/`
PDFs, supplementary materials, and other files attached to references.

## Obsidian Integration

### Recommended Plugins

1. **Citations Plugin** (hans/obsidian-citation-plugin)
   - Reads BibTeX and CSL JSON
   - Creates literature notes automatically
   - Insert citations and references

2. **Zotero Integration** (mgmeyers/obsidian-zotero-integration)
   - Requires Better BibTeX for Zotero
   - Direct Zotero connection
   - Import bibliographies and PDF annotations

### Literature Notes

Literature notes are stored in `/research/literature/` with YAML frontmatter containing citation metadata. This format is:
- Human-readable
- Git-friendly
- Queryable with Dataview plugin
- Convertible to BibTeX/CSL JSON

**Format**:
```yaml
---
title: "Paper Title"
authors: [Author1, Author2]
year: 2024
venue: Journal Name
type: journal-article
doi: 10.xxxx/xxxxx
url: https://doi.org/...
tags: [topic1, topic2]
---
```

## Citation Workflow

### Option 1: Zotero + Better BibTeX (Recommended)

1. **Install Zotero**: https://www.zotero.org/
2. **Install Better BibTeX**: https://retorque.re/zotero-better-bibtex/
3. **Configure Export**:
   - Right-click library → "Export Library"
   - Format: "Better BibLaTeX"
   - Check "Keep updated"
   - Save to: `/Users/byron/projects/research/citations/library.bib`
4. **Add References**: Zotero automatically updates the .bib file

### Option 2: Manual BibTeX Editing

Edit .bib files directly in VS Code or any text editor. Use standard BibTeX entry types:

```bibtex
@article{key2024,
  author = {Last, First and Last, First},
  title = {Paper Title},
  journal = {Journal Name},
  year = {2024},
  volume = {52},
  pages = {1--10},
  doi = {10.xxxx/xxxxx}
}
```

### Option 3: Obsidian YAML → BibTeX

Create literature notes with YAML frontmatter, then convert to BibTeX using a script (to be created) or manually.

## Export Compatibility

### Supported Tools

| Tool | BibTeX | CSL JSON | Notes |
|------|--------|----------|-------|
| **Zotero** | ✅ Import/Export | ✅ Export | Best with Better BibTeX |
| **Mendeley** | ✅ Import/Export | ❌ No | Cannot customize citation keys |
| **EndNote** | ✅ Import/Export | ❌ No | Poor LaTeX compatibility |
| **LaTeX** | ✅ Native | ❌ No | BibTeX/BibLaTeX only |
| **Pandoc** | ✅ Yes | ✅ Yes | Both formats supported |
| **Word** | ✅ Via Zotero | ✅ Via Zotero | Requires plugin |

**Recommendation**: Use **BibTeX** as primary format for maximum compatibility.

## Adding New References

### From DOI

```bash
# Use doi2bib.org or curl
curl -LH "Accept: application/x-bibtex" https://doi.org/10.1093/nar/gkad927
```

### From PubMed

Use Zotero browser connector or PubMed export feature (select "MEDLINE" format, convert to BibTeX).

### Manual Entry

Add to appropriate .bib file following BibTeX syntax.

## Quality Control

### Citation Keys

Use consistent format: `authorYEARkeyword`
- Example: `demagalhaes2024hagr`
- Lowercase, no spaces
- Include first author, year, and keyword

### Required Fields by Entry Type

**@article**:
```bibtex
author, title, journal, year, volume, pages, doi (recommended)
```

**@book**:
```bibtex
author/editor, title, year, publisher, isbn (recommended)
```

**@inproceedings**:
```bibtex
author, title, booktitle, year, pages (optional), doi (recommended)
```

## Related Documentation

- **Literature Notes**: `/research/literature/` - Markdown notes on papers
- **Datasets**: `/research/datasets/` - Data source documentation
- **Experiments**: `/research/experiments/` - Citing datasets and papers

---

**Last Updated**: 2025-12-02
**Total References**: 3 (GenAge-related)
**Format**: BibTeX (primary), CSL JSON (optional)

## Directory Map
```
citations/
├── cellage-references.bib
├── genage-references.bib
└── immunos-references.bib
```
