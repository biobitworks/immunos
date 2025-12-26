# immunOS Preprint - Dashboard

## Project Overview

**Type**: Preprint manuscript
**Status**: Draft
**Last Updated**: 2025-12-24

Main immunOS preprint for publication.

---

## Quick Navigation

### Core Files
- [[immunos-preprint|Main Manuscript]] - Primary document
- [[immunos-preprint-v1|Replication Draft]] - SciFact baseline study
- [[README|Project Overview]] - Build instructions and setup
- [[LINKS|Links]] - Obsidian-friendly navigation
- [[bibliography|References]] - BibTeX bibliography

### Directories
- `figures/` - Figure assets for manuscript
- `supplementary/` - Supplementary materials
- `data/` - Analysis data and results
- `scripts/` - Analysis and processing scripts
- `journal/` - Project log entries

---

## Manuscript Sections

1. Abstract
2. Introduction
3. Methods
4. Results
5. Discussion
6. References

---

## Build Commands

### PDF Output
```bash
pandoc immunos-preprint.md \
  --bibliography=bibliography.bib \
  --citeproc \
  -o immunos-preprint.pdf
```

### DOCX Output
```bash
pandoc immunos-preprint.md \
  --bibliography=bibliography.bib \
  --citeproc \
  -o immunos-preprint.docx
```

---

## Related Projects

- [[../immunos-mcp/INDEX|IMMUNOS-MCP]] - MCP server
- [[../immunos81/INDEX|IMMUNOS81]] - Medical diagnosis system
- [[../HOME|Projects Hub]] - Main navigation

---

## Status Tracking

**Current Phase**: Template setup complete

**Next Steps**:
- [ ] Define title and author information
- [ ] Draft abstract
- [ ] Write introduction
- [ ] Document methods
- [ ] Add results
- [ ] Write discussion
- [ ] Add references to bibliography.bib
- [ ] Create figures
- [ ] Compile first draft

---

**Dashboard Created**: 2025-12-24
