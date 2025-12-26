# immunOS Preprint

## Summary
Preprint workspace with manuscript, replication draft, and supporting data/scripts.


## Overview

Preprint manuscript for immunOS research.

## Project Structure

```
immunos-preprint/
├── README.md              # This file
├── LINKS.md               # Obsidian-friendly links to project resources
├── immunos-preprint.md    # Main manuscript
├── immunos-preprint-v1.md # Replication draft (SciFact baseline)
├── figures/               # Figure assets
├── supplementary/         # Supplementary materials
├── data/                  # Analysis data
├── scripts/               # Analysis scripts
├── journal/               # Project log entries
└── bibliography.bib       # BibTeX references
```

## Build Instructions

### Prerequisites

- LaTeX distribution (for PDF compilation)
- Pandoc (for format conversion)
- BibTeX (for reference management)

### Compile to PDF

```bash
pandoc immunos-preprint.md \
  --bibliography=bibliography.bib \
  --citeproc \
  -o immunos-preprint.pdf
```

### Compile to DOCX

```bash
pandoc immunos-preprint.md \
  --bibliography=bibliography.bib \
  --citeproc \
  -o immunos-preprint.docx
```

## Development Workflow

1. Edit `immunos-preprint.md` for manuscript content
2. Add figures to `figures/`
3. Add citations to `bibliography.bib`
4. Run analysis scripts in `scripts/`
5. Store results in `data/`
6. Add supplementary materials to `supplementary/`

## Status

**Type**: Preprint
**Status**: Draft
**Last Updated**: 2025-12-24

## Related Projects

- [IMMUNOS-MCP](../immunos-mcp/) - MCP server implementation
- [IMMUNOS81](../immunos81/) - Medical diagnosis system

## Links

See [LINKS.md](LINKS.md) for Obsidian-friendly navigation.

## Directory Map
```
immunos-preprint/
├── data/
├── figures/
├── journal/
├── scripts/
├── supplementary/
├── INDEX.md
├── LINKS.md
├── bibliography.bib
├── immunos-preprint-v1.md
└── immunos-preprint.md
```
