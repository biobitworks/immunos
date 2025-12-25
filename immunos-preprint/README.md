# immunOS Preprint

## Overview

Preprint manuscript for immunOS research.

## Project Structure

```
immunos-preprint/
├── README.md              # This file
├── LINKS.md               # Obsidian-friendly links to project resources
├── immunos-preprint.md    # Main manuscript
├── figures/               # Figure assets
├── supplementary/         # Supplementary materials
├── data/                  # Analysis data
├── scripts/               # Analysis scripts
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

- [immunOS Replication Preprint](../immunos-replication-preprint/) - SciFact baseline replication study
- [IMMUNOS-MCP](../immunos-mcp/) - MCP server implementation
- [IMMUNOS81](../immunos81/) - Medical diagnosis system

## Links

See [LINKS.md](LINKS.md) for Obsidian-friendly navigation.
