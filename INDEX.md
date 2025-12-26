# Projects Index

Comprehensive index of all projects in the workspace.

**Last Updated**: 2025-12-24

---

## Quick Start

- [[HOME|Projects Hub]] - Main navigation dashboard
- [[CLAUDE|Claude Context]] - AI assistant instructions

---

## Active Projects

### Research & Publications

- [[immunos-preprint/INDEX|immunOS Preprint]] - Main immunOS preprint manuscript (Draft)
- [[immunos-preprint/README|immunOS Preprint]] - Main manuscript + replication draft
- [[papers/INDEX|Papers Archive]] - Academic papers with figures (32 papers tracked)
- [[research/INDEX|Aging Biology Research]] - Longevity science research hub

### Software Projects

- [[immunos-mcp/INDEX|IMMUNOS-MCP]] - MCP server for code security (MVP Complete)
- [[immunos81/INDEX|IMMUNOS81]] - Medical diagnosis system (Active Development)

### Professional Development

- [[portfolio/README|Portfolio & Professional Development]] - Career planning and strategic goals
- [[biobitworks-logistics/README|Biobitworks Logistics]] - Community lab network
- [[research/projects/bioviztech/README|BioViz Research Project]] - ECM research and lab setup

---

## Project Categories

### By Type
- **Preprints**: immunos-preprint
- **Research**: Aging Biology Research, Papers Archive, BioViz Research
- **Software**: IMMUNOS-MCP, IMMUNOS81
- **Professional**: Portfolio, Biobitworks Logistics

### By Status
- **Active Development**: IMMUNOS81, BioViz Research, Biobitworks Logistics
- **MVP Complete**: IMMUNOS-MCP
- **Draft**: immunos-preprint (includes replication v1)
- **Ongoing**: Papers Archive, Aging Biology Research, Portfolio

---

## IMMUNOS System

### Core Scripts
Located in `/scripts/`:
- `immunos_recover.py` - Context recovery
- `immunos_journal.py` - Daily journaling
- `immunos_memory.py` - Memory management
- `immunos_snapshot.py` - Session snapshots
- `immunos_dashboard.py` - Web dashboard
- `immunos_cite_verify.py` - Citation verification
- `immunos_nk_scan.py` - Anomaly detection
- `immunos_token_analyzer.py` - Token optimization

### System Data
Located in `.immunos/`:
- `memory/` - Conversations, decisions, snapshots
- `journal/` - Daily journals
- `recovery/` - Recovery summaries
- `db/` - SQLite databases

---

## Quick Actions

### Start Session
```bash
python3 scripts/immunos_recover.py && cat .immunos/recovery/CONTEXT_RECOVERY.md
```

### Create Snapshot
```bash
python3 scripts/immunos_snapshot.py create --trigger manual --summary "Current work"
```

### Check Health
```bash
python3 scripts/immunos_memory.py stats
python3 scripts/immunos_todo.py list --overdue
```

---

## Workspace Structure

```
/Users/byron/projects/
├── .immunos/                    # IMMUNOS system data
├── biobitworks-logistics/       # Community lab network
├── immunos-mcp/                 # MCP server project
├── immunos-preprint/            # Main immunOS preprint
├── immunos-preprint/# Main preprint + SciFact replication study
├── immunos81/                   # Medical diagnosis system
├── papers/                      # Papers archive
├── portfolio/                   # Professional development
├── research/                    # Aging biology research
│   ├── experiments/            # Research experiments
│   └── projects/bioviztech/    # BioViz project
├── scripts/                     # IMMUNOS system scripts
├── CLAUDE.md                    # Claude context instructions
├── HOME.md                      # Projects hub dashboard
└── INDEX.md                     # This file
```

---

## External Links

### Resources
- [bioRxiv](https://www.biorxiv.org/) - Preprint server
- [PubMed](https://pubmed.ncbi.nlm.nih.gov/) - Literature search
- [Zotero](https://www.zotero.org/) - Reference management

### Tools
- [Pandoc](https://pandoc.org/) - Document conversion
- [Ollama](https://ollama.ai/) - Local AI models
- [Claude Code](https://claude.com/claude-code) - AI assistant

---

**Workspace**: `/Users/byron/projects/`
**Platform**: macOS (Darwin 25.1.0)
