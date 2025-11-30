---
type: vault-documentation
tags: [obsidian, research, documentation]
created: 2025-11-30
---

# Obsidian Research Vault

This directory (`/Users/byron/projects/`) is configured as an **Obsidian vault** for research documentation and knowledge management.

## 🎯 Purpose

Track daily work, experiments, code changes, and research progress in a structured, interconnected knowledge base optimized for:
- Writing research papers
- Reproducing experiments
- Documenting technical decisions
- Maintaining code documentation
- Building research literature database

## 📂 Vault Structure

```
projects/
├── daily/                    # Daily journal entries
│   └── 2025-11-30.md        # One file per day
├── projects/                 # Project-specific content
│   ├── immunos-mcp/
│   │   ├── journal/         # Daily project logs
│   │   ├── code-mirror/     # Python → Markdown mirrors
│   │   ├── api/             # API documentation
│   │   └── diagrams/        # Mermaid diagrams
│   ├── immunos81/
│   └── experiments/
├── research/
│   ├── experiments/         # Experiment logs
│   ├── literature/          # Paper summaries
│   ├── ideas/              # Hypotheses
│   └── metrics/            # Performance data
├── templates/               # Note templates
└── scripts/                 # Automation
```

## 🚀 Quick Start

### Open Vault in Obsidian

1. Launch Obsidian
2. "Open folder as vault"
3. Select: `/Users/byron/projects/`

### Daily Workflow

```bash
# Create today's note
./scripts/create-daily-note.sh

# After coding, update vault
./scripts/update-vault.sh all
```

## 📊 Current Content

**As of 2025-11-30**:

| Type | Count | Location |
|------|-------|----------|
| Daily Notes | 1 | `daily/` |
| Project Journals | 1 | `projects/*/journal/` |
| Experiments | 1 | `research/experiments/` |
| Literature Notes | 3 | `research/literature/` |
| Code Mirrors | 38 | `projects/*/code-mirror/` |
| API Docs | 18 | `projects/*/api/` |
| Diagrams | 2 | `projects/*/diagrams/` |
| Templates | 2 | `templates/` |
| Scripts | 5 | `scripts/` |

**Total**: ~70 markdown files

## 🔗 Key Entry Points

### Today's Work
- [[daily/2025-11-30]] - Daily journal
- [[projects/immunos-mcp/journal/2025-11-30]] - Project details

### QML-AiNet Research
- [[research/experiments/qml-ainet-validation-2025-11-30]] - Validation results (92.9% accuracy)
- [[research/literature/pang-coghill-2015-qml-ainet]] - Paper summary
- [[projects/immunos-mcp/code-mirror/src/algorithms/qml_ainet]] - Implementation
- [[projects/immunos-mcp/diagrams/qml-ainet-algorithm]] - Algorithm flow

### System Architecture
- [[projects/immunos-mcp/diagrams/system-overview]] - Full system
- [[projects/immunos-mcp/api/README]] - API index

### Automation
- [[scripts/README]] - All automation scripts

## 🛠️ Automation Scripts

### 1. Update Entire Vault
```bash
./scripts/update-vault.sh all        # All projects
./scripts/update-vault.sh immunos-mcp # Single project
```

**Updates**:
- Code mirrors (Python → Markdown)
- API documentation (from docstrings)
- Vault statistics

### 2. Code Documentation
```bash
# Mirror code
python3 scripts/sync-code-to-obsidian.py immunos-mcp

# Generate API docs
python3 scripts/generate-api-docs.py immunos-mcp
```

### 3. Create Notes
```bash
# Daily note
./scripts/create-daily-note.sh

# Experiment log
./scripts/create-experiment-log.sh experiment-name
```

## 📝 Templates

### Daily Note (`templates/daily-note.md`)

**Auto-fills**:
- Today's date
- Yesterday/tomorrow links
- Project journal links

**Sections**:
- Summary & key achievements
- Project journals (with links)
- Files created/modified
- Experiments conducted
- Technical decisions
- Performance metrics
- Challenges & solutions
- Next steps
- Reflections

### Experiment Log (`templates/experiment-log.md`)

**Structure**:
- Hypothesis
- Method (parameters, datasets, metrics)
- Results (tables + detailed findings)
- Analysis (performance, key findings, comparison)
- Reproducibility (code, commands, environment)
- Conclusions (validation, contributions, limitations)
- Next experiments

## 🎓 How to Use

### For Research Papers

The vault provides everything needed for publication:

1. **Experiment Logs** → Results section
   - Reproducible methodology
   - Performance tables
   - Statistical analysis

2. **Literature Notes** → Related Work section
   - Paper summaries
   - Comparisons
   - Evolution of ideas

3. **Daily Notes** → Research timeline
   - What was tried when
   - Why decisions were made
   - How challenges were solved

4. **Code + API Docs** → Implementation section
   - Algorithm pseudocode
   - Architecture diagrams
   - Example usage

### For Code Documentation

- **Code Mirrors**: Source code with syntax highlighting
- **API Docs**: Auto-generated from docstrings (classes, functions, parameters)
- **Diagrams**: Mermaid visualizations (flowcharts, sequences, graphs)
- **Examples**: Usage demonstrations

### For Daily Work Tracking

1. **Morning**: Open today's daily note
2. **During work**: Document decisions, experiments, code changes
3. **Evening**: Link to project journals, experiments, new files
4. **Weekly**: Run `update-vault.sh` to sync everything

## 🔍 Obsidian Features to Use

### Graph View
- Daily notes form timeline backbone
- Projects as hubs
- Literature creates research lineage
- Experiments connect to code

### Search
- Cmd/Ctrl+O: Quick file switcher
- Cmd/Ctrl+Shift+F: Full-text search
- Tag search: `tag:#qml-ainet`

### Links
- Internal: `[[path/to/file|Display Name]]`
- With line numbers: `file.md#heading`
- Block references: `[[file^block-id]]`

### Tags
Vault uses hierarchical tags:
- `#immunos` / `#qml-ainet` / `#opt-ainet` - Algorithms
- `#experiment` / `#validation` - Research
- `#api` / `#documentation` - Code
- `#daily-log` / `#research` - Organization

## 📈 Vault Statistics Script

The `update-vault.sh` script reports:

```
╔═══════════════════════════════════════╗
║   Vault Statistics                    ║
╚═══════════════════════════════════════╝

  Daily Notes:       1
  Project Journals:  1
  Experiments:       1
  Literature Notes:  3
  Code Mirrors:      38
  API Docs:          18

  Total Documents:   62
```

## 🎯 Research Goals Tracked

### Completed ✅
- QML-AiNet implementation (first public code)
- Validation on 7 datasets (92.9% accuracy)
- LLM integration (6 agents, 3 models)
- Obsidian vault setup

### In Progress 🔄
- API documentation generation
- Literature database expansion
- Experiment scaling tests

### Planned 📋
- Network suppression integration
- Real-world dataset validation
- Research paper writing
- Offline deployment bundle

See: [[daily/2025-11-30#Next Steps]]

## 💡 Tips

### Linking Strategy
- Link experiments to daily notes for timeline
- Link code to API docs for details
- Link literature to experiments for context
- Link diagrams to implementations

### Organization
- One daily note per day (no exceptions)
- One experiment log per experiment
- One literature note per paper
- Project journals capture technical depth

### Maintenance
- **Daily**: Create daily note, document work
- **Weekly**: Run `update-vault.sh all`
- **Monthly**: Review graph, consolidate notes

## 🔄 Integration with Git

The vault is git-friendly:

```bash
cd /Users/byron/projects
git add daily/ research/ projects/
git commit -m "Daily log and experiment results"
git push
```

**Tip**: `.gitignore` code mirrors (they're auto-generated from source)

## 📚 Documentation

- [[scripts/README|Automation Scripts Documentation]]
- [[templates/daily-note|Daily Note Template]]
- [[templates/experiment-log|Experiment Log Template]]
- [[research/literature/README|Literature Index]]
- [[projects/immunos-mcp/api/README|API Documentation Index]]

## ❓ Troubleshooting

### Scripts Not Executable
```bash
chmod +x scripts/*.sh scripts/*.py
```

### Code Mirrors Out of Date
```bash
python3 scripts/sync-code-to-obsidian.py immunos-mcp
```

### Broken Links
- Use Obsidian's "Detect broken links" (Cmd+P)
- Check file paths are relative to vault root
- Ensure referenced files exist

### Obsidian Not Finding Vault
- Vault root must be `/Users/byron/projects/`
- Check Settings → Files & Links → "Detect all file extensions"
- Reload vault (Cmd+R)

## 🌐 Export & Sharing

### Static Site
- Obsidian Publish (official)
- Jekyll/Hugo with markdown files
- MkDocs with Material theme

### PDF Export
- Per-note: Right-click → Export to PDF
- Batch: Use Pandoc with scripts

### Markdown Archive
Files are plain markdown - share directly or zip folder

## 🏆 Success Criteria

This vault succeeds when:
- ✅ You can write a research paper entirely from vault notes
- ✅ Any experiment can be reproduced from documentation
- ✅ Every decision has traceable context
- ✅ Code changes are automatically documented
- ✅ Daily progress is consistently captured

## 🙏 Tools Used

- **[Obsidian](https://obsidian.md/)** - Knowledge base
- **Python 3** - Automation scripts
- **Bash** - Utility scripts
- **Mermaid** - Diagram generation
- **Markdown** - Universal format
- **Git** - Version control

---

**Created**: 2025-11-30
**Last Updated**: 2025-11-30
**Vault Version**: 1.0
**Total Notes**: ~70

**Start here**: [[daily/2025-11-30|Today's Work]]
