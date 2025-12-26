---
type: documentation
tags: [scripts, automation, maintenance]
created: 2025-11-30
---

# Vault Automation Scripts

Scripts for maintaining and updating the Obsidian research vault.

## Available Scripts

### 1. `update-vault.sh` - Full Vault Update

Updates all generated content in the vault.

**Usage**:
```bash
# Update all projects
./scripts/update-vault.sh all

# Update specific project
./scripts/update-vault.sh immunos-mcp
```

**What it does**:
- Mirrors code files to markdown
- Regenerates API documentation
- Updates vault statistics
- Reports changes

**When to run**:
- After significant code changes
- Before writing research papers
- Weekly maintenance

---

### 2. `sync-code-to-obsidian.py` - Code Mirror

Mirrors Python source code as markdown files for Obsidian viewing.

**Usage**:
```bash
python3 scripts/sync-code-to-obsidian.py <project-name>

# Example
python3 scripts/sync-code-to-obsidian.py immunos-mcp
```

**Features**:
- Hash-based change detection (only updates modified files)
- Preserves directory structure
- Adds YAML frontmatter with metadata
- Syntax highlighting with code blocks
- Creates index files

**Output**: `projects/<project>/code-mirror/`

---

### 3. `generate-api-docs.py` - API Documentation

Extracts docstrings from Python files and generates markdown documentation.

**Usage**:
```bash
python3 scripts/generate-api-docs.py <project-name>

# Example
python3 scripts/generate-api-docs.py immunos-mcp
```

**Features**:
- Parses Python AST for docstrings
- Documents classes, methods, functions
- Extracts type annotations
- Generates navigable index
- Links to source code

**Output**: `projects/<project>/api/`

---

### 4. `create-daily-note.sh` - Daily Note Creation

Creates a new daily note from template.

**Usage**:
```bash
# Create for today
./scripts/create-daily-note.sh

# Create for specific date
./scripts/create-daily-note.sh 2025-12-01
```

**Features**:
- Uses template: `templates/daily-note.md`
- Auto-fills date placeholders
- Calculates yesterday/tomorrow links
- Optional: Opens in Obsidian

**Output**: `daily/YYYY-MM-DD.md`

---

### 5. `create-experiment-log.sh` - Experiment Log Creation

Creates a new experiment log from template.

**Usage**:
```bash
./scripts/create-experiment-log.sh <experiment-name>

# Example
./scripts/create-experiment-log.sh qml-scaling-test
```

**Features**:
- Uses template: `templates/experiment-log.md`
- Auto-fills experiment name and date
- Creates with standard structure
- Optional: Opens in editor

**Output**: `research/experiments/<name>-YYYY-MM-DD.md`

---

## Automation Workflow

### Daily Workflow
```bash
# 1. Create today's daily note
./scripts/create-daily-note.sh

# 2. Work on code...

# 3. Update vault before end of day
./scripts/update-vault.sh all
```

### Experiment Workflow
```bash
# 1. Create experiment log
./scripts/create-experiment-log.sh my-experiment

# 2. Fill in hypothesis and method

# 3. Run experiment

# 4. Document results in log

# 5. Link from daily note
```

### Code Documentation Workflow
```bash
# After adding/updating code with docstrings:
python3 scripts/sync-code-to-obsidian.py immunos-mcp
python3 scripts/generate-api-docs.py immunos-mcp
```

## Script Requirements

### Python Scripts
- **Python**: 3.8+
- **Dependencies**: None (uses only standard library)

### Bash Scripts
- **Shell**: bash 4.0+
- **OS**: macOS or Linux
- **Tools**: sed, date, find

## Making Scripts Executable

```bash
chmod +x scripts/*.sh
chmod +x scripts/*.py
```

## Scheduled Automation (Optional)

### Using cron (Linux/macOS)

Add to crontab (`crontab -e`):

```cron
# Update vault daily at 11:55 PM
55 23 * * * cd /Users/byron/projects && ./scripts/update-vault.sh all >> logs/vault-update.log 2>&1

# Create tomorrow's daily note at midnight
0 0 * * * cd /Users/byron/projects && ./scripts/create-daily-note.sh >> logs/daily-note.log 2>&1
```

### Using launchd (macOS)

Create `~/Library/LaunchAgents/com.vault.update.plist`:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>Label</key>
    <string>com.vault.update</string>
    <key>ProgramArguments</key>
    <array>
        <string>/Users/byron/projects/scripts/update-vault.sh</string>
        <string>all</string>
    </array>
    <key>StartCalendarInterval</key>
    <dict>
        <key>Hour</key>
        <integer>23</integer>
        <key>Minute</key>
        <integer>55</integer>
    </dict>
</dict>
</plist>
```

Load with: `launchctl load ~/Library/LaunchAgents/com.vault.update.plist`

## Troubleshooting

### Permission Denied
```bash
chmod +x scripts/<script-name>
```

### Script Not Found
```bash
# Run from projects root
cd /Users/byron/projects
./scripts/update-vault.sh
```

### Python Import Errors
```bash
# Use python3 explicitly
python3 scripts/generate-api-docs.py immunos-mcp
```

### Date Format Errors (create-daily-note.sh)
```bash
# Use YYYY-MM-DD format
./scripts/create-daily-note.sh 2025-12-01  # Correct
./scripts/create-daily-note.sh 12/01/2025  # Wrong
```

## Adding New Scripts

1. Create script in `scripts/` directory
2. Add shebang line (`#!/bin/bash` or `#!/usr/bin/env python3`)
3. Add description and usage comments
4. Make executable: `chmod +x scripts/<script>`
5. Document in this README
6. Test thoroughly

## Links

- [[../templates/|Templates]]
- [[../daily/|Daily Notes]]
- [[../research/experiments/|Experiments]]
- [[../projects/immunos-mcp/api/|API Documentation]]

---

**Last Updated**: 2025-11-30
**Total Scripts**: 5
**Language**: Python 3, Bash

## Summary
Contains 4 subdirectories and 50 files. Key subfolders: archive/, logs/, static/, templates/.

## Directory Map
```
scripts/
├── archive/
├── logs/
├── static/
├── templates/
├── DASHBOARD_FEATURE_ROADMAP.md
├── README_SPATIALLIBD_QUERIES.md
├── analyze-genage.py
├── create-daily-note.sh
├── create-experiment-log.sh
├── download-genage.py
├── export_code_snapshots.py
├── generate-api-docs.py
├── immunos_api.py
├── immunos_assist.py
├── immunos_baseline.py
├── immunos_chat.py
├── immunos_cite_verify.py
├── immunos_config.py
├── immunos_dashboard.py
├── immunos_data_logger.py
├── immunos_download_datasets.sh
├── immunos_features_code.py
├── immunos_features_emotion.py
├── immunos_features_hallucination.py
├── immunos_features_network.py
├── immunos_features_research.py
├── immunos_github_validator.py
├── immunos_gui.py
├── immunos_handoff.py
├── immunos_inbox.py
├── immunos_journal.py
├── immunos_log_session.py
├── immunos_memory.py
├── immunos_model_manager.py
├── immunos_negsel.py
├── immunos_nk_scan.py
├── immunos_recover.py
├── immunos_review.py
├── immunos_routing.py
└── immunos_sentinel.py
└── ... (14 more)
```
