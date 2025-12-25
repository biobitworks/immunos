# IMMUNOS Installation Architecture

## Design Principle

IMMUNOS must maintain strict separation between:
1. **System files** - Installed software (read-only for users)
2. **User data** - Projects, memories, configurations (user-writable)

This enables standard software installation patterns (pip, brew, apt).

---

## Target Directory Structure

### System Installation (Read-Only)

```
$IMMUNOS_HOME/                    # e.g., /usr/local/lib/immunos
├── bin/
│   ├── immunos                   # Main CLI entry point
│   ├── immunos-dashboard         # Dashboard launcher
│   └── immunos-recover           # Recovery tool
├── lib/
│   ├── immunos_api.py
│   ├── immunos_dashboard.py
│   ├── immunos_memory.py
│   ├── immunos_snapshot.py
│   ├── immunos_recover.py
│   ├── immunos_nk_scan.py
│   └── ...
├── static/
│   ├── js/
│   │   └── monitor.js
│   └── css/
├── templates/
│   └── monitor.html
├── kb/system/                    # System KB (shipped with software)
│   ├── immunos-dashboard.md
│   ├── immunos-models.md
│   ├── immunos-publications.md
│   └── index.md
└── config/
    └── defaults.json             # Default configuration
```

### User Data Directory (User-Writable)

```
$IMMUNOS_USER_DATA/               # e.g., ~/.immunos
├── config/
│   ├── settings.json             # User settings
│   ├── orchestrator.json         # Backend config
│   └── projects.json             # Project definitions
├── memory/
│   ├── snapshots/
│   │   └── snap_*.json
│   ├── conversations/
│   │   └── mem_*.json
│   └── index.json
├── recovery/
│   └── CONTEXT_RECOVERY.md
├── kb/                           # User KB (project-specific)
│   ├── my-project-notes.md
│   └── custom-references.md
├── cache/
│   └── models.json               # Cached model list
└── db/
    └── dashboard.db              # SQLite database
```

### Project Working Directories (User-Defined)

```
/path/to/user/projects/
├── project-alpha/
│   ├── .immunos-project.json     # Project-specific config
│   ├── src/
│   └── docs/
├── project-beta/
│   ├── .immunos-project.json
│   └── ...
└── ...
```

---

## Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `IMMUNOS_HOME` | `/usr/local/lib/immunos` | System installation directory |
| `IMMUNOS_USER_DATA` | `~/.immunos` | User data directory |
| `IMMUNOS_CONFIG` | `~/.immunos/config` | User configuration |
| `IMMUNOS_KB_PATH` | System + User KB | Colon-separated KB paths |

---

## Configuration Hierarchy

Settings are loaded in order (later overrides earlier):

1. **System defaults**: `$IMMUNOS_HOME/config/defaults.json`
2. **User settings**: `$IMMUNOS_USER_DATA/config/settings.json`
3. **Project settings**: `./immunos-project.json`
4. **Environment variables**: `IMMUNOS_*`
5. **CLI arguments**: `--port`, `--config`, etc.

---

## Knowledge Base Paths

KB files are searched in order:

1. **Project KB**: `./.immunos/kb/` (current directory)
2. **User KB**: `$IMMUNOS_USER_DATA/kb/`
3. **System KB**: `$IMMUNOS_HOME/kb/system/`

### API Behavior

`/api/kb/index` returns merged list:
```json
{
  "pages": [
    {"name": "my-notes", "source": "user", "path": "~/.immunos/kb/my-notes.md"},
    {"name": "immunos-dashboard", "source": "system", "path": "/usr/local/lib/immunos/kb/system/immunos-dashboard.md"}
  ]
}
```

---

## Installation Methods

### pip (Recommended)

```bash
pip install immunos

# Creates:
# - Entry points in PATH (immunos, immunos-dashboard)
# - Package in site-packages
# - User data in ~/.immunos (on first run)
```

### Homebrew (macOS)

```bash
brew install immunos

# Creates:
# - Binaries in /usr/local/bin/
# - Lib in /usr/local/lib/immunos/
# - User data in ~/.immunos (on first run)
```

### Manual/Development

```bash
git clone https://github.com/biobitworks/immunos.git
cd immunos
pip install -e .

# Development mode - uses repo directly
# User data still in ~/.immunos
```

---

## Project Configuration

### .immunos-project.json

```json
{
  "name": "My Project",
  "version": "1.0.0",
  "working_directory": "/path/to/project",
  "models": {
    "default": "deepseek-r1:14b",
    "code": "qwen2.5-coder:7b",
    "quick": "qwen2.5:1.5b"
  },
  "kb": {
    "include": ["immunos-models", "project-notes"],
    "exclude": []
  },
  "context": {
    "auto_snapshot": true,
    "snapshot_interval_hours": 4
  }
}
```

### Per-Project KB

Projects can have local KB in `.immunos/kb/`:

```
my-project/
├── .immunos/
│   └── kb/
│       ├── project-architecture.md
│       └── api-reference.md
├── .immunos-project.json
└── src/
```

---

## Migration Path

### Current → Target

```bash
# Current (development)
/Users/byron/projects/
├── scripts/           → $IMMUNOS_HOME/lib/
├── docs/kb/           → $IMMUNOS_HOME/kb/system/
├── .immunos/          → ~/.immunos/ (already correct)
└── data/              → User project data (stays)

# Migration script (future)
immunos migrate --from /Users/byron/projects --to-system
```

---

## API Changes Required

### Current (development)

```python
base_path = Path('/Users/byron/projects')
kb_path = base_path / 'docs' / 'kb'
```

### Target (installed)

```python
import os
from pathlib import Path

IMMUNOS_HOME = Path(os.getenv('IMMUNOS_HOME', '/usr/local/lib/immunos'))
IMMUNOS_USER = Path(os.getenv('IMMUNOS_USER_DATA', Path.home() / '.immunos'))

def get_kb_paths():
    """Return KB search paths in priority order."""
    paths = []

    # Project KB (if in a project)
    project_kb = Path.cwd() / '.immunos' / 'kb'
    if project_kb.exists():
        paths.append(project_kb)

    # User KB
    user_kb = IMMUNOS_USER / 'kb'
    if user_kb.exists():
        paths.append(user_kb)

    # System KB
    system_kb = IMMUNOS_HOME / 'kb' / 'system'
    paths.append(system_kb)

    return paths
```

---

## Dashboard Considerations

### Current State

Dashboard serves files from:
- `scripts/templates/monitor.html`
- `scripts/static/js/monitor.js`
- `docs/kb/*.md`

### Target State

Dashboard serves from:
- `$IMMUNOS_HOME/templates/` (system templates)
- `$IMMUNOS_HOME/static/` (system static files)
- KB from merged paths (system + user + project)

User customization via:
- `$IMMUNOS_USER_DATA/config/dashboard.json` (theme, defaults)
- Custom KB in user/project directories

---

## Security Considerations

1. **System files**: Read-only, owned by root/admin
2. **User files**: User-owned, writable
3. **Project files**: Project-specific permissions
4. **No secrets in system KB**: API keys, paths in user config only
5. **Sandboxed execution**: Dashboard cannot write to system paths

---

## Development vs Production

| Aspect | Development | Production |
|--------|-------------|------------|
| System path | Git repo | Installed package |
| User data | `~/.immunos` | `~/.immunos` |
| KB source | `docs/kb/` | `$IMMUNOS_HOME/kb/system/` |
| Config | Local files | Layered config |
| Updates | `git pull` | `pip upgrade immunos` |

---

**Last Updated**: 2025-12-25
**Status**: Planning/Development
**Target**: v1.0.0 release
