# IMMUNOS Knowledge Base

## Overview

The IMMUNOS Knowledge Base contains documentation for the AI-assisted research system. Files are stored in `docs/kb/*.md` and accessible via the dashboard.

---

## Available Documentation

### System Documentation

| Document | Description |
|----------|-------------|
| [immunos-architecture](immunos-architecture.md) | Installation structure, system vs user separation |
| [immunos-packaging](immunos-packaging.md) | pip/brew/npm/apt/rpm packaging standards, OS & arch support |
| [immunos-dashboard](immunos-dashboard.md) | Dashboard buttons, panels, and features |
| [immunos-models](immunos-models.md) | Model architecture and cell type mappings |
| [immunos-publications](immunos-publications.md) | Publications, datasets, and references (IEEE format) |

### Project Templates

| Project Type | Working Directory | Recommended Models |
|--------------|-------------------|-------------------|
| **IMMUNOS Core** | `/Users/byron/projects` | Dendritic (research), B Cell (code) |
| **Prion Clock** | `/Users/byron/projects/prion-clock` | Dendritic (research) |
| **Papers Archive** | `/Users/byron/projects/papers` | NK Cell (quick), Dendritic (analysis) |

---

## Quick Reference

### Dashboard URL
```
http://localhost:5001/monitor
```

### Start Dashboard
```bash
python3 scripts/immunos_dashboard.py --port 5001
```

### Recovery Commands
```bash
# Restore context
python3 scripts/immunos_recover.py

# View recovery
cat .immunos/recovery/CONTEXT_RECOVERY.md

# Create snapshot
python3 scripts/immunos_snapshot.py create --trigger manual --summary "Description"
```

### Model Commands
```bash
# List models
ollama list

# Start Ollama
ollama serve

# Run model interactively
ollama run qwen2.5-coder:7b
```

---

## Project Configuration

Projects are saved in browser localStorage and include:

| Field | Description |
|-------|-------------|
| `name` | Project identifier |
| `workingDirectory` | Local path or server URL |
| `models` | Selected Ollama models |
| `kb` | Selected KB files |
| `uploadedFiles` | Session files (content preserved) |
| `activeModel` | Current chat model |
| `createdAt` | Creation timestamp |
| `lastUsed` | Last access timestamp |

### Example Project Structure
```json
{
  "IMMUNOS Core": {
    "workingDirectory": "/Users/byron/projects",
    "models": ["qwen2.5-coder:7b", "deepseek-r1:14b"],
    "kb": ["immunos-dashboard", "immunos-models"],
    "activeModel": "deepseek-r1:14b",
    "createdAt": "2025-12-25T14:30:00.000Z",
    "lastUsed": "2025-12-25T14:35:00.000Z"
  }
}
```

---

## Cell Type Reference

| Cell Type | Role | Model Examples | Best For |
|-----------|------|----------------|----------|
| **NK Cell** | Quick Scanner | qwen2.5:1.5b | Fast queries, anomaly detection |
| **B Cell** | Code Verifier | qwen2.5-coder:7b | Implementation, debugging |
| **Dendritic** | Research Analyst | deepseek-r1:14b | Analysis, planning, synthesis |
| **T Cell** | Memory Agent | qwen2.5:7b | Context, learning, retrieval |
| **Thymus** | Central Validator | Claude Opus 4.5 | Complex tasks, validation |

---

## API Endpoints

### Recovery & Snapshots
| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/recovery/run` | POST | Run context recovery |
| `/api/recovery/status` | GET | Get snapshot status |
| `/api/snapshot/create` | POST | Create new snapshot |

### Knowledge Base
| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/kb/index` | GET | List all KB pages |
| `/api/kb/page?name=X` | GET | Get page content |

### Models & Chat
| Endpoint | Method | Description |
|----------|--------|-------------|
| `/api/models/list` | GET | List Ollama models |
| `/api/chat` | POST | Send chat message |
| `/api/ollama/start` | POST | Start Ollama server |
| `/api/ollama/stop` | POST | Stop Ollama server |

---

## File Locations

### IMMUNOS System
```
/Users/byron/projects/
├── .immunos/
│   ├── memory/
│   │   ├── snapshots/     # Session snapshots
│   │   ├── conversations/ # Stored conversations
│   │   └── index.json     # Memory index
│   ├── recovery/
│   │   └── CONTEXT_RECOVERY.md
│   └── config/
├── scripts/
│   ├── immunos_*.py       # All IMMUNOS tools
│   ├── static/js/         # Dashboard JS
│   └── templates/         # Dashboard HTML
└── docs/
    └── kb/                # Knowledge Base files
```

---

## Adding New KB Pages

1. Create markdown file in `docs/kb/`:
   ```bash
   touch docs/kb/my-topic.md
   ```

2. Add title as H1:
   ```markdown
   # My Topic Title
   ```

3. Page will appear in dashboard KB list automatically

4. Reference in other docs:
   ```markdown
   See [my-topic](my-topic.md) for details.
   ```

---

**Last Updated**: 2025-12-25
**Version**: v20251225-kb5
