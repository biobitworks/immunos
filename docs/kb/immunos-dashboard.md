# IMMUNOS Dashboard Guide

## Overview

The IMMUNOS Dashboard is a web-based interface for managing the IMMUNOS AI-assisted research system. It provides real-time monitoring, model management, and context persistence tools.

**URL**: http://localhost:5001/monitor
**Version**: v20251225-kb5

---

## Header Buttons

### Recover Context (Purple)
**Purpose**: Restore context from the most recent snapshot for session continuity.

**What it does**:
1. Runs `python3 scripts/immunos_recover.py`
2. Reads the latest snapshot from `.immunos/memory/snapshots/`
3. Generates `.immunos/recovery/CONTEXT_RECOVERY.md`
4. Displays recovery summary in chat area

**When to use**:
- Starting a new Claude session
- After a context reset
- When switching between AI models
- To verify current context state

**Output includes**:
- Snapshot name and timestamp
- Number of memories recovered
- Recent file modifications
- Active todos
- User preferences

---

### Create Snapshot (Teal)
**Purpose**: Save current session state for future recovery.

**What it does**:
1. Prompts for a summary description
2. Runs `python3 scripts/immunos_snapshot.py create --trigger manual --summary "..."`
3. Creates snapshot JSON in `.immunos/memory/snapshots/snap_YYYY-MM-DD_HHMMSS.json`
4. Updates `latest.json` symlink
5. Regenerates recovery context

**Snapshot contains**:
- Timestamp and trigger type
- Summary description
- Recent file modifications (last 24h)
- Active tasks from todo system
- Git repository state
- User preferences

**When to use**:
- Before ending a work session
- After completing a significant task
- Before switching to a different AI model
- Every 2-4 hours during active work

**Valid triggers**: `manual`, `auto_threshold`, `reset_imminent`, `task_complete`

---

## Status Bar

### Connection Status (Left)
- **Green dot**: WebSocket connected to Flask server
- **Red dot**: Disconnected - refresh page or restart server

### Ollama Status (Middle)
- **Green dot**: Ollama server running at localhost:11434
- **Red dot**: Ollama stopped - use Start/Stop buttons

### Start/Stop Ollama (Right)
- **Start Ollama**: Launches `ollama serve` in background
- **Stop Ollama**: Terminates Ollama process

---

## Models Panel

Displays available Ollama models with IMMUNOS cell type labels.

### Model Labels
| Label | Cell Type | Models | Use Case |
|-------|-----------|--------|----------|
| `[NK Cell]` | Natural Killer | qwen2.5:1.5b, phi | Quick scans, simple queries |
| `[B Cell]` | B Cell | qwen2.5-coder:7b | Code implementation, debugging |
| `[Dendritic]` | Dendritic Cell | deepseek-r1:14b | Research analysis, planning |
| `[T Cell]` | T Cell | qwen2.5:7b | Memory, context, adaptive learning |
| `[Embed]` | Embedding | bge-m3, mxbai-embed-large | Vector embeddings (not for chat) |

### Model Checkboxes
- Check/uncheck to include models in project configurations
- All models checked by default

### Active Model Dropdown
- Select which model to use for chat
- Shows cell type label and size

---

## T Cell Memory Panel

Context management system using biological immune system metaphor.

### Panel Header
- **DNA icon**: Indicates memory/context functionality
- **File count**: Shows total selected context files
- **Include in context**: Master toggle for context injection

### Project Management
| Button | Function |
|--------|----------|
| **Select Project** | Dropdown of saved projects |
| **Save** | Save current selections as new project |
| **Load** | Restore project configuration |
| **New** | Clear all selections, start fresh |

**Project saves**:
- Selected models (checkboxes)
- Selected KB files
- Active chat model
- Uploaded files (content preserved)

**Storage**: `localStorage` key `immunos_projects`

### File Upload
- **Drag and drop** or **click to upload**
- Accepted formats: `.md`, `.txt`, `.json`, `.py`, `.js`
- Files stored in browser memory (session only)
- Each file has checkbox for inclusion

### Knowledge Base Files
- Loaded from `/api/kb/index` (server-side)
- Files from `docs/kb/*.md`
- **Select All**: Check all KB files
- **Clear**: Uncheck all KB files

### How Context Works
1. When you send a message, `getSelectedKBContext()` is called
2. Selected KB files are fetched from server
3. Selected uploaded files are included from memory
4. Context is prepended to your message:
   ```
   ### Knowledge Base Context:
   --- File Title ---
   [file content]

   ### Uploaded Files:
   --- filename.md ---
   [file content]

   ### User Question:
   [your message]
   ```

### Auto-Highlight
When T Cell model (qwen2.5:7b) is selected:
- Panel gets blue ring highlight
- Header changes to "T Cell Memory (Active)"

---

## Chat Area

### Message Types
| Type | Appearance | Source |
|------|------------|--------|
| User | Blue bubble, right-aligned | Your messages |
| Assistant | Gray bubble, left-aligned | Model responses |
| System | Gray bubble, "System" label | Status messages |
| Error | Red bubble, centered | Error messages |

### Context Indicator
When context files are selected, shows:
```
[Using N context file(s): X KB + Y uploaded]
```

### Sending Messages
1. Type in input field
2. Press Enter or click Send
3. Message + context sent to `/api/chat`
4. Response streamed via WebSocket

---

## Startup Behavior

On page load, the dashboard:
1. Initializes WebSocket connection
2. Sets up event listeners
3. Checks Ollama status (direct API call)
4. Loads available models
5. Loads KB index
6. Loads saved projects from localStorage
7. **Checks recovery status** - warns if snapshot > 4 hours old

---

## API Endpoints Used

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/api/recovery/run` | POST | Run context recovery |
| `/api/recovery/status` | GET | Get snapshot age |
| `/api/snapshot/create` | POST | Create new snapshot |
| `/api/kb/index` | GET | List KB files |
| `/api/kb/page?name=X` | GET | Get KB file content |
| `/api/chat` | POST | Send chat message |
| `/api/ollama/start` | POST | Start Ollama server |
| `/api/ollama/stop` | POST | Stop Ollama server |
| `http://localhost:11434/api/tags` | GET | List Ollama models (direct) |

---

## Keyboard Shortcuts

| Key | Action |
|-----|--------|
| Enter | Send message |

---

## Troubleshooting

### Changes not appearing after refresh
1. Check version tag matches expected (e.g., `v20251225-kb5`)
2. If not, restart Flask server: `pkill -f immunos_dashboard && python3 scripts/immunos_dashboard.py --port 5001`
3. Hard refresh: Cmd+Shift+R (Mac) or Ctrl+Shift+R (Windows)

### Models not loading
1. Check Ollama status indicator
2. If red, click "Start Ollama"
3. Verify: `curl http://localhost:11434/api/tags`

### Context not being included
1. Check "Include in context" is checked
2. Verify files have checkmarks
3. Check file count shows > 0

### Snapshot creation fails
1. Check Flask server logs: `cat /tmp/dashboard.log`
2. Valid triggers: manual, auto_threshold, reset_imminent, task_complete

---

## File Locations

| Path | Purpose |
|------|---------|
| `scripts/immunos_dashboard.py` | Flask server |
| `scripts/immunos_api.py` | API endpoints |
| `scripts/templates/monitor.html` | Dashboard HTML |
| `scripts/static/js/monitor.js` | Dashboard JavaScript |
| `docs/kb/*.md` | Knowledge Base files |
| `.immunos/memory/snapshots/` | Snapshot storage |
| `.immunos/recovery/` | Recovery context |

---

**Last Updated**: 2025-12-25
**Maintained by**: IMMUNOS System
