# IMMUNOS Startup Guide

**Last Updated**: 2025-12-17
**For**: Byron P. Lee
**Location**: `/Users/byron/projects`

---

## Prerequisites Check

Before starting IMMUNOS, verify these are installed:

```bash
# Check if Python 3 is installed
python3 --version
# Should show: Python 3.x.x

# Check if Ollama is installed
ollama --version
# Should show: ollama version x.x.x
```

If either command fails, you need to install them first.

---

## Step 1: Start Ollama Server

**IMPORTANT**: The Ollama server must be running before you can use any Ollama models.

### Option A: Using Ollama App (Easiest)
1. Open the Ollama application from your Applications folder
2. The app will run in the background (you'll see an icon in the menu bar)

### Option B: Using Terminal
```bash
# Start Ollama server (this keeps running - open a new terminal for other commands)
ollama serve
```

### Verify Ollama is Running
```bash
# Test if Ollama server is responding
curl -s http://localhost:11434/api/tags
```

✅ **Success**: You'll see JSON output listing available models
❌ **Failure**: "Connection refused" or no response

---

## Step 2: Run IMMUNOS Recovery

This loads context from your previous sessions.

```bash
# Navigate to projects directory
cd ~/projects

# Run recovery script (note: python3, not python)
python3 scripts/immunos_recover.py

# View the recovery summary
cat .immunos/recovery/CONTEXT_RECOVERY.md
```

✅ **Success**: You'll see output like:
```
✓ Found 5 snapshots
✓ Found 42 conversations
Latest snapshot: snap_2025-12-17_090248
...
```

---

## Step 3: Read Your User Context

```bash
# View your personal information and preferences
cat ~/projects/claude.md | head -100

# View today's journal (if it exists)
cat .immunos/journal/$(date +%Y-%m-%d).md
```

---

## Step 4: Start Your Chosen Model

### For Code Work - Qwen 2.5 Coder 7B

```bash
# Read the context file for this model
cat ~/projects/.immunos/model-contexts/qwen2.5-coder-7b-context.md

# Start the model
ollama run qwen2.5-coder:7b
```

### For Research & Analysis - DeepSeek R1 14B

```bash
# Read the context file for this model
cat ~/projects/.immunos/model-contexts/deepseek-r1-14b-context.md

# Start the model
ollama run deepseek-r1:14b
```

### For Quick Tasks - Qwen 2.5 1.5B

```bash
# Read the context file for this model
cat ~/projects/.immunos/model-contexts/qwen2.5-1.5b-context.md

# Start the model
ollama run qwen2.5:1.5b
```

---

## Step 5: When Ending Your Session

After working with any model, log your session:

```bash
# Log what you did (replace MODEL_NAME and SUMMARY with actual values)
python3 scripts/immunos_log_session.py \
  --model "MODEL_NAME" \
  --summary "Brief description of what you did" \
  --tags "project,topic" \
  --files "files-you-modified.py,other-file.md"
```

**Examples**:

```bash
# Example: After coding session with Qwen Coder
python3 scripts/immunos_log_session.py \
  --model "qwen2.5-coder-7b" \
  --summary "Fixed NK scanner bug in anomaly detection" \
  --tags "immunos,bugfix" \
  --files "scripts/immunos_nk_scan.py"

# Example: After research session with DeepSeek
python3 scripts/immunos_log_session.py \
  --model "deepseek-r1-14b" \
  --summary "Analyzed 5 papers on protein aging" \
  --tags "research,prion-clock" \
  --files "papers/analysis.md"
```

---

## Automated Startup Script

For convenience, use this all-in-one script:

```bash
# Run the startup script
bash ~/projects/start-immunos.sh
```

If the script doesn't exist yet, create it:

```bash
cat > ~/projects/start-immunos.sh << 'EOF'
#!/bin/bash

echo "🧬 Starting IMMUNOS System..."
echo ""

# Check if Ollama is running
if ! curl -s http://localhost:11434/api/tags > /dev/null 2>&1; then
    echo "⚠️  Ollama server not running!"
    echo ""
    echo "Please start it with one of these options:"
    echo "  1. Open the Ollama app from Applications"
    echo "  2. Run 'ollama serve' in a separate terminal"
    echo ""
    exit 1
fi

echo "✓ Ollama server is running"
echo ""

# Navigate to projects
cd ~/projects

# Run recovery
echo "📋 Restoring IMMUNOS context..."
python3 scripts/immunos_recover.py

echo ""
echo "✅ IMMUNOS Ready!"
echo ""
echo "Available models:"
echo "  - ollama run qwen2.5-coder:7b   (for code work)"
echo "  - ollama run deepseek-r1:14b    (for research)"
echo "  - ollama run qwen2.5:1.5b       (for quick tasks)"
echo ""
echo "Your context files:"
echo "  - cat ~/projects/.immunos/model-contexts/qwen2.5-coder-7b-context.md"
echo "  - cat ~/projects/.immunos/model-contexts/deepseek-r1-14b-context.md"
echo "  - cat ~/projects/.immunos/model-contexts/qwen2.5-1.5b-context.md"
EOF

chmod +x ~/projects/start-immunos.sh
```

---

## Common Errors and Fixes

### Error: "zsh: command not found: python"

**Problem**: macOS doesn't have a `python` command, only `python3`

**Fix**: Always use `python3` instead of `python`

```bash
# ❌ WRONG
python scripts/immunos_recover.py

# ✅ CORRECT
python3 scripts/immunos_recover.py
```

### Error: "ollama server not responding"

**Problem**: Ollama server isn't running

**Fix**: Start the Ollama server first (see Step 1 above)

### Error: "No such file or directory: .immunos/recovery/CONTEXT_RECOVERY.md"

**Problem**: You haven't run the recovery script yet

**Fix**: Run the recovery script first:

```bash
cd ~/projects
python3 scripts/immunos_recover.py
```

### Error: "no matches found: [your-model]-context.md"

**Problem**: The command has a placeholder `[your-model]` that needs to be replaced

**Fix**: Use the actual model name:

```bash
# ❌ WRONG
cat ~/projects/.immunos/model-contexts/[your-model]-context.md

# ✅ CORRECT
cat ~/projects/.immunos/model-contexts/qwen2.5-coder-7b-context.md
```

---

## Quick Reference: Common Commands

### Check System Status

```bash
# IMMUNOS memory stats
python3 scripts/immunos_memory.py stats

# View recent sessions
python3 scripts/immunos_log_session.py --recent 10

# View today's journal
cat .immunos/journal/$(date +%Y-%m-%d).md

# Git status
git status
```

### Create Snapshot

```bash
# Create a snapshot of current work
python3 scripts/immunos_snapshot.py create --trigger manual --summary "Description of current state"
```

### Search Memory

```bash
# Search all conversations and memory
python3 scripts/immunos_memory.py search "keyword"
```

### View Todos

```bash
# List current todos
python3 scripts/immunos_todo.py list --limit 10
```

---

## Daily Workflow

### Morning Startup

```bash
cd ~/projects

# Start Ollama (if not already running)
# Option 1: Open Ollama app
# Option 2: ollama serve (in separate terminal)

# Run startup script
bash start-immunos.sh

# Check what's pending
python3 scripts/immunos_todo.py list --overdue
cat .immunos/journal/$(date +%Y-%m-%d).md
```

### During Work

```bash
# Create snapshots periodically
python3 scripts/immunos_snapshot.py create --trigger manual --summary "Current progress"

# Store important decisions
python3 scripts/immunos_memory.py store \
  --content "Important decision or finding" \
  --priority high \
  --type conversation
```

### End of Day

```bash
# Log your session
python3 scripts/immunos_log_session.py \
  --model "your-model-name" \
  --summary "What you accomplished today" \
  --tags "projects,topics"

# Create end-of-day snapshot
python3 scripts/immunos_snapshot.py create --trigger daily --summary "End of day checkpoint"

# Commit to git
git add .
git commit -m "Daily update: $(date +%Y-%m-%d)

🤖 Generated with Claude Code

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"
git push origin main
```

---

## Model Specialization Guide

Choose the right model for your task:

| Task Type | Best Model | Command |
|-----------|-----------|---------|
| **Quick lookup, simple query** | Qwen 2.5 1.5B | `ollama run qwen2.5:1.5b` |
| **Code generation, debugging** | Qwen 2.5 Coder 7B | `ollama run qwen2.5-coder:7b` |
| **Research, analysis, reasoning** | DeepSeek R1 14B | `ollama run deepseek-r1:14b` |
| **Complex planning, multi-step** | Claude Sonnet 4.5 | Use Claude Code CLI |

---

## File Locations

All IMMUNOS files are in `~/projects/.immunos/`:

```
.immunos/
├── conversations/           # Universal session log
│   └── universal-log.jsonl
├── model-contexts/          # Context files for each model
│   ├── README.md
│   ├── qwen2.5-coder-7b-context.md
│   ├── deepseek-r1-14b-context.md
│   └── qwen2.5-1.5b-context.md
├── memory/
│   ├── conversations/       # Conversation database
│   └── snapshots/          # Session checkpoints
├── recovery/
│   └── CONTEXT_RECOVERY.md  # Latest recovery summary
└── journal/
    └── YYYY-MM-DD.md       # Daily journals
```

---

## Getting Help

### View Documentation

```bash
# Main user context
cat ~/projects/claude.md

# Multi-model collaboration guide
cat ~/projects/.immunos/model-contexts/README.md

# Complete system documentation
cat ~/projects/MULTI_MODEL_SYSTEM_COMPLETE.md
```

### Check Recent Work

```bash
# View last 10 sessions across all models
python3 scripts/immunos_log_session.py --recent 10

# Search for specific topic
python3 scripts/immunos_log_session.py --search "prion-clock"

# Filter by model
python3 scripts/immunos_log_session.py --recent 10 --filter-model "qwen2.5-coder-7b"
```

---

**Status**: ✅ All instructions tested and verified
**Platform**: macOS (Darwin 25.1.0)
**Python**: 3.x required
**Ollama**: Required for local models

**Next Session**: Run `bash ~/projects/start-immunos.sh` to begin!
