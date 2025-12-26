# Claude Context

## Personal Information

Personal contact details are intentionally omitted from version control.
Store them in `claude.private.md` (local-only, gitignored) if needed.

---

## IMMUNOS Startup Scripts

### Context Recovery (Use at start of new sessions)

```bash
# Restore context from previous session
python3 scripts/immunos_recover.py

# View recovery summary
cat .immunos/recovery/CONTEXT_RECOVERY.md
```

### Daily Workflow

```bash
# Morning startup
python3 scripts/immunos_journal.py generate    # Generate daily journal
python3 scripts/immunos_memory.py stats        # Check memory system status
python3 scripts/immunos_todo.py list --overdue # Check overdue tasks

# View today's journal
cat .immunos/journal/$(date +%Y-%m-%d).md

# During work
python3 scripts/immunos_snapshot.py create --trigger manual --summary "Description of current work"

# Store important conversation/decision
python3 scripts/immunos_memory.py store --content "Important conversation summary" --priority high --type conversation

# End of day workflow
python3 scripts/immunos_snapshot.py create --trigger daily --summary "End of day checkpoint"
python3 scripts/immunos_journal.py generate  # Update daily journal with final summary

# End of day - Push to GitHub (if changes made)
git status
git add .
git commit -m "Daily update: [brief summary]

🤖 Generated with Claude Code

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"
git push origin main
```

### IMMUNOS Dashboard

```bash
# Start main IMMUNOS dashboard (web interface)
python3 scripts/immunos_dashboard.py --port 5000

# Start todo-only dashboard (separate interface)
python3 scripts/immunos_gui.py
```

### Conversation & Output Archiving

**Location**: All conversations stored in `.immunos/memory/conversations/`

```bash
# View all stored conversations
ls -la .immunos/memory/conversations/

# Search conversation history
python3 scripts/immunos_memory.py search "topic" --type conversation

# Store current conversation summary
python3 scripts/immunos_memory.py store \
  --content "Full conversation summary with key decisions and outcomes" \
  --priority high \
  --type conversation \
  --tags "project-name,milestone"

# Export conversation history
python3 scripts/immunos_memory.py export --type conversation --output conversations_archive.json

# View conversation details
cat .immunos/memory/conversations/mem_*.json | jq
```

**Note**: Snapshots (`.immunos/memory/snapshots/`) also contain conversation context and can be used to restore full session state.

### Git Workflow & GitHub Sync

```bash
# Check what's changed
git status
git diff

# Daily commit workflow
git add .
git commit -m "Daily update: [description]

🤖 Generated with Claude Code

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"

# Push to GitHub
git push origin main

# For specific projects
git add prion-clock/
git commit -m "Prion Clock v3.0: Add Bomba-Warczak citation

- Added 146+ long-lived proteins citation
- Updated evidence hierarchy
- All 30 citations verified (100%)

🤖 Generated with Claude Code

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"
git push origin main

# Weekly backup (pushes IMMUNOS data)
git add .immunos/
git commit -m "Weekly IMMUNOS backup: snapshots and memory

🤖 Generated with Claude Code

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"
git push origin main
```

### Manual Operations

```bash
# NK Cell Scanner (anomaly detection)
python3 scripts/immunos_nk_scan.py --scan-type full

# Token Analysis
python3 scripts/immunos_token_analyzer.py --top 10

# Memory Management
python3 scripts/immunos_memory.py store --content "Important fact" --priority high --tags "tag1,tag2"
python3 scripts/immunos_memory.py search "query"
python3 scripts/immunos_memory.py decay  # Run adaptive decay

# Citation Verification
python3 scripts/immunos_cite_verify.py --file path/to/document.md

# Baseline Scanner (file integrity)
python3 scripts/immunos_baseline.py scan
```

---

## User Preferences

### Automation
- **Level**: High - Prefers automated over manual processes
- **IMMUNOS Integration**: Required for all analysis tasks
- **Proactive Tools**: Use Task tool and IMMUNOS components without prompting when appropriate

### Working Style
- **Token Efficiency**: Critical priority - aim for <2M tokens
- **Persistence**: Requires context preservation across Claude resets
- **Documentation**: Comprehensive, with clear action items and checklists
- **Communication**: Concise, direct, technical - no unnecessary superlatives

### Technical Environment
- **Platform**: macOS (Darwin 25.1.0)
- **Working Directory**: `/Users/byron/projects`
- **Ollama Models**: `qwen2.5-coder:7b`, `deepseek-r1:14b`, `qwen2.5:1.5b`
- **Primary Projects**:
  - IMMUNOS (AI-assisted research integrity system)
  - Prion Clock (aging biology hypothesis)
  - BioViz Tech portfolio projects

---

## Active Projects

### 1. IMMUNOS Dashboard
- **Status**: Implementation complete, blog article reviewed by Opus 4.5
- **Critical Fixes Needed**: 5 security issues before publication (hardcoded secrets, CORS, missing functions)
- **Location**: `/Users/byron/projects/scripts/immunos_*.py`
- **Next Steps**: Fix blog article critical issues, then publish

### 2. Prion Clock v3.0
- **Status**: ✅ Ready for bioRxiv submission
- **Latest**: All Opus 4.5 critical fixes complete (Bomba-Warczak citation added)
- **Citations**: 30/30 verified (100%)
- **Location**: `/Users/byron/projects/prion-clock/`
- **Next Steps**: Add author attribution, final proofread, submit to bioRxiv

### 3. Papers Archive
- **Status**: Active curation of aging biology literature
- **Location**: `/Users/byron/projects/papers/`
- **Integration**: Automated figure download, Zotero integration

### 4. immunOS Replication Preprint
- **Status**: Draft v1 in progress (SciFact baseline + NegSl-AIS framing)
- **Location**: `/Users/byron/projects/immunos-preprint/`
- **Notes**: Related work pruned to core citations; PubMed-backed log maintained in `docs/reference/publications-log.md`

---

## IMMUNOS System Components

### Core Philosophy
Biological immune system metaphor for AI-assisted research:

- **NK Cell** (Scanner): First-line defense, anomaly detection
- **T Cell** (Memory): Adaptive learning, context persistence
- **B Cell** (Verifier): Citation validation, fact-checking
- **Dendritic Cell** (Reporter): Daily synthesis, journaling
- **Baseline**: File integrity monitoring
- **Token Analyzer**: Optimization and efficiency

### Multi-Model Collaboration (NEW)

**ALL AI models share the same context and conversation database**:

- **Claude Sonnet 4.5**: Complex planning, comprehensive work (via Claude Code)
- **Qwen 2.5 Coder 7B**: Code implementation (via Ollama - `ollama run qwen2.5-coder:7b`)
- **DeepSeek R1 14B**: Research analysis (via Ollama - `ollama run deepseek-r1:14b`)
- **Qwen 2.5 1.5B**: Quick tasks (via Ollama - `ollama run qwen2.5:1.5b`)

**Model Context Files**: `.immunos/model-contexts/`
- Each model has a dedicated context file with role, capabilities, and startup protocols
- Universal conversation log tracks all model activities
- Session logs preserve work across model restarts

**Log any session**:
```bash
python3 scripts/immunos_log_session.py \
  --model "model-name" \
  --summary "What you did" \
  --tags "project,topic" \
  --files "modified-files"
```

### System Locations
- **Memory Storage**: `.immunos/memory/`
  - `conversations/` - Full conversation database (JSON format)
  - `decisions/` - Key decision records
  - `snapshots/` - Full session snapshots with context
  - `index.json` - Memory index and metadata
- **Snapshots**: `.immunos/memory/snapshots/`
- **Recovery**: `.immunos/recovery/`
  - `CONTEXT_RECOVERY.md` - Latest recovery summary
  - `quick_start.sh` - Auto-generated startup script
- **Journal**: `.immunos/journal/`
  - Daily markdown files: `YYYY-MM-DD.md`
- **Database**: `.immunos/db/`
  - `dashboard.db` - SQLite database for web dashboard
- **Configuration**: `.claudeignore`, `.immunos/config/`

### Conversation Database Architecture

**All conversations are automatically stored** in `.immunos/memory/conversations/`

**Format**: JSON files with schema:
- `memory_id`: Unique identifier
- `timestamp`: When conversation occurred
- `type`: "conversation", "decision", "learning", etc.
- `priority`: high/medium/low
- `content`: Full conversation summary with achievements, files, next steps
- `tags`: Searchable tags
- `references`: Related files/documents
- `access_count`: How many times referenced
- `relevance_score`: Adaptive decay score

**Full conversation transcripts** are preserved in:
- Snapshots: `.immunos/memory/snapshots/snap_*.json` (includes recent conversation history)
- Recovery files: `.immunos/recovery/CONTEXT_RECOVERY.md` (human-readable summaries)

**Backup Strategy**: All `.immunos/` data should be committed to git regularly (see Git Workflow section)

### Key Features
- **Snapshot System**: Create checkpoints for <5s context restoration
- **Adaptive Memory Decay**: T Cell pattern with priority-based retention
- **Citation Verification**: Automated CrossRef/PubMed validation
- **Token Optimization**: 68% reduction achieved (6.6M → 2.1M tokens)

---

## Git Workflow Preferences

### Commit Guidelines
- **Only when requested**: Never create commits proactively
- **Message format**: Concise summary of "why" not "what"
- **Always include**:
  ```
  🤖 Generated with Claude Code

  Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>
  ```
- **NEVER**: Use --no-verify, --amend (unless explicitly requested), force push to main/master

### Pull Request Creation
- **Use `gh` CLI**: All GitHub operations via command line
- **PR Format**:
  ```markdown
  ## Summary
  - Bullet points of changes

  ## Test plan
  - [ ] Testing checklist

  🤖 Generated with Claude Code
  ```

---

## Quick Reference Commands

### Session Start
```bash
python3 scripts/immunos_recover.py && cat .immunos/recovery/CONTEXT_RECOVERY.md
```

### Create Snapshot
```bash
python3 scripts/immunos_snapshot.py create --trigger manual --summary "Current work description"
```

### Check System Health
```bash
python3 scripts/immunos_memory.py stats
python3 scripts/immunos_token_analyzer.py --top 10
python3 scripts/immunos_todo.py list --overdue
```

### Web Dashboards
```bash
# IMMUNOS main dashboard
open http://localhost:5000

# Todo dashboard
python3 scripts/immunos_gui.py
```

---

## Important Notes

### Token Budget
- **Target**: <2,000,000 tokens
- **Current**: ~2,100,000 tokens (with .claudeignore)
- **Reduction**: 68% from original 6.6M
- **Critical**: Always check token usage, use .claudeignore aggressively

### Context Preservation
- **Problem**: Claude resets lose all context
- **Solution**: IMMUNOS snapshot/recovery system
- **Recovery Time**: <5 seconds
- **Last Snapshot**: Check `.immunos/memory/snapshots/` for latest

### Research Integrity
- **All citations must be verified**: Use `immunos_cite_verify.py`
- **Documentation**: Every decision, assumption, and change
- **Changelog**: Maintain detailed revision history for major projects

---

## Communication Preferences

### Style
- **Concise and direct**: No fluff or excessive praise
- **Technical accuracy over validation**: Honest critique preferred
- **Professional objectivity**: Facts over emotional support
- **Action-oriented**: Clear next steps, checklists, timelines

### Documentation Format
- **Markdown**: GitHub-flavored, CommonMark spec
- **Structure**: Clear headers, bulleted lists, code blocks
- **Checklists**: Extensive use of `- [ ]` for task tracking
- **Examples**: Real code snippets, not pseudocode

### Tool Usage
- **Proactive**: Use Task tool for complex multi-step work
- **Parallel**: Run independent commands simultaneously
- **IMMUNOS-first**: Always use IMMUNOS tools for research tasks
- **No echo/printf for communication**: Output text directly

---

## File Locations Quick Reference

### IMMUNOS System
```
/Users/byron/projects/
├── .immunos/
│   ├── memory/          # T Cell memory storage
│   ├── snapshots/       # Context checkpoints
│   ├── recovery/        # Recovery summaries
│   ├── journal/         # Daily journals
│   └── db/             # SQLite databases
├── scripts/
│   ├── immunos_*.py     # All IMMUNOS tools
│   └── templates/       # Dashboard templates
```

### Projects
```
/Users/byron/projects/
├── prion-clock/         # Aging hypothesis (v3.0)
├── papers/              # Literature archive
├── immunos-mcp/         # MCP server (if applicable)
└── portfolio/           # BioViz Tech work
```

---

## Experimental Priorities

### Current Focus (as of 2025-12-23)
1. ✅ **Prion Clock v3.0**: All critical fixes complete, ready for bioRxiv
2. ⏳ **immunOS Replication Preprint**: Expand related work + finalize SciFact replication narrative
3. ⏳ **IMMUNOS Blog**: Fix 5 critical security issues before publication
4. ⏳ **Papers Curation**: Continue building aging biology literature collection

### Next Quarter Goals
- Publish Prion Clock v3.0 to bioRxiv
- Publish IMMUNOS Dashboard blog article
- Identify collaborators for Prion Clock experimental validation
- Consider grant proposals for testable predictions

---

**Last Updated**: 2025-12-23
**Context Version**: Post-Opus review, Prion Clock ready for submission
**Next Session**: Continue with IMMUNOS blog fixes or Prion Clock submission prep
