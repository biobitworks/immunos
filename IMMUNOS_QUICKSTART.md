# IMMUNOS Quick Start Guide
**Artificial Immune System Orchestrator for Local + Airgapped AI**

Last Updated: December 22, 2025

---

## 🎯 What is IMMUNOS?

IMMUNOS is a **multi-agent artificial immune system** that provides:
- **Scientific Validation** - Verify citations, catch hallucinations
- **Security Layer** - Prevent errors from propagating to research
- **Offline-first Routing** - Use local models when airgapped
- **Token Savings** - Lower costs via intelligent routing
- **RAIT Compliance** - Rigorous, Accurate, Interpretable, Transparent

**Naming note**: IMMUNOS is the product name (not an operating system). The core Python package still uses
`immunos_mcp` for compatibility.

---

## 🚀 Quick Start

### 1. Chat with IMMUNOS (Terminal)
```bash
python3 scripts/immunos_chat.py
```

Commands in chat:
- `/help` - Show commands
- `/status` - Token usage
- `/handoff` - Save context
- `/clear` - Clear history
- `exit` - Quit

### 2. Use IMMUNOS for Tasks (FREE via Ollama)
```bash
# Code review
python3 scripts/immunos_assist.py review myfile.py

# Generate code
python3 scripts/immunos_assist.py code "sort a dictionary by value"

# Reasoning
python3 scripts/immunos_assist.py reason "explain P vs NP"

# Analyze file
python3 scripts/immunos_assist.py analyze myfile.txt "summarize this"
```

### 3. Dashboard & Documentation
```bash
# Real-time monitoring
http://localhost:5001/monitor

# Full documentation
http://localhost:5001/docs

# Check handoff status
python3 scripts/immunos_handoff.py check
```

---

## 💡 How It Works

### Automatic Model Routing

**Offline (Ollama)**:
- `qwen2.5-coder:7b` → Code tasks, reviews, generation
- `deepseek-r1:14b` → Reasoning, planning, problem solving
- `qwen2.5:1.5b` → Simple queries, file analysis

**Online (Configurable)**:
- Claude Code
- ChatGPT
- OpenRouter
- Local server endpoints

### Context Handoff System

When online models approach limits (150k warning, 180k critical):

1. **Warning Phase** (>150k tokens)
   - Routine tasks automatically route to local models
   - Dashboard shows warning

2. **Critical Phase** (>180k tokens)
   - All tasks route to local models
   - Context saved automatically
   - Resume when online models are available

**Manual Handoff**:
```bash
python3 scripts/immunos_handoff.py save --task "Current work" --steps "Next step 1" "Next step 2"
```

**Claude Code + ChatGPT handoff prompts**:
See `docs/kb/handoff.md` for prompts that include file locations and safe context.

---

## 📊 Token Savings

### Current Session Example:
```
Claude:  1,500 tokens  ($0.01)
Ollama:  1,477 tokens  ($0.00 - FREE)
SAVED:   $0.01        (56% reduction)
```

### Projected Savings (with full use):
- **Before IMMUNOS**: ~$1.50/day (200k Claude tokens)
- **With IMMUNOS**: ~$0.45/day (70% routed to Ollama)
- **Savings**: ~$1.05/day (**70% reduction**)

View in dashboard: http://localhost:5001/monitor → "Tokens & Savings" tab

---

## 🔬 Scientific Credibility (RAIT Principles)

### Rigorous
- 100% citation verification
- Multi-detector consensus
- Biological immune system algorithms

### Accurate
- Negative selection (trained on "self", detects "non-self")
- Cross-validation across domains
- Primary source verification

### Interpretable
- Every detection includes confidence score
- Reasoning provided for all decisions
- Human-readable outputs

### Transparent
- All routing decisions logged
- Open source code
- Auditable validation trails

---

## 🧬 Multi-Modal Analysis

### Current Support:
- **Text**: Hallucination detection, citation verification
- **Code**: Vulnerability scanning, style review
- **Data**: Anomaly detection, pattern recognition

### Future Support:
- **Images**: Figure validation, manipulation detection
- **Video**: Temporal analysis, protocol verification

---

## 🎓 Use Cases

### Research Paper Writing
```bash
# Verify citations
python3 scripts/immunos_chat.py
> Verify that Li et al 2025 found 146+ long-lived proteins in oocytes
```

### Code Development
```bash
# Review for security
python3 scripts/immunos_assist.py review authentication.py
```

### Data Analysis
```bash
# Reason about statistics
python3 scripts/immunos_assist.py reason "Is my p-value significant at 0.05?"
```

---

## 🎛️ Dashboard Features

Navigate to: http://localhost:5001/monitor

**Tabs**:
1. **Chat** - Orchestrator chat (primary interface)
2. **Overview** - System status, recent events, issues, Spleen Summary
3. **KB** - Knowledge base docs in-app
4. **Thymus** - Training intake queue and status
5. **Domains** - Per-domain monitoring
6. **Orchestrator Map** - Agent network visualization
7. **System** - Models + tokens + savings

**Chat Interface**:
- Live chat with IMMUNOS in browser
- Orchestrator/router modes
- Online/offline backend controls

---

## 🔧 Advanced Usage

### Check Token Status
```bash
python3 scripts/immunos_handoff.py check
```

Output:
```
Session ID: xxx-xxx
Token Usage: 122,000 / 200,000 (61.0%)
Should Handoff: NO
Reason: normal
```

### List Available Handoffs
```bash
python3 scripts/immunos_handoff.py list
```

### Load Saved Handoff
```bash
python3 scripts/immunos_handoff.py load --handoff-id 20251219_123000
```

---

## 📚 Documentation

**Full Documentation**: http://localhost:5001/docs

**GitHub**: https://github.com/biobitworks/immunos

**Topics Covered**:
- RAIT principles in detail
- Negative selection model
- Multi-modal analysis capabilities
- Terminal usage examples
- Dashboard guide
- Scientific credibility

---

## 💬 Chat Commands Reference

While in `python3 scripts/immunos_chat.py`:

| Command | Description |
|---------|-------------|
| `/help` | Show all commands |
| `/status` | Token usage & session info |
| `/handoff` | Create context handoff |
| `/clear` | Clear conversation history |
| `exit` or `quit` | Exit chat |

---

## 🔄 Workflow Example

### Typical Research Session:

1. **Start Chat**
   ```bash
   python3 scripts/immunos_chat.py
   ```

2. **Ask Questions** (Routes automatically to best model)
   ```
   You: Summarize the UBE2V1 study by Li et al 2025
   IMMUNOS: [Uses qwen2.5:1.5b for simple summary]
   ```

3. **Complex Reasoning** (Routes to DeepSeek)
   ```
   You: Why would maternal protein aggregation affect offspring lifespan?
   IMMUNOS: [Uses deepseek-r1:14b for reasoning]
   ```

4. **Check Savings**
   - Open http://localhost:5001/monitor
   - See token savings in real-time

5. **Context Handoff** (If approaching limits)
   ```
   /handoff
   ```

6. **Resume Later** (When Claude returns)
   ```bash
   python3 scripts/immunos_handoff.py load
   ```

---

## ⚡ Performance

### Model Response Times (Approximate):
- `qwen2.5:1.5b` - ~500ms (lightweight)
- `qwen2.5-coder:7b` - ~1-2s (code tasks)
- `deepseek-r1:14b` - ~3-5s (deep reasoning)
- `claude-sonnet-4.5` - ~2-4s (complex tasks)

### Token Processing:
- Ollama: **FREE**, unlimited local usage
- Claude: Costs apply, but 70-80% reduction via routing

---

## 🎯 Next Steps

1. ✅ **Try the chat**: `python3 scripts/immunos_chat.py`
2. ✅ **Review the docs**: http://localhost:5001/docs
3. ✅ **Monitor savings**: http://localhost:5001/monitor
4. 🔄 **Train detectors**: (Phase 3 - Coming next)

---

## 📧 Support

**Issues**: https://github.com/biobitworks/immunos/issues
**Email**: byron@biobitworks.com
**ORCID**: 0000-0002-4925-4795

---

**Remember**: IMMUNOS is your **negative selection layer** to Claude. It validates outputs, catches hallucinations, and maintains scientific rigor (RAIT) while saving you money through intelligent routing!
