---
title: VS Code LLM Quick Reference
tags: [vscode, cheatsheet, reference]
created: 2025-11-30
---

# VS Code LLM Quick Reference

## ⚡ Setup

```bash
# One-command setup
./scripts/setup-vscode-llm.sh
```

## ⌨️ Keyboard Shortcuts

| Action | Shortcut | Description |
|--------|----------|-------------|
| **Chat** | `Cmd+L` | Open Continue chat panel |
| **Inline Edit** | `Cmd+I` | Edit selected code with AI |
| **Accept** | `Tab` | Accept autocomplete suggestion |
| **Reject** | `Esc` | Dismiss suggestion |
| **New Chat** | `Cmd+Shift+L` | Start fresh conversation |

## 💬 Custom Commands

In chat, type these slash commands:

| Command | Description | Example |
|---------|-------------|---------|
| `/test` | Generate pytest tests | Select function → `/test` |
| `/docstring` | Add Google-style docs | Select class → `/docstring` |
| `/security` | Security scan | Select code → `/security` |
| `/optimize` | Performance tips | Select function → `/optimize` |
| `/commit` | Generate commit msg | `/commit` |
| `/cmd` | Shell command | "find all .py files" |
| `/edit` | Edit selected code | Select → `/edit make async` |
| `/comment` | Add comments | Select → `/comment` |

## 🤖 Model Selection

Click model dropdown in Continue chat:

| Model | Best For | Speed | Use When |
|-------|----------|-------|----------|
| **Qwen 2.5 Coder (7B)** | Code generation, refactoring | ⚡⚡ Medium | Daily coding (default) |
| **DeepSeek R1 (14B)** | Debugging, reasoning | ⚡ Slow | Complex problems |
| **Qwen Fast (1.5B)** | Quick questions | ⚡⚡⚡ Fast | Simple queries |

## 🎯 Common Use Cases

### 1. Generate Tests
```
1. Select a function
2. Cmd+L to open chat
3. Type: /test
4. Review and accept
```

### 2. Add Documentation
```
1. Select class or function
2. Cmd+L
3. Type: /docstring
4. Get Google-style docstrings
```

### 3. Debug Code
```
1. Select problematic code
2. Cmd+L
3. Ask: "Why is this not working?"
4. Use DeepSeek R1 for complex issues
```

### 4. Refactor Code
```
1. Select code
2. Cmd+I (inline edit)
3. Type: "Refactor using dataclasses"
4. Review changes, press Cmd+Y to accept
```

### 5. Security Review
```
1. Select code
2. Cmd+L
3. Type: /security
4. Get vulnerability report
```

### 6. Quick Completions
```
1. Just start typing
2. AI suggests next code
3. Press Tab to accept
4. (Uses fast qwen2.5:1.5b model)
```

## 📝 Pro Tips

### Context Awareness
Open related files before asking questions:
```
# Continue automatically sees:
- Current file
- Open files in editor
- Recent edits
```

### Multi-file Changes
```
1. Open all files you want to change
2. In chat: "Update all API calls to use async/await"
3. Continue suggests changes across files
```

### Codebase Questions
```
"How does QML-AiNet differ from Opt-AiNet in this project?"
→ Continue searches entire codebase
```

### Temperature Control
Already optimized in config:
- Code generation: Low temp (0.2) = deterministic
- Debugging: Medium temp (0.3) = balanced
- Autocomplete: Very low (0.1) = consistent

### Custom Model Per Task
In chat:
```
1. Switch to DeepSeek R1 for debugging
2. Switch back to Qwen Coder for implementation
3. Models remember conversation context
```

## 🔧 Configuration Files

### Continue Config
`.vscode/continue-config.json`
```json
{
  "models": [...],           // Available models
  "customCommands": [...],   // Your /commands
  "tabAutocompleteModel": {} // Fast model for Tab
}
```

### VS Code Settings
`.vscode/settings.json`
```json
{
  "ollama.defaultModel": "qwen2.5-coder:7b",
  "editor.inlineSuggest.enabled": true
}
```

## 🚨 Troubleshooting

### No Suggestions?
```bash
# 1. Check Ollama is running
ollama list

# 2. Reload VS Code
Cmd+Shift+P → "Reload Window"
```

### Slow Responses?
```
1. Switch to qwen2.5:1.5b (faster model)
2. Close other applications
3. Check: top | grep ollama
```

### Models Not Loading?
```bash
# Restart Ollama
brew services restart ollama

# Check models installed
ollama list
```

### Chat Not Opening?
```
1. Check extension installed:
   Cmd+Shift+X → Search "Continue"

2. Reload window:
   Cmd+Shift+P → "Reload Window"
```

## 📊 Performance

### Autocomplete Speed
- **qwen2.5:1.5b**: ~100ms ✅
- **qwen2.5-coder:7b**: ~500ms ⚠️
- **deepseek-r1:14b**: ~2s ❌

**Default**: Fast model (qwen2.5:1.5b) for autocomplete

### Chat Response Time
- **Simple question**: 2-5 seconds
- **Code generation**: 5-10 seconds
- **Complex debugging**: 10-30 seconds

### Memory Usage
All 3 models loaded: ~15GB RAM
One model at a time: ~5-10GB RAM

## 🎓 Learning Path

### Day 1: Basic Chat
1. Install: `./scripts/setup-vscode-llm.sh`
2. Try: `Cmd+L` → Ask "What does this function do?"
3. Experiment with `/test` and `/docstring`

### Day 2: Inline Edits
1. Try: `Cmd+I` → "Add error handling"
2. Select code → `Cmd+I` → "Make this async"
3. Practice accepting/rejecting changes

### Day 3: Autocomplete
1. Start typing a function
2. Press `Tab` to accept suggestions
3. Adjust to AI's coding style

### Week 1: Custom Commands
1. Add your own commands to config
2. Create project-specific prompts
3. Build personal workflow

## 🔗 Links

- **Full Guide**: [[VS-Code-LLM-Integration]]
- **Continue Docs**: https://continue.dev/docs
- **Ollama Models**: https://ollama.com/library
- **Config**: `.vscode/continue-config.json`

## 📱 Quick Commands Summary

```bash
# Chat
Cmd+L                    # Open chat
/test                    # Generate tests
/docstring              # Add docstrings
/security               # Security scan
/optimize               # Performance tips

# Inline
Cmd+I                    # Edit selection
Cmd+Y                    # Accept changes
Cmd+N                    # Reject changes

# Autocomplete
Tab                      # Accept suggestion
Esc                      # Dismiss suggestion

# Model switching
Click dropdown in chat   # Change model
```

---

**Last Updated**: 2025-11-30
**Models**: qwen2.5-coder:7b, deepseek-r1:14b, qwen2.5:1.5b
**Total Size**: ~15GB
**Status**: ✅ Ready to use

**Setup**: `./scripts/setup-vscode-llm.sh`
