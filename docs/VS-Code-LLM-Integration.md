---
title: VS Code Local LLM Integration Guide
tags: [vscode, ollama, llm, development]
created: 2025-11-30
---

# VS Code Local LLM Integration

Guide to using your local Ollama models (qwen2.5-coder, deepseek-r1, qwen2.5) as coding assistants in VS Code.

## 🎯 Available Models

| Model | Size | Best For | Speed | Context |
|-------|------|----------|-------|---------|
| **qwen2.5-coder:7b** | 4.7GB | Code generation, refactoring | Medium | 32K tokens |
| **deepseek-r1:14b** | 9.0GB | Complex reasoning, debugging | Slow | 16K tokens |
| **qwen2.5:1.5b** | 1.0GB | Fast completions, simple tasks | Fast | 8K tokens |

## 🚀 Option 1: Continue.dev (Recommended)

### Features
- ✅ **Tab autocomplete** (like GitHub Copilot)
- ✅ **Chat interface** (Cmd+L)
- ✅ **Inline edits** (Cmd+I)
- ✅ **Context awareness** (reads open files)
- ✅ **Custom slash commands**
- ✅ **Codebase indexing**

### Installation

```bash
# 1. Install VS Code extension
code --install-extension Continue.continue

# 2. Configuration is already set up at:
#    .vscode/continue-config.json
```

### Usage

#### Chat (Cmd+L)
```
You: "How does QML-AiNet's mutation operator work?"
Qwen 2.5 Coder: [Explains from your codebase]
```

#### Inline Edit (Cmd+I)
1. Select code
2. Press Cmd+I
3. Type instruction: "Add type hints and docstrings"
4. Accept changes

#### Tab Autocomplete
- Just start typing
- Press Tab to accept suggestion
- Uses fast qwen2.5:1.5b model

#### Custom Commands

We've configured these commands:

```bash
/test        # Generate pytest tests
/docstring   # Add Google-style docstrings
/optimize    # Performance improvements
/security    # Security vulnerability scan
/commit      # Generate commit message
```

**Example**:
1. Select a function
2. Type `/test` in chat
3. Get comprehensive pytest tests

### Model Selection

Continue lets you switch models mid-conversation:

1. Open chat (Cmd+L)
2. Click model dropdown
3. Choose:
   - **Qwen 2.5 Coder** - Default, best for code
   - **DeepSeek R1** - Complex debugging
   - **Qwen Fast** - Quick questions

### Configuration Location

Your configuration is in `.vscode/continue-config.json`:

```json
{
  "models": [
    {
      "title": "Qwen 2.5 Coder (7B)",
      "provider": "ollama",
      "model": "qwen2.5-coder:7b",
      "apiBase": "http://localhost:11434"
    }
    // ... more models
  ],
  "customCommands": [
    {
      "name": "test",
      "prompt": "Write comprehensive tests..."
    }
    // ... more commands
  ]
}
```

## 🔧 Option 2: Ollama Extension

Simpler, chat-only alternative.

### Installation

```bash
code --install-extension ollama.ollama-vscode
```

### Usage

1. Open Command Palette (Cmd+Shift+P)
2. Type "Ollama: Chat"
3. Select model
4. Ask questions

**Pros**: Simple, lightweight
**Cons**: No autocomplete, no inline edits

## 💬 Option 3: MCP Server Integration

Use IMMUNOS-MCP as a VS Code tool via Model Context Protocol.

### Setup (Future)

```bash
# In your MCP client config
{
  "mcpServers": {
    "immunos": {
      "command": "python",
      "args": ["-m", "immunos_mcp.server"],
      "env": {
        "OLLAMA_HOST": "http://localhost:11434"
      }
    }
  }
}
```

**Capabilities**:
- Code security analysis via immune agents
- Pattern matching with B Cell
- Anomaly detection with NK Cell
- Multi-agent threat assessment

## 🎨 Use Cases

### 1. Code Generation
```
Prompt: "Create a new immune agent for feature extraction based on the existing B Cell agent"
Model: Qwen 2.5 Coder
```

### 2. Debugging
```
Prompt: "Why is network suppression not working in my Opt-AiNet implementation?"
Model: DeepSeek R1 (better reasoning)
```

### 3. Documentation
```
Command: Select function → /docstring
Model: Qwen 2.5 Coder
Result: Google-style docstrings added
```

### 4. Testing
```
Command: Select class → /test
Model: Qwen 2.5 Coder
Result: Comprehensive pytest suite
```

### 5. Security Review
```
Command: Select code → /security
Model: DeepSeek R1
Result: Vulnerability analysis
```

### 6. Refactoring
```
Inline Edit: "Extract this into smaller functions with proper error handling"
Model: Qwen 2.5 Coder
```

## 📊 Performance Comparison

| Task | qwen2.5-coder:7b | deepseek-r1:14b | qwen2.5:1.5b |
|------|------------------|-----------------|--------------|
| **Autocomplete** | ⚠️ Slow | ❌ Too slow | ✅ Perfect |
| **Code generation** | ✅ Excellent | ✅ Excellent | ⚠️ Basic |
| **Debugging** | ✅ Good | ✅✅ Best | ❌ Weak |
| **Documentation** | ✅ Excellent | ✅ Good | ⚠️ Basic |
| **Security analysis** | ✅ Good | ✅✅ Best | ❌ Weak |
| **Simple questions** | ✅ Fast | ⚠️ Slow | ✅✅ Fastest |

**Recommended Strategy**:
- **Autocomplete**: qwen2.5:1.5b (configured in Continue)
- **Daily coding**: qwen2.5-coder:7b
- **Complex problems**: deepseek-r1:14b

## 🔥 Pro Tips

### 1. Context Awareness

Continue automatically includes:
- Current file
- Related imports
- Recently edited files

**Tip**: Open relevant files before asking questions

### 2. Codebase Indexing

Continue can index your entire codebase:

```json
{
  "docs": [
    {
      "title": "IMMUNOS-MCP Codebase",
      "startUrl": "file:///Users/byron/projects/immunos-mcp"
    }
  ]
}
```

Then ask: "How is QML-AiNet different from Opt-AiNet in this codebase?"

### 3. Temperature Settings

Our config optimizes for code:
- **qwen2.5-coder**: temp=0.2 (deterministic)
- **deepseek-r1**: temp=0.3 (slightly creative)
- **qwen2.5**: temp=0.1 (very deterministic)

Lower = more consistent, higher = more creative

### 4. Custom Prompts

Add your own commands to `.vscode/continue-config.json`:

```json
{
  "customCommands": [
    {
      "name": "ainet",
      "prompt": "Explain this code in terms of artificial immune system concepts (antibodies, affinity, clonal selection, etc.)",
      "description": "Explain using AIS terminology"
    }
  ]
}
```

### 5. Keyboard Shortcuts

| Action | Shortcut | Description |
|--------|----------|-------------|
| Open chat | Cmd+L | Start conversation |
| Inline edit | Cmd+I | Edit selected code |
| Accept suggestion | Tab | Accept autocomplete |
| Reject suggestion | Esc | Dismiss autocomplete |
| New chat | Cmd+Shift+L | Fresh conversation |

## 🚨 Troubleshooting

### Models Not Showing Up

```bash
# 1. Check Ollama is running
ollama list

# 2. Restart Ollama
brew services restart ollama

# 3. Reload VS Code window
Cmd+Shift+P → "Reload Window"
```

### Autocomplete Not Working

```json
// In .vscode/settings.json
{
  "editor.inlineSuggest.enabled": true,
  "editor.quickSuggestions": {
    "other": true
  }
}
```

### Slow Responses

- Switch to faster model (qwen2.5:1.5b)
- Close other applications
- Check CPU usage: `top | grep ollama`

### Out of Memory

```bash
# Check memory usage
ollama ps

# If multiple models loaded, restart Ollama
brew services restart ollama

# Reduce context length in continue-config.json
```

## 📚 Resources

- **Continue.dev**: https://continue.dev/docs
- **Ollama Models**: https://ollama.com/library
- **MCP Protocol**: https://modelcontextprotocol.io
- **Qwen Documentation**: https://github.com/QwenLM/Qwen2.5-Coder

## 🔗 Related

- [[../projects/immunos-mcp/docs/Model-Selection-By-Agent-Role|Model Selection Guide]]
- [[../projects/immunos-mcp/docs/Offline-Deployment-Architecture|Offline Deployment]]
- [[../scripts/README|Automation Scripts]]

## 🎯 Next Steps

1. **Install Continue**: `code --install-extension Continue.continue`
2. **Try chat**: Cmd+L → Ask about your code
3. **Try autocomplete**: Start typing, press Tab
4. **Customize**: Edit `.vscode/continue-config.json`
5. **Explore**: Try custom commands (/test, /docstring, etc.)

---

**Created**: 2025-11-30
**Models Available**: 3 (qwen2.5-coder:7b, deepseek-r1:14b, qwen2.5:1.5b)
**Total Size**: ~15GB
**Status**: ✅ Ready to use
