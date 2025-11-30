---
title: VS Code LLM Manual Setup (No CLI Required)
tags: [vscode, manual-setup, guide]
created: 2025-11-30
---

# Manual VS Code LLM Setup

If you can't use the automated script, follow these steps manually.

## 1. Install Continue.dev Extension

### Via VS Code UI
1. Open VS Code
2. Click Extensions icon (or Cmd+Shift+X)
3. Search for "Continue"
4. Click "Install" on "Continue - Codestral, Claude, and more"

### Via URL
Visit: https://marketplace.visualstudio.com/items?itemName=Continue.continue

## 2. Configure Continue for Ollama

1. After installing Continue, you'll see a Continue icon in the left sidebar
2. Click it to open the Continue panel
3. Click the gear icon (⚙️) in Continue panel
4. This opens `~/.continue/config.json`

**Replace the entire file with**:

```json
{
  "models": [
    {
      "title": "Qwen 2.5 Coder (7B)",
      "provider": "ollama",
      "model": "qwen2.5-coder:7b",
      "apiBase": "http://localhost:11434",
      "contextLength": 32768,
      "completionOptions": {
        "temperature": 0.2,
        "topP": 0.95,
        "maxTokens": 2048
      }
    },
    {
      "title": "DeepSeek R1 (14B)",
      "provider": "ollama",
      "model": "deepseek-r1:14b",
      "apiBase": "http://localhost:11434",
      "contextLength": 16384,
      "completionOptions": {
        "temperature": 0.3,
        "topP": 0.9,
        "maxTokens": 4096
      }
    },
    {
      "title": "Qwen 2.5 Fast (1.5B)",
      "provider": "ollama",
      "model": "qwen2.5:1.5b",
      "apiBase": "http://localhost:11434",
      "contextLength": 8192,
      "completionOptions": {
        "temperature": 0.1,
        "topP": 0.95,
        "maxTokens": 512
      }
    }
  ],
  "tabAutocompleteModel": {
    "title": "Fast Autocomplete",
    "provider": "ollama",
    "model": "qwen2.5:1.5b",
    "apiBase": "http://localhost:11434"
  },
  "slashCommands": [
    {
      "name": "edit",
      "description": "Edit selected code"
    },
    {
      "name": "comment",
      "description": "Add comments to code"
    },
    {
      "name": "share",
      "description": "Export conversation"
    },
    {
      "name": "cmd",
      "description": "Generate shell command"
    },
    {
      "name": "commit",
      "description": "Generate commit message"
    }
  ],
  "customCommands": [
    {
      "name": "test",
      "prompt": "Write comprehensive tests for the selected code using pytest. Include edge cases and docstrings.",
      "description": "Generate pytest tests"
    },
    {
      "name": "docstring",
      "prompt": "Add detailed docstrings to the selected code following Google style guide. Include Args, Returns, Raises, and Examples.",
      "description": "Add Google-style docstrings"
    },
    {
      "name": "optimize",
      "prompt": "Analyze the selected code for performance improvements. Suggest optimizations for speed and memory usage.",
      "description": "Optimize code performance"
    },
    {
      "name": "security",
      "prompt": "Review the selected code for security vulnerabilities. Check for common issues like SQL injection, XSS, command injection, and unsafe deserialization.",
      "description": "Security code review"
    }
  ],
  "allowAnonymousTelemetry": false
}
```

3. Save the file
4. Reload VS Code window (Cmd+Shift+P → "Reload Window")

## 3. Verify Setup

1. **Check Ollama is running**:
   ```bash
   curl http://localhost:11434/api/tags
   ```

   If this fails, start Ollama:
   ```bash
   ollama serve
   # Or if installed via Homebrew:
   brew services start ollama
   ```

2. **Test Continue chat**:
   - Press `Cmd+L` in VS Code
   - You should see the Continue chat panel
   - Try asking: "What is this codebase about?"

3. **Test autocomplete**:
   - Open a Python file
   - Start typing a function
   - You should see AI suggestions (may take a few seconds first time)

## 4. Optional: Enable Inline Suggestions

Add to your VS Code settings (Cmd+, to open settings):

**File** → **Preferences** → **Settings** → Search for "settings.json"

Add these lines:

```json
{
  "editor.inlineSuggest.enabled": true,
  "editor.quickSuggestions": {
    "comments": true,
    "strings": true,
    "other": true
  }
}
```

## 5. Test Your Setup

### Test Chat (Cmd+L)
```
You: "Explain how QML-AiNet works"
AI: [Should respond with explanation]
```

### Test Inline Edit (Cmd+I)
1. Select some code
2. Press Cmd+I
3. Type: "Add type hints"
4. AI should suggest changes

### Test Custom Commands
In chat:
```
Select a function → /test
→ Should generate pytest tests
```

## 6. Troubleshooting

### Continue not connecting to Ollama?

**Check Ollama is running**:
```bash
ollama list
```

If not running:
```bash
# Start manually
ollama serve

# Or with Homebrew
brew services start ollama
```

**Check config path**:
- Continue config: `~/.continue/config.json`
- Should have `"apiBase": "http://localhost:11434"`

### Models not found?

**Install missing models**:
```bash
ollama pull qwen2.5-coder:7b
ollama pull deepseek-r1:14b
ollama pull qwen2.5:1.5b
```

### Autocomplete not working?

1. Check settings.json has `"editor.inlineSuggest.enabled": true`
2. Reload window: Cmd+Shift+P → "Reload Window"
3. Try typing again (first suggestion may be slow)

### Chat opens but no response?

1. Check Continue panel for error messages
2. Verify Ollama: `curl http://localhost:11434/api/tags`
3. Check model is downloaded: `ollama list`
4. Restart VS Code

## 7. Next Steps

Once working:

1. **Read the guides**:
   - [[VS-Code-LLM-Integration]] - Full guide
   - [[VS-Code-LLM-Quick-Reference]] - Cheatsheet

2. **Try the commands**:
   - `/test` - Generate tests
   - `/docstring` - Add docs
   - `/security` - Security scan

3. **Customize**:
   - Add your own custom commands to config
   - Adjust temperature/tokens for your style
   - Create project-specific prompts

## 8. Quick Reference

**Keyboard Shortcuts**:
- `Cmd+L` - Open chat
- `Cmd+I` - Inline edit
- `Tab` - Accept suggestion
- `Esc` - Dismiss suggestion

**Chat Commands**:
- `/test` - Generate tests
- `/docstring` - Add docstrings
- `/security` - Security scan
- `/optimize` - Performance tips
- `/commit` - Commit message

**Config Locations**:
- Continue: `~/.continue/config.json`
- VS Code: `.vscode/settings.json` (workspace)
- VS Code: `~/Library/Application Support/Code/User/settings.json` (global)

---

**Setup Time**: 5-10 minutes
**No CLI Required**: ✅
**Works Offline**: ✅ (after models downloaded)

**Support**: [[VS-Code-LLM-Integration]]
