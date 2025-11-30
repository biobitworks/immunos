#!/bin/bash
# Setup VS Code LLM Integration
#
# Installs Continue.dev extension and configures local Ollama models
#
# Usage:
#   ./scripts/setup-vscode-llm.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECTS_ROOT="$(dirname "$SCRIPT_DIR")"

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${BLUE}╔═══════════════════════════════════════╗${NC}"
echo -e "${BLUE}║  VS Code LLM Integration Setup       ║${NC}"
echo -e "${BLUE}╚═══════════════════════════════════════╝${NC}"
echo ""

# Check if code command is available
if ! command -v code &> /dev/null; then
    echo -e "${RED}✗ VS Code 'code' command not found${NC}"
    echo ""
    echo "Please install VS Code command line tools:"
    echo "  1. Open VS Code"
    echo "  2. Cmd+Shift+P"
    echo "  3. Type: 'Shell Command: Install code command in PATH'"
    echo ""
    exit 1
fi

echo -e "${GREEN}✓ VS Code found${NC}"
echo ""

# Check if Ollama is running
echo -e "${BLUE}► Checking Ollama...${NC}"
if ! curl -s http://localhost:11434/api/tags > /dev/null 2>&1; then
    echo -e "${RED}✗ Ollama is not running${NC}"
    echo ""
    echo "Starting Ollama..."
    if command -v ollama &> /dev/null; then
        ollama serve &
        sleep 3
        echo -e "${GREEN}✓ Ollama started${NC}"
    else
        echo -e "${RED}✗ Ollama not installed${NC}"
        echo "Run: brew install ollama"
        exit 1
    fi
else
    echo -e "${GREEN}✓ Ollama is running${NC}"
fi
echo ""

# Check which models are installed
echo -e "${BLUE}► Checking installed models...${NC}"
INSTALLED_MODELS=$(ollama list | tail -n +2 | awk '{print $1}')

REQUIRED_MODELS=(
    "qwen2.5-coder:7b"
    "deepseek-r1:14b"
    "qwen2.5:1.5b"
)

MISSING_MODELS=()
for model in "${REQUIRED_MODELS[@]}"; do
    if echo "$INSTALLED_MODELS" | grep -q "^$model"; then
        echo -e "${GREEN}  ✓ $model${NC}"
    else
        echo -e "${YELLOW}  ⚠ $model (not installed)${NC}"
        MISSING_MODELS+=("$model")
    fi
done
echo ""

if [ ${#MISSING_MODELS[@]} -gt 0 ]; then
    echo -e "${YELLOW}Some models are missing. Install them?${NC}"
    echo "This may take a while and require ~15GB of space."
    echo ""
    echo "Missing models:"
    for model in "${MISSING_MODELS[@]}"; do
        echo "  - $model"
    done
    echo ""
    read -p "Install missing models? (y/n) " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        for model in "${MISSING_MODELS[@]}"; do
            echo -e "${BLUE}Downloading $model...${NC}"
            ollama pull "$model"
        done
        echo -e "${GREEN}✓ All models installed${NC}"
    else
        echo -e "${YELLOW}⚠ Skipping model installation${NC}"
        echo "  You can install later with: ollama pull <model-name>"
    fi
    echo ""
fi

# Install Continue.dev extension
echo -e "${BLUE}► Installing Continue.dev extension...${NC}"
if code --list-extensions | grep -q "Continue.continue"; then
    echo -e "${GREEN}✓ Continue.dev already installed${NC}"
else
    code --install-extension Continue.continue
    echo -e "${GREEN}✓ Continue.dev installed${NC}"
fi
echo ""

# Optional: Install Ollama extension
echo -e "${BLUE}► Install Ollama extension? (simpler, chat-only)${NC}"
read -p "Install Ollama extension? (y/n) " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    if code --list-extensions | grep -q "ollama.ollama-vscode"; then
        echo -e "${GREEN}✓ Ollama extension already installed${NC}"
    else
        code --install-extension ollama.ollama-vscode
        echo -e "${GREEN}✓ Ollama extension installed${NC}"
    fi
fi
echo ""

# Check configuration files
echo -e "${BLUE}► Checking configuration files...${NC}"

VSCODE_DIR="$PROJECTS_ROOT/.vscode"
if [ ! -d "$VSCODE_DIR" ]; then
    mkdir -p "$VSCODE_DIR"
    echo -e "${GREEN}✓ Created .vscode directory${NC}"
fi

if [ -f "$VSCODE_DIR/continue-config.json" ]; then
    echo -e "${GREEN}✓ Continue configuration found${NC}"
else
    echo -e "${YELLOW}⚠ Continue configuration not found${NC}"
    echo "  Expected at: $VSCODE_DIR/continue-config.json"
fi

if [ -f "$VSCODE_DIR/settings.json" ]; then
    echo -e "${GREEN}✓ VS Code settings found${NC}"
else
    echo -e "${YELLOW}⚠ VS Code settings not found${NC}"
    echo "  Expected at: $VSCODE_DIR/settings.json"
fi
echo ""

# Optional: Install embedding model for codebase search
echo -e "${BLUE}► Install embedding model for codebase search?${NC}"
echo "This enables semantic search across your codebase in Continue."
read -p "Install nomic-embed-text? (~300MB) (y/n) " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    if ollama list | grep -q "nomic-embed-text"; then
        echo -e "${GREEN}✓ nomic-embed-text already installed${NC}"
    else
        ollama pull nomic-embed-text
        echo -e "${GREEN}✓ nomic-embed-text installed${NC}"
    fi
fi
echo ""

# Summary
echo -e "${BLUE}╔═══════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   Setup Complete!                     ║${NC}"
echo -e "${BLUE}╚═══════════════════════════════════════╝${NC}"
echo ""
echo "✅ VS Code LLM integration is ready!"
echo ""
echo "📖 Quick Start:"
echo ""
echo "  1. Restart VS Code (Cmd+Q, reopen)"
echo "  2. Open Command Palette (Cmd+Shift+P)"
echo "  3. Type 'Continue' to see available commands"
echo ""
echo "⌨️  Keyboard Shortcuts:"
echo ""
echo "  Cmd+L       - Open Continue chat"
echo "  Cmd+I       - Inline code edit"
echo "  Tab         - Accept autocomplete"
echo ""
echo "💡 Try These Commands:"
echo ""
echo "  /test       - Generate pytest tests"
echo "  /docstring  - Add Google-style docstrings"
echo "  /security   - Security vulnerability scan"
echo "  /optimize   - Performance improvements"
echo ""
echo "📚 Documentation:"
echo "  $PROJECTS_ROOT/docs/VS-Code-LLM-Integration.md"
echo ""
echo "🎯 Configured Models:"
echo "  • qwen2.5-coder:7b  - Code generation (default)"
echo "  • deepseek-r1:14b   - Complex reasoning"
echo "  • qwen2.5:1.5b      - Fast autocomplete"
echo ""
echo "Enjoy coding with local AI! 🚀"
echo ""
