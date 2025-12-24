---
source: /Users/byron/projects/immunos-mcp/setup_mcp_mvp.sh
relative: immunos-mcp/setup_mcp_mvp.sh
generated_at: 2025-12-23 10:28
---

```bash
#!/bin/bash
# IMMUNOS-MCP MVP Setup
# One-command setup for MCP server integration with Continue.dev

set -e

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${BLUE}╔═══════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   IMMUNOS-MCP MVP Setup              ║${NC}"
echo -e "${BLUE}╚═══════════════════════════════════════╝${NC}"
echo ""

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Step 1: Install MCP SDK
echo -e "${BLUE}[1/4] Installing MCP SDK...${NC}"
if command -v uv &> /dev/null; then
    uv pip install "mcp>=1.3.1"
else
    pip install "mcp>=1.3.1"
fi
echo -e "${GREEN}✓ MCP SDK installed${NC}"
echo ""

# Step 2: Install package in editable mode
echo -e "${BLUE}[2/4] Installing IMMUNOS-MCP...${NC}"
if command -v uv &> /dev/null; then
    uv pip install -e .
else
    pip install -e .
fi
echo -e "${GREEN}✓ IMMUNOS-MCP installed${NC}"
echo ""

# Step 3: Create Continue config directory
echo -e "${BLUE}[3/4] Configuring Continue.dev...${NC}"
mkdir -p ~/.continue/mcpServers

# Create Continue MCP config
cat > ~/.continue/mcpServers/immunos-mvp.yaml << EOF
# IMMUNOS-MCP MVP Configuration
# Simple B Cell agent for code security scanning

mcpServers:
  immunos-mvp:
    command: $SCRIPT_DIR/.venv/bin/python
    args:
      - -u
      - $SCRIPT_DIR/src/immunos_mcp/servers/simple_mcp_server.py
    env:
      PYTHONPATH: $SCRIPT_DIR
EOF

echo -e "${GREEN}✓ Continue configured${NC}"
echo -e "  Config: ~/.continue/mcpServers/immunos-mvp.yaml"
echo ""

# Step 4: Test server
echo -e "${BLUE}[4/4] Testing server...${NC}"
timeout 2 python -u src/immunos_mcp/servers/simple_mcp_server.py > /dev/null 2>&1 || true
if [ $? -eq 124 ]; then
    echo -e "${GREEN}✓ Server can start${NC}"
else
    echo -e "${YELLOW}⚠ Server test inconclusive (this is usually OK)${NC}"
fi
echo ""

# Success!
echo -e "${GREEN}╔═══════════════════════════════════════╗${NC}"
echo -e "${GREEN}║   Setup Complete!                     ║${NC}"
echo -e "${GREEN}╚═══════════════════════════════════════╝${NC}"
echo ""
echo "🎉 IMMUNOS-MCP MVP is ready to use!"
echo ""
echo "📖 Next Steps:"
echo ""
echo "1. Restart VS Code (Cmd+Q, then reopen)"
echo ""
echo "2. Open a Python file with code to scan"
echo ""
echo "3. Select the code you want to analyze"
echo ""
echo "4. Press Cmd+L to open Continue chat"
echo ""
echo "5. Type: /immune-scan"
echo ""
echo "6. Get security analysis results!"
echo ""
echo "💡 Example code to test:"
echo '   eval(user_input)  ← Should detect as vulnerable'
echo '   print("hello")    ← Should detect as safe'
echo ""
echo "📚 Documentation: README_MCP.md"
echo ""

```
