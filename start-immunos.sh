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
