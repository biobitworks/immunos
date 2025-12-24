---
source: /Users/byron/projects/scripts/create-experiment-log.sh
relative: scripts/create-experiment-log.sh
generated_at: 2025-12-23 10:28
---

```bash
#!/bin/bash
# Create Experiment Log
#
# Creates a new experiment log from template
#
# Usage:
#   ./scripts/create-experiment-log.sh <experiment-name>
#
# Example:
#   ./scripts/create-experiment-log.sh qml-scaling-test

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECTS_ROOT="$(dirname "$SCRIPT_DIR")"
TEMPLATE="$PROJECTS_ROOT/templates/experiment-log.md"
EXPERIMENTS_DIR="$PROJECTS_ROOT/research/experiments"

# Check arguments
if [ $# -eq 0 ]; then
    echo "Usage: ./scripts/create-experiment-log.sh <experiment-name>"
    echo ""
    echo "Example:"
    echo "  ./scripts/create-experiment-log.sh qml-scaling-test"
    exit 1
fi

EXPERIMENT_NAME="$1"
DATE=$(date +%Y-%m-%d)
OUTPUT_FILE="$EXPERIMENTS_DIR/${EXPERIMENT_NAME}-${DATE}.md"

# Check if file already exists
if [ -f "$OUTPUT_FILE" ]; then
    echo "Experiment log already exists: $OUTPUT_FILE"
    exit 1
fi

# Create experiments directory if needed
mkdir -p "$EXPERIMENTS_DIR"

# Check if template exists
if [ ! -f "$TEMPLATE" ]; then
    echo "Error: Template not found at $TEMPLATE"
    exit 1
fi

# Replace template variables
sed "s/experiment-name/$EXPERIMENT_NAME/g" "$TEMPLATE" | \
sed "s/{{date:YYYY-MM-DD}}/$DATE/g" > "$OUTPUT_FILE"

echo "✓ Created experiment log: $OUTPUT_FILE"
echo ""
echo "Next steps:"
echo "  1. Fill in the experiment details"
echo "  2. Link from daily note: [[daily/$DATE]]"
echo "  3. Run your experiment"
echo "  4. Document results in the log"
echo ""
echo "Open in editor? (y/n)"
read -r response
if [ "$response" = "y" ]; then
    if command -v code &> /dev/null; then
        code "$OUTPUT_FILE"
    else
        open "$OUTPUT_FILE"
    fi
fi

```
