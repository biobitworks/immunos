---
source: /Users/byron/projects/scripts/create-daily-note.sh
relative: scripts/create-daily-note.sh
generated_at: 2025-12-23 10:28
---

```bash
#!/bin/bash
# Create Daily Note
#
# Creates a new daily note from template for today (or specified date)
#
# Usage:
#   ./scripts/create-daily-note.sh           # Create for today
#   ./scripts/create-daily-note.sh 2025-12-01  # Create for specific date

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECTS_ROOT="$(dirname "$SCRIPT_DIR")"
TEMPLATE="$PROJECTS_ROOT/templates/daily-note.md"
DAILY_DIR="$PROJECTS_ROOT/daily"

# Get date (default to today)
if [ $# -eq 0 ]; then
    DATE=$(date +%Y-%m-%d)
else
    DATE="$1"
fi

# Validate date format
if ! [[ "$DATE" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
    echo "Error: Invalid date format. Use YYYY-MM-DD"
    exit 1
fi

OUTPUT_FILE="$DAILY_DIR/$DATE.md"

# Check if file already exists
if [ -f "$OUTPUT_FILE" ]; then
    echo "Daily note already exists: $OUTPUT_FILE"
    echo "Open in Obsidian? (y/n)"
    read -r response
    if [ "$response" = "y" ]; then
        open "obsidian://open?vault=projects&file=daily/$DATE"
    fi
    exit 0
fi

# Create daily directory if needed
mkdir -p "$DAILY_DIR"

# Check if template exists
if [ ! -f "$TEMPLATE" ]; then
    echo "Error: Template not found at $TEMPLATE"
    exit 1
fi

# Replace template variables
sed "s/{{date:YYYY-MM-DD}}/$DATE/g" "$TEMPLATE" > "$OUTPUT_FILE"

# Calculate yesterday and tomorrow
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    YESTERDAY=$(date -v-1d -j -f "%Y-%m-%d" "$DATE" +%Y-%m-%d)
    TOMORROW=$(date -v+1d -j -f "%Y-%m-%d" "$DATE" +%Y-%m-%d)
else
    # Linux
    YESTERDAY=$(date -d "$DATE - 1 day" +%Y-%m-%d)
    TOMORROW=$(date -d "$DATE + 1 day" +%Y-%m-%d)
fi

# Replace yesterday and tomorrow placeholders
sed -i.bak "s/{{yesterday}}/$YESTERDAY/g" "$OUTPUT_FILE"
sed -i.bak "s/{{tomorrow}}/$TOMORROW/g" "$OUTPUT_FILE"
rm "$OUTPUT_FILE.bak"

echo "✓ Created daily note: $OUTPUT_FILE"
echo ""
echo "Open in Obsidian? (y/n)"
read -r response
if [ "$response" = "y" ]; then
    open "obsidian://open?vault=projects&file=daily/$DATE"
fi

```
