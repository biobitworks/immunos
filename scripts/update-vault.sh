#!/bin/bash
# Update Obsidian Vault
#
# Runs all maintenance scripts to keep the vault up-to-date:
# - Sync code mirrors
# - Regenerate API documentation
# - Update statistics
# - Check for broken links (future)
#
# Usage:
#   ./scripts/update-vault.sh [project-name]
#   ./scripts/update-vault.sh all

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECTS_ROOT="$(dirname "$SCRIPT_DIR")"

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}╔═══════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   Obsidian Vault Update Script       ║${NC}"
echo -e "${BLUE}╚═══════════════════════════════════════╝${NC}"
echo ""

# Determine which projects to update
if [ $# -eq 0 ] || [ "$1" = "all" ]; then
    PROJECTS=("immunos-mcp" "immunos81" "experiments")
else
    PROJECTS=("$1")
fi

# Update each project
for project in "${PROJECTS[@]}"; do
    echo -e "${GREEN}► Updating project: ${project}${NC}"
    echo ""

    # Check if project exists
    if [ ! -d "$PROJECTS_ROOT/$project" ]; then
        echo -e "${YELLOW}⚠ Project directory not found: $project${NC}"
        echo -e "${YELLOW}  Skipping...${NC}"
        echo ""
        continue
    fi

    # 1. Sync code to markdown
    echo -e "${BLUE}  [1/2] Syncing code to markdown...${NC}"
    if python3 "$SCRIPT_DIR/sync-code-to-obsidian.py" "$project"; then
        echo -e "${GREEN}  ✓ Code mirrored successfully${NC}"
    else
        echo -e "${YELLOW}  ⚠ Code mirror failed${NC}"
    fi
    echo ""

    # 2. Generate API documentation
    echo -e "${BLUE}  [2/2] Generating API documentation...${NC}"
    if python3 "$SCRIPT_DIR/generate-api-docs.py" "$project"; then
        echo -e "${GREEN}  ✓ API docs generated${NC}"
    else
        echo -e "${YELLOW}  ⚠ API docs generation failed${NC}"
    fi
    echo ""

    echo -e "${GREEN}✓ Project ${project} updated${NC}"
    echo ""
done

# Generate vault statistics
echo -e "${BLUE}► Generating vault statistics...${NC}"

# Count files by type
DAILY_NOTES=$(find "$PROJECTS_ROOT/daily" -name "*.md" 2>/dev/null | wc -l || echo 0)
PROJECT_JOURNALS=$(find "$PROJECTS_ROOT/projects" -path "*/journal/*.md" 2>/dev/null | wc -l || echo 0)
EXPERIMENTS=$(find "$PROJECTS_ROOT/research/experiments" -name "*.md" 2>/dev/null | wc -l || echo 0)
LITERATURE=$(find "$PROJECTS_ROOT/research/literature" -name "*.md" -not -name "README.md" 2>/dev/null | wc -l || echo 0)
CODE_MIRRORS=$(find "$PROJECTS_ROOT/projects" -path "*/code-mirror/*.md" 2>/dev/null | wc -l || echo 0)
API_DOCS=$(find "$PROJECTS_ROOT/projects" -path "*/api/*.md" 2>/dev/null | wc -l || echo 0)

# Display statistics
echo ""
echo -e "${BLUE}╔═══════════════════════════════════════╗${NC}"
echo -e "${BLUE}║   Vault Statistics                    ║${NC}"
echo -e "${BLUE}╚═══════════════════════════════════════╝${NC}"
echo ""
echo "  Daily Notes:       $DAILY_NOTES"
echo "  Project Journals:  $PROJECT_JOURNALS"
echo "  Experiments:       $EXPERIMENTS"
echo "  Literature Notes:  $LITERATURE"
echo "  Code Mirrors:      $CODE_MIRRORS"
echo "  API Docs:          $API_DOCS"
echo ""
TOTAL=$((DAILY_NOTES + PROJECT_JOURNALS + EXPERIMENTS + LITERATURE + CODE_MIRRORS + API_DOCS))
echo -e "${GREEN}  Total Documents:   $TOTAL${NC}"
echo ""

echo -e "${GREEN}✓ Vault update complete!${NC}"
echo ""
echo -e "Next steps:"
echo -e "  - Review changes in Obsidian"
echo -e "  - Update [[daily/$(date +%Y-%m-%d)]] with today's work"
echo -e "  - Check for any broken links"
echo ""
