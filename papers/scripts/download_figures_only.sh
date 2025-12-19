#!/bin/bash
# Download figures for a Nature paper using curl
# Usage: ./download_figures_only.sh DOI

set -e

DOI="$1"
if [ -z "$DOI" ]; then
    echo "Usage: $0 DOI"
    echo "Example: $0 10.1038/s41467-025-64275-4"
    exit 1
fi

# Convert DOI to directory name
DIR=$(echo "$DOI" | tr '/' '_')

# Extract components from DOI
# Example: 10.1038/s43587-025-01021-x
DOI_SUFFIX=$(echo "$DOI" | sed 's|^10.1038/||')
JOURNAL_CODE=$(echo "$DOI_SUFFIX" | sed 's/^s//' | cut -d'-' -f1)
YEAR_SHORT=$(echo "$DOI_SUFFIX" | cut -d'-' -f2)
ARTICLE_FULL=$(echo "$DOI_SUFFIX" | cut -d'-' -f3)

# Handle year (e.g., "025" -> "2025")
if [ "${#YEAR_SHORT}" -eq 3 ]; then
    YEAR="2${YEAR_SHORT}"
else
    YEAR="$YEAR_SHORT"
fi

# Remove leading zeros and trailing letters from article ID
ARTICLE_ID=$(echo "$ARTICLE_FULL" | sed 's/^0*//' | sed 's/[a-z]*$//')

echo "Downloading figures for $DOI"
echo "  Journal: $JOURNAL_CODE"
echo "  Year: $YEAR"
echo "  Article: $ARTICLE_ID"
echo "  Directory: $DIR"
echo ""

# Create directory if needed
mkdir -p "$DIR"

# Download figures
COUNT=0
CONSECUTIVE_FAILS=0
for i in {1..20}; do
    FIG_URL="https://media.springernature.com/lw1200/springer-static/image/art%3A${DOI}/MediaObjects/${JOURNAL_CODE}_${YEAR}_${ARTICLE_ID}_Fig${i}_HTML.png"
    FIG_PATH="$DIR/figure_${i}.png"

    # Add delay between requests (except first)
    if [ $i -gt 1 ]; then
        sleep 2
    fi

    # Download
    HTTP_CODE=$(curl -s -w "%{http_code}" -o "$FIG_PATH" "$FIG_URL")

    # Check if successful and is real PNG
    if [ "$HTTP_CODE" = "200" ] && [ -f "$FIG_PATH" ]; then
        # Check if it's a real PNG (starts with PNG signature)
        if file "$FIG_PATH" | grep -q "PNG image"; then
            SIZE=$(du -h "$FIG_PATH" | cut -f1)
            echo "  ✓ Figure $i downloaded ($SIZE)"
            COUNT=$((COUNT + 1))
            CONSECUTIVE_FAILS=0
        else
            rm -f "$FIG_PATH"
            CONSECUTIVE_FAILS=$((CONSECUTIVE_FAILS + 1))
        fi
    else
        rm -f "$FIG_PATH"
        CONSECUTIVE_FAILS=$((CONSECUTIVE_FAILS + 1))
    fi

    # Stop after 3 consecutive failures
    if [ $CONSECUTIVE_FAILS -ge 3 ]; then
        break
    fi
done

echo ""
echo "Downloaded $COUNT figures to $DIR/"
