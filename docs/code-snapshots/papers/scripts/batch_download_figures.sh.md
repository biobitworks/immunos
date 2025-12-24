---
source: /Users/byron/projects/papers/scripts/batch_download_figures.sh
relative: papers/scripts/batch_download_figures.sh
generated_at: 2025-12-23 10:28
---

```bash
#!/bin/bash
# Batch download figures for failed papers

set -e

cd /Users/byron/projects/papers

# List of failed papers as DOI|FigureID pairs
PAPERS=(
    "10.1038/s43587-025-01014-w|43587_2025_1014"
    "10.1038/s42003-025-09219-w|42003_2025_9219"
    "10.1038/s41587-025-02830-6|41587_2025_2830"
    "10.1038/s41514-025-00290-5|41514_2025_290"
    "10.1038/s41598-025-22076-1|41598_2025_22076"
    "10.1038/s41467-025-65664-5|41467_2025_65664"
    "10.1038/s43587-025-01011-z|43587_2025_1011"
    "10.1038/s41591-025-03999-8|41591_2025_3999"
    "10.1038/s43587-025-00961-8|43587_2025_961"
    "10.1038/s43587-025-00986-z|43587_2025_986"
    "10.1038/s41467-018-03770-3|41467_2018_3770"
    "10.1038/s41598-025-26154-2|41598_2025_26154"
    "10.1038/s41467-025-64652-z|41467_2025_64652"
    "10.1038/s43587-025-00943-w|43587_2025_943"
    "10.1038/s43856-025-01239-1|43856_2025_1239"
    "10.1038/s41467-025-64275-4|41467_2025_64275"
    "10.1038/s41467-025-63429-8|41467_2025_63429"
    "10.1038/s41598-025-25545-9|41598_2025_25545"
    "10.1038/s41467-025-65297-8|41467_2025_65297"
    "10.1038/s43856-025-01222-w|43856_2025_1222"
    "10.1038/s41514-025-00287-0|41514_2025_287"
    "10.1038/s41467-025-63229-0|41467_2025_63229"
    "10.1038/s41598-025-27042-5|41598_2025_27042"
    "10.1038/s41467-025-64835-8|41467_2025_64835"
    "10.1038/s43587-025-01021-x|43587_2025_1021"
)

echo "Downloading figures for ${#PAPERS[@]} papers..."
echo "=================================================================="

PAPER_COUNT=0
for PAPER in "${PAPERS[@]}"; do
    PAPER_COUNT=$((PAPER_COUNT + 1))
    DOI=$(echo "$PAPER" | cut -d'|' -f1)
    FIG_ID=$(echo "$PAPER" | cut -d'|' -f2)
    DIR=$(echo "$DOI" | tr '/' '_')

    echo ""
    echo "[$PAPER_COUNT/${#PAPERS[@]}] $DOI"

    # Add delay between papers
    if [ $PAPER_COUNT -gt 1 ]; then
        sleep 5
    fi

    COUNT=0
    CONSECUTIVE_FAILS=0
    for i in {1..20}; do
        FIG_URL="https://media.springernature.com/lw1200/springer-static/image/art%3A${DOI}/MediaObjects/${FIG_ID}_Fig${i}_HTML.png"
        FIG_PATH="$DIR/figure_${i}.png"

        # Add delay between figures
        if [ $i -gt 1 ]; then
            sleep 2
        fi

        # Download
        HTTP_CODE=$(curl -s -w "%{http_code}" -o "$FIG_PATH" "$FIG_URL")

        # Check if successful
        if [ "$HTTP_CODE" = "200" ] && [ -f "$FIG_PATH" ]; then
            if file "$FIG_PATH" | grep -q "PNG image"; then
                SIZE=$(du -h "$FIG_PATH" | cut -f1)
                echo "  ✓ Figure $i ($SIZE)"
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

        if [ $CONSECUTIVE_FAILS -ge 3 ]; then
            break
        fi
    done

    echo "  Total: $COUNT figures"
done

echo ""
echo "=================================================================="
echo "Download complete!"
echo ""
echo "Summary:"
find . -maxdepth 2 -name "figure_*.png" -size +10k | sed 's|/figure.*||' | sort -u | wc -l | xargs echo "Papers with figures:"

```
