#!/bin/bash
# Retry downloading failed papers with rate limiting

set -e

# List of failed papers (DOIs)
DOIS=(
    "10.1038/s43587-025-01014-w"
    "10.1038/s42003-025-09219-w"
    "10.1038/s41587-025-02830-6"
    "10.1038/s41514-025-00290-5"
    "10.1038/s41598-025-22076-1"
    "10.1038/s41467-025-65664-5"
    "10.1038/s43587-025-01011-z"
    "10.1038/s41591-025-03999-8"
    "10.1038/s43587-025-00961-8"
    "10.1038/s43587-025-00986-z"
    "10.1038/s41467-018-03770-3"
    "10.1038/s41598-025-26154-2"
    "10.1038/s41467-025-64652-z"
    "10.1038/s43587-025-00943-w"
    "10.1038/s43856-025-01239-1"
    "10.1038/s41467-025-64275-4"
    "10.1038/s41467-025-63429-8"
    "10.1038/s41598-025-25545-9"
    "10.1038/s41467-025-65297-8"
    "10.1038/s43856-025-01222-w"
    "10.1038/s41514-025-00287-0"
    "10.1038/s41467-025-63229-0"
    "10.1038/s41598-025-27042-5"
    "10.1038/s41467-025-64835-8"
)

echo "Starting batch retry of ${#DOIS[@]} failed papers with rate limiting..."
echo "========================================================================"

# Download each paper with 5 second delay between papers
for i in "${!DOIS[@]}"; do
    doi="${DOIS[$i]}"
    echo ""
    echo "[$((i+1))/${#DOIS[@]}] Downloading $doi..."

    # Add delay before each download (except first)
    if [ $i -gt 0 ]; then
        echo "  Waiting 5 seconds to avoid rate limiting..."
        sleep 5
    fi

    # Run download script
    if ./scripts/download_paper.py --skip-license-check "$doi"; then
        echo "  ✓ Success"
    else
        echo "  ⚠️  Failed (exit code $?)"
    fi
done

echo ""
echo "========================================================================"
echo "Batch retry complete!"
echo ""
echo "Checking results..."
find . -maxdepth 2 -name "paper.pdf" -type f | wc -l | xargs echo "Papers with PDFs:"
