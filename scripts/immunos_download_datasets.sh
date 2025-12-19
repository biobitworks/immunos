#!/bin/bash
"""
IMMUNOS Dataset Download Script
================================

Downloads datasets for all 5 domains.

Strategy: Start with small datasets for quick testing
- Emotion: CK+ (requires manual download - license)
- Hallucination: TruthfulQA (git clone)
- Network: NSL-KDD (wget)
- Code: DiverseVul (git clone)
- Research: SciFact (git clone)

Total size for small datasets: ~600MB
Estimated download time: 5-10 minutes
"""

set -e  # Exit on error

# Create data directory
DATA_DIR="${HOME}/immunos_data"
echo "Creating data directory: ${DATA_DIR}"
mkdir -p "${DATA_DIR}"/{emotion,hallucination,network,code,research}

# Color output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}IMMUNOS Dataset Download${NC}"
echo "================================"
echo ""

# ============================================================================
# 1. HALLUCINATION DETECTION (easiest - git clone)
# ============================================================================
echo -e "${GREEN}[1/5] Downloading Hallucination Detection Datasets${NC}"

if [ ! -d "${DATA_DIR}/hallucination/truthfulqa" ]; then
    echo "Downloading TruthfulQA (817 questions, ~5MB)..."
    git clone --depth 1 https://github.com/sylinrl/TruthfulQA "${DATA_DIR}/hallucination/truthfulqa"
    echo "✓ TruthfulQA downloaded"
else
    echo "✓ TruthfulQA already exists"
fi

if [ ! -d "${DATA_DIR}/hallucination/halueval" ]; then
    echo "Downloading HaluEval (35K samples, ~500MB)..."
    git clone --depth 1 https://github.com/RUCAIBox/HaluEval "${DATA_DIR}/hallucination/halueval"
    echo "✓ HaluEval downloaded"
else
    echo "✓ HaluEval already exists"
fi

echo ""

# ============================================================================
# 2. NETWORK INTRUSION DETECTION
# ============================================================================
echo -e "${GREEN}[2/5] Downloading Network Intrusion Datasets${NC}"

if [ ! -f "${DATA_DIR}/network/KDDTrain+.txt" ]; then
    echo "Downloading NSL-KDD (148K records, ~25MB)..."
    cd "${DATA_DIR}/network"

    # Download NSL-KDD dataset (use curl for macOS compatibility)
    curl -sL -o "KDDTrain+.txt" "https://raw.githubusercontent.com/defcom17/NSL_KDD/master/KDDTrain%2B.txt"
    curl -sL -o "KDDTest+.txt" "https://raw.githubusercontent.com/defcom17/NSL_KDD/master/KDDTest%2B.txt"
    curl -sL -o "KDDTrain+_20Percent.txt" "https://raw.githubusercontent.com/defcom17/NSL_KDD/master/KDDTrain%2B_20Percent.txt"

    echo "✓ NSL-KDD downloaded"
    cd - > /dev/null
else
    echo "✓ NSL-KDD already exists"
fi

echo ""

# ============================================================================
# 3. CODE VULNERABILITY DETECTION
# ============================================================================
echo -e "${GREEN}[3/5] Downloading Code Vulnerability Datasets${NC}"

if [ ! -d "${DATA_DIR}/code/diversevul" ]; then
    echo "Downloading DiverseVul (18K functions, ~1GB)..."
    git clone --depth 1 https://github.com/wagner-group/diversevul "${DATA_DIR}/code/diversevul"
    echo "✓ DiverseVul downloaded"
else
    echo "✓ DiverseVul already exists"
fi

echo ""

# ============================================================================
# 4. RESEARCH PAPER VERIFICATION
# ============================================================================
echo -e "${GREEN}[4/5] Downloading Research Verification Datasets${NC}"

if [ ! -d "${DATA_DIR}/research/scifact" ]; then
    echo "Downloading SciFact (1.4K claims, ~50MB)..."
    git clone --depth 1 https://github.com/allenai/scifact "${DATA_DIR}/research/scifact"
    echo "✓ SciFact downloaded"
else
    echo "✓ SciFact already exists"
fi

echo ""

# ============================================================================
# 5. EMOTION RECOGNITION (manual download required)
# ============================================================================
echo -e "${GREEN}[5/5] Emotion Recognition Datasets${NC}"
echo -e "${YELLOW}⚠ CK+ dataset requires manual download (license agreement)${NC}"
echo "Download from: http://www.jeffcohn.net/Resources/"
echo "After downloading, extract to: ${DATA_DIR}/emotion/ckplus"
echo ""
echo -e "${YELLOW}⚠ FER2013 dataset requires Kaggle account${NC}"
echo "Download from: https://www.kaggle.com/datasets/msambare/fer2013"
echo "After downloading, extract to: ${DATA_DIR}/emotion/fer2013"
echo ""

# ============================================================================
# SUMMARY
# ============================================================================
echo "================================"
echo -e "${GREEN}Download Summary${NC}"
echo "================================"
echo ""

# Check what's downloaded
downloaded=0
total=5

if [ -d "${DATA_DIR}/hallucination/truthfulqa" ]; then
    echo "✓ Hallucination: TruthfulQA"
    ((downloaded++))
fi

if [ -f "${DATA_DIR}/network/KDDTrain+.txt" ]; then
    echo "✓ Network: NSL-KDD"
    ((downloaded++))
fi

if [ -d "${DATA_DIR}/code/diversevul" ]; then
    echo "✓ Code: DiverseVul"
    ((downloaded++))
fi

if [ -d "${DATA_DIR}/research/scifact" ]; then
    echo "✓ Research: SciFact"
    ((downloaded++))
fi

if [ -d "${DATA_DIR}/emotion/ckplus" ]; then
    echo "✓ Emotion: CK+"
    ((downloaded++))
else
    echo "⚠ Emotion: CK+ (manual download required)"
fi

echo ""
echo "Downloaded: ${downloaded}/${total} datasets"
echo ""

# Disk usage
echo "Disk usage:"
du -sh "${DATA_DIR}" 2>/dev/null || echo "No data yet"

echo ""
echo -e "${GREEN}Next steps:${NC}"
echo "1. Verify datasets: ls -R ${DATA_DIR}"
echo "2. Install dependencies: pip install -r requirements.txt"
echo "3. Test extractors: python scripts/test_all_extractors.py"
echo "4. Train models: python scripts/immunos_train.py --quick"
