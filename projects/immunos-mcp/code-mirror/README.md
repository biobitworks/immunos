---
project: immunos-mcp
source: README.md
type: code-mirror
language: md
size: 6676
modified: 2025-11-17T10:30:15.481588
hash: f756209b10cbb6b21d0ab2e804a82112
tags: [code, source, auto-generated]
---

> [!info] Source File
> **Path**: `README.md`
> **Size**: 6676 bytes
> **Modified**: 2025-11-17
> **Project**: [[projects/immunos-mcp/|immunos-mcp]]
```markdown
# IMMUNOS-MCP

**Artificial Immune System Multi-Agent LLM Server via Model Context Protocol**

An MCP server implementing biological immune system principles for adaptive pattern recognition, anomaly detection, and content validation using multiple specialized AI agents.

## Overview

IMMUNOS-MCP creates a multi-agent system where each agent plays a specific immune cell role:

| Immune Component | AI Agent Role | Capabilities |
|-----------------|---------------|--------------|
| **T Helper Cells** | Orchestrator | Coordinates multi-agent workflows, routes tasks |
| **B Cells** | Pattern Matcher | Recognizes specific patterns, learns from examples |
| **NK Cells** | Anomaly Detector | Detects novel threats without prior training (Negative Selection) |
| **Dendritic Cells** | Context Processor | Extracts features, processes signals |
| **Memory Cells** | Knowledge Cache | Stores and retrieves learned patterns |
| **T Killer Cells** | Validator | Quality control, adversarial detection |
| **T Regulatory Cells** | Calibrator | Confidence calibration, threshold adjustment |

## Architecture

```
┌─────────────────────────────────────────────────────┐
│         MCP Orchestrator (T Helper Agent)            │
│       Coordinates immune response workflow            │
└──────────────────┬──────────────────────────────────┘
                   │ MCP Protocol (JSON-RPC)
        ┌──────────┼──────────┬──────────────┐
        │          │          │              │
┌───────▼──┐ ┌────▼────┐ ┌───▼──────┐ ┌────▼─────┐
│ B Cell   │ │ NK Cell │ │Dendritic │ │ Memory   │
│ Agent    │ │ Agent   │ │   Cell   │ │  Agent   │
│          │ │         │ │  Agent   │ │          │
│Pattern   │ │Anomaly  │ │Feature   │ │Cached    │
│Matching  │ │Detection│ │Extract   │ │Results   │
└──────────┘ └─────────┘ └──────────┘ └──────────┘
```

## Key Features

### 1. Negative Selection Algorithm (NK Cells)
- Trained on "self" (normal patterns)
- Detects "non-self" (anomalies, threats) without explicit training
- Zero-shot anomaly detection using LLM reasoning

### 2. Adaptive Pattern Recognition (B Cells)
- Hybrid approach: Traditional Immunos-81 + LLM embeddings
- Online learning from new examples
- Affinity maturation for improved accuracy

### 3. Multi-Signal Processing (Dendritic Cells)
- PAMP (Pathogen-Associated Molecular Patterns): Known threats
- Danger signals: Context-based warnings
- Safe signals: Verified benign patterns

### 4. Semantic Memory (Memory Cells)
- Vector database for fast pattern retrieval
- Few-shot learning from cached examples
- Memory consolidation and pruning

## Installation

```bash
cd immunos-mcp
uv pip install -e .

# With development dependencies
uv pip install -e ".[dev]"
```

## Usage

### As MCP Server

```bash
# Start the orchestrator server
uv run python src/servers/orchestrator_server.py
```

### Programmatic API

```python
from immunos_mcp.agents import BCellAgent, NKCellAgent, DendriticAgent
from immunos_mcp.core import Antigen

# Create agents
bcell = BCellAgent()
nk_cell = NKCellAgent()
dendritic = DendriticAgent()

# Process input
antigen = Antigen.from_text("Suspicious email content here...")
features = dendritic.extract_features(antigen)

# Pattern matching
affinity = bcell.calculate_affinity(features)

# Anomaly detection
is_anomaly, confidence = nk_cell.detect_novelty(features)

print(f"Pattern match: {affinity:.2f}")
print(f"Is anomaly: {is_anomaly} (confidence: {confidence:.2f})")
```

### Example: Code Review

```python
from immunos_mcp import ImmunosMCP

# Initialize system
immunos = ImmunosMCP()

# Train on your codebase (self)
immunos.train_on_normal_code("path/to/your/repo")

# Review new code (detect non-self)
result = immunos.review_code("""
import os
eval(os.environ.get('MALICIOUS_CODE'))
""")

print(result.classification)  # "anomalous"
print(result.confidence)      # 0.95
print(result.explanation)     # "NK Cell detected: eval() with environment variable"
```

## Applications

### 1. Self/Non-Self Code Security
- Train on organization's codebase
- Detect malicious or suspicious patterns
- Automated security review for PRs

### 2. Conversation Safety
- Detect prompt injection attempts
- Flag adversarial inputs
- Adaptive learning from reports

### 3. Multi-Agent Validation
- Ensemble reasoning for high-stakes decisions
- Consensus building across agents
- Identify uncertain cases for human review

### 4. Online Learning
- Continuously improve from feedback
- No full retraining required
- Adaptive to new patterns

## Project Structure

```
immunos-mcp/
├── src/
│   ├── servers/
│   │   └── orchestrator_server.py    # Main MCP server
│   ├── agents/
│   │   ├── bcell_agent.py            # Pattern matching
│   │   ├── nk_cell_agent.py          # Anomaly detection
│   │   ├── dendritic_agent.py        # Feature extraction
│   │   └── memory_agent.py           # Knowledge cache
│   ├── core/
│   │   ├── antigen.py                # Data representations
│   │   ├── affinity.py               # Similarity calculations
│   │   └── protocols.py              # Shared structures
│   └── config/
│       └── mcp_config.json           # Configuration
├── tests/
├── examples/
└── docs/
```

## Development

### Running Tests

```bash
pytest tests/
```

### Code Formatting

```bash
black src/ tests/
ruff check src/ tests/
```

## Related Projects

- **immunos81**: Original Immunos-81 implementation (Hunt & Cooke, 2000)
- **MCP SDK**: Model Context Protocol framework
- **ChromaDB**: Vector database for memory agent

## References

- Hunt, J. E., & Cooke, D. E. (2000). An immune system-based pattern recognition approach to the classification of outcome from clinical procedures. *JAMIA*, 7(1), 28-41.
- Model Context Protocol: https://github.com/modelcontextprotocol
- Negative Selection Algorithm: de Castro & Von Zuben (2002)
- Dendritic Cell Algorithm: Greensmith et al. (2005)

## License

MIT License (consistent with immunos81 project)

```
