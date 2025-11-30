---
project: immunos-mcp
source: GETTING_STARTED.md
type: code-mirror
language: md
size: 10998
modified: 2025-11-25T13:27:28.545551
hash: 79ee1d21ef2bfe20c2547b3145b5e52b
tags: [code, source, auto-generated]
---

> [!info] Source File
> **Path**: `GETTING_STARTED.md`
> **Size**: 10998 bytes
> **Modified**: 2025-11-25
> **Project**: [[projects/immunos-mcp/|immunos-mcp]]
```markdown
# IMMUNOS-MCP: Getting Started Guide

## What We've Built

An Artificial Immune System (AIS) based multi-agent framework where AI agents play biological immune cell roles:

### Core Components Implemented

✅ **Core Modules** (`src/core/`):
- `antigen.py`: Data representation (text, code, structured data)
- `affinity.py`: Similarity calculations (traditional + embedding-based)
- `protocols.py`: Shared data structures and communication protocols

✅ **B Cell Agent** (`src/agents/bcell_agent.py`):
- Pattern matching inspired by antibody-producing B cells
- Learns from labeled examples
- Calculates affinity (similarity) to learned patterns
- Groups patterns into clones for collective recognition
- Supports SHA and RHA recognition strategies
- Hybrid approach: traditional Immunos-81 + LLM embeddings

✅ **NK Cell Agent** (`src/agents/nk_cell_agent.py`):
- Anomaly detection using Negative Selection Algorithm
- Trains ONLY on "self" (normal) data
- Generates detectors that don't match self
- Detects "non-self" (anomalies) without explicit training
- Zero-shot threat detection capability
- Danger theory integration for context-aware detection

✅ **Demo** (`examples/simple_demo.py`):
- Email spam classification example
- Shows B Cell and NK Cell working together
- Demonstrates supervised + unsupervised detection

## Project Structure

```
immunos-mcp/
├── src/
│   ├── core/                  ✅ Core data structures
│   │   ├── antigen.py        → Data representation
│   │   ├── affinity.py       → Similarity calculations
│   │   └── protocols.py      → Communication protocols
│   │
│   ├── agents/                ✅ Immune cell agents
│   │   ├── bcell_agent.py    → Pattern matching (done)
│   │   ├── nk_cell_agent.py  → Anomaly detection (done)
│   │   ├── dendritic_agent.py  → Feature extraction (TODO)
│   │   └── memory_agent.py   → Knowledge cache (TODO)
│   │
│   └── servers/               🔜 MCP servers (TODO)
│       └── orchestrator_server.py → Main MCP endpoint
│
├── examples/
│   └── simple_demo.py        ✅ Working demo
│
├── tests/                     🔜 Unit tests (TODO)
├── README.md                  ✅ Project overview
└── pyproject.toml            ✅ Dependencies
```

## Installation

### Dependencies

```bash
cd immunos-mcp

# Core dependencies (numpy, pydantic)
uv pip install -e .

# Optional: MCP support
uv pip install -e ".[mcp]"

# Optional: LLM integration
uv pip install -e ".[llm]"

# Optional: Vector database
uv pip install -e ".[vector]"

# Development dependencies
uv pip install -e ".[dev]"
```

**Note**: ChromaDB has compatibility issues with Python 3.14. Use Python 3.10-3.12 for vector database features.

## Quick Start

### 1. Basic Usage (Python API)

```python
from src.core.antigen import Antigen
from src.agents.bcell_agent import BCellAgent
from src.agents.nk_cell_agent import NKCellAgent
import numpy as np

# Create training data
train_antigens = [
    Antigen.from_text("Normal email text", class_label="normal"),
    Antigen.from_text("SPAM WIN MONEY!!!", class_label="spam"),
]

# Train B Cell (pattern matching)
bcell = BCellAgent(agent_name="bcell_001")
bcell.train(train_antigens, embeddings=None)  # Add embeddings for better accuracy

# Train NK Cell (anomaly detection) - only on normal data
normal_antigens = [a for a in train_antigens if a.class_label == "normal"]
nk_cell = NKCellAgent(agent_name="nk_001")
nk_cell.train_on_self(normal_antigens, embeddings=None)

# Test new input
test_antigen = Antigen.from_text("Free prizes click now!")

# B Cell classification
result = bcell.recognize(test_antigen, strategy="sha")
print(f"Class: {result.predicted_class}, Confidence: {result.confidence:.2f}")

# NK Cell anomaly detection
anomaly = nk_cell.detect_novelty(test_antigen)
print(f"Is anomaly: {anomaly.is_anomaly}, Score: {anomaly.anomaly_score:.2f}")
```

### 2. Run the Demo

```bash
# With proper Python environment (3.10-3.12 + numpy)
python examples/simple_demo.py
```

**Demo Output**:
- Trains B Cell on 5 normal + 4 spam emails
- Trains NK Cell on 5 normal emails only
- Tests on 4 new emails
- Shows both agents' predictions and combined decision

## Key Concepts

### 1. B Cell Agent (Pattern Matching)

**Biological Inspiration**: B cells produce antibodies that bind to specific antigens.

**Implementation**:
- Learns patterns from labeled training examples
- Calculates **affinity** (similarity) to learned patterns
- Groups patterns into **clones** (populations)
- Calculates **avidity** (collective affinity): `sum(affinities) * log(1 + clone_size)`
- Recogn ition strategies:
  - **SHA**: Simple Highest Avidity (winner takes all)
  - **RHA**: Relative Highest Avidity (5% threshold for uncertainty)

**Use Cases**:
- Text classification (spam, sentiment)
- Code classification (bug, feature, security)
- Multi-class pattern recognition

### 2. NK Cell Agent (Anomaly Detection)

**Biological Inspiration**: NK cells detect abnormal cells without prior exposure. T cells undergo negative selection in the thymus to eliminate self-reactive cells.

**Implementation**:
- Trains ONLY on "self" (normal) patterns
- Generates **negative selection detectors** that don't match self
- Detects "non-self" by checking if input:
  1. Doesn't match self patterns
  2. Matches non-self detectors
- Provides zero-shot anomaly detection

**Use Cases**:
- Intrusion detection (train on normal network traffic)
- Code security (train on safe code, detect malicious)
- Content moderation (train on acceptable content, flag violations)

### 3. Hybrid Affinity Calculation

Combines two methods for robustness:

1. **Traditional** (Immunos-81):
   - Numeric: `exp(-|v1 - v2| / sigma)`
   - Nominal: Exact match (1.0 or 0.0)
   - String: Character overlap

2. **Embedding-based**:
   - Cosine similarity of LLM embeddings
   - Semantic understanding

3. **Hybrid** (default):
   - Weighted combination (30% traditional + 70% embedding)
   - Best of both: structured + semantic

## Immune System Mapping

| Biological Component | Computational Role | IMMUNOS-MCP Implementation |
|---------------------|-------------------|---------------------------|
| **B Cells** | Antibody production, pattern recognition | `BCellAgent`: Learns specific patterns, calculates affinity |
| **Plasma Cells** | Antibody factories | `Clone`: Groups of patterns with collective avidity |
| **NK Cells** | Innate immunity, autonomous killing | `NKCellAgent`: Negative selection, zero-shot detection |
| **T Helper Cells** | Orchestration, coordination | `OrchestratorServer` (TODO): Coordinates all agents |
| **Dendritic Cells** | Antigen presentation, signal processing | `DendriticAgent` (TODO): Feature extraction |
| **Memory Cells** | Rapid recall of past infections | `MemoryAgent` (TODO): Cached pattern retrieval |
| **T Killer Cells** | Eliminate infected cells | `TKillerAgent` (TODO): Validation, filtering |
| **T Regulatory Cells** | Balance, prevent autoimmunity | `TRegulatoryAgent` (TODO): Calibration, threshold tuning |

## Applications

### 1. Code Security Review

```python
# Train on your organization's codebase (self)
normal_code_files = load_codebase("path/to/repo")
nk_cell.train_on_self(normal_code_files)

# Review new PR
pr_code = Antigen.from_code(pr_content, language="python")
result = nk_cell.detect_novelty(pr_code)

if result.is_anomaly:
    print(f"⚠️ Suspicious code detected: {result.explanation}")
```

### 2. Multi-Agent Validation

```python
# Multiple B Cell agents with different perspectives
bcell1 = BCellAgent("conservative_classifier")
bcell2 = BCellAgent("aggressive_classifier")

# Get predictions from both
result1 = bcell1.recognize(test_input)
result2 = bcell2.recognize(test_input)

# Ensemble decision
if result1.predicted_class == result2.predicted_class:
    print(f"High confidence: {result1.predicted_class}")
else:
    print(f"Uncertain: {result1.predicted_class} vs {result2.predicted_class}")
```

### 3. Online Learning

```python
# Add new pattern without retraining
new_example = Antigen.from_text("New spam pattern", class_label="spam")
bcell.add_pattern(new_example, embedding=new_embedding)

# System adapts immediately
```

## Next Steps (TODO)

### Phase 1: Complete Core Agents (Priority)

- [ ] **Dendritic Cell Agent**: Feature extraction, signal processing
- [ ] **Memory Agent**: Vector database integration (ChromaDB/FAISS)

### Phase 2: MCP Server Integration

- [ ] **Orchestrator Server**: Main MCP endpoint that coordinates agents
- [ ] Expose tools: `classify_input`, `detect_anomaly`, `explain_decision`
- [ ] Expose resources: `recognition_strategies`, `system_metrics`
- [ ] Handle JSON-RPC communication

### Phase 3: LLM Integration

- [ ] Real LLM embeddings (OpenAI, Anthropic, local models)
- [ ] LLM-as-judge for validation (T Killer Cell)
- [ ] LLM-based affinity maturation
- [ ] Zero-shot reasoning for NK Cell

### Phase 4: Advanced Features

- [ ] T Killer Agent: Output validation, adversarial detection
- [ ] T Regulatory Agent: Confidence calibration
- [ ] Affinity maturation: LLM-guided pattern improvement
- [ ] Danger theory: Multi-signal integration

### Phase 5: Production Ready

- [ ] Comprehensive unit tests
- [ ] Performance benchmarks
- [ ] Docker deployment
- [ ] API documentation
- [ ] Example applications

## Research Background

This implementation is based on research from:

1. **Immunos-81** (Hunt & Cooke, 2000):
   - Original AIS for medical diagnosis
   - SHA/RHA recognition strategies
   - Affinity/avidity calculations

2. **Negative Selection Algorithm** (Forrest et al., 1994):
   - Self/non-self discrimination
   - Anomaly detection without labeled anomalies
   - Thymic selection analogy

3. **Danger Theory** (Matzinger, 2002):
   - Context-aware immune response
   - Multiple signal types (PAMP, danger, safe)
   - Applied to dendritic cell algorithm

4. **Recent AIS-LLM Hybrids** (2024-2025):
   - Bio AI Agent: Multi-agent for CAR-T therapy
   - Ab-Affinity: LLM + genetic algorithms
   - Elloe AI: Immune system for LLM safety

## Contributing

This project is in early development. Key areas for contribution:

1. **Dendritic Cell Implementation**: Feature extraction algorithms
2. **Memory Agent**: Vector database integration and caching
3. **MCP Server**: Full server implementation with tool exposure
4. **LLM Integration**: Embeddings, reasoning, validation
5. **Benchmarks**: Compare against standard ML algorithms
6. **Applications**: Real-world use cases and demos

## License

MIT License (consistent with immunos81 parent project)

## References

- Hunt, J. E., & Cooke, D. E. (2000). Immunos-81 paper. *JAMIA*, 7(1), 28-41.
- Model Context Protocol: https://github.com/modelcontextprotocol
- Your research summary: `/Users/byron/projects/immunos-mcp/RESEARCH_SUMMARY.md` (see agent research output)

```
