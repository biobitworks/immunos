# IMMUNOS Model Architecture

## Overview

IMMUNOS uses a biological immune system metaphor for AI-assisted research. Each model maps to an immune cell type based on its function.

## Model-to-Cell Mappings

### Chat Models (Interactive Agents)

| Model | Cell Type | Function | Size | Use Case |
|-------|-----------|----------|------|----------|
| `qwen2.5:1.5b` | **NK Cell** | Quick Scanner | 986 MB | Fast anomaly detection, quick lookups, simple queries |
| `qwen2.5-coder:7b` | **B Cell** | Code Verifier | 4.7 GB | Code generation, debugging, implementation, validation |
| `deepseek-r1:14b` | **Dendritic Cell** | Research Analyst | 9.0 GB | Research analysis, reasoning, strategic planning, reporting |
| `qwen2.5:7b` | **T Cell** | Memory Agent | 4.7 GB | Context persistence, adaptive learning, memory retrieval |

### Embedding Models (Memory Infrastructure)

| Model | Cell Type | Function | Size | Use Case |
|-------|-----------|----------|------|----------|
| `bge-m3` | **T Cell Memory** | Multilingual Embeddings | 1.2 GB | Cross-language semantic search, international content |
| `mxbai-embed-large` | **T Cell Memory** | High-Quality Embeddings | 669 MB | Precise semantic matching, research retrieval |
| `nomic-embed-text` | **T Cell Memory** | Fast Embeddings | 274 MB | Quick lookups, real-time search |

### Cloud Models (Orchestration)

| Model | Cell Type | Function | Access | Use Case |
|-------|-----------|----------|--------|----------|
| Claude Opus 4.5 | **Thymus** | Central Validator | Cloud API | Complex planning, multi-step tasks, validation |
| Claude Sonnet 4.5 | **Orchestrator** | Task Router | Cloud API | Task routing, coordination, comprehensive work |

## Cell Type Descriptions

### NK Cell (Natural Killer)
- **Role**: First-line defense, anomaly detection
- **Behavior**: Fast, aggressive scanning without prior sensitization
- **Model**: Small, fast models for quick decisions
- **Tasks**: Quick scans, anomaly flags, simple lookups

### B Cell
- **Role**: Verification and validation
- **Behavior**: Produces specific "antibodies" (code, fixes) for identified problems
- **Model**: Code-specialized models
- **Tasks**: Code generation, bug fixes, implementation verification

### Dendritic Cell
- **Role**: Information processing and reporting
- **Behavior**: Captures, processes, and presents information
- **Model**: Large reasoning models
- **Tasks**: Research analysis, daily synthesis, journaling, strategic planning

### T Cell
- **Role**: Adaptive memory and learning
- **Behavior**: Remembers past encounters, provides long-term immunity
- **Model**: Context-aware models + embedding infrastructure
- **Tasks**: Context persistence, memory retrieval, pattern recognition

### Thymus
- **Role**: Central training and validation
- **Behavior**: Trains T cells, negative selection (removing self-reactive cells)
- **Model**: Most capable cloud models
- **Tasks**: Complex validation, training oversight, quality control

## Ollama Commands

```bash
# List installed models
ollama list

# Pull missing models
ollama pull qwen2.5:7b          # T Cell chat model
ollama pull qwen2.5:1.5b        # NK Cell (if missing)
ollama pull qwen2.5-coder:7b    # B Cell
ollama pull deepseek-r1:14b     # Dendritic Cell
ollama pull bge-m3              # T Cell Memory (multilingual)
ollama pull mxbai-embed-large   # T Cell Memory (high quality)
ollama pull nomic-embed-text    # T Cell Memory (fast)

# Run a model interactively
ollama run qwen2.5-coder:7b

# Check model info
ollama show qwen2.5-coder:7b
```

## Dashboard Integration

The IMMUNOS dashboard displays models with their cell type labels:
- `[NK Cell] qwen2.5:1.5b (1.0GB)` - Quick Tasks
- `[B Cell] qwen2.5-coder:7b (4.7GB)` - Code Implementation
- `[Dendritic] deepseek-r1:14b (9.0GB)` - Research Analysis
- `[T Cell] qwen2.5:7b (4.7GB)` - Memory & Context

Embedding models are used internally for RAG and don't appear in chat dropdown.

## Recommended Additions

For enhanced capabilities, consider adding:

| Model | Cell Type | Size | Benefit |
|-------|-----------|------|---------|
| `qwen2.5:7b` | T Cell | 4.7 GB | Better context handling than 1.5b |
| `llama3.2:3b` | NK Cell Alt | 2.0 GB | Alternative fast scanner |
| `codellama:7b` | B Cell Alt | 3.8 GB | Alternative code model |

## Memory Requirements

| Configuration | Total VRAM | Models Loaded |
|---------------|------------|---------------|
| Minimal | 4 GB | 1 small model |
| Standard | 8 GB | 1-2 medium models |
| Full | 16+ GB | Multiple models, fast switching |

Your M1 Pro (11.8 GB) can run most configurations with 1 model loaded at a time.

---

**Last Updated**: 2025-12-25
**Maintained by**: IMMUNOS System
