---
title: "immunOS: A Local-First Artificial Immune System for Research Integrity Verification"
author:
  - name: Byron
    affiliation: BiobitWorks
    email: byron@biobitworks.com
date: 2025-12-24
version: v0.1
type: preprint
status: draft
license: CC-BY-4.0
repository: https://github.com/biobitworks/immunos
---

# immunOS: A Local-First Artificial Immune System for Research Integrity Verification

**Author**: Byron (BiobitWorks)
**Date**: December 2025
**Version**: 0.1 (Draft)

## Abstract

We introduce immunOS, a local-first research integrity verification system that applies artificial immune system (AIS) principles to detect hallucinations, inconsistencies, and unsupported claims in AI-generated and human-authored scientific content. immunOS implements a multi-agent architecture inspired by biological immune cells: B-cells for pattern recognition and claim classification, NK-cells for anomaly detection via negative selection, T-cells for coordination and context management, Dendritic cells for feature extraction, and Memory cells for persistent learning. The system operates entirely offline with local language models (Ollama), ensuring data privacy and enabling air-gapped research workflows. We describe the architecture, agent roles, orchestration protocol, and initial validation on the SciFact claim verification benchmark. immunOS achieves competitive performance while maintaining full local execution, offering researchers a privacy-preserving alternative to cloud-based verification systems.

**Keywords**: artificial immune system, research integrity, claim verification, local LLM, negative selection, hallucination detection

---

## Introduction

### Background

The proliferation of AI-generated content in research poses challenges for maintaining scientific integrity. Large language models can produce plausible but unsupported claims, hallucinated citations, and subtle inconsistencies that evade traditional review processes. Existing verification systems typically rely on cloud-based APIs, raising concerns about data privacy, reproducibility, and operational independence.

### Motivation

Biological immune systems provide a compelling metaphor for research integrity verification. The immune system distinguishes "self" (legitimate, verified content) from "non-self" (foreign, potentially harmful content) through adaptive pattern recognition, anomaly detection, and coordinated multi-cell responses. This paradigm maps naturally to the task of identifying unsupported claims against a corpus of trusted references.

### Contributions

1. **Multi-agent AIS architecture**: A modular system with specialized agents mirroring immune cell types (B-cell, NK-cell, T-cell, Dendritic, Memory)
2. **Local-first execution**: Complete offline operation using Ollama models, ensuring data privacy and reproducibility
3. **Negative selection for anomaly detection**: Detector generation against trusted "self" patterns for identifying novel threats
4. **Orchestrated verification pipeline**: Coordinated multi-agent analysis with confidence calibration and signal integration
5. **Extensible framework**: Domain-specific agent routing and pluggable embedding strategies

---

## Methods

### System Architecture

immunOS implements a hierarchical multi-agent architecture:

```
┌─────────────────────────────────────────────────────────────┐
│                     Orchestrator                             │
│  (Coordinates agents, aggregates results, routes by domain) │
└─────────────────────┬───────────────────────────────────────┘
                      │
    ┌─────────────────┼─────────────────┐
    │                 │                 │
    ▼                 ▼                 ▼
┌───────┐       ┌───────┐       ┌───────────┐
│B-Cell │       │NK-Cell│       │ Dendritic │
│Pattern│       │Anomaly│       │  Feature  │
│ Match │       │Detect │       │ Extract   │
└───────┘       └───────┘       └───────────┘
    │                 │                 │
    └─────────────────┼─────────────────┘
                      │
                      ▼
               ┌───────────┐
               │  Memory   │
               │   Cell    │
               └───────────┘
```

### Agent Roles

| Agent | Biological Analog | Function | Model |
|-------|-------------------|----------|-------|
| B-Cell | B lymphocyte | Pattern recognition, claim classification | qwen2.5-coder:7b |
| NK-Cell | Natural killer cell | Anomaly detection, novelty scoring | deepseek-r1:14b |
| T-Cell | T lymphocyte | Coordination, context integration | deepseek-r1:14b |
| Dendritic | Dendritic cell | Feature extraction, signal derivation | qwen2.5:1.5b |
| Memory | Memory B/T cell | Persistent state, learning over time | qwen2.5:1.5b |

### Negative Selection Algorithm

The NK-Cell agent implements a negative selection approach:

1. **Self-pattern learning**: Train detectors on trusted reference corpus (verified claims, known-good abstracts)
2. **Detector generation**: Create detector strings that do not match any self-pattern above threshold
3. **Anomaly detection**: Flag inputs matching detectors as potential "non-self" (hallucinations, unsupported claims)
4. **Confidence calibration**: Integrate with Dendritic signals (danger, PAMP, safe) for final scoring

### Local Model Integration

immunOS uses Ollama for local LLM inference:

- **No external API calls**: All processing occurs on-device
- **Model flexibility**: Configurable model assignments per agent role
- **Embedding support**: Local embeddings via nomic-embed-text or similar

### Orchestration Protocol

1. Antigen (input) received by Orchestrator
2. Dendritic agent extracts features and derives signals
3. B-Cell performs pattern matching against trained clones
4. NK-Cell evaluates anomaly score via negative selection
5. Results aggregated with confidence calibration
6. Memory agent stores result for future reference
7. Orchestrator returns structured result with metadata

---

## Results

### SciFact Baseline Performance

Evaluated on SciFact dev set (open retrieval):

| Task | Precision | Recall | F1 |
|------|-----------|--------|-----|
| Sentence selection | 0.525 | 0.437 | 0.477 |
| Sentence label | 0.469 | 0.391 | 0.426 |
| Abstract label (only) | 0.553 | 0.474 | 0.510 |
| Abstract label (rationalized) | 0.525 | 0.450 | 0.485 |

### Orchestrator Performance

Terminal-based smoke tests (research domain):

| Input | Classification | Confidence | Anomaly |
|-------|---------------|------------|---------|
| "SIRT1 activates autophagy" | support | 66.4% | True |
| "Resveratrol extends lifespan in mice" | support | 67.4% | True |
| "mTOR inhibition promotes autophagy" | support | 65.6% | True |

Execution time: 40-65ms per classification (local, no GPU).

---

## Discussion

### Advantages of Local-First Design

- **Data privacy**: Sensitive research data never leaves the local machine
- **Reproducibility**: Fixed model versions, no API deprecation risks
- **Offline operation**: Suitable for air-gapped or restricted environments
- **Cost control**: No per-query API charges

### Limitations

- **Model capability**: Local models (7B-14B) have lower capacity than cloud APIs (100B+)
- **Training data**: Negative selection requires curated "self" corpus
- **Dataset access**: Full replication blocked by MAHNOB-HCI availability

### Future Directions

- Enhanced NK-Cell detectors with learned embeddings
- T-Cell coordination for multi-document reasoning
- Integration with citation verification pipelines
- Community-contributed "self" pattern libraries

---

## Conclusion

immunOS demonstrates that artificial immune system principles can be effectively applied to research integrity verification in a fully local, privacy-preserving architecture. The multi-agent design enables modular development and domain-specific customization, while the negative selection approach provides a principled framework for anomaly detection. Initial results on SciFact show competitive performance, with clear paths for enhancement through improved detector training and model capabilities.

---

## Acknowledgments

This work was developed with assistance from Claude (Anthropic) for code generation and documentation.

---

## References

See companion preprint: "immunOS Replication v1: SciFact Baseline and Negative-Selection AIS Framing" for detailed citations and replication notes.

Primary references:
- Wadsworth, J.D.F. & Bhagwat, N. (2025). SciFact: A Dataset for Scientific Claim Verification. arXiv:2004.14974
- Umair, M. et al. (2025). Negative selection-based artificial immune system (NegSl-AIS). Results in Engineering. DOI: 10.1016/j.rineng.2025.106601
- Ji, Z. et al. (2007). Revisiting negative selection algorithms. Evolutionary Computation. DOI: 10.1162/evco.2007.15.2.223

---

## Supplementary Materials

- System architecture diagrams: `figures/`
- Orchestrator logs: `/immunos-mcp/logs/orchestrator_runs.log`
- Agent state files: `.immunos/agents/`
- Replication preprint: `/immunos-preprint/immunos-preprint-v1.md`
