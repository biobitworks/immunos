---
project: immunos-mcp
type: api-documentation
source: src/agents/nk_cell_agent.py
generated: 2025-11-30
tags: [api, documentation, immunos-mcp]
---

# nk_cell_agent

**Source**: [[../code-mirror/src/agents/nk_cell_agent.py|src/agents/nk_cell_agent.py]]

NK Cell Agent for IMMUNOS-MCP

Anomaly detection agent using Negative Selection Algorithm.
Inspired by Natural Killer cells and T-cell negative selection in the thymus.

## Contents

### Classes
- [NKCellAgent](#nkcellagent)

## Classes

### NKCellAgent

**Source**: `src/agents/nk_cell_agent.py:18`

NK Cell Agent - Anomaly Detector via Negative Selection

Mimics NK cell and thymic selection:
- Trains on "self" (normal) data
- Generates detectors that DON'T match self
- Detects "non-self" (anomalies) without explicit training
- Provides zero-shot threat detection

#### Methods

##### `__init__(`self, agent_name: str, detection_threshold: float, num_detectors: int)`

Initialize NK Cell agent.

Args:
    agent_name: Unique identifier
    detection_threshold: Threshold for anomaly detection (0-1)
    num_detectors: Number of negative selection detectors to generate

##### `train_on_self(`self, normal_data: List[Antigen], embeddings: Optional[List[np.ndarray]])`

Train on normal ("self") data.

The negative selection algorithm learns what is normal,
then detects deviations from normal as anomalies.

Args:
    normal_data: List of normal (non-anomalous) antigens
    embeddings: Optional embeddings for normal data

##### `_generate_detectors(`self)`

Generate detectors using negative selection.

Detectors are patterns that DON'T match self.
They are generated randomly and filtered to ensure they don't match normal data.

##### `detect_novelty(`self, antigen: Antigen, antigen_embedding: Optional[np.ndarray])` → `AnomalyResult`

Detect if antigen is novel/anomalous using negative selection.

An antigen is anomalous if it:
1. Doesn't match self patterns (non-self)
2. Matches one or more negative selection detectors

Args:
    antigen: Antigen to test
    antigen_embedding: Optional embedding

Returns:
    AnomalyResult with detection outcome

##### `_generate_explanation(`self, is_anomaly: bool, max_self_sim: float, max_detector_sim: float, best_pattern: Optional[Pattern], detector_matches: List[str])` → `str`

Generate human-readable explanation

##### `assess_danger(`self, antigen: Antigen, context_signals: Optional[Dict[(str, float)]], antigen_embedding: Optional[np.ndarray])` → `AnomalyResult`

Assess danger level using context signals (danger theory).

Combines negative selection with contextual danger signals.

Args:
    antigen: Antigen to assess
    context_signals: Dict of signal_type -> strength (0-1)
    antigen_embedding: Optional embedding

Returns:
    AnomalyResult with danger assessment

##### `zero_shot_detect(`self, antigen: Antigen, description: str)` → `AnomalyResult`

Zero-shot anomaly detection using LLM reasoning.

For cases where embeddings aren't available or additional reasoning is needed.

Args:
    antigen: Antigen to assess
    description: Text description for LLM analysis

Returns:
    AnomalyResult from LLM reasoning

##### `get_statistics(`self)` → `Dict[(str, Any)]`

Get agent statistics

##### `to_agent_response(`self, result: AnomalyResult, success: bool, execution_time: float)` → `AgentResponse`

Convert anomaly result to agent response

---

## Links

- [[../code-mirror/src/agents/nk_cell_agent.py|Source Code]]
- [[../../projects/immunos-mcp/|Project Root]]
- [[./README|API Index]]

