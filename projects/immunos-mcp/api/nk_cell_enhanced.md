---
project: immunos-mcp
type: api-documentation
source: src/agents/nk_cell_enhanced.py
generated: 2025-11-30
tags: [api, documentation, immunos-mcp]
---

# nk_cell_enhanced

**Source**: [[../code-mirror/src/agents/nk_cell_enhanced.py|src/agents/nk_cell_enhanced.py]]

Enhanced NK Cell Agent for IMMUNOS-MCP

Implements NegSl-AIS methodology for improved anomaly detection:
- Adaptive threshold calculation (min_distance, mean_distance, percentile)
- Per-class detector generation (one-vs-rest strategy)
- Enhanced validation rule: d(detector, nearest_self) > τ
- Quality metrics tracking

## Contents

### Classes
- [ThresholdMethod](#thresholdmethod)
- [ClassDetectorSet](#classdetectorset)
- [EnhancedNKCellAgent](#enhancednkcellagent)

## Classes

### ThresholdMethod

**Inherits from**: `Enum`

**Source**: `src/agents/nk_cell_enhanced.py:22`

Methods for calculating adaptive thresholds

---

### ClassDetectorSet

**Source**: `src/agents/nk_cell_enhanced.py:29`

Detectors for a specific class (one-vs-rest)

#### Methods

##### `__init__(`self, class_label: str, num_detectors: int)`

*No documentation available.*

##### `add_detector(`self, detector: np.ndarray, quality_score: float)`

Add detector with quality score

##### `get_statistics(`self)` → `Dict[(str, Any)]`

Get detector statistics

---

### EnhancedNKCellAgent

**Source**: `src/agents/nk_cell_enhanced.py:56`

Enhanced NK Cell Agent with NegSl-AIS methodology

Improvements over base NK Cell:
- Adaptive threshold calculation based on self-pattern distribution
- Per-class detector generation (one-vs-rest strategy)
- Enhanced validation rule: d(detector, nearest_self) > τ
- Quality metrics for detector effectiveness

#### Methods

##### `__init__(`self, agent_name: str, threshold_method: str, detectors_per_class: int, max_generation_attempts: int)`

Initialize Enhanced NK Cell agent.

Args:
    agent_name: Unique identifier
    threshold_method: "min_distance", "mean_distance", or "percentile"
    detectors_per_class: Number of detectors to generate per class (15-25 recommended)
    max_generation_attempts: Maximum attempts to generate each detector

##### `train_on_self(`self, normal_data: List[Antigen], embeddings: Optional[List[np.ndarray]], class_labels: Optional[List[str]])`

Train on normal ("self") data with optional class labels.

If class_labels provided, generates per-class detectors (one-vs-rest).
Otherwise, treats all data as single "self" class.

Args:
    normal_data: List of normal (non-anomalous) antigens
    embeddings: Optional embeddings for normal data
    class_labels: Optional class labels for one-vs-rest strategy

##### `_calculate_adaptive_threshold(`self, self_embeddings: List[np.ndarray])` → `float`

Calculate adaptive threshold based on self-pattern distribution.

Args:
    self_embeddings: Embeddings of self patterns

Returns:
    Adaptive threshold value

##### `_generate_class_detectors(`self, class_label: str, self_embeddings: List[np.ndarray], threshold: float)` → `ClassDetectorSet`

Generate detectors for a specific class using enhanced validation.

Implements NegSl-AIS validation rule:
d(detector, nearest_self) > τ

Args:
    class_label: Class label for this detector set
    self_embeddings: Embeddings of self patterns
    threshold: Adaptive threshold

Returns:
    ClassDetectorSet with generated detectors

##### `detect_novelty(`self, antigen: Antigen, antigen_embedding: Optional[np.ndarray], target_class: Optional[str])` → `AnomalyResult`

Detect if antigen is novel/anomalous using enhanced negative selection.

Args:
    antigen: Antigen to test
    antigen_embedding: Optional embedding
    target_class: Optional target class for one-vs-rest detection

Returns:
    AnomalyResult with detection outcome

##### `_generate_explanation(`self, is_anomaly: bool, anomaly_score: float, best_class: Optional[str], detector_matches: List[Tuple[(str, int)]])` → `str`

Generate human-readable explanation

##### `get_statistics(`self)` → `Dict[(str, Any)]`

Get comprehensive agent statistics

##### `to_agent_response(`self, result: AnomalyResult, success: bool, execution_time: float)` → `AgentResponse`

Convert anomaly result to agent response

---

## Links

- [[../code-mirror/src/agents/nk_cell_enhanced.py|Source Code]]
- [[../../projects/immunos-mcp/|Project Root]]
- [[./README|API Index]]

