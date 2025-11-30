---
project: immunos-mcp
type: api-documentation
source: src/core/protocols.py
generated: 2025-11-30
tags: [api, documentation, immunos-mcp]
---

# protocols

**Source**: [[../code-mirror/src/core/protocols.py|src/core/protocols.py]]

Protocol definitions for IMMUNOS-MCP

Defines shared data structures and communication protocols
between immune system agents.

## Contents

### Classes
- [RecognitionStrategy](#recognitionstrategy)
- [SignalType](#signaltype)
- [RecognitionResult](#recognitionresult)
- [AnomalyResult](#anomalyresult)
- [Signal](#signal)
- [ExtractedFeatures](#extractedfeatures)
- [AgentRequest](#agentrequest)
- [AgentResponse](#agentresponse)
- [Pattern](#pattern)
- [Detector](#detector)
- [Clone](#clone)

## Classes

### RecognitionStrategy

**Inherits from**: `Enum`

**Source**: `src/core/protocols.py:13`

Recognition strategies for classification

---

### SignalType

**Inherits from**: `Enum`

**Source**: `src/core/protocols.py:21`

Signal types for dendritic cell processing

---

### RecognitionResult

**Source**: `src/core/protocols.py:30`

Result from pattern recognition

#### Methods

##### `to_dict(`self)` → `Dict[(str, Any)]`

Convert to dictionary for JSON serialization

---

### AnomalyResult

**Source**: `src/core/protocols.py:54`

Result from anomaly detection

#### Methods

##### `to_dict(`self)` → `Dict[(str, Any)]`

Convert to dictionary

---

### Signal

**Source**: `src/core/protocols.py:78`

Signal processed by dendritic cells

#### Methods

##### `to_dict(`self)` → `Dict[(str, Any)]`

Convert to dictionary

---

### ExtractedFeatures

**Source**: `src/core/protocols.py:98`

Features extracted by dendritic cells

#### Methods

##### `get_overall_signal(`self)` → `float`

Calculate overall signal strength

##### `to_dict(`self)` → `Dict[(str, Any)]`

Convert to dictionary

---

### AgentRequest

**Source**: `src/core/protocols.py:132`

Request sent to an agent

#### Methods

##### `to_dict(`self)` → `Dict[(str, Any)]`

Convert to dictionary

---

### AgentResponse

**Source**: `src/core/protocols.py:150`

Response from an agent

#### Methods

##### `to_dict(`self)` → `Dict[(str, Any)]`

Convert to dictionary

---

### Pattern

**Source**: `src/core/protocols.py:174`

Learned pattern (like B cell antibody)

#### Methods

##### `to_dict(`self)` → `Dict[(str, Any)]`

Convert to dictionary

---

### Detector

**Source**: `src/core/protocols.py:202`

Detector for negative selection (NK cell)

#### Methods

##### `matches_self(`self, antigen: Any)` → `bool`

Check if antigen matches self patterns (returns True if normal)

##### `is_non_self(`self, antigen: Any)` → `bool`

Check if antigen is non-self (anomaly)

##### `to_dict(`self)` → `Dict[(str, Any)]`

Convert to dictionary

---

### Clone

**Source**: `src/core/protocols.py:229`

Clone population (group of B cells with similar patterns)

#### Methods

##### `__post_init__(`self)`

*No documentation available.*

##### `calculate_avidity(`self, affinity_scores: List[float])` → `float`

Calculate avidity using Immunos-81 formula

##### `to_dict(`self)` → `Dict[(str, Any)]`

Convert to dictionary

---

## Links

- [[../code-mirror/src/core/protocols.py|Source Code]]
- [[../../projects/immunos-mcp/|Project Root]]
- [[./README|API Index]]

