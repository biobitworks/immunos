---
project: immunos-mcp
type: api-documentation
source: src/agents/bcell_agent.py
generated: 2025-11-30
tags: [api, documentation, immunos-mcp]
---

# bcell_agent

**Source**: [[../code-mirror/src/agents/bcell_agent.py|src/agents/bcell_agent.py]]

B Cell Agent for IMMUNOS-MCP

Pattern matching agent inspired by B cells and antibody production.
Learns specific patterns and calculates affinity to new antigens.

## Contents

### Classes
- [BCellAgent](#bcellagent)

## Classes

### BCellAgent

**Source**: `src/agents/bcell_agent.py:17`

B Cell Agent - Pattern Matcher

Mimics B cell behavior:
- Learns patterns from training examples (antibody production)
- Calculates affinity to new antigens (binding strength)
- Groups into clones for collective recognition (plasma cells)
- Supports affinity maturation (improving pattern matching)

#### Methods

##### `__init__(`self, agent_name: str, affinity_method: str)`

Initialize B Cell agent.

Args:
    agent_name: Unique identifier for this agent
    affinity_method: "traditional", "embedding", or "hybrid"

##### `train(`self, training_data: List[Antigen], embeddings: Optional[List[np.ndarray]])`

Train on antigens (create patterns/antibodies).

Args:
    training_data: List of training antigens with labels
    embeddings: Optional pre-computed embeddings

##### `_create_clones(`self)`

Group patterns into clones by class label

##### `calculate_affinity(`self, antigen: Antigen, pattern: Pattern, antigen_embedding: Optional[np.ndarray])` → `float`

Calculate affinity between antigen and pattern.

Args:
    antigen: Test antigen
    pattern: Learned pattern
    antigen_embedding: Optional embedding for antigen

Returns:
    Affinity score (0-1)

##### `recognize(`self, antigen: Antigen, antigen_embedding: Optional[np.ndarray], strategy: str)` → `RecognitionResult`

Recognize antigen using learned patterns.

Args:
    antigen: Antigen to classify
    antigen_embedding: Optional embedding
    strategy: "sha" (Simple Highest Avidity) or "rha" (Relative Highest Avidity)

Returns:
    RecognitionResult with prediction and confidence

##### `_apply_sha(`self, avidity_scores: Dict[(str, float)])` → `RecognitionResult`

Simple Highest Avidity strategy

##### `_apply_rha(`self, avidity_scores: Dict[(str, float)], threshold: float)` → `RecognitionResult`

Relative Highest Avidity strategy with 5% threshold

##### `add_pattern(`self, antigen: Antigen, embedding: Optional[np.ndarray])`

Add new pattern online (incremental learning).

Args:
    antigen: New training example
    embedding: Optional embedding

##### `mature_affinity(`self, antigen: Antigen, embedding: Optional[np.ndarray], iterations: int)` → `Pattern`

Affinity maturation: improve pattern matching through mutation and selection.

Args:
    antigen: Target antigen
    embedding: Optional embedding
    iterations: Number of maturation iterations

Returns:
    Best matured pattern

##### `get_statistics(`self)` → `Dict[(str, Any)]`

Get agent statistics

##### `to_agent_response(`self, result: RecognitionResult, success: bool, execution_time: float)` → `AgentResponse`

Convert recognition result to agent response

---

## Links

- [[../code-mirror/src/agents/bcell_agent.py|Source Code]]
- [[../../projects/immunos-mcp/|Project Root]]
- [[./README|API Index]]

