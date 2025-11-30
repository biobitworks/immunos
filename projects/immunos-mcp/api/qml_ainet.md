---
project: immunos-mcp
type: api-documentation
source: src/algorithms/qml_ainet.py
generated: 2025-11-30
tags: [api, documentation, immunos-mcp]
---

# qml_ainet

**Source**: [[../code-mirror/src/algorithms/qml_ainet.py|src/algorithms/qml_ainet.py]]

QML-AiNet: Qualitative Model Learning using Artificial Immune Network

Implementation based on:
- Pang & Coghill (2015) - "QML-AiNet: an immune network approach to learning
  qualitative differential equation models"

Key differences from Opt-AiNet:
1. Modified mutation operator: uniform random selection from constraint space
2. Integer encoding for constraint selection
3. Fitness based on model-data consistency
4. Application to qualitative differential equation learning

## Contents

### Classes
- [QualitativeValue](#qualitativevalue)
- [QualitativeChange](#qualitativechange)
- [QualitativeRelation](#qualitativerelation)
- [QualitativeConstraint](#qualitativeconstraint)
- [QDEModel](#qdemodel)
- [QMLAiNet](#qmlainet)

## Classes

### QualitativeValue

**Inherits from**: `Enum`

**Source**: `src/algorithms/qml_ainet.py:32`

Qualitative values for variables.

---

### QualitativeChange

**Inherits from**: `Enum`

**Source**: `src/algorithms/qml_ainet.py:42`

Qualitative changes (derivatives).

---

### QualitativeRelation

**Inherits from**: `Enum`

**Source**: `src/algorithms/qml_ainet.py:51`

Relations between variables.

---

### QualitativeConstraint

**Source**: `src/algorithms/qml_ainet.py:61`

A qualitative constraint in a QDE model.

Examples:
- d(requests)/dt = INCREASING
- latency = HIGH when requests > threshold
- d(errors)/dt ∝+ requests  (errors increase with requests)

#### Methods

##### `__str__(`self)`

*No documentation available.*

##### `__repr__(`self)`

*No documentation available.*

---

### QDEModel

**Source**: `src/algorithms/qml_ainet.py:85`

Qualitative Differential Equation Model.

Represents system behavior using qualitative constraints.

#### Methods

##### `__str__(`self)`

*No documentation available.*

##### `satisfies(`self, observation: Dict[(str, Any)])` → `bool`

Check if observation satisfies this QDE model.

Args:
    observation: Dict mapping variable names to observed values/changes

Returns:
    True if observation is consistent with model

---

### QMLAiNet

**Inherits from**: `OptAiNet`

**Source**: `src/algorithms/qml_ainet.py:133`

QML-AiNet: Qualitative Model Learning using Immune Network.

Extends Opt-AiNet with:
- Integer encoding for constraint selection
- Modified mutation operator (uniform random selection)
- Qualitative model representation
- Model-data consistency evaluation

#### Methods

##### `__init__(`self, constraint_space: List[List[QualitativeConstraint]], observations: List[Dict[(str, Any)]], variables: List[str], population_size: int, clone_multiplier: int, affinity_threshold: float, mutation_beta: float, num_random: int, max_generations: int)`

Initialize QML-AiNet.

Args:
    constraint_space: List of constraint options for each "slot"
                    (position in model). Each slot contains possible
                    constraints for that position.
    observations: Training data - list of observation dicts
    variables: Variable names in the system
    population_size: Initial population size
    clone_multiplier: Clones per antibody
    affinity_threshold: Network suppression threshold
    mutation_beta: Mutation decay parameter
    num_random: Random antibodies per generation
    max_generations: Max iterations

##### `initialize_population(`self)`

Initialize population with random constraint selections.

##### `_evaluate_qde_model(`self, position: np.ndarray)` → `float`

Evaluate QDE model fitness.

Fitness = consistency with observations.

Args:
    position: Integer array selecting constraints

Returns:
    Fitness score (higher = more consistent with data)

##### `_mutate(`self, position: np.ndarray, mutation_rate: float)` → `np.ndarray`

QML-AiNet modified mutation operator.

Key difference from Opt-AiNet: Uses UNIFORM RANDOM selection
from constraint space instead of Gaussian noise.

This is because adjacent constraints don't have neighborhood
relationships in qualitative constraint space.

Args:
    position: Original constraint selection (integer array)
    mutation_rate: Probability of mutating each position

Returns:
    Mutated position

##### `_evaluate(`self, position: np.ndarray)` → `float`

Override base evaluation to use QDE model fitness.

##### `decode_model(`self, position: np.ndarray)` → `QDEModel`

Decode antibody position to QDE model.

Args:
    position: Integer array selecting constraints

Returns:
    QDEModel object

##### `get_best_model(`self)` → `QDEModel`

Get the best QDE model found.

Returns:
    Best QDEModel

##### `get_top_models(`self, k: int)` → `List[QDEModel]`

Get top k QDE models.

Args:
    k: Number of models to return

Returns:
    List of top QDEModel objects

---

## Links

- [[../code-mirror/src/algorithms/qml_ainet.py|Source Code]]
- [[../../projects/immunos-mcp/|Project Root]]
- [[./README|API Index]]

