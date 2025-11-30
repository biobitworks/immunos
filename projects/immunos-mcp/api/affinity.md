---
project: immunos-mcp
type: api-documentation
source: src/core/affinity.py
generated: 2025-11-30
tags: [api, documentation, immunos-mcp]
---

# affinity

**Source**: [[../code-mirror/src/core/affinity.py|src/core/affinity.py]]

Affinity calculation for IMMUNOS-MCP

Combines traditional Immunos-81 affinity calculations with
LLM embedding-based semantic similarity.

## Contents

### Classes
- [AffinityResult](#affinityresult)
- [AffinityCalculator](#affinitycalculator)
- [DistanceMetric](#distancemetric)

## Classes

### AffinityResult

**Source**: `src/core/affinity.py:15`

Result of affinity calculation

#### Methods

##### `__float__(`self)` → `float`

*No documentation available.*

---

### AffinityCalculator

**Source**: `src/core/affinity.py:26`

Calculate affinity (similarity) between antigens and patterns.

Supports multiple calculation methods:
- Traditional: Exponential decay for numeric, exact match for nominal
- Embedding: Cosine similarity of LLM embeddings
- Hybrid: Weighted combination of traditional + embedding

#### Methods

##### `__init__(`self, method: str, embedding_weight: float)`

Initialize affinity calculator.

Args:
    method: "traditional", "embedding", or "hybrid"
    embedding_weight: Weight for embedding similarity (0-1) in hybrid mode

##### `calculate(`self, antigen1: Any, antigen2: Any, embeddings1: Optional[np.ndarray], embeddings2: Optional[np.ndarray])` → `AffinityResult`

Calculate affinity between two antigens.

Args:
    antigen1: First antigen
    antigen2: Second antigen
    embeddings1: Optional embedding vector for antigen1
    embeddings2: Optional embedding vector for antigen2

Returns:
    AffinityResult with score and metadata

##### `_traditional_affinity(`self, val1: Any, val2: Any)` → `AffinityResult`

Traditional Immunos-81 affinity calculation.

For numeric: exp(-|v1 - v2| / sigma)
For nominal: 1.0 if equal, 0.0 otherwise
For text: Character-level distance

##### `_embedding_affinity(`self, emb1: np.ndarray, emb2: np.ndarray)` → `AffinityResult`

Embedding-based affinity using cosine similarity.

Args:
    emb1: Embedding vector 1
    emb2: Embedding vector 2

Returns:
    AffinityResult with cosine similarity score

##### `calculate_avidity(`self, affinities: List[float], clone_size: int)` → `float`

Calculate avidity (collective affinity) for a clone.

Uses Immunos-81 formula: sum(affinities) * log(1 + clone_size)

Args:
    affinities: List of affinity scores
    clone_size: Number of members in clone

Returns:
    Avidity score

---

### DistanceMetric

**Source**: `src/core/affinity.py:186`

Distance metrics for pattern matching

#### Methods

##### `euclidean(`vec1: np.ndarray, vec2: np.ndarray)` → `float`

Euclidean distance

##### `manhattan(`vec1: np.ndarray, vec2: np.ndarray)` → `float`

Manhattan (L1) distance

##### `cosine_distance(`vec1: np.ndarray, vec2: np.ndarray)` → `float`

Cosine distance (1 - cosine similarity)

##### `hamming(`vec1: np.ndarray, vec2: np.ndarray)` → `float`

Hamming distance (for binary vectors)

---

## Links

- [[../code-mirror/src/core/affinity.py|Source Code]]
- [[../../projects/immunos-mcp/|Project Root]]
- [[./README|API Index]]

