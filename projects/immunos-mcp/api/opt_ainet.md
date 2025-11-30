---
project: immunos-mcp
type: api-documentation
source: src/algorithms/opt_ainet.py
generated: 2025-11-30
tags: [api, documentation, immunos-mcp]
---

# opt_ainet

**Source**: [[../code-mirror/src/algorithms/opt_ainet.py|src/algorithms/opt_ainet.py]]

Opt-AiNet: Optimization Artificial Immune Network

Implementation based on:
- de Castro & Timmis (2002) - "An artificial immune network for multimodal function optimization"
- Clever Algorithms documentation

Core principles:
- Network suppression for diversity maintenance
- Clonal selection with fitness-based mutation
- Multi-modal optimization through network interactions

## Contents

### Classes
- [Antibody](#antibody)
- [OptAiNet](#optainet)

## Classes

### Antibody

**Source**: `src/algorithms/opt_ainet.py:21`

Antibody (candidate solution) in the immune network.

Attributes:
    position: Solution vector in search space
    fitness: Fitness score (higher is better)
    age: Generation count since creation
    clone_count: Number of clones to generate

#### Methods

##### `__hash__(`self)`

*No documentation available.*

---

### OptAiNet

**Source**: `src/algorithms/opt_ainet.py:40`

Optimization Artificial Immune Network (Opt-AiNet).

Multi-modal optimization using immune network principles:
- Clonal selection for exploitation
- Network suppression for diversity
- Random insertion for exploration

#### Methods

##### `__init__(`self, fitness_function: Callable[([np.ndarray], float)], dimensions: int, bounds: Tuple[(float, float)], population_size: int, clone_multiplier: int, affinity_threshold: float, mutation_beta: float, num_random: int, max_generations: int, minimize: bool)`

Initialize Opt-AiNet.

Args:
    fitness_function: Function to maximize (or minimize if minimize=True)
    dimensions: Dimensionality of search space
    bounds: (min, max) bounds for each dimension
    population_size: Initial population size
    clone_multiplier: Number of clones per antibody
    affinity_threshold: Suppression threshold (as fraction of search space)
    mutation_beta: Controls mutation decay (higher = less mutation)
    num_random: Number of random antibodies to insert each generation
    max_generations: Maximum iterations
    minimize: If True, minimizes fitness; if False, maximizes

##### `initialize_population(`self)`

Initialize population with random antibodies.

##### `_evaluate(`self, position: np.ndarray)` → `float`

Evaluate fitness (handle minimization).

##### `_update_best(`self)`

Update best antibody found so far.

##### `clone_and_mutate(`self, antibody: Antibody)` → `List[Antibody]`

Clone antibody and mutate clones.

Better antibodies (higher fitness) are:
- Cloned more
- Mutated less

Args:
    antibody: Parent antibody to clone

Returns:
    List of mutated clones

##### `_mutate(`self, position: np.ndarray, mutation_rate: float)` → `np.ndarray`

Mutate position with Gaussian noise.

Args:
    position: Original position
    mutation_rate: Mutation strength (0-1)

Returns:
    Mutated position (clamped to bounds)

##### `network_suppression(`self, population: List[Antibody])` → `List[Antibody]`

Network suppression: eliminate similar, lower-fitness antibodies.

This maintains diversity by preventing redundant solutions.

Args:
    population: Current population

Returns:
    Suppressed population (without redundant antibodies)

##### `add_random_antibodies(`self, num_random: int)` → `List[Antibody]`

Add random antibodies for exploration.

Args:
    num_random: Number of random antibodies to generate

Returns:
    List of random antibodies

##### `optimize(`self, verbose: bool)` → `Tuple[(np.ndarray, float)]`

Run Opt-AiNet optimization.

Args:
    verbose: Print progress

Returns:
    Tuple of (best_position, best_fitness)

##### `get_all_optima(`self, top_k: int)` → `List[Tuple[(np.ndarray, float)]]`

Get multiple optimal solutions (for multi-modal optimization).

Args:
    top_k: Number of top solutions to return

Returns:
    List of (position, fitness) tuples

---

## Links

- [[../code-mirror/src/algorithms/opt_ainet.py|Source Code]]
- [[../../projects/immunos-mcp/|Project Root]]
- [[./README|API Index]]

