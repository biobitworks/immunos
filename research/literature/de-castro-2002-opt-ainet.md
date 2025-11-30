---
title: "An Artificial Immune Network for Multimodal Function Optimization"
authors: [Leandro N. de Castro, Jonathan Timmis]
year: 2002
venue: IEEE Congress on Evolutionary Computation
type: conference-paper
tags: [opt-ainet, artificial-immune-systems, optimization, clonal-selection, network-suppression]
relevance: foundational
status: referenced
project: immunos-mcp
---

# Opt-AiNet: Artificial Immune Network for Optimization

**Full Citation**: de Castro, L. N., & Timmis, J. (2002). An Artificial Immune Network for Multimodal Function Optimization. *IEEE Congress on Evolutionary Computation (CEC)*, 1, 699-704.

## Summary

This foundational paper introduces **Opt-AiNet**, an artificial immune network algorithm for continuous multimodal optimization. It combines clonal selection principles with network suppression to find multiple optima simultaneously.

## Key Contributions

1. **Network Suppression Mechanism**
   - Eliminates similar solutions to maintain diversity
   - Prevents population collapse to single optimum
   - Enables multimodal optimization (finding multiple peaks)

2. **Adaptive Clonal Expansion**
   - Higher fitness → more clones
   - Fitness-proportional mutation rate
   - Balances exploitation vs exploration

3. **Continuous Optimization Framework**
   - Real-valued position vectors
   - Gaussian mutation operator
   - Affinity-based distance metric

## Algorithm Overview

```
INPUT: Objective function f(x), search bounds
OUTPUT: Population of high-fitness solutions (multiple optima)

1. Initialize random population
2. Evaluate fitness for all antibodies
3. Select top n antibodies
4. Clone proportionally to fitness
5. Mutate clones (Gaussian, fitness-dependent rate)
6. Evaluate clone fitness
7. Network suppression (remove similar solutions)
8. Add d random antibodies
9. Repeat until convergence
```

## Key Mechanisms

### Clonal Expansion
```
n_clones = round(β × N × fitness_normalized)

Where:
  β = clone multiplier (typically 10)
  N = population size
  fitness_normalized ∈ [0, 1]
```

### Gaussian Mutation
```
mutation_rate = α × e^(-ρ × fitness_normalized)
x_mutated = x + N(0, mutation_rate)

Where:
  α = decay parameter (controls exploration)
  ρ = controls exploitation-exploration tradeoff
```

### Network Suppression
```
For each antibody (sorted by fitness):
  If distance to any kept antibody < threshold:
    Suppress (remove)
  Else:
    Keep

Result: Diverse population covering multiple optima
```

## Comparison with Other Algorithms

| Algorithm | Multi-Modal | Diversity | Convergence Speed |
|-----------|-------------|-----------|-------------------|
| **Opt-AiNet** | ✅ Yes | High (suppression) | Medium |
| Genetic Algorithm | ⚠️ Depends | Medium (niching) | Fast |
| CLONALG | ❌ No | Low | Medium |
| Particle Swarm | ⚠️ Depends | Medium | Fast |

**Advantage**: Native multimodal optimization without niching hacks

## Results from Paper

**Benchmark Functions**:
1. Schaffer F6 (2D, multiple peaks)
2. Shubert (2D, 760 local optima)
3. Rastrigin (10D, highly multimodal)

**Performance**:
- Successfully found multiple global optima
- Maintained diverse population (network suppression)
- Comparable accuracy to genetic algorithms
- 10-100x slower than single-optimum methods (acceptable tradeoff)

## Relevance to QML-AiNet

**QML-AiNet** (Pang & Coghill 2015) adapts Opt-AiNet for discrete spaces:

| Aspect | Opt-AiNet (2002) | QML-AiNet (2015) |
|--------|------------------|------------------|
| **Space** | Continuous | Discrete |
| **Mutation** | Gaussian | Uniform random |
| **Fitness** | f(x) objective | Model-data consistency |
| **Encoding** | Real vectors | Integer indices |
| **Suppression** | Same | Same (inherited) |

**Key Insight**: Network suppression is domain-agnostic and works for both continuous and discrete optimization.

## Our Implementation

**Status**: ✅ Base class for QML-AiNet

**Code**: [[../../projects/immunos-mcp/code-mirror/src/algorithms/opt_ainet|opt_ainet.py]]
- 550 lines
- Fully documented
- Parameterized for extension

**Usage**:
```python
from algorithms.opt_ainet import OptAiNet

# Continuous optimization
opt = OptAiNet(
    pop_size=30,
    clone_mult=10,
    affinity_threshold=0.1
)
best_solutions = opt.optimize(objective_fn, bounds)
```

## Critical Insights

### Why Network Suppression?

**Problem**: Without suppression, population converges to single best optimum
**Solution**: Suppress similar antibodies → maintain diversity → find multiple optima

**Analogy**: Immune system maintains repertoire of different antibodies for different pathogens (not just one super-antibody)

### Fitness-Proportional Cloning

**High fitness** → More clones → More mutations → Exploit promising regions
**Low fitness** → Few clones → Less computation → Efficient resource use

**Result**: Natural balance between exploration and exploitation

## Limitations

1. **Slow convergence** (more iterations than GA/PSO)
2. **Parameter sensitivity** (α, ρ, β, threshold)
3. **Continuous spaces only** (fixed by QML-AiNet)
4. **No theoretical convergence guarantee**

## Applications

### From Paper
- Function optimization (benchmarks)
- Engineering design
- Parameter tuning

### Our Extension
- **Code behavior analysis** (via QML-AiNet)
- Pattern matching with clonal selection
- Multi-modal threat detection

## Future Work

### From Paper
- Adaptive parameter tuning
- Hybrid with other evolutionary algorithms
- Theoretical convergence analysis

### Our Extensions
- [ ] Apply network suppression to B Cell agent
- [ ] Multi-modal malware detection (different attack vectors)
- [ ] Adaptive population sizing based on threat level

## Related Papers

- [[pang-coghill-2015-qml-ainet|Pang & Coghill (2015)]] - QML-AiNet (our implementation)
- [[de-castro-2001-clonal-selection|de Castro & Von Zuben (2002)]] - CLONALG (predecessor)
- [[dasgupta-1999-immune-algorithms|Dasgupta (1999)]] - Immune algorithm survey
- [[timmis-2000-ais-review|Timmis (2000)]] - AIS theoretical foundations

## Links

- [[../../projects/immunos-mcp/journal/2025-11-30|Implementation Journal]]
- [[../../projects/immunos-mcp/diagrams/system-overview|System Architecture]]
- [[../experiments/qml-ainet-validation-2025-11-30|Validation Results]]
- [[../ideas/network-suppression-integration|Network Suppression for B Cells]]

## Tags for Search

#opt-ainet #artificial-immune-systems #multimodal-optimization #clonal-selection #network-suppression #continuous-optimization #evolutionary-computation

---

**Relevance**: Foundation for QML-AiNet implementation
**Implementation Status**: ✅ Base class completed
**Paper Availability**: IEEE Xplore
**Citations**: 1000+ (highly influential)
