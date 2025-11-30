---
project: immunos-mcp
type: algorithm-diagram
tags: [diagram, mermaid, qml-ainet, algorithm]
created: 2025-11-30
---

# QML-AiNet Algorithm Flow

Detailed algorithm flow for Qualitative Model Learning using Artificial Immune Networks.

## Algorithm Overview

```mermaid
flowchart TD
    Start([Start: Training Data + Constraint Space]) --> Init

    Init[Initialize Population<br/>Random constraint selections]
    Init --> Eval1[Evaluate All Antibodies<br/>Model-data consistency]

    Eval1 --> Sort1[Sort by Fitness<br/>Best to worst]
    Sort1 --> Select[Select Top n Antibodies]

    Select --> Clone[Clonal Expansion<br/>Fitness-proportional cloning]
    Clone --> Mutate[Modified Mutation<br/>Uniform random selection]

    Mutate --> Eval2[Evaluate Clones]
    Eval2 --> Combine[Combine Population + Clones]

    Combine --> Suppress[Network Suppression<br/>Remove similar solutions]
    Suppress --> Random[Add d Random Antibodies<br/>Maintain diversity]

    Random --> Update[Update Best Solution]
    Update --> Check{Converged?}

    Check -->|No| Eval1
    Check -->|Yes| Decode[Decode Best Position<br/>to QDE Model]

    Decode --> End([Return: QDE Model])

    style Init fill:#e1f5fe
    style Mutate fill:#fff3e0
    style Suppress fill:#fce4ec
    style Decode fill:#e8f5e9
```

## Key Differences from Opt-AiNet

### 1. Encoding

```mermaid
graph LR
    subgraph "Opt-AiNet (Continuous)"
        O1[Position: Real Numbers]
        O2["[0.5, -1.2, 0.8]"]
        O1 --> O2
    end

    subgraph "QML-AiNet (Discrete)"
        Q1[Position: Constraint Indices]
        Q2["[2, 0, 1]"]
        Q3[Maps to Constraints]
        Q1 --> Q2 --> Q3
    end

    style O2 fill:#e3f2fd
    style Q2 fill:#fff3e0
    style Q3 fill:#fff3e0
```

### 2. Mutation Operator

**Opt-AiNet (Gaussian)**:
```python
# Continuous space mutation
mutation_rate = alpha * e^(-beta * fitness)
mutated[i] = position[i] + N(0, mutation_rate)
```

**QML-AiNet (Uniform Random)**:
```python
# Discrete space mutation
for i in range(len(position)):
    if random() < mutation_rate:
        position[i] = random_choice(constraint_space[i])
```

### 3. Fitness Function

**Opt-AiNet**:
- Numeric objective function (e.g., minimize f(x))

**QML-AiNet**:
- Model-data consistency
- Fraction of observations satisfied by model

```mermaid
graph LR
    Model[QDE Model] --> Obs1[Observation 1]
    Model --> Obs2[Observation 2]
    Model --> Obs3[Observation 3]

    Obs1 -->|Satisfied| S1[✓]
    Obs2 -->|Satisfied| S2[✓]
    Obs3 -->|Not Satisfied| S3[✗]

    S1 --> Fitness["Fitness = 2/3 = 0.67"]
    S2 --> Fitness
    S3 --> Fitness

    style S1 fill:#c8e6c9
    style S2 fill:#c8e6c9
    style S3 fill:#ffcdd2
```

## Network Suppression

Eliminates antibodies that are too similar to higher-fitness solutions.

```mermaid
graph TD
    Pop[Population<br/>Sorted by Fitness] --> A1

    A1{For each antibody}
    A1 --> Comp[Compare to kept antibodies]

    Comp --> Dist[Calculate distance]
    Dist --> Thr{Distance < threshold?}

    Thr -->|Yes| Suppress[Suppress<br/>Too similar]
    Thr -->|No| Keep[Keep<br/>Sufficiently different]

    Suppress --> A1
    Keep --> Add[Add to suppressed population]
    Add --> A1

    A1 -->|Done| Result[Return: Diverse population]

    style Suppress fill:#ffcdd2
    style Keep fill:#c8e6c9
```

**Purpose**: Maintains diversity, prevents premature convergence

**Example**:
```
Original Population (10 antibodies):
  [0,1,2] fitness=0.9  ← Best
  [0,1,2] fitness=0.85 ← Suppressed (too similar to best)
  [0,1,3] fitness=0.8  ← Kept (different)
  [1,2,0] fitness=0.7  ← Kept (different)
  ...

Suppressed Population (6 antibodies):
  Only diverse solutions kept
```

## Constraint Space Example

```mermaid
graph TD
    subgraph "Slot 0: Request Trend"
        S0_0[STEADY]
        S0_1[INCREASING]
        S0_2[RAPID_INCREASE]
    end

    subgraph "Slot 1: Latency Relation"
        S1_0["latency ∝+ requests"]
        S1_1["latency = LOW"]
        S1_2["latency = HIGH"]
    end

    subgraph "Slot 2: Error Relation"
        S2_0["errors ∝+ latency"]
        S2_1["errors = ZERO"]
        S2_2["errors ∝+ requests"]
    end

    Position["Position [1, 0, 2]"] --> S0_1
    Position --> S1_0
    Position --> S2_2

    S0_1 --> Model
    S1_0 --> Model
    S2_2 --> Model

    Model[QDE Model] --> Eval[Evaluate against<br/>observations]
    Eval --> Fit[Fitness Score]

    style Position fill:#fff3e0
    style Model fill:#e8f5e9
    style Fit fill:#c8e6c9
```

**Search Space Size**: 3 × 3 × 3 = 27 combinations

## Evolution Process

```mermaid
graph TB
    subgraph "Generation 1"
        G1_1["[0,0,0] f=0.2"]
        G1_2["[1,1,1] f=0.5"]
        G1_3["[2,2,2] f=0.3"]
        G1_Best["Best: [1,1,1]"]

        G1_1 -.-> G1_Best
        G1_2 -.-> G1_Best
        G1_3 -.-> G1_Best
    end

    subgraph "Cloning & Mutation"
        C1["Clone [1,1,1] × 5"]
        M1["Mutate → [1,0,1]"]
        M2["Mutate → [1,1,2]"]
        M3["Mutate → [0,1,1]"]
    end

    subgraph "Generation 2"
        G2_1["[1,0,1] f=0.7"]
        G2_2["[1,1,2] f=0.9"]
        G2_3["[0,1,1] f=0.6"]
        G2_Best["Best: [1,1,2]"]

        G2_1 -.-> G2_Best
        G2_2 -.-> G2_Best
        G2_3 -.-> G2_Best
    end

    G1_Best --> C1
    C1 --> M1
    C1 --> M2
    C1 --> M3
    M1 --> G2_1
    M2 --> G2_2
    M3 --> G2_3

    style G1_Best fill:#fff3e0
    style G2_Best fill:#c8e6c9
```

Fitness improves: 0.5 → 0.7 → 0.9

## Complete QML-AiNet Pseudocode

```python
def qml_ainet(constraint_space, observations, parameters):
    # Initialize
    population = initialize_random(constraint_space, parameters.pop_size)
    evaluate_fitness(population, observations)

    for generation in range(parameters.max_gen):
        # Sort by fitness
        population.sort(key=lambda ab: ab.fitness, reverse=True)

        # Clone best antibodies
        clones = []
        for antibody in population[:parameters.n_select]:
            n_clones = round(parameters.clone_mult * antibody.fitness)
            for _ in range(n_clones):
                clone = antibody.clone()
                mutation_rate = exp(-parameters.beta * antibody.fitness)
                mutate_uniform(clone, constraint_space, mutation_rate)
                clones.append(clone)

        # Evaluate clones
        evaluate_fitness(clones, observations)

        # Network suppression
        combined = population + clones
        population = network_suppression(combined, parameters.threshold)

        # Add diversity
        for _ in range(parameters.d_random):
            population.append(initialize_random(constraint_space, 1)[0])

        # Update best
        best = max(population, key=lambda ab: ab.fitness)
        if converged(best, previous_best):
            break

    # Decode best solution
    model = decode_to_qde_model(best.position, constraint_space)
    return model
```

## Example Run

### Input
```python
variables = ["requests", "latency", "errors"]

constraint_space = [
    [STEADY, INCREASING, RAPID_INCREASE],
    [latency∝+requests, latency=LOW, latency=HIGH],
    [errors∝+latency, errors=ZERO, errors∝+requests]
]

observations = [
    {requests: STEADY, latency: LOW, errors: ZERO},
    {requests: RAPID_INCREASE, latency: HIGH, errors: HIGH}
]
```

### Evolution
```
Gen 0: Best=[0,1,1] Fitness=0.50 (1/2 observations)
Gen 5: Best=[1,0,0] Fitness=0.50
Gen 10: Best=[1,0,2] Fitness=1.00 (2/2 observations) ← FOUND!
```

### Output Model
```
QDE Model:
  Variables: requests, latency, errors
  Constraints:
    - d(requests)/dt = INCREASING
    - latency ∝+ requests
    - errors ∝+ requests
  Fitness: 1.0000
```

## Links

- [[../code-mirror/src/algorithms/qml_ainet|QML-AiNet Implementation]]
- [[../code-mirror/src/algorithms/opt_ainet|Opt-AiNet Base Class]]
- [[../../research/experiments/qml-ainet-validation-2025-11-30|Validation Results]]
- [[../../research/literature/pang-coghill-2015-qml-ainet|Original Paper]]

---

*Generated: 2025-11-30*
