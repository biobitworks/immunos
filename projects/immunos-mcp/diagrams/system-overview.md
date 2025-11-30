---
project: immunos-mcp
type: architecture-diagram
tags: [diagram, mermaid, system-overview, architecture]
created: 2025-11-30
---

# IMMUNOS-MCP System Overview

High-level architecture showing all components and their relationships.

## System Architecture

```mermaid
graph TB
    subgraph "Input"
        Code[Code Samples]
        Patterns[Known Patterns]
    end

    subgraph "Core Framework"
        Antigen[Antigen<br/>Data Representation]
        Affinity[Affinity Calculator<br/>Similarity Metrics]
    end

    subgraph "6 Immune Agents"
        Macro[Macrophage<br/>Triage]
        Dend[Dendritic Cell<br/>Feature Extraction]
        BCell[B Cell<br/>Pattern Matching]
        NKCell[NK Cell<br/>Anomaly Detection]
        QML[QML Agent<br/>Behavior Learning]
        TCell[T Cell<br/>Coordination]
    end

    subgraph "Algorithms"
        OptAiNet[Opt-AiNet<br/>Multi-modal Optimization]
        QMLAiNet[QML-AiNet<br/>Qualitative Learning]
    end

    subgraph "LLM Enhancement"
        Ollama[Ollama Runtime]
        Qwen[qwen2.5-coder:7b<br/>Code Analysis]
        DeepSeek[deepseek-r1:14b<br/>Reasoning]
        QwenSmall[qwen2.5:1.5b<br/>Fast Triage]
    end

    subgraph "Output"
        Result[Classification Result]
        Explanation[Detailed Explanation]
        Confidence[Confidence Score]
    end

    Code --> Antigen
    Patterns --> BCell

    Antigen --> Macro
    Macro --> Dend
    Dend --> BCell
    Dend --> NKCell
    Dend --> QML

    BCell --> Affinity
    BCell --> TCell
    NKCell --> TCell
    QML --> TCell

    TCell --> Result
    TCell --> Explanation
    TCell --> Confidence

    BCell -.uses.-> OptAiNet
    QML -.uses.-> QMLAiNet

    Macro -.LLM.-> QwenSmall
    Dend -.LLM.-> Qwen
    BCell -.LLM.-> Qwen
    NKCell -.LLM.-> DeepSeek
    TCell -.LLM.-> DeepSeek
    QML -.LLM.-> DeepSeek

    Qwen -.runs on.-> Ollama
    DeepSeek -.runs on.-> Ollama
    QwenSmall -.runs on.-> Ollama

    style Antigen fill:#e1f5fe
    style Affinity fill:#e1f5fe
    style Macro fill:#fff3e0
    style Dend fill:#fff3e0
    style BCell fill:#e8f5e9
    style NKCell fill:#fce4ec
    style QML fill:#f3e5f5
    style TCell fill:#fff9c4
    style OptAiNet fill:#ede7f6
    style QMLAiNet fill:#ede7f6
    style Ollama fill:#e0f2f1
    style Result fill:#c8e6c9
```

## Component Descriptions

### Core Framework
- **Antigen**: Represents code samples to be analyzed
- **Affinity Calculator**: Measures similarity between antigens and patterns

### Immune Agents
1. **Macrophage** (qwen2.5:1.5b)
   - First responder
   - Quick triage: LOW/MEDIUM/HIGH priority
   - Routes to appropriate agents

2. **Dendritic Cell** (qwen2.5-coder:7b)
   - Feature extraction from code
   - Identifies API calls, patterns, structures
   - Presents antigens to other agents

3. **B Cell** (qwen2.5-coder:7b)
   - Pattern matching with clonal selection
   - Learns from training examples
   - Uses semantic similarity via LLM

4. **NK Cell** (deepseek-r1:14b)
   - Anomaly detection
   - Self vs non-self discrimination
   - Deep reasoning about suspicious behavior

5. **QML Agent** (deepseek-r1:14b)
   - Qualitative model learning
   - Behavioral pattern recognition
   - Constraint-based reasoning

6. **T Cell** (deepseek-r1:14b)
   - Coordinator of all agents
   - Synthesizes findings
   - Makes final classification decision

### Algorithms
- **Opt-AiNet**: Multi-modal optimization with network suppression
- **QML-AiNet**: Qualitative differential equation learning

### LLM Runtime
- **Ollama**: Local LLM inference (offline-capable)
- **Models**: 3 models totaling ~15GB

## Data Flow

```mermaid
sequenceDiagram
    participant User
    participant Macro as Macrophage
    participant Dend as Dendritic
    participant BCell as B Cell
    participant NK as NK Cell
    participant QML as QML Agent
    participant TCell as T Cell

    User->>Macro: Submit Code
    Macro->>Macro: Triage (priority)

    Macro->>Dend: Forward (HIGH priority)
    Dend->>Dend: Extract Features
    Note over Dend: API calls, patterns,<br/>structure analysis

    par Parallel Analysis
        Dend->>BCell: Present Antigen
        BCell->>BCell: Pattern Match
        BCell-->>TCell: Result + Confidence
    and
        Dend->>NK: Present Antigen
        NK->>NK: Anomaly Detection
        NK-->>TCell: Result + Reasoning
    and
        Dend->>QML: Present Antigen
        QML->>QML: Qualitative Model
        QML-->>TCell: Behavior Assessment
    end

    TCell->>TCell: Coordinate Findings
    TCell->>TCell: Synthesize Decision
    TCell->>User: Classification + Explanation
```

## Links

- [[../code-mirror/src/core/antigen|Antigen Implementation]]
- [[../code-mirror/src/core/affinity|Affinity Calculator]]
- [[../code-mirror/src/agents/bcell_agent|B Cell Agent]]
- [[../code-mirror/src/agents/nk_cell_enhanced|NK Cell Agent]]
- [[../code-mirror/src/algorithms/opt_ainet|Opt-AiNet]]
- [[../code-mirror/src/algorithms/qml_ainet|QML-AiNet]]
- [[../code-mirror/examples/llm_enhanced_agents|LLM-Enhanced Agents]]

## References

- [[../../research/literature/de-castro-2002-opt-ainet|de Castro & Timmis (2002)]] - Opt-AiNet
- [[../../research/literature/pang-coghill-2015-qml-ainet|Pang & Coghill (2015)]] - QML-AiNet
- [[../docs/Model-Selection-By-Agent-Role|Model Selection Guide]]

---

*Generated: 2025-11-30*
