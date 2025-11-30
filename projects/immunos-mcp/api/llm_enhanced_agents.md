---
project: immunos-mcp
type: api-documentation
source: immunos-mcp/examples/llm_enhanced_agents.py
generated: 2025-11-30
tags: [api, documentation, immunos-mcp]
---

# llm_enhanced_agents

**Source**: [[../code-mirror/immunos-mcp/examples/llm_enhanced_agents.py|immunos-mcp/examples/llm_enhanced_agents.py]]

LLM-Enhanced Immune System Agents for IMMUNOS-MCP

Wraps base agents with Ollama LLM capabilities for:
- Semantic code understanding
- Natural language explanations
- Improved anomaly detection
- Qualitative reasoning

Model assignments (16GB configuration):
- B Cell + Dendritic: qwen2.5-coder:7b (code analysis)
- NK Cell + T Cell + QML: deepseek-r1:14b (deep reasoning)
- Macrophage: qwen2.5:1.5b (fast triage)

## Contents

### Classes
- [LLMEnhancedBCellAgent](#llmenhancedbcellagent)
- [LLMEnhancedNKCellAgent](#llmenhancednkcellagent)
- [LLMEnhancedTCellAgent](#llmenhancedtcellagent)
- [LLMEnhancedDendriticAgent](#llmenhanceddendriticagent)
- [LLMEnhancedMacrophageAgent](#llmenhancedmacrophageagent)
- [LLMEnhancedQMLAgent](#llmenhancedqmlagent)

## Classes

### LLMEnhancedBCellAgent

**Inherits from**: `BCellAgent`

**Source**: `immunos-mcp/examples/llm_enhanced_agents.py:85`

B Cell agent with LLM-powered semantic understanding.

Model: qwen2.5-coder:7b (4.7GB)
Purpose: Understand code semantics for better pattern matching

#### Methods

##### `__init__(`self, agent_name: str, affinity_method: str, llm_model: str)`

*No documentation available.*

##### `calculate_semantic_similarity(`self, code1: str, code2: str)` → `float`

Calculate semantic similarity using LLM.

Args:
    code1: First code snippet
    code2: Second code snippet

Returns:
    Similarity score (0.0 to 1.0)

##### `explain_pattern_match(`self, antigen: Antigen, result: RecognitionResult)` → `str`

Generate natural language explanation of pattern match.

Args:
    antigen: Test antigen
    result: Recognition result

Returns:
    Human-readable explanation

##### `recognize(`self, antigen: Antigen, antigen_embedding: Optional[np.ndarray], strategy: str, explain: bool)` → `RecognitionResult`

Recognize with optional LLM explanation.

Args:
    antigen: Antigen to classify
    antigen_embedding: Optional embedding
    strategy: Recognition strategy
    explain: Generate LLM explanation

Returns:
    Recognition result with optional LLM explanation

---

### LLMEnhancedNKCellAgent

**Source**: `immunos-mcp/examples/llm_enhanced_agents.py:207`

NK Cell agent with LLM-powered anomaly reasoning.

Model: deepseek-r1:14b (8GB)
Purpose: Deep reasoning about anomalies and threats

#### Methods

##### `__init__(`self, agent_name: str, llm_model: str)`

*No documentation available.*

##### `learn_baseline(`self, normal_samples: List[str])`

Learn baseline behavior from normal samples.

Args:
    normal_samples: List of normal code samples

##### `detect_anomaly(`self, code: str, explain: bool)` → `Tuple[(bool, float, str)]`

Detect if code is anomalous.

Args:
    code: Code to analyze
    explain: Generate detailed explanation

Returns:
    (is_anomalous, confidence, explanation)

##### `recognize(`self, antigen: Antigen, explain: bool)` → `RecognitionResult`

Recognize (detect anomalies) in antigen.

Args:
    antigen: Code antigen
    explain: Generate detailed explanation

Returns:
    Recognition result

---

### LLMEnhancedTCellAgent

**Source**: `immunos-mcp/examples/llm_enhanced_agents.py:351`

T Cell coordinator with LLM-powered decision making.

Model: deepseek-r1:14b (8GB)
Purpose: Coordinate multiple agents and synthesize findings

#### Methods

##### `__init__(`self, agent_name: str, llm_model: str)`

*No documentation available.*

##### `coordinate(`self, antigen: Antigen, agent_results: List[RecognitionResult], explain: bool)` → `RecognitionResult`

Coordinate multiple agent responses.

Args:
    antigen: Test antigen
    agent_results: Results from other agents
    explain: Generate detailed coordination explanation

Returns:
    Coordinated result

---

### LLMEnhancedDendriticAgent

**Source**: `immunos-mcp/examples/llm_enhanced_agents.py:446`

Dendritic cell with LLM-powered feature extraction.

Model: qwen2.5-coder:7b (4.7GB)
Purpose: Extract semantic features from code

#### Methods

##### `__init__(`self, agent_name: str, llm_model: str)`

*No documentation available.*

##### `extract_features(`self, code: str)` → `Dict[(str, Any)]`

Extract features from code using LLM.

Args:
    code: Source code

Returns:
    Feature dictionary

##### `process_antigen(`self, antigen: Antigen)` → `Antigen`

Process antigen and add extracted features.

Args:
    antigen: Input antigen

Returns:
    Antigen with enhanced features

---

### LLMEnhancedMacrophageAgent

**Source**: `immunos-mcp/examples/llm_enhanced_agents.py:533`

Macrophage with LLM-powered fast triage.

Model: qwen2.5:1.5b (1.0GB)
Purpose: Quick initial assessment

#### Methods

##### `__init__(`self, agent_name: str, llm_model: str)`

*No documentation available.*

##### `triage(`self, code: str)` → `Dict[(str, Any)]`

Quick triage of code sample.

Args:
    code: Source code

Returns:
    Triage result with priority and category

---

### LLMEnhancedQMLAgent

**Source**: `immunos-mcp/examples/llm_enhanced_agents.py:606`

QML agent with LLM-powered qualitative reasoning.

Model: deepseek-r1:14b (8GB)
Purpose: Learn qualitative behavioral models

#### Methods

##### `__init__(`self, agent_name: str, llm_model: str)`

*No documentation available.*

##### `learn_behavior_model(`self, code_samples: List[str], labels: List[str])` → `Dict[(str, Any)]`

Learn qualitative behavioral model.

Args:
    code_samples: Training code samples
    labels: Labels (e.g., "safe", "malicious")

Returns:
    Learned model description

##### `recognize(`self, antigen: Antigen, model: Dict[(str, Any)])` → `RecognitionResult`

Recognize using qualitative model.

Args:
    antigen: Test antigen
    model: Learned model

Returns:
    Recognition result

---

## Links

- [[../code-mirror/immunos-mcp/examples/llm_enhanced_agents.py|Source Code]]
- [[../../projects/immunos-mcp/|Project Root]]
- [[./README|API Index]]

