---
project: immunos-mcp
source: __init__.py
type: code-mirror
language: py
size: 733
modified: 2025-11-17T10:31:59.035382
hash: 9cd73382a2d48dc04fdf94c24ac996af
description: "Core modules for IMMUNOS-MCP"
tags: [code, source, auto-generated]
---

> [!info] Source File
> **Path**: `__init__.py`
> **Size**: 733 bytes
> **Modified**: 2025-11-17
> **Project**: [[projects/immunos-mcp/|immunos-mcp]]
```python
"""Core modules for IMMUNOS-MCP"""

from .antigen import Antigen, AntigenBatch, DataType
from .affinity import AffinityCalculator, AffinityResult, DistanceMetric
from .protocols import (
    RecognitionResult,
    AnomalyResult,
    RecognitionStrategy,
    Signal,
    SignalType,
    ExtractedFeatures,
    Pattern,
    Detector,
    Clone,
    AgentRequest,
    AgentResponse,
)

__all__ = [
    "Antigen",
    "AntigenBatch",
    "DataType",
    "AffinityCalculator",
    "AffinityResult",
    "DistanceMetric",
    "RecognitionResult",
    "AnomalyResult",
    "RecognitionStrategy",
    "Signal",
    "SignalType",
    "ExtractedFeatures",
    "Pattern",
    "Detector",
    "Clone",
    "AgentRequest",
    "AgentResponse",
]

```
