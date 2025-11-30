---
project: immunos-mcp
source: __init__.py
type: code-mirror
language: py
size: 697
modified: 2025-11-17T10:31:51.477331
hash: 050dbcabc374543a349be2c368583048
description: "IMMUNOS-MCP: Artificial Immune System Multi-Agent LLM Server  A Model Context Protocol server implementing biological immune system principles for adaptive pattern recognition and anomaly detection."
tags: [code, source, auto-generated]
---

> [!info] Source File
> **Path**: `__init__.py`
> **Size**: 697 bytes
> **Modified**: 2025-11-17
> **Project**: [[projects/immunos-mcp/|immunos-mcp]]
```python
"""
IMMUNOS-MCP: Artificial Immune System Multi-Agent LLM Server

A Model Context Protocol server implementing biological immune system
principles for adaptive pattern recognition and anomaly detection.
"""

__version__ = "0.1.0"

from .core.antigen import Antigen, AntigenBatch, DataType
from .core.affinity import AffinityCalculator, AffinityResult
from .core.protocols import (
    RecognitionResult,
    AnomalyResult,
    RecognitionStrategy,
    Signal,
    SignalType,
)

__all__ = [
    "Antigen",
    "AntigenBatch",
    "DataType",
    "AffinityCalculator",
    "AffinityResult",
    "RecognitionResult",
    "AnomalyResult",
    "RecognitionStrategy",
    "Signal",
    "SignalType",
]

```
