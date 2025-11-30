---
project: immunos-mcp
source: __init__.py
type: code-mirror
language: py
size: 519
modified: 2025-11-26T12:33:16.778252
hash: 117933927c6dba21d9f28026ec79d0ea
description: "Immune-Inspired Optimization Algorithms  Implementations of Opt-AiNet and QML-AiNet for optimization and qualitative model learning."
tags: [code, source, auto-generated]
---

> [!info] Source File
> **Path**: `__init__.py`
> **Size**: 519 bytes
> **Modified**: 2025-11-26
> **Project**: [[projects/immunos-mcp/|immunos-mcp]]
```python
"""
Immune-Inspired Optimization Algorithms

Implementations of Opt-AiNet and QML-AiNet for optimization and
qualitative model learning.
"""

from .opt_ainet import OptAiNet, Antibody
from .qml_ainet import (
    QMLAiNet,
    QDEModel,
    QualitativeConstraint,
    QualitativeValue,
    QualitativeChange,
    QualitativeRelation
)

__all__ = [
    "OptAiNet",
    "Antibody",
    "QMLAiNet",
    "QDEModel",
    "QualitativeConstraint",
    "QualitativeValue",
    "QualitativeChange",
    "QualitativeRelation",
]

```
