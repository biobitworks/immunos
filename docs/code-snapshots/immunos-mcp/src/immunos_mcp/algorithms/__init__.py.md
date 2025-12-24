---
source: /Users/byron/projects/immunos-mcp/src/immunos_mcp/algorithms/__init__.py
relative: immunos-mcp/src/immunos_mcp/algorithms/__init__.py
generated_at: 2025-12-23 10:28
---

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
