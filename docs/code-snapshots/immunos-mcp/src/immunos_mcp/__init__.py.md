---
source: /Users/byron/projects/immunos-mcp/src/immunos_mcp/__init__.py
relative: immunos-mcp/src/immunos_mcp/__init__.py
generated_at: 2025-12-23 10:28
---

```python
"""
IMMUNOS: Artificial Immune System Multi-Agent System

Standalone AIS implementation. MCP is optional for packaging or integration.
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
