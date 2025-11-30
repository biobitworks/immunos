---
project: immunos-mcp
source: __init__.py
type: code-mirror
language: py
size: 289
modified: 2025-11-25T14:16:08.400073
hash: 125d0e43eb6e8c0f27bece299118aa82
description: "Code Security Scanner Datasets  Contains curated examples of safe and vulnerable code patterns for training the immune system agents."
tags: [code, source, auto-generated]
---

> [!info] Source File
> **Path**: `__init__.py`
> **Size**: 289 bytes
> **Modified**: 2025-11-25
> **Project**: [[projects/immunos-mcp/|immunos-mcp]]
```python
"""
Code Security Scanner Datasets

Contains curated examples of safe and vulnerable code patterns
for training the immune system agents.
"""

from .safe_patterns import SAFE_PATTERNS
from .vulnerable_patterns import VULNERABLE_PATTERNS

__all__ = ["SAFE_PATTERNS", "VULNERABLE_PATTERNS"]

```
