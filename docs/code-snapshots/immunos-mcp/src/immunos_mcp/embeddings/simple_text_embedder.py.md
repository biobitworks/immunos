---
source: /Users/byron/projects/immunos-mcp/src/immunos_mcp/embeddings/simple_text_embedder.py
relative: immunos-mcp/src/immunos_mcp/embeddings/simple_text_embedder.py
generated_at: 2025-12-23 10:28
---

```python
"""
Lightweight, deterministic text embedder for offline use.
"""

from __future__ import annotations

import hashlib
import re
from typing import Iterable

import numpy as np


class SimpleTextEmbedder:
    """
    Hash-based embedder for short text snippets.

    Produces deterministic fixed-size vectors without external dependencies.
    """

    def __init__(self, dim: int = 128):
        self.dim = dim

    def _tokenize(self, text: str) -> Iterable[str]:
        return re.findall(r"[a-zA-Z0-9_]+", text.lower())

    def embed(self, text: str) -> np.ndarray:
        vector = np.zeros(self.dim, dtype=np.float32)
        for token in self._tokenize(text):
            digest = hashlib.md5(token.encode("utf-8")).hexdigest()
            idx = int(digest[:8], 16) % self.dim
            sign = 1.0 if int(digest[8], 16) < 8 else -1.0
            vector[idx] += sign

        norm = np.linalg.norm(vector)
        if norm > 0:
            vector = vector / norm
        return vector

```
