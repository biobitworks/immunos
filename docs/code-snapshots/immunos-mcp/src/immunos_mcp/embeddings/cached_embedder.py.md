---
source: /Users/byron/projects/immunos-mcp/src/immunos_mcp/embeddings/cached_embedder.py
relative: immunos-mcp/src/immunos_mcp/embeddings/cached_embedder.py
generated_at: 2025-12-23 10:28
---

```python
"""
Cached embedder placeholder.

Provides a minimal interface for future embedding caching.
"""

from typing import Dict, Optional, List


class CachedEmbedder:
    """In-memory cache for embeddings."""

    def __init__(self):
        self._cache: Dict[str, List[float]] = {}

    def get(self, key: str) -> Optional[List[float]]:
        """Get cached embedding by key."""
        return self._cache.get(key)

    def set(self, key: str, vector: List[float]) -> None:
        """Store embedding by key."""
        self._cache[key] = vector

```
