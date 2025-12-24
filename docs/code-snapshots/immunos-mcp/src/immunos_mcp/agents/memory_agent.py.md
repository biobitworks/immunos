---
source: /Users/byron/projects/immunos-mcp/src/immunos_mcp/agents/memory_agent.py
relative: immunos-mcp/src/immunos_mcp/agents/memory_agent.py
generated_at: 2025-12-23 10:28
---

```python
"""
Memory Agent for IMMUNOS

Stores and retrieves prior observations (file-backed).
"""

import json
from pathlib import Path
from typing import Dict, Any, Optional


class MemoryAgent:
    """In-memory placeholder for memory retrieval."""

    def __init__(self, agent_name: str = "memory_001", path: Optional[Path] = None):
        self.agent_name = agent_name
        self.path = path or Path(".immunos") / "memory" / "immunos_memory.json"
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._store: Dict[str, Any] = self._load()

    def store(self, key: str, value: Any) -> None:
        """Store a value by key."""
        self._store[key] = value
        self._save()

    def retrieve(self, key: str) -> Optional[Any]:
        """Retrieve a stored value by key."""
        return self._store.get(key)

    def clear(self) -> None:
        """Clear all stored memory."""
        self._store = {}
        self._save()

    def _load(self) -> Dict[str, Any]:
        """Load memory from disk if available."""
        if not self.path.exists():
            return {}
        try:
            with self.path.open("r", encoding="utf-8") as handle:
                return json.load(handle)
        except (json.JSONDecodeError, OSError):
            return {}

    def _save(self) -> None:
        """Persist memory to disk."""
        with self.path.open("w", encoding="utf-8") as handle:
            json.dump(self._store, handle, indent=2)

```
