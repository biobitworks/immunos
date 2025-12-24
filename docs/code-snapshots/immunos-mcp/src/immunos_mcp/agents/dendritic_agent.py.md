---
source: /Users/byron/projects/immunos-mcp/src/immunos_mcp/agents/dendritic_agent.py
relative: immunos-mcp/src/immunos_mcp/agents/dendritic_agent.py
generated_at: 2025-12-23 10:28
---

```python
"""
Dendritic Agent for IMMUNOS

Extracts features and signals from antigens for downstream agents.
"""

from typing import Dict, Any
from ..core.antigen import Antigen


class DendriticAgent:
    """Simple feature extractor and signal aggregator."""

    def __init__(self, agent_name: str = "dendritic_001"):
        self.agent_name = agent_name

    def extract_features(self, antigen: Antigen) -> Dict[str, Any]:
        """Return basic text features for downstream processing."""
        text = antigen.get_text_content()
        length = len(text)
        lines = text.count("\n") + 1
        tokens = len(text.split())

        return {
            "length": length,
            "lines": lines,
            "tokens": tokens,
        }

```
