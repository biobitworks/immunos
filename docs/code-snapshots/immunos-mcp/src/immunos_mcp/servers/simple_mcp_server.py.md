---
source: /Users/byron/projects/immunos-mcp/src/immunos_mcp/servers/simple_mcp_server.py
relative: immunos-mcp/src/immunos_mcp/servers/simple_mcp_server.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
IMMUNOS-MCP Simple Server (MVP)

Minimal MCP server exposing B Cell pattern matching agent to Continue.dev.
This is a simplified version for quick setup and testing.

Usage:
    python src/immunos_mcp/servers/simple_mcp_server.py

Or install and run:
    uv pip install -e .
    immunos-mcp-mvp
"""

import asyncio
import sys
from typing import Any, Dict

try:
    from mcp.server import Server
    from mcp.server.stdio import stdio_server
    from mcp.types import Tool, TextContent
except ImportError:
    print("ERROR: MCP SDK not installed. Run: uv pip install mcp>=1.3.1")
    sys.exit(1)

from immunos_mcp.core.antigen import Antigen, DataType
from immunos_mcp.orchestrator import ImmunosOrchestrator


class SimpleMCPServer:
    """Minimal MCP server wrapper for IMMUNOS orchestrator"""

    def __init__(self):
        self.server = Server("immunos-mcp-mvp")
        self.orchestrator = ImmunosOrchestrator()
        self._register_tools()

    def _setup_agent(self):
        """Kept for API compatibility (no-op)."""
        return None

    def _register_tools(self):
        """Register MCP tool for Continue"""

        @self.server.list_tools()
        async def list_tools() -> list[Tool]:
            return [
                Tool(
                    name="immune_scan",
                    description="Scan code for security vulnerabilities using IMMUNOS B Cell pattern matching. Returns classification (safe/vulnerable) with confidence score.",
                    inputSchema={
                        "type": "object",
                        "properties": {
                            "code": {
                                "type": "string",
                                "description": "The code snippet to analyze for security issues"
                            }
                        },
                        "required": ["code"]
                    }
                )
            ]

        @self.server.call_tool()
        async def call_tool(name: str, arguments: Dict[str, Any]) -> list[TextContent]:
            if name == "immune_scan":
                return await self._scan_code(arguments)
            else:
                raise ValueError(f"Unknown tool: {name}")

    async def _scan_code(self, args: Dict[str, Any]) -> list[TextContent]:
        """Run B Cell pattern matching on code"""
        code = args['code']

        # Create antigen
        antigen = Antigen(
            data=code,
            data_type=DataType.CODE,
            identifier="user_code"
        )

        # Run recognition
        result = self.orchestrator.analyze(antigen)

        # Determine risk level
        if result.anomaly:
            if result.confidence > 0.8:
                risk = "🔴 HIGH RISK"
                recommendation = "Do NOT use this code. Security review required."
            elif result.confidence > 0.5:
                risk = "🟡 MEDIUM RISK"
                recommendation = "Proceed with caution. Manual review recommended."
            else:
                risk = "🟢 LOW RISK (uncertain)"
                recommendation = "Possibly safe, but low confidence. Review recommended."
        else:  # safe
            if result.confidence > 0.8:
                risk = "🟢 SAFE"
                recommendation = "Code appears safe. Standard review process."
            else:
                risk = "🟡 UNCLEAR"
                recommendation = "Classification uncertain. Manual review recommended."

        # Format output
        output = f"""## 🧬 IMMUNOS Security Scan

**Classification:** {result.classification.upper() if result.classification else 'UNKNOWN'}
**Risk Level:** {risk}
**Confidence:** {result.confidence:.1%}

### Pattern Matching Results
**Avidity Scores:**
"""

        avidity_scores = result.metadata.get("bcell", {}).get("avidity_scores", {})
        if avidity_scores:
            for label, score in sorted(avidity_scores.items(), key=lambda x: x[1], reverse=True):
                bar = '█' * int(score * 20)
                output += f"- {label}: {score:.3f} {bar}\n"
        else:
            output += "No pattern matches found\n"

        output += f"\n### Recommendation\n{recommendation}\n"

        output += f"\n### Analysis Details\n"
        output += f"- Patterns matched: {len(avidity_scores)}\n"
        output += f"- Strategy: Simple Highest Avidity (SHA)\n"
        output += f"- Uncertain: {'Yes' if result.metadata.get('bcell', {}).get('is_uncertain') else 'No'}\n"
        output += f"- Agents: {', '.join(result.agents)}\n"
        output += f"- Signals: {result.signals}\n"
        output += f"- Model Roles: {result.model_roles}\n"

        output += "\n---\n*Powered by IMMUNOS-MCP B Cell Agent*"

        return [TextContent(type="text", text=output)]

    async def run(self):
        """Run the MCP server"""
        print("IMMUNOS-MCP MVP Server starting...", file=sys.stderr)
        print("Waiting for MCP client connection via stdio...", file=sys.stderr)

        async with stdio_server() as (read_stream, write_stream):
            await self.server.run(
                read_stream,
                write_stream,
                self.server.create_initialization_options()
            )


async def main():
    """Entry point"""
    server = SimpleMCPServer()
    await server.run()


if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("\nServer stopped", file=sys.stderr)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

```
