---
project: immunos-mcp
source: pyproject.toml
type: code-mirror
language: toml
size: 801
modified: 2025-11-25T14:21:46.410401
hash: c536b09d5d89ffde6f6d136d81d9788a
tags: [code, source, auto-generated]
---

> [!info] Source File
> **Path**: `pyproject.toml`
> **Size**: 801 bytes
> **Modified**: 2025-11-25
> **Project**: [[projects/immunos-mcp/|immunos-mcp]]
```toml
[project]
name = "immunos-mcp"
version = "0.1.0"
description = "Artificial Immune System MCP Server - Multi-agent LLM system with immune-inspired roles"
authors = [{name = "IMMUNOS-MCP Team"}]
readme = "README.md"
requires-python = ">=3.10"
dependencies = [
    "numpy>=1.26.0",
    "pydantic>=2.0.0",
]

[project.optional-dependencies]
mcp = [
    "mcp>=1.3.1",
]
llm = [
    "anthropic>=0.39.0",
    "openai>=1.0.0",
]
vector = [
    "chromadb>=0.5.0",
]
dev = [
    "pytest>=8.0.0",
    "pytest-asyncio>=0.24.0",
    "black>=24.0.0",
    "ruff>=0.7.0",
]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.hatch.build.targets.wheel]
packages = ["src"]

[tool.black]
line-length = 100
target-version = ["py310"]

[tool.ruff]
line-length = 100
target-version = "py310"

```
