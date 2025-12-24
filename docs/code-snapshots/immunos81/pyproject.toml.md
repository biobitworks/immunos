---
source: /Users/byron/projects/immunos81/pyproject.toml
relative: immunos81/pyproject.toml
generated_at: 2025-12-23 10:28
---

```toml
[project]
name = "immunos81"
version = "1.0.0"
description = "Classical Artificial Immune System for Pattern Recognition (Hunt & Cooke, 2000)"
authors = [{name = "Byron"}]
readme = "README.md"
requires-python = ">=3.10"
license = {text = "MIT"}
keywords = ["artificial immune system", "machine learning", "pattern recognition", "classification"]

dependencies = [
    "numpy>=1.24.0",
]

[project.optional-dependencies]
viz = [
    "matplotlib>=3.7.0",
    "plotly>=5.14.0",
]
dashboard = [
    "dash>=2.9.0",
    "plotly>=5.14.0",
]
dev = [
    "numpy>=1.24.0",
    "matplotlib>=3.7.0",
    "plotly>=5.14.0",
    "dash>=2.9.0",
    "scikit-learn>=1.2.0",
    "pytest>=7.3.0",
]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.black]
line-length = 100
target-version = ["py310"]

[tool.ruff]
line-length = 100
target-version = "py310"

```
