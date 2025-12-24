---
source: /Users/byron/projects/rockbeatspaper/environments/my_env/pyproject.toml
relative: rockbeatspaper/environments/my_env/pyproject.toml
generated_at: 2025-12-23 10:28
---

```toml
[project]
name = "my-env"
description = "Your environment description here"
tags = ["placeholder-tag", "train", "eval"]
version = "0.1.0"
requires-python = ">=3.10"
dependencies = [
    "verifiers>=0.1.7.post0",
]

[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[tool.hatch.build]
include = ["my_env.py", "pyproject.toml"] 

[tool.verifiers.eval]
num_examples = 5
rollouts_per_example = 3

```
