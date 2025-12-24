---
source: /Users/byron/projects/scripts/export_code_snapshots.py
relative: scripts/export_code_snapshots.py
generated_at: 2025-12-23 10:28
---

````python
#!/usr/bin/env python3
"""
Export code files to markdown snapshots for Obsidian/VSCode/Xcode viewing.
"""

from __future__ import annotations

import argparse
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List

CODE_EXTENSIONS: Dict[str, str] = {
    ".py": "python",
    ".js": "javascript",
    ".jsx": "jsx",
    ".ts": "typescript",
    ".tsx": "tsx",
    ".mjs": "javascript",
    ".cjs": "javascript",
    ".c": "c",
    ".h": "c",
    ".cpp": "cpp",
    ".hpp": "cpp",
    ".cc": "cpp",
    ".go": "go",
    ".rs": "rust",
    ".java": "java",
    ".kt": "kotlin",
    ".swift": "swift",
    ".cs": "csharp",
    ".rb": "ruby",
    ".php": "php",
    ".sh": "bash",
    ".zsh": "bash",
    ".bash": "bash",
    ".ps1": "powershell",
    ".sql": "sql",
    ".proto": "proto",
    ".yaml": "yaml",
    ".yml": "yaml",
    ".toml": "toml",
    ".ini": "ini",
    ".cfg": "ini",
    ".json": "json",
    ".xml": "xml",
    ".html": "html",
    ".css": "css",
    ".scss": "scss",
    ".mdx": "mdx",
    ".dockerfile": "docker",
}

SPECIAL_FILENAMES: Dict[str, str] = {
    "makefile": "make",
    "dockerfile": "docker",
    "docker-compose.yml": "yaml",
    "docker-compose.yaml": "yaml",
}

SKIP_DIRS = {
    ".git",
    "node_modules",
    ".venv",
    "venv",
    ".immunos",
    ".pytest_cache",
    "__pycache__",
    ".cache",
    ".tox",
    "dist",
    "build",
    ".egg-info",
    ".mypy_cache",
    ".idea",
    ".vscode",
    ".obsidian",
}


def iter_code_files(root: Path, output_root: Path) -> Iterable[Path]:
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        current = Path(dirpath)

        # Skip output root and known dirs
        if output_root in current.parents or current == output_root:
            dirnames[:] = []
            continue

        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS and not d.startswith(".")]

        for name in filenames:
            lower = name.lower()
            if lower in SPECIAL_FILENAMES:
                yield current / name
                continue

            suffix = Path(name).suffix.lower()
            if suffix in CODE_EXTENSIONS:
                yield current / name


def language_for(path: Path) -> str:
    lower = path.name.lower()
    if lower in SPECIAL_FILENAMES:
        return SPECIAL_FILENAMES[lower]
    suffix = path.suffix.lower()
    return CODE_EXTENSIONS.get(suffix, "")


def write_snapshot(src_path: Path, root: Path, output_root: Path) -> None:
    rel = src_path.relative_to(root)
    out_path = output_root / rel
    out_path = out_path.with_name(out_path.name + ".md")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    raw = src_path.read_bytes()
    text = raw.decode("utf-8", errors="replace")
    fence = "```"
    if fence in text:
        fence = "````"

    lang = language_for(src_path)
    header = (
        "---\n"
        f"source: {src_path}\n"
        f"relative: {rel.as_posix()}\n"
        f"generated_at: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n"
        "---\n\n"
    )

    content = f"{header}{fence}{lang}\n{text}\n{fence}\n"
    out_path.write_text(content, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Export code files to markdown snapshots")
    parser.add_argument("--root", default="/Users/byron/projects", help="Root folder to scan")
    parser.add_argument("--output", default="/Users/byron/projects/docs/code-snapshots", help="Output folder")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    output_root = Path(args.output).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    count = 0
    for src in iter_code_files(root, output_root):
        write_snapshot(src, root, output_root)
        count += 1

    print(f"Exported {count} code files to {output_root}")


if __name__ == "__main__":
    main()

````
