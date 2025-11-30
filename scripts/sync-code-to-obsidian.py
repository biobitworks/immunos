#!/usr/bin/env python3
"""
Sync Code to Obsidian

Mirrors Python source code as Markdown files for viewing in Obsidian.
Preserves directory structure and adds metadata frontmatter.

Usage:
    python sync-code-to-obsidian.py [project_name]
    python sync-code-to-obsidian.py immunos-mcp
    python sync-code-to-obsidian.py --all  # All projects

Features:
- Converts .py → .md with syntax highlighting
- Adds YAML frontmatter with metadata
- Preserves directory structure
- Links to API docs and related files
- Skips __pycache__, .venv, etc.
"""

import os
import sys
from pathlib import Path
from typing import List, Dict, Any
from datetime import datetime
import hashlib

# Projects root
PROJECTS_ROOT = Path("/Users/byron/projects")
VAULT_ROOT = PROJECTS_ROOT  # Obsidian vault is at projects root

# Projects to sync
PROJECTS = {
    "immunos-mcp": PROJECTS_ROOT / "immunos-mcp",
    "immunos81": PROJECTS_ROOT / "immunos81",
    "rockbeatspaper": PROJECTS_ROOT / "rockbeatspaper",
    "experiments": PROJECTS_ROOT / "experiments",
}

# Directories to skip
SKIP_DIRS = {
    "__pycache__",
    ".venv",
    "venv",
    ".git",
    ".pytest_cache",
    ".mypy_cache",
    "node_modules",
    ".obsidian",
    "code-mirror",  # Don't mirror the mirror!
}

# File extensions to mirror
MIRROR_EXTENSIONS = {".py", ".sh", ".yml", ".yaml", ".toml", ".md"}


def get_file_metadata(filepath: Path, project_root: Path) -> Dict[str, Any]:
    """Extract metadata from source file."""
    stat = filepath.stat()

    # Relative path from project root
    rel_path = filepath.relative_to(project_root)

    # Calculate file hash for change detection
    with open(filepath, 'rb') as f:
        file_hash = hashlib.md5(f.read()).hexdigest()

    # Extract first docstring if Python file
    docstring = None
    if filepath.suffix == ".py":
        try:
            with open(filepath, 'r') as f:
                content = f.read()
                # Simple docstring extraction (first triple-quoted string)
                if '"""' in content:
                    start = content.find('"""') + 3
                    end = content.find('"""', start)
                    if end != -1:
                        docstring = content[start:end].strip()
        except:
            pass

    return {
        "source_path": str(rel_path),
        "file_type": filepath.suffix[1:],  # Remove dot
        "size_bytes": stat.st_size,
        "modified": datetime.fromtimestamp(stat.st_mtime).isoformat(),
        "file_hash": file_hash,
        "docstring": docstring,
    }


def create_frontmatter(metadata: Dict[str, Any], project_name: str) -> str:
    """Generate YAML frontmatter for markdown file."""
    lines = ["---"]
    lines.append(f"project: {project_name}")
    lines.append(f"source: {metadata['source_path']}")
    lines.append(f"type: code-mirror")
    lines.append(f"language: {metadata['file_type']}")
    lines.append(f"size: {metadata['size_bytes']}")
    lines.append(f"modified: {metadata['modified']}")
    lines.append(f"hash: {metadata['file_hash']}")

    if metadata.get("docstring"):
        # Escape and truncate docstring
        doc = metadata["docstring"].replace("\n", " ")[:200]
        lines.append(f"description: \"{doc}\"")

    lines.append("tags: [code, source, auto-generated]")
    lines.append("---")
    lines.append("")
    return "\n".join(lines)


def get_syntax_highlight_lang(extension: str) -> str:
    """Get syntax highlighting language for file extension."""
    lang_map = {
        ".py": "python",
        ".sh": "bash",
        ".yml": "yaml",
        ".yaml": "yaml",
        ".toml": "toml",
        ".md": "markdown",
        ".json": "json",
        ".txt": "text",
    }
    return lang_map.get(extension, "text")


def mirror_file(source_path: Path, mirror_path: Path, project_name: str) -> bool:
    """
    Mirror a single source file to markdown.

    Returns:
        True if file was updated, False if skipped (no changes)
    """
    # Check if mirror already exists and is up to date
    if mirror_path.exists():
        # Check hash in frontmatter
        try:
            with open(mirror_path, 'r') as f:
                content = f.read()
                if "hash:" in content:
                    existing_hash = content.split("hash:")[1].split("\n")[0].strip()
                    current_hash = get_file_metadata(source_path, source_path.parent).get("file_hash")
                    if existing_hash == current_hash:
                        return False  # No changes, skip
        except:
            pass

    # Get metadata
    metadata = get_file_metadata(source_path, source_path.parent)

    # Read source content
    try:
        with open(source_path, 'r') as f:
            source_content = f.read()
    except UnicodeDecodeError:
        # Binary file, skip
        return False

    # Create markdown content
    frontmatter = create_frontmatter(metadata, project_name)

    # Add navigation links
    nav_links = [
        f"\n> [!info] Source File",
        f"> **Path**: `{metadata['source_path']}`",
        f"> **Size**: {metadata['size_bytes']} bytes",
        f"> **Modified**: {metadata['modified'][:10]}",
        f"> **Project**: [[projects/{project_name}/|{project_name}]]",
        "",
    ]

    # Add syntax-highlighted code block
    lang = get_syntax_highlight_lang(source_path.suffix)
    code_block = [
        f"```{lang}",
        source_content,
        "```",
        "",
    ]

    # Combine all parts
    markdown_content = frontmatter + "\n".join(nav_links) + "\n".join(code_block)

    # Write mirror file
    mirror_path.parent.mkdir(parents=True, exist_ok=True)
    with open(mirror_path, 'w') as f:
        f.write(markdown_content)

    return True


def mirror_project(project_name: str, project_root: Path, verbose: bool = True) -> Dict[str, int]:
    """
    Mirror all source files from a project.

    Returns:
        Statistics dict with counts
    """
    mirror_root = VAULT_ROOT / "projects" / project_name / "code-mirror"

    stats = {
        "scanned": 0,
        "mirrored": 0,
        "skipped": 0,
        "errors": 0,
    }

    if verbose:
        print(f"\n🔄 Mirroring {project_name}...")
        print(f"   Source: {project_root}")
        print(f"   Mirror: {mirror_root}")

    # Walk through project directory
    for root, dirs, files in os.walk(project_root):
        # Skip excluded directories
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]

        for filename in files:
            source_path = Path(root) / filename

            # Check extension
            if source_path.suffix not in MIRROR_EXTENSIONS:
                continue

            stats["scanned"] += 1

            # Calculate mirror path
            rel_path = source_path.relative_to(project_root)
            mirror_path = mirror_root / rel_path.with_suffix(".md")

            # Mirror the file
            try:
                updated = mirror_file(source_path, mirror_path, project_name)
                if updated:
                    stats["mirrored"] += 1
                    if verbose:
                        print(f"   ✓ {rel_path}")
                else:
                    stats["skipped"] += 1
            except Exception as e:
                stats["errors"] += 1
                if verbose:
                    print(f"   ✗ {rel_path}: {e}")

    if verbose:
        print(f"\n   📊 Stats:")
        print(f"      Scanned:  {stats['scanned']}")
        print(f"      Mirrored: {stats['mirrored']}")
        print(f"      Skipped:  {stats['skipped']}")
        print(f"      Errors:   {stats['errors']}")

    return stats


def create_index_file(project_name: str, project_root: Path, stats: Dict[str, int]):
    """Create an index file for the project code mirror."""
    mirror_root = VAULT_ROOT / "projects" / project_name / "code-mirror"
    index_path = mirror_root / "README.md"

    # Collect all mirrored files organized by directory
    files_by_dir = {}
    for root, dirs, files in os.walk(mirror_root):
        if files:
            rel_dir = Path(root).relative_to(mirror_root)
            files_by_dir[str(rel_dir)] = sorted([f for f in files if f.endswith(".md") and f != "README.md"])

    # Generate index content
    content = [
        "---",
        f"project: {project_name}",
        "type: code-mirror-index",
        f"generated: {datetime.now().isoformat()}",
        "tags: [index, code-mirror, auto-generated]",
        "---",
        "",
        f"# {project_name} - Code Mirror",
        "",
        f"Mirrored source code from [[projects/{project_name}/|{project_name}]] project.",
        "",
        "## Statistics",
        "",
        f"- **Files Scanned**: {stats['scanned']}",
        f"- **Files Mirrored**: {stats['mirrored']}",
        f"- **Files Skipped**: {stats['skipped']} (no changes)",
        f"- **Errors**: {stats['errors']}",
        f"- **Last Updated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Directory Structure",
        "",
    ]

    # Add directory tree
    for dir_path in sorted(files_by_dir.keys()):
        files = files_by_dir[dir_path]
        if dir_path == ".":
            content.append("### Root")
        else:
            content.append(f"### {dir_path}/")
        content.append("")
        for filename in files:
            file_path = Path(dir_path) / filename if dir_path != "." else filename
            display_name = filename.replace(".md", "")
            content.append(f"- [[{file_path.as_posix().replace('.md', '')}|{display_name}]]")
        content.append("")

    # Add links
    content.extend([
        "## Quick Links",
        "",
        f"- [[projects/{project_name}/|Project Home]]",
        f"- [[projects/{project_name}/journal/|Project Journal]]",
        f"- [[projects/{project_name}/api/|API Documentation]]",
        f"- [[daily/{datetime.now().strftime('%Y-%m-%d')}|Today's Log]]",
        "",
        "---",
        "",
        "*Auto-generated by sync-code-to-obsidian.py*",
    ])

    # Write index file
    with open(index_path, 'w') as f:
        f.write("\n".join(content))

    print(f"   📑 Created index: {index_path.name}")


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Sync source code to Obsidian vault")
    parser.add_argument("project", nargs="?", help="Project name (or --all)")
    parser.add_argument("--all", action="store_true", help="Sync all projects")
    parser.add_argument("--quiet", action="store_true", help="Minimal output")

    args = parser.parse_args()

    verbose = not args.quiet

    if verbose:
        print("=" * 70)
        print("Code Mirror - Sync to Obsidian")
        print("=" * 70)

    # Determine which projects to sync
    if args.all or (args.project and args.project == "--all"):
        projects_to_sync = PROJECTS.items()
    elif args.project:
        if args.project not in PROJECTS:
            print(f"❌ Unknown project: {args.project}")
            print(f"Available: {', '.join(PROJECTS.keys())}")
            return 1
        projects_to_sync = [(args.project, PROJECTS[args.project])]
    else:
        # Default: sync all projects
        projects_to_sync = PROJECTS.items()

    # Sync each project
    total_stats = {"scanned": 0, "mirrored": 0, "skipped": 0, "errors": 0}

    for project_name, project_root in projects_to_sync:
        if not project_root.exists():
            if verbose:
                print(f"⚠️  Skipping {project_name}: {project_root} does not exist")
            continue

        stats = mirror_project(project_name, project_root, verbose)

        # Create index
        create_index_file(project_name, project_root, stats)

        # Accumulate stats
        for key in total_stats:
            total_stats[key] += stats[key]

    # Final summary
    if verbose:
        print("\n" + "=" * 70)
        print("Summary")
        print("=" * 70)
        print(f"Total Files Scanned:  {total_stats['scanned']}")
        print(f"Total Files Mirrored: {total_stats['mirrored']}")
        print(f"Total Files Skipped:  {total_stats['skipped']}")
        print(f"Total Errors:         {total_stats['errors']}")
        print("")
        print("✅ Code mirror sync complete!")
        print(f"📁 View in Obsidian: {VAULT_ROOT}")

    return 0 if total_stats["errors"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
