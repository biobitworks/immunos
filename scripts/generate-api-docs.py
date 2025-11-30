#!/usr/bin/env python3
"""
Generate API Documentation from Docstrings

Extracts docstrings from Python source files and generates
markdown documentation for Obsidian.

Usage:
    python scripts/generate-api-docs.py <project-name>

Example:
    python scripts/generate-api-docs.py immunos-mcp
"""

import ast
import sys
from pathlib import Path
from typing import List, Dict, Any, Optional
from datetime import datetime


class DocstringExtractor(ast.NodeVisitor):
    """Extract docstrings from Python AST."""

    def __init__(self, module_path: Path):
        self.module_path = module_path
        self.module_doc = None
        self.classes = []
        self.functions = []

    def visit_Module(self, node: ast.Module):
        """Extract module-level docstring."""
        self.module_doc = ast.get_docstring(node)
        self.generic_visit(node)

    def visit_ClassDef(self, node: ast.ClassDef):
        """Extract class docstring and methods."""
        class_info = {
            'name': node.name,
            'docstring': ast.get_docstring(node),
            'methods': [],
            'line': node.lineno
        }

        # Extract base classes
        bases = [self._get_name(base) for base in node.bases]
        class_info['bases'] = bases

        # Extract methods
        for item in node.body:
            if isinstance(item, ast.FunctionDef):
                method_info = self._extract_function(item)
                method_info['is_method'] = True
                class_info['methods'].append(method_info)

        self.classes.append(class_info)

    def visit_FunctionDef(self, node: ast.FunctionDef):
        """Extract function docstring."""
        # Only extract top-level functions (not methods)
        if not isinstance(node, ast.AsyncFunctionDef):
            # Check if this is a top-level function
            for parent in ast.walk(ast.parse("")):
                if hasattr(parent, 'body') and node in getattr(parent, 'body', []):
                    func_info = self._extract_function(node)
                    func_info['is_method'] = False
                    self.functions.append(func_info)
                    break

    def _extract_function(self, node: ast.FunctionDef) -> Dict[str, Any]:
        """Extract function/method information."""
        func_info = {
            'name': node.name,
            'docstring': ast.get_docstring(node),
            'line': node.lineno,
            'args': [],
            'returns': None
        }

        # Extract arguments
        for arg in node.args.args:
            arg_name = arg.arg
            arg_type = self._get_annotation(arg.annotation)
            func_info['args'].append({'name': arg_name, 'type': arg_type})

        # Extract return type
        if node.returns:
            func_info['returns'] = self._get_annotation(node.returns)

        return func_info

    def _get_annotation(self, node) -> Optional[str]:
        """Get type annotation as string."""
        if node is None:
            return None
        if isinstance(node, ast.Name):
            return node.id
        if isinstance(node, ast.Constant):
            return str(node.value)
        if isinstance(node, ast.Subscript):
            value = self._get_annotation(node.value)
            slice_val = self._get_annotation(node.slice)
            return f"{value}[{slice_val}]"
        return ast.unparse(node)

    def _get_name(self, node) -> str:
        """Get name from node."""
        if isinstance(node, ast.Name):
            return node.id
        if isinstance(node, ast.Attribute):
            return ast.unparse(node)
        return str(node)


def extract_docstrings(file_path: Path) -> DocstringExtractor:
    """Extract all docstrings from a Python file."""
    with open(file_path, 'r') as f:
        source = f.read()

    try:
        tree = ast.parse(source)
        extractor = DocstringExtractor(file_path)
        extractor.visit(tree)
        return extractor
    except SyntaxError as e:
        print(f"Error parsing {file_path}: {e}")
        return None


def generate_markdown(extractor: DocstringExtractor, project_name: str,
                     source_path: Path, api_path: Path) -> str:
    """Generate markdown documentation from extracted docstrings."""
    relative_path = source_path.relative_to(source_path.parent.parent.parent)

    # Create frontmatter
    md = f"""---
project: {project_name}
type: api-documentation
source: {relative_path}
generated: {datetime.now().strftime('%Y-%m-%d')}
tags: [api, documentation, {project_name}]
---

# {source_path.stem}

**Source**: [[../code-mirror/{relative_path}|{relative_path}]]

"""

    # Add module docstring
    if extractor.module_doc:
        md += f"{extractor.module_doc}\n\n"
    else:
        md += f"*Module: `{source_path.stem}`*\n\n"

    # Table of contents
    md += "## Contents\n\n"
    if extractor.classes:
        md += "### Classes\n"
        for cls in extractor.classes:
            md += f"- [{cls['name']}](#{cls['name'].lower()})\n"
        md += "\n"
    if extractor.functions:
        md += "### Functions\n"
        for func in extractor.functions:
            md += f"- [{func['name']}](#{func['name'].lower()})\n"
        md += "\n"

    # Document classes
    if extractor.classes:
        md += "## Classes\n\n"
        for cls in extractor.classes:
            md += f"### {cls['name']}\n\n"

            # Inheritance
            if cls['bases']:
                bases_str = ', '.join(cls['bases'])
                md += f"**Inherits from**: `{bases_str}`\n\n"

            # Source location
            md += f"**Source**: `{relative_path}:{cls['line']}`\n\n"

            # Docstring
            if cls['docstring']:
                md += f"{cls['docstring']}\n\n"

            # Methods
            if cls['methods']:
                md += "#### Methods\n\n"
                for method in cls['methods']:
                    md += f"##### `{method['name']}(`"

                    # Arguments
                    args = []
                    for arg in method['args']:
                        if arg['type']:
                            args.append(f"{arg['name']}: {arg['type']}")
                        else:
                            args.append(arg['name'])
                    md += ", ".join(args)
                    md += ")`"

                    # Return type
                    if method['returns']:
                        md += f" → `{method['returns']}`"

                    md += "\n\n"

                    # Method docstring
                    if method['docstring']:
                        md += f"{method['docstring']}\n\n"
                    else:
                        md += "*No documentation available.*\n\n"

            md += "---\n\n"

    # Document functions
    if extractor.functions:
        md += "## Functions\n\n"
        for func in extractor.functions:
            md += f"### `{func['name']}(`"

            # Arguments
            args = []
            for arg in func['args']:
                if arg['type']:
                    args.append(f"{arg['name']}: {arg['type']}")
                else:
                    args.append(arg['name'])
            md += ", ".join(args)
            md += ")`"

            # Return type
            if func['returns']:
                md += f" → `{func['returns']}`"

            md += "\n\n"

            # Source location
            md += f"**Source**: `{relative_path}:{func['line']}`\n\n"

            # Docstring
            if func['docstring']:
                md += f"{func['docstring']}\n\n"
            else:
                md += "*No documentation available.*\n\n"

            md += "---\n\n"

    # Links
    md += "## Links\n\n"
    md += f"- [[../code-mirror/{relative_path}|Source Code]]\n"
    md += f"- [[../../projects/{project_name}/|Project Root]]\n"
    md += f"- [[./README|API Index]]\n\n"

    return md


def generate_api_index(project_name: str, api_path: Path,
                       documented_files: List[Path]) -> str:
    """Generate API documentation index."""
    md = f"""---
project: {project_name}
type: api-index
generated: {datetime.now().strftime('%Y-%m-%d')}
tags: [api, documentation, index]
---

# {project_name.upper()} API Documentation

Auto-generated API documentation from Python docstrings.

**Last Updated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## Modules

"""

    # Group by category
    categories = {
        'algorithms': [],
        'agents': [],
        'core': [],
        'examples': [],
        'other': []
    }

    for file_path in documented_files:
        relative = file_path.relative_to(file_path.parent.parent.parent)
        parts = relative.parts

        if 'algorithms' in parts:
            categories['algorithms'].append(file_path)
        elif 'agents' in parts:
            categories['agents'].append(file_path)
        elif 'core' in parts:
            categories['core'].append(file_path)
        elif 'examples' in parts:
            categories['examples'].append(file_path)
        else:
            categories['other'].append(file_path)

    # Generate sections
    for category, files in categories.items():
        if not files:
            continue

        md += f"### {category.title()}\n\n"
        for file_path in sorted(files):
            module_name = file_path.stem
            relative = file_path.relative_to(file_path.parent.parent.parent)
            md += f"- [[./{module_name}|{module_name}]] - `{relative}`\n"
        md += "\n"

    # Statistics
    md += f"""## Statistics

- **Total Modules**: {len(documented_files)}
- **Algorithms**: {len(categories['algorithms'])}
- **Agents**: {len(categories['agents'])}
- **Core**: {len(categories['core'])}
- **Examples**: {len(categories['examples'])}

## Links

- [[../../projects/{project_name}/|Project Root]]
- [[../code-mirror/|Source Code Mirror]]
- [[../../daily/{datetime.now().strftime('%Y-%m-%d')}|Today's Log]]

---

*Auto-generated by `scripts/generate-api-docs.py`*
"""

    return md


def main():
    if len(sys.argv) < 2:
        print("Usage: python generate-api-docs.py <project-name>")
        sys.exit(1)

    project_name = sys.argv[1]

    # Paths
    script_dir = Path(__file__).parent
    projects_root = script_dir.parent
    project_path = projects_root / "projects" / project_name
    source_path = project_path.parent.parent / project_name / "src"
    examples_path = project_path.parent.parent / project_name / "examples"
    api_path = project_path / "api"

    # Create API directory
    api_path.mkdir(parents=True, exist_ok=True)

    print(f"Generating API documentation for {project_name}...")
    print(f"Source: {source_path}")
    print(f"Output: {api_path}")

    documented_files = []

    # Process source files
    if source_path.exists():
        for py_file in source_path.rglob("*.py"):
            if py_file.name.startswith("_") and py_file.name != "__init__.py":
                continue

            print(f"Processing {py_file.relative_to(source_path)}...")

            extractor = extract_docstrings(py_file)
            if extractor:
                markdown = generate_markdown(extractor, project_name, py_file, api_path)

                # Write markdown file
                md_filename = py_file.stem + ".md"
                md_path = api_path / md_filename
                with open(md_path, 'w') as f:
                    f.write(markdown)

                documented_files.append(py_file)
                print(f"  → {md_path.relative_to(projects_root)}")

    # Process examples
    if examples_path.exists():
        for py_file in examples_path.glob("*.py"):
            if py_file.name.startswith("_"):
                continue

            print(f"Processing {py_file.name}...")

            extractor = extract_docstrings(py_file)
            if extractor:
                markdown = generate_markdown(extractor, project_name, py_file, api_path)

                # Write markdown file
                md_filename = py_file.stem + ".md"
                md_path = api_path / md_filename
                with open(md_path, 'w') as f:
                    f.write(markdown)

                documented_files.append(py_file)
                print(f"  → {md_path.relative_to(projects_root)}")

    # Generate index
    index_md = generate_api_index(project_name, api_path, documented_files)
    index_path = api_path / "README.md"
    with open(index_path, 'w') as f:
        f.write(index_md)

    print(f"\n✓ Generated API documentation for {len(documented_files)} modules")
    print(f"✓ Index: {index_path.relative_to(projects_root)}")


if __name__ == "__main__":
    main()
