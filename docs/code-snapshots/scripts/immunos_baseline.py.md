---
source: /Users/byron/projects/scripts/immunos_baseline.py
relative: scripts/immunos_baseline.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
IMMUNOS Baseline Scanner

Create a baseline snapshot of ~/projects/ for change detection.
Generates file tree with metadata, hashes, and AST summaries for code files.

Usage:
    python immunos_baseline.py [--full-scan] [--output baseline.json]
"""

import os
import sys
import json
import hashlib
import ast
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Set
import mimetypes

# Directories to skip
SKIP_DIRS = {
    '.git', '.svn', '.hg',
    'node_modules', '__pycache__', '.pytest_cache',
    '.venv', 'venv', 'env',
    '.DS_Store', '.obsidian',
    '.immunos',  # Skip our own data directory
}

# File extensions to skip (binaries)
SKIP_EXTENSIONS = {
    '.pyc', '.pyo', '.so', '.dylib', '.dll',
    '.exe', '.bin', '.obj', '.o',
    '.zip', '.tar', '.gz', '.bz2', '.7z',
    '.jpg', '.jpeg', '.png', '.gif', '.bmp', '.ico',
    '.mp3', '.mp4', '.avi', '.mov',
    '.db-wal', '.db-shm',  # SQLite temp files
}

# Large files to skip (>10MB)
MAX_FILE_SIZE = 10 * 1024 * 1024


class BaselineScanner:
    def __init__(self, projects_root: Path):
        self.projects_root = projects_root
        self.file_count = 0
        self.dir_count = 0
        self.skipped_count = 0
        self.total_size = 0
        self.file_signatures = {}
        self.file_tree = {}
        self.known_patterns = {
            'paper_directories': [],
            'metadata_files': [],
            'protocols': [],
            'scripts': [],
            'readme_files': [],
        }

    def should_skip_file(self, file_path: Path) -> bool:
        """Check if file should be skipped"""
        # Skip by extension
        if file_path.suffix in SKIP_EXTENSIONS:
            return True

        # Skip large files
        try:
            if file_path.stat().st_size > MAX_FILE_SIZE:
                return True
        except OSError:
            return True

        # Skip hidden files
        if file_path.name.startswith('.') and file_path.name not in {'.gitignore', '.env.example'}:
            return True

        return False

    def compute_file_hash(self, file_path: Path) -> str:
        """Compute SHA256 hash of file"""
        sha256 = hashlib.sha256()
        try:
            with open(file_path, 'rb') as f:
                for chunk in iter(lambda: f.read(8192), b''):
                    sha256.update(chunk)
            return sha256.hexdigest()
        except Exception as e:
            return f"error:{str(e)}"

    def extract_python_ast(self, file_path: Path) -> Optional[Dict]:
        """Extract AST summary from Python file"""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()

            tree = ast.parse(content)

            classes = [node.name for node in ast.walk(tree) if isinstance(node, ast.ClassDef)]
            functions = [node.name for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)]
            imports = []

            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imports.extend([alias.name for alias in node.names])
                elif isinstance(node, ast.ImportFrom):
                    if node.module:
                        imports.append(node.module)

            return {
                'classes': classes[:20],  # Limit to first 20
                'functions': functions[:20],
                'imports': list(set(imports[:30])),  # Unique, limit to 30
                'line_count': len(content.splitlines())
            }
        except Exception as e:
            return {'error': str(e)}

    def detect_language(self, file_path: Path) -> str:
        """Detect file language/type"""
        suffix = file_path.suffix.lower()
        language_map = {
            '.py': 'python',
            '.js': 'javascript',
            '.ts': 'typescript',
            '.sh': 'shell',
            '.bash': 'shell',
            '.md': 'markdown',
            '.json': 'json',
            '.yaml': 'yaml',
            '.yml': 'yaml',
            '.sql': 'sql',
            '.html': 'html',
            '.css': 'css',
            '.txt': 'text',
        }
        return language_map.get(suffix, 'unknown')

    def scan_file(self, file_path: Path, relative_path: str) -> Optional[Dict]:
        """Scan a single file and extract metadata"""
        try:
            stat = file_path.stat()

            # Basic metadata
            file_info = {
                'hash': self.compute_file_hash(file_path),
                'size': stat.st_size,
                'modified': datetime.fromtimestamp(stat.st_mtime).isoformat(),
                'permissions': oct(stat.st_mode)[-3:],
                'language': self.detect_language(file_path),
            }

            # Extract AST for Python files
            if file_path.suffix == '.py':
                file_info['ast_summary'] = self.extract_python_ast(file_path)

            self.total_size += stat.st_size
            self.file_count += 1

            return file_info

        except Exception as e:
            print(f"  Error scanning {relative_path}: {e}", file=sys.stderr)
            self.skipped_count += 1
            return None

    def detect_patterns(self):
        """Detect known patterns in file structure"""
        for rel_path in self.file_signatures.keys():
            # Paper directories (10.1038_*)
            if rel_path.startswith('papers/10.'):
                dir_name = rel_path.split('/')[1]
                if dir_name not in self.known_patterns['paper_directories']:
                    self.known_patterns['paper_directories'].append(dir_name)

            # Metadata files
            if rel_path.endswith('/metadata.json'):
                self.known_patterns['metadata_files'].append(rel_path)

            # Protocol files
            if '/protocols/' in rel_path and rel_path.endswith('.md'):
                self.known_patterns['protocols'].append(rel_path)

            # Scripts
            if '/scripts/' in rel_path and rel_path.endswith('.py'):
                self.known_patterns['scripts'].append(rel_path)

            # README files
            if 'README.md' in rel_path:
                self.known_patterns['readme_files'].append(rel_path)

    def build_tree(self):
        """Build hierarchical file tree"""
        tree = {}

        for rel_path in self.file_signatures.keys():
            parts = rel_path.split('/')
            current = tree

            for i, part in enumerate(parts[:-1]):
                if part not in current:
                    current[part] = {}
                current = current[part]

            # Leaf file
            current[parts[-1]] = {
                'type': 'file',
                'size': self.file_signatures[rel_path]['size'],
                'language': self.file_signatures[rel_path]['language']
            }

        self.file_tree = tree

    def scan(self) -> Dict:
        """Perform full baseline scan"""
        print(f"Starting baseline scan of {self.projects_root}")
        print(f"Skip dirs: {', '.join(SKIP_DIRS)}")
        print()

        # Walk through projects directory
        for root, dirs, files in os.walk(self.projects_root):
            # Filter out skip directories
            dirs[:] = [d for d in dirs if d not in SKIP_DIRS]

            self.dir_count += len(dirs)
            root_path = Path(root)
            relative_root = root_path.relative_to(self.projects_root)

            # Show progress
            if self.file_count % 100 == 0:
                print(f"  Scanned {self.file_count} files, {self.dir_count} dirs...", end='\r')

            # Scan files in current directory
            for file_name in files:
                file_path = root_path / file_name
                relative_path = str(relative_root / file_name)

                # Skip unwanted files
                if self.should_skip_file(file_path):
                    self.skipped_count += 1
                    continue

                # Scan file
                file_info = self.scan_file(file_path, relative_path)
                if file_info:
                    self.file_signatures[relative_path] = file_info

        print(f"\n\nScan complete!")
        print(f"  Files scanned: {self.file_count}")
        print(f"  Directories: {self.dir_count}")
        print(f"  Files skipped: {self.skipped_count}")
        print(f"  Total size: {self.total_size / (1024*1024):.1f} MB")

        # Detect patterns
        print("\nDetecting patterns...")
        self.detect_patterns()

        print(f"  Paper directories: {len(self.known_patterns['paper_directories'])}")
        print(f"  Protocol files: {len(self.known_patterns['protocols'])}")
        print(f"  Scripts: {len(self.known_patterns['scripts'])}")

        # Build tree
        print("\nBuilding file tree...")
        self.build_tree()

        # Create baseline data structure
        baseline = {
            'scan_date': datetime.now().isoformat(),
            'projects_root': str(self.projects_root),
            'file_count': self.file_count,
            'dir_count': self.dir_count,
            'total_size_mb': round(self.total_size / (1024*1024), 2),
            'file_signatures': self.file_signatures,
            'file_tree': self.file_tree,
            'known_patterns': self.known_patterns,
        }

        return baseline


def main():
    parser = argparse.ArgumentParser(description='Create baseline snapshot of projects')
    parser.add_argument('--full-scan', action='store_true',
                        help='Full scan including large files')
    parser.add_argument('--output', type=str,
                        default='/Users/byron/projects/.immunos/baseline.json',
                        help='Output file path')
    parser.add_argument('--projects-root', type=str,
                        default='/Users/byron/projects',
                        help='Projects root directory')

    args = parser.parse_args()

    projects_root = Path(args.projects_root)
    if not projects_root.exists():
        print(f"Error: Projects root not found: {projects_root}")
        sys.exit(1)

    # Create scanner and run
    scanner = BaselineScanner(projects_root)
    baseline = scanner.scan()

    # Save baseline
    output_path = Path(args.output)
    print(f"\nSaving baseline to {output_path}...")

    with open(output_path, 'w') as f:
        json.dump(baseline, f, indent=2)

    print(f"✓ Baseline saved ({output_path.stat().st_size / 1024:.1f} KB)")

    # Print summary
    print("\n" + "="*60)
    print("BASELINE SUMMARY")
    print("="*60)
    print(f"Scan Date: {baseline['scan_date']}")
    print(f"Files: {baseline['file_count']:,}")
    print(f"Directories: {baseline['dir_count']:,}")
    print(f"Total Size: {baseline['total_size_mb']:.1f} MB")
    print(f"Output: {output_path}")
    print("="*60)


if __name__ == '__main__':
    main()

```
