#!/usr/bin/env python3
"""
IMMUNOS Quick Document Review

Quickly review any file or directory with AI-powered summaries.

Usage:
    python immunos_review.py <file_or_directory>
    python immunos_review.py --quick papers/
    python immunos_review.py biobitworks-logistics/README.md
"""

import os
import sys
import json
import argparse
import requests
from pathlib import Path
from typing import Optional, Dict


class QuickReviewer:
    """Quick document review with Ollama"""

    def __init__(self, ollama_url: str = "http://localhost:11434"):
        self.ollama_url = ollama_url

    def call_ollama(self, model: str, prompt: str) -> Optional[str]:
        """Call Ollama API"""
        try:
            response = requests.post(
                f"{self.ollama_url}/api/generate",
                json={
                    "model": model,
                    "prompt": prompt,
                    "stream": False,
                },
                timeout=30
            )
            if response.status_code == 200:
                return response.json().get('response', '').strip()
            else:
                return f"Error: API returned {response.status_code}"
        except Exception as e:
            return f"Error: {str(e)}"

    def get_file_stats(self, file_path: Path) -> Dict:
        """Get basic file statistics"""
        try:
            stat = file_path.stat()
            content = file_path.read_text(encoding='utf-8', errors='ignore')
            lines = content.splitlines()

            return {
                'size': stat.st_size,
                'lines': len(lines),
                'language': self.detect_language(file_path),
                'modified': stat.st_mtime,
            }
        except Exception as e:
            return {'error': str(e)}

    def detect_language(self, file_path: Path) -> str:
        """Detect file language"""
        suffix = file_path.suffix.lower()
        language_map = {
            '.py': 'Python',
            '.js': 'JavaScript',
            '.ts': 'TypeScript',
            '.sh': 'Shell',
            '.md': 'Markdown',
            '.json': 'JSON',
            '.yaml': 'YAML',
            '.yml': 'YAML',
        }
        return language_map.get(suffix, 'Unknown')

    def review_file(self, file_path: Path, quick: bool = False) -> str:
        """Review a single file"""
        # Get file stats
        stats = self.get_file_stats(file_path)

        # Read content
        try:
            content = file_path.read_text(encoding='utf-8', errors='ignore')
        except Exception as e:
            return f"Error reading file: {e}"

        # Truncate for review
        preview_content = content[:3000] if len(content) > 3000 else content

        # Generate summary with Ollama
        model = "qwen2.5:1.5b" if quick else "qwen2.5-coder:7b"

        prompt = f"""Summarize this {stats.get('language', 'code')} file in 2-3 sentences. Focus on:
1. What it does
2. Key functionality
3. Notable features or patterns

File: {file_path.name}

Content:
{preview_content}

Summary:"""

        summary = self.call_ollama(model, prompt)

        # Format output
        output = f"""# Quick Review: {file_path.name}

## Summary
{summary}

## Key Stats
- **Language**: {stats.get('language', 'Unknown')}
- **Lines**: {stats.get('lines', 0):,}
- **Size**: {stats.get('size', 0) / 1024:.1f} KB
- **Last Modified**: {self.format_timestamp(stats.get('modified', 0))}

## File Type
{stats.get('language', 'Unknown')}
"""

        # Check for common patterns
        if stats.get('language') == 'Python':
            # Count classes and functions
            import re
            classes = len(re.findall(r'^\s*class\s+\w+', content, re.MULTILINE))
            functions = len(re.findall(r'^\s*def\s+\w+', content, re.MULTILINE))
            imports = len(re.findall(r'^\s*import\s+', content, re.MULTILINE))

            output += f"""
## Code Structure
- **Classes**: {classes}
- **Functions**: {functions}
- **Import Statements**: {imports}
"""

        return output

    def review_directory(self, dir_path: Path, quick: bool = False) -> str:
        """Review a directory"""
        # Count files by type
        file_counts = {}
        total_size = 0
        file_list = []

        for item in dir_path.rglob('*'):
            if item.is_file() and not item.name.startswith('.'):
                suffix = item.suffix or 'no-extension'
                file_counts[suffix] = file_counts.get(suffix, 0) + 1
                total_size += item.stat().st_size

                # Collect interesting files
                if suffix in ['.py', '.js', '.md', '.json', '.yaml', '.sh']:
                    file_list.append(item.relative_to(dir_path))

        # Generate summary with Ollama
        model = "qwen2.5:1.5b"

        # List first 10 files
        file_sample = '\n'.join(str(f) for f in file_list[:10])

        prompt = f"""Describe what this directory contains based on the file list. In 2-3 sentences, what is this project/directory for?

Directory: {dir_path.name}
Total files: {sum(file_counts.values())}

Sample files:
{file_sample}

Description:"""

        summary = self.call_ollama(model, prompt)

        # Format output
        output = f"""# Quick Review: {dir_path.name}/

## Summary
{summary}

## Directory Stats
- **Total Files**: {sum(file_counts.values()):,}
- **Total Size**: {total_size / (1024*1024):.1f} MB
- **File Types**: {len(file_counts)}

## Files by Type
"""

        for suffix, count in sorted(file_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
            output += f"- **{suffix}**: {count} files\n"

        if len(file_list) > 0:
            output += f"""
## Key Files
"""
            for f in file_list[:10]:
                output += f"- `{f}`\n"

        return output

    def format_timestamp(self, timestamp: float) -> str:
        """Format Unix timestamp"""
        from datetime import datetime
        dt = datetime.fromtimestamp(timestamp)
        return dt.strftime("%Y-%m-%d %H:%M:%S")


def main():
    parser = argparse.ArgumentParser(description='Quick document review')
    parser.add_argument('path', type=str, help='File or directory to review')
    parser.add_argument('--quick', action='store_true',
                        help='Quick mode (faster, uses smaller model)')
    parser.add_argument('--ollama-url', type=str,
                        default='http://localhost:11434',
                        help='Ollama API URL')
    parser.add_argument('--output', type=str,
                        help='Output file (default: print to stdout)')

    args = parser.parse_args()

    path = Path(args.path)
    if not path.exists():
        print(f"Error: Path not found: {path}", file=sys.stderr)
        sys.exit(1)

    # Create reviewer
    reviewer = QuickReviewer(ollama_url=args.ollama_url)

    # Review
    print(f"Reviewing {path}...")
    if path.is_file():
        output = reviewer.review_file(path, quick=args.quick)
    else:
        output = reviewer.review_directory(path, quick=args.quick)

    # Output
    if args.output:
        output_path = Path(args.output)
        output_path.write_text(output)
        print(f"✓ Review saved to {output_path}")
    else:
        print()
        print(output)


if __name__ == '__main__':
    main()
