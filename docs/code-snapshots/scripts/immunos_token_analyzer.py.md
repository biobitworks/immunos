---
source: /Users/byron/projects/scripts/immunos_token_analyzer.py
relative: scripts/immunos_token_analyzer.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
IMMUNOS Token Analyzer

Analyze token usage across ~/projects/ using tiktoken (Claude's tokenizer).
Identify token hotspots, calculate information density, and suggest optimizations.

Usage:
    python immunos_token_analyzer.py [--top 50] [--output token_report.json]
"""

import os
import sys
import json
import argparse
import tiktoken
from pathlib import Path
from typing import Dict, List, Tuple
from datetime import datetime
from collections import defaultdict


# Directories to skip
SKIP_DIRS = {
    '.git', '.svn', '.hg',
    'node_modules', '__pycache__', '.pytest_cache',
    '.venv', 'venv', 'env',
    '.immunos',
}

# File extensions to analyze
ANALYZE_EXTENSIONS = {
    '.py', '.js', '.ts', '.jsx', '.tsx',
    '.md', '.txt', '.rst',
    '.json', '.yaml', '.yml',
    '.sh', '.bash',
    '.sql', '.html', '.css',
    '.java', '.go', '.rs', '.cpp', '.c', '.h',
}


class TokenAnalyzer:
    """Analyze token usage using tiktoken"""

    def __init__(self):
        # Claude uses cl100k_base encoding (same as GPT-4)
        try:
            self.encoder = tiktoken.get_encoding("cl100k_base")
        except Exception as e:
            print(f"Error loading tiktoken encoder: {e}")
            print("Falling back to approximation (1 token ≈ 4 chars)")
            self.encoder = None

        self.file_tokens = {}
        self.total_tokens = 0
        self.total_files = 0
        self.by_extension = defaultdict(lambda: {'count': 0, 'tokens': 0})
        self.by_directory = defaultdict(lambda: {'count': 0, 'tokens': 0})

    def count_tokens(self, text: str) -> int:
        """Count tokens in text"""
        if self.encoder:
            try:
                return len(self.encoder.encode(text))
            except Exception:
                # Fallback to approximation
                return len(text) // 4
        else:
            # Approximation: 1 token ≈ 4 characters
            return len(text) // 4

    def analyze_file(self, file_path: Path, relative_path: str) -> Dict:
        """Analyze token usage in a single file"""
        try:
            # Read file content
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()

            # Count tokens
            tokens = self.count_tokens(content)
            chars = len(content)
            lines = len(content.splitlines())

            # Calculate density (tokens per character)
            density = tokens / chars if chars > 0 else 0

            # Analyze patterns
            issues = []

            # Check for excessive whitespace
            blank_lines = content.count('\n\n\n')  # 3+ consecutive newlines
            if blank_lines > 10:
                issues.append({
                    'type': 'whitespace',
                    'severity': 'low',
                    'description': f'{blank_lines} instances of excessive blank lines',
                    'savings_estimate': blank_lines * 2,  # Rough estimate
                })

            # Check for embedded data (Python)
            if file_path.suffix == '.py':
                # Count dict/list literals >20 lines
                import re
                data_blocks = re.findall(r'[\[\{][^\[\{\]\}]{500,}[\]\}]', content, re.DOTALL)
                if data_blocks:
                    data_tokens = sum(self.count_tokens(block) for block in data_blocks)
                    issues.append({
                        'type': 'embedded_data',
                        'severity': 'medium',
                        'description': f'{len(data_blocks)} large data structures embedded in code',
                        'savings_estimate': int(data_tokens * 0.9),  # 90% could be externalized
                    })

            # Check for verbose documentation
            if file_path.suffix == '.md':
                # If file is >2000 tokens, suggest adding summary
                if tokens > 2000:
                    issues.append({
                        'type': 'verbose_doc',
                        'severity': 'medium',
                        'description': f'Large documentation file ({tokens} tokens)',
                        'savings_estimate': int(tokens * 0.6),  # Could save 60% with summary
                    })

            return {
                'file': relative_path,
                'tokens': tokens,
                'chars': chars,
                'lines': lines,
                'density': density,
                'extension': file_path.suffix,
                'issues': issues,
                'optimization_potential': sum(i['savings_estimate'] for i in issues),
            }

        except Exception as e:
            return {
                'file': relative_path,
                'error': str(e),
                'tokens': 0,
            }

    def should_analyze(self, file_path: Path) -> bool:
        """Check if file should be analyzed"""
        # Skip by extension
        if file_path.suffix not in ANALYZE_EXTENSIONS:
            return False

        # Skip hidden files
        if file_path.name.startswith('.'):
            return False

        # Skip very large files (>1MB)
        try:
            if file_path.stat().st_size > 1024 * 1024:
                return False
        except OSError:
            return False

        return True

    def analyze_directory(self, projects_root: Path):
        """Analyze all files in directory tree"""
        print(f"Starting token analysis of {projects_root}")
        print(f"Using tiktoken encoder: {self.encoder is not None}")
        print()

        file_count = 0

        for root, dirs, files in os.walk(projects_root):
            # Filter out skip directories
            dirs[:] = [d for d in dirs if d not in SKIP_DIRS]

            root_path = Path(root)
            relative_root = root_path.relative_to(projects_root)

            for file_name in files:
                file_path = root_path / file_name

                # Skip unwanted files
                if not self.should_analyze(file_path):
                    continue

                file_count += 1
                if file_count % 50 == 0:
                    print(f"  Analyzed {file_count} files, {self.total_tokens:,} tokens...", end='\r')

                # Analyze file
                relative_path = str(relative_root / file_name)
                result = self.analyze_file(file_path, relative_path)

                if 'error' not in result:
                    self.file_tokens[relative_path] = result
                    self.total_tokens += result['tokens']
                    self.total_files += 1

                    # Update stats by extension
                    ext = result['extension']
                    self.by_extension[ext]['count'] += 1
                    self.by_extension[ext]['tokens'] += result['tokens']

                    # Update stats by directory
                    dir_name = str(relative_root).split('/')[0] if relative_root.parts else 'root'
                    self.by_directory[dir_name]['count'] += 1
                    self.by_directory[dir_name]['tokens'] += result['tokens']

        print(f"\n\nToken analysis complete!")
        print(f"  Files analyzed: {self.total_files:,}")
        print(f"  Total tokens: {self.total_tokens:,}")
        print(f"  Average tokens/file: {self.total_tokens // self.total_files if self.total_files else 0:,}")
        print()

    def generate_report(self, top_n: int = 50) -> Dict:
        """Generate comprehensive token report"""
        # Sort files by token count
        sorted_files = sorted(
            self.file_tokens.items(),
            key=lambda x: x[1]['tokens'],
            reverse=True
        )

        # Top token consumers
        top_files = [
            {
                'file': f,
                'tokens': data['tokens'],
                'lines': data['lines'],
                'density': data['density'],
                'optimization_potential': data.get('optimization_potential', 0),
                'issues': data.get('issues', []),
            }
            for f, data in sorted_files[:top_n]
        ]

        # Files with highest optimization potential
        optimizable = sorted(
            [
                {
                    'file': f,
                    'tokens': data['tokens'],
                    'savings': data.get('optimization_potential', 0),
                    'savings_pct': (data.get('optimization_potential', 0) / data['tokens'] * 100) if data['tokens'] else 0,
                }
                for f, data in self.file_tokens.items()
                if data.get('optimization_potential', 0) > 100
            ],
            key=lambda x: x['savings'],
            reverse=True
        )

        # Token usage by extension
        extension_stats = [
            {
                'extension': ext,
                'files': stats['count'],
                'tokens': stats['tokens'],
                'avg_tokens': stats['tokens'] // stats['count'] if stats['count'] else 0,
                'pct_of_total': (stats['tokens'] / self.total_tokens * 100) if self.total_tokens else 0,
            }
            for ext, stats in sorted(
                self.by_extension.items(),
                key=lambda x: x[1]['tokens'],
                reverse=True
            )
        ]

        # Token usage by directory
        directory_stats = [
            {
                'directory': dir_name,
                'files': stats['count'],
                'tokens': stats['tokens'],
                'avg_tokens': stats['tokens'] // stats['count'] if stats['count'] else 0,
                'pct_of_total': (stats['tokens'] / self.total_tokens * 100) if self.total_tokens else 0,
            }
            for dir_name, stats in sorted(
                self.by_directory.items(),
                key=lambda x: x[1]['tokens'],
                reverse=True
            )
        ]

        # Calculate total optimization potential
        total_savings = sum(
            data.get('optimization_potential', 0)
            for data in self.file_tokens.values()
        )

        report = {
            'scan_date': datetime.now().isoformat(),
            'total_files': self.total_files,
            'total_tokens': self.total_tokens,
            'avg_tokens_per_file': self.total_tokens // self.total_files if self.total_files else 0,
            'total_optimization_potential': total_savings,
            'potential_savings_pct': (total_savings / self.total_tokens * 100) if self.total_tokens else 0,
            'top_token_consumers': top_files,
            'high_optimization_potential': optimizable[:30],
            'by_extension': extension_stats,
            'by_directory': directory_stats,
        }

        return report

    def print_summary(self, report: Dict):
        """Print summary to console"""
        print("="*70)
        print("TOKEN ANALYSIS SUMMARY")
        print("="*70)

        print(f"\n📊 Overall Statistics")
        print(f"  Total Files: {report['total_files']:,}")
        print(f"  Total Tokens: {report['total_tokens']:,}")
        print(f"  Avg Tokens/File: {report['avg_tokens_per_file']:,}")
        print(f"  Optimization Potential: {report['total_optimization_potential']:,} tokens ({report['potential_savings_pct']:.1f}%)")

        print(f"\n🔥 Top 10 Token Consumers")
        for i, item in enumerate(report['top_token_consumers'][:10], 1):
            opt = item.get('optimization_potential', 0)
            opt_str = f" [can save ~{opt:,} tokens]" if opt > 100 else ""
            print(f"  {i:2}. {item['file']}")
            print(f"      {item['tokens']:,} tokens ({item['lines']:,} lines){opt_str}")

        print(f"\n💡 Top 10 Optimization Opportunities")
        for i, item in enumerate(report['high_optimization_potential'][:10], 1):
            print(f"  {i:2}. {item['file']}")
            print(f"      {item['tokens']:,} tokens → save ~{item['savings']:,} tokens ({item['savings_pct']:.0f}%)")

        print(f"\n📁 Token Usage by File Type")
        for item in report['by_extension'][:10]:
            print(f"  {item['extension']:8} {item['files']:4} files  {item['tokens']:10,} tokens  ({item['pct_of_total']:4.1f}%)")

        print(f"\n📂 Token Usage by Directory")
        for item in report['by_directory'][:10]:
            print(f"  {item['directory']:30} {item['files']:4} files  {item['tokens']:10,} tokens  ({item['pct_of_total']:4.1f}%)")

        print("\n" + "="*70)


def main():
    parser = argparse.ArgumentParser(description='Analyze token usage across projects')
    parser.add_argument('--projects-root', type=str,
                        default='/Users/byron/projects',
                        help='Projects root directory')
    parser.add_argument('--top', type=int, default=50,
                        help='Number of top files to include in report')
    parser.add_argument('--output', type=str,
                        default='/Users/byron/projects/.immunos/token_baseline.json',
                        help='Output report file')

    args = parser.parse_args()

    projects_root = Path(args.projects_root)
    if not projects_root.exists():
        print(f"Error: Projects root not found: {projects_root}")
        sys.exit(1)

    # Create analyzer and run
    analyzer = TokenAnalyzer()
    analyzer.analyze_directory(projects_root)

    # Generate report
    report = analyzer.generate_report(top_n=args.top)

    # Save report
    output_path = Path(args.output)
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)

    print(f"✓ Report saved to {output_path} ({output_path.stat().st_size / 1024:.1f} KB)")
    print()

    # Print summary
    analyzer.print_summary(report)


if __name__ == '__main__':
    main()

```
