#!/usr/bin/env python3
"""
IMMUNOS NK Cell Scanner

Use Negative Selection algorithm with Ollama LLMs to detect anomalies:
- Security vulnerabilities (hardcoded secrets, permissions)
- Code quality issues (complexity, duplicates)
- Structural anomalies (naming, organization)
- Documentation issues (missing docs, outdated info)

Usage:
    python immunos_nk_scan.py [--baseline baseline.json] [--ollama-url http://localhost:11434]
"""

import os
import sys
import json
import re
import argparse
import requests
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional
from dataclasses import dataclass, asdict
from enum import Enum

# Import negative selection components
try:
    from immunos_negsel import (
        NegativeSelectionClassifier,
        CodeFeatureExtractor,
        ScannerFusion
    )
    NEGSEL_AVAILABLE = True
except ImportError:
    NEGSEL_AVAILABLE = False
    print("Warning: immunos_negsel module not found. Negative selection disabled.", file=sys.stderr)


class Severity(Enum):
    """Anomaly severity levels"""
    HIGH = "HIGH"
    MEDIUM = "MEDIUM"
    LOW = "LOW"


class Category(Enum):
    """Anomaly categories"""
    SECURITY = "security"
    CODE_QUALITY = "code_quality"
    STRUCTURAL = "structural"
    DOCUMENTATION = "documentation"
    TASK_MANAGEMENT = "task_management"


@dataclass
class Anomaly:
    """Detected anomaly"""
    category: str
    severity: str
    file_path: str
    line_number: Optional[int]
    pattern: str
    description: str
    recommendation: str
    ollama_analysis: Optional[str] = None


class NKCellScanner:
    """NK Cell scanner using negative selection with Ollama"""

    # Detection mode
    MODE_PATTERN = "pattern"      # Traditional regex patterns
    MODE_NEGSEL = "negsel"        # Negative selection ML
    MODE_HYBRID = "hybrid"        # Both methods combined

    # Security patterns to detect
    SECURITY_PATTERNS = {
        'api_key': r'(?i)(api[_-]?key|apikey)["\']?\s*[:=]\s*["\']([a-zA-Z0-9_\-]{20,})["\']',
        'password': r'(?i)(password|passwd|pwd)["\']?\s*[:=]\s*["\']([^"\']{8,})["\']',
        'secret': r'(?i)(secret|token)["\']?\s*[:=]\s*["\']([a-zA-Z0-9_\-]{20,})["\']',
        'aws_key': r'(?i)(aws[_-]?access|aws[_-]?secret)["\']?\s*[:=]\s*["\']([A-Z0-9]{20,})["\']',
        'private_key': r'-----BEGIN (RSA |DSA |EC )?PRIVATE KEY-----',
    }

    # Code quality patterns
    QUALITY_MARKERS = {
        'todo': r'(?i)#\s*(TODO|FIXME|XXX|HACK|BUG)',
        'noqa': r'#\s*noqa',
        'type_ignore': r'#\s*type:\s*ignore',
    }

    def __init__(self, ollama_url: str = "http://localhost:11434",
                 mode: str = "pattern", enable_negsel: bool = False):
        self.ollama_url = ollama_url
        self.anomalies: List[Anomaly] = []
        self.mode = mode if NEGSEL_AVAILABLE else self.MODE_PATTERN

        # Negative selection components
        self.negsel_classifier: Optional[NegativeSelectionClassifier] = None
        self.feature_extractor: Optional[CodeFeatureExtractor] = None
        self.scanner_fusion: Optional[ScannerFusion] = None

        # Initialize negative selection if enabled
        if enable_negsel and NEGSEL_AVAILABLE:
            self._initialize_negsel()

    def _initialize_negsel(self):
        """Initialize negative selection classifier"""
        try:
            # Try to load trained detector
            self.negsel_classifier = NegativeSelectionClassifier.load_from_db("SAFE")
            self.feature_extractor = CodeFeatureExtractor()
            self.scanner_fusion = ScannerFusion()
            print("✓ Loaded trained negative selection detector from database")
        except ValueError:
            # No trained detector found
            print("⚠ No trained detector found. Run detector training first.")
            print("  Use: python scripts/immunos_negsel.py train <safe_code_dir>")
            self.negsel_classifier = None
        except Exception as e:
            print(f"⚠ Failed to initialize negative selection: {e}", file=sys.stderr)
            self.negsel_classifier = None

    def call_ollama(self, model: str, prompt: str) -> Optional[str]:
        """Call Ollama API for analysis"""
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
                return response.json().get('response', '')
            else:
                print(f"  Warning: Ollama API error {response.status_code}", file=sys.stderr)
                return None
        except Exception as e:
            print(f"  Warning: Ollama call failed: {e}", file=sys.stderr)
            return None

    def check_security_patterns(self, file_path: str, content: str):
        """Check for security vulnerabilities using static patterns"""
        lines = content.splitlines()

        for pattern_name, pattern in self.SECURITY_PATTERNS.items():
            matches = re.finditer(pattern, content)
            for match in matches:
                # Find line number
                line_num = content[:match.start()].count('\n') + 1

                # Extract context
                context_line = lines[line_num - 1] if line_num <= len(lines) else ""

                # Don't flag examples or comments
                if '#' in context_line and context_line.strip().startswith('#'):
                    continue
                if 'example' in context_line.lower() or 'placeholder' in context_line.lower():
                    continue

                self.anomalies.append(Anomaly(
                    category=Category.SECURITY.value,
                    severity=Severity.HIGH.value,
                    file_path=file_path,
                    line_number=line_num,
                    pattern=pattern_name,
                    description=f"Hardcoded {pattern_name.replace('_', ' ')} detected",
                    recommendation=f"Move to environment variable or secrets manager",
                ))

    def check_file_permissions(self, file_path: Path, relative_path: str):
        """Check for suspicious file permissions"""
        try:
            stat = file_path.stat()
            perms = oct(stat.st_mode)[-3:]

            # World-writable files
            if perms[-1] in ['2', '3', '6', '7']:
                self.anomalies.append(Anomaly(
                    category=Category.SECURITY.value,
                    severity=Severity.MEDIUM.value,
                    file_path=relative_path,
                    line_number=None,
                    pattern="world_writable",
                    description=f"World-writable file: {perms}",
                    recommendation=f"Change permissions to 600 or 644: chmod 600 {file_path}",
                ))

            # Executable configs or data files
            if file_path.suffix in ['.json', '.yaml', '.yml', '.txt', '.md']:
                if perms[0] in ['7', '5', '3', '1']:  # Owner executable
                    self.anomalies.append(Anomaly(
                        category=Category.SECURITY.value,
                        severity=Severity.LOW.value,
                        file_path=relative_path,
                        line_number=None,
                        pattern="unexpected_executable",
                        description=f"Data file marked executable: {perms}",
                        recommendation=f"Remove execute permission: chmod -x {file_path}",
                    ))
        except Exception:
            pass

    def check_code_quality(self, file_path: str, content: str):
        """Check for code quality issues"""
        lines = content.splitlines()

        # Check function length (Python)
        if file_path.endswith('.py'):
            in_function = False
            function_start = 0
            function_name = ""

            for i, line in enumerate(lines, 1):
                # Detect function start
                if re.match(r'\s*def\s+(\w+)', line):
                    if in_function and (i - function_start) > 100:
                        # Previous function was too long
                        self.anomalies.append(Anomaly(
                            category=Category.CODE_QUALITY.value,
                            severity=Severity.MEDIUM.value,
                            file_path=file_path,
                            line_number=function_start,
                            pattern="long_function",
                            description=f"Function '{function_name}' is {i - function_start} lines long",
                            recommendation="Consider refactoring into smaller functions",
                        ))

                    in_function = True
                    function_start = i
                    match = re.match(r'\s*def\s+(\w+)', line)
                    function_name = match.group(1) if match else "unknown"

                # Detect class or end of file
                elif line.startswith('class ') or i == len(lines):
                    if in_function and (i - function_start) > 100:
                        self.anomalies.append(Anomaly(
                            category=Category.CODE_QUALITY.value,
                            severity=Severity.MEDIUM.value,
                            file_path=file_path,
                            line_number=function_start,
                            pattern="long_function",
                            description=f"Function '{function_name}' is {i - function_start} lines long",
                            recommendation="Consider refactoring into smaller functions",
                        ))
                    in_function = False

        # Check for TODO/FIXME markers
        for pattern_name, pattern in self.QUALITY_MARKERS.items():
            matches = re.finditer(pattern, content)
            for match in matches:
                line_num = content[:match.start()].count('\n') + 1
                context = lines[line_num - 1] if line_num <= len(lines) else ""

                if pattern_name == 'todo':
                    self.anomalies.append(Anomaly(
                        category=Category.CODE_QUALITY.value,
                        severity=Severity.LOW.value,
                        file_path=file_path,
                        line_number=line_num,
                        pattern="todo_marker",
                        description=f"TODO/FIXME marker: {context.strip()}",
                        recommendation="Address or track as issue",
                    ))

    def check_documentation(self, file_path: str, content: str):
        """Check for documentation issues"""
        # Check for Python files without docstrings
        if file_path.endswith('.py'):
            # Check if file has module docstring
            lines = content.strip().splitlines()
            if lines and not (lines[0].startswith('"""') or lines[0].startswith("'''")):
                # Check if second or third line has docstring (after shebang/encoding)
                has_docstring = False
                for i in range(min(5, len(lines))):
                    if lines[i].startswith('"""') or lines[i].startswith("'''"):
                        has_docstring = True
                        break

                if not has_docstring and len(content) > 200:  # Only flag non-trivial files
                    self.anomalies.append(Anomaly(
                        category=Category.DOCUMENTATION.value,
                        severity=Severity.LOW.value,
                        file_path=file_path,
                        line_number=1,
                        pattern="missing_docstring",
                        description="Python file missing module docstring",
                        recommendation="Add docstring describing module purpose",
                    ))

    def check_structural(self, file_path: str, baseline: Dict):
        """Check for structural anomalies"""
        # Check naming conventions
        path_parts = file_path.split('/')

        # Check for mixed case in directory names (should be lowercase with hyphens)
        for part in path_parts[:-1]:  # All except filename
            if part and re.search(r'[A-Z]', part) and '_' not in part:
                # Has uppercase but not using snake_case
                self.anomalies.append(Anomaly(
                    category=Category.STRUCTURAL.value,
                    severity=Severity.LOW.value,
                    file_path=file_path,
                    line_number=None,
                    pattern="inconsistent_naming",
                    description=f"Directory '{part}' uses mixed case (prefer lowercase-with-hyphens)",
                    recommendation=f"Rename to: {part.lower().replace('_', '-')}",
                ))

    def check_todo_anomalies(self, baseline: Dict):
        """Check for todo management anomalies"""
        from datetime import datetime

        todo_path = Path('/Users/byron/projects/todo/.index/todos.json')

        if not todo_path.exists():
            return  # Todo system not initialized

        try:
            with open(todo_path, 'r') as f:
                todo_index = json.load(f)

            todos = todo_index.get('todos', {})
            today = datetime.now().date()

            # Detect overdue todos
            for todo_id, todo_data in todos.items():
                if todo_data['status'] == 'done':
                    continue

                if todo_data.get('due'):
                    due_date_str = todo_data['due']
                    # Handle ISO format
                    if 'T' in due_date_str:
                        due_date = datetime.fromisoformat(due_date_str.replace('Z', '+00:00')).date()
                    else:
                        due_date = datetime.strptime(due_date_str, '%Y-%m-%d').date()

                    days_overdue = (today - due_date).days

                    if days_overdue > 7:
                        # Severely overdue (HIGH)
                        self.anomalies.append(Anomaly(
                            category=Category.TASK_MANAGEMENT.value,
                            severity=Severity.HIGH.value,
                            file_path=f"todo/{todo_data.get('file', 'unknown')}",
                            line_number=None,
                            pattern="overdue_task",
                            description=f"Todo '{todo_data['title']}' is {days_overdue} days overdue",
                            recommendation=f"Complete or reschedule: python scripts/immunos_todo.py complete {todo_id}",
                        ))
                    elif days_overdue > 0:
                        # Recently overdue (MEDIUM)
                        self.anomalies.append(Anomaly(
                            category=Category.TASK_MANAGEMENT.value,
                            severity=Severity.MEDIUM.value,
                            file_path=f"todo/{todo_data.get('file', 'unknown')}",
                            line_number=None,
                            pattern="overdue_task",
                            description=f"Todo '{todo_data['title']}' is {days_overdue} days overdue",
                            recommendation=f"Address urgently: python scripts/immunos_todo.py complete {todo_id}",
                        ))

            # Detect todos without due dates
            no_due_date = [t for t in todos.values()
                           if t['status'] in ['next', 'waiting'] and not t.get('due')]

            if len(no_due_date) > 10:
                self.anomalies.append(Anomaly(
                    category=Category.TASK_MANAGEMENT.value,
                    severity=Severity.LOW.value,
                    file_path="todo/.index/todos.json",
                    line_number=None,
                    pattern="no_due_date",
                    description=f"{len(no_due_date)} active todos without due dates",
                    recommendation="Add due dates: python scripts/immunos_todo.py move TODO-ID next --due YYYY-MM-DD",
                ))

            # Detect stale waiting todos (>30 days)
            for todo_id, todo_data in todos.items():
                if todo_data['status'] == 'waiting':
                    created_str = todo_data['created']
                    if 'T' in created_str:
                        created = datetime.fromisoformat(created_str.replace('Z', '+00:00')).date()
                    else:
                        created = datetime.strptime(created_str, '%Y-%m-%d').date()

                    days_waiting = (today - created).days

                    if days_waiting > 30:
                        self.anomalies.append(Anomaly(
                            category=Category.TASK_MANAGEMENT.value,
                            severity=Severity.MEDIUM.value,
                            file_path=f"todo/{todo_data.get('file', 'unknown')}",
                            line_number=None,
                            pattern="stale_waiting",
                            description=f"Todo '{todo_data['title']}' waiting {days_waiting} days",
                            recommendation=f"Follow up or move to someday: python scripts/immunos_todo.py move {todo_id} someday",
                        ))

        except Exception as e:
            # Silently skip if todo system has issues
            pass

    def analyze_with_ollama(self, file_path: str, content: str, model: str) -> Optional[str]:
        """Use Ollama for semantic analysis"""
        # Only analyze code files, skip if too large
        if len(content) > 5000 or not file_path.endswith(('.py', '.js', '.ts', '.sh')):
            return None

        prompt = f"""Analyze this code file for potential issues. Focus on:
1. Security vulnerabilities (injection, unsafe operations)
2. Code quality problems (complexity, anti-patterns)
3. Best practice violations

File: {file_path}

Code:
{content[:2000]}...

Reply with a brief JSON array of issues found, each with: {{"severity": "high|medium|low", "issue": "description", "recommendation": "fix"}}
If no issues, reply with empty array []."""

        return self.call_ollama(model, prompt)

    def check_negsel(self, file_path: str) -> Optional[tuple]:
        """
        Check file using negative selection algorithm.

        Returns:
            (classification, confidence, severity) or None
        """
        if not self.negsel_classifier or not self.feature_extractor:
            return None

        try:
            # Extract features
            features = self.feature_extractor.extract_from_file(file_path)

            # Classify
            classification, confidence = self.negsel_classifier.predict_single(features)

            # Determine severity based on classification and confidence
            if classification == "non-self":
                if confidence > 0.8:
                    severity = Severity.HIGH
                elif confidence > 0.5:
                    severity = Severity.MEDIUM
                else:
                    severity = Severity.LOW

                return (classification, confidence, severity)
            else:
                # It's "self" (safe) - no anomaly
                return None

        except Exception as e:
            print(f"  NegSel error on {file_path}: {e}", file=sys.stderr)
            return None

    def scan_file(self, file_path: Path, relative_path: str, baseline: Dict):
        """Scan a single file for anomalies"""
        # Skip binary files
        if file_path.suffix in ['.pyc', '.so', '.dylib', '.db', '.sqlite', '.pdf', '.png', '.jpg']:
            return

        # Skip large files
        try:
            if file_path.stat().st_size > 100000:  # 100KB limit for detailed analysis
                return
        except:
            return

        # Read file content
        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
        except Exception:
            return

        # Pattern-based checks (always run)
        if self.mode in [self.MODE_PATTERN, self.MODE_HYBRID]:
            self.check_security_patterns(relative_path, content)
            self.check_file_permissions(file_path, relative_path)
            self.check_code_quality(relative_path, content)
            self.check_documentation(relative_path, content)
            self.check_structural(relative_path, baseline)

        # Negative selection check (if enabled)
        if self.mode in [self.MODE_NEGSEL, self.MODE_HYBRID] and self.negsel_classifier:
            negsel_result = self.check_negsel(str(file_path))

            if negsel_result:
                classification, confidence, severity = negsel_result

                # Add anomaly
                self.anomalies.append(Anomaly(
                    category=Category.SECURITY.value,
                    severity=severity.value,
                    file_path=relative_path,
                    line_number=None,
                    pattern="negsel_detection",
                    description=f"Negative selection classified as '{classification}' (confidence: {confidence:.2f})",
                    recommendation=f"Review file for anomalous patterns. ML detector flagged as non-self.",
                ))

    def scan_directory(self, projects_root: Path, baseline: Dict, use_ollama: bool = False):
        """Scan entire directory tree"""
        print(f"NK Cell scan starting...")
        print(f"Detection mode: {self.mode}")
        print(f"Ollama integration: {'enabled' if use_ollama else 'disabled (use --ollama to enable)'}")
        if self.negsel_classifier:
            print(f"Negative selection: ✓ Active ({len(self.negsel_classifier.valid_detectors)} detectors loaded)")
        else:
            print(f"Negative selection: ✗ Not available")
        print()

        # Check todo system first
        print("📋 Scanning todo management...")
        self.check_todo_anomalies(baseline)

        file_count = 0
        skip_dirs = {'.git', 'node_modules', '__pycache__', '.venv', '.immunos'}

        for root, dirs, files in os.walk(projects_root):
            dirs[:] = [d for d in dirs if d not in skip_dirs]

            for file_name in files:
                file_path = Path(root) / file_name
                relative_path = str(file_path.relative_to(projects_root))

                # Skip hidden files
                if file_name.startswith('.'):
                    continue

                file_count += 1
                if file_count % 50 == 0:
                    print(f"  Scanned {file_count} files, found {len(self.anomalies)} anomalies...", end='\r')

                self.scan_file(file_path, relative_path, baseline)

        print(f"\n\nNK Cell scan complete!")
        print(f"  Files scanned: {file_count}")
        print(f"  Anomalies detected: {len(self.anomalies)}")
        print()

    def generate_report(self) -> Dict:
        """Generate anomaly report"""
        # Categorize anomalies
        by_category = {}
        by_severity = {}

        for anomaly in self.anomalies:
            by_category.setdefault(anomaly.category, []).append(anomaly)
            by_severity.setdefault(anomaly.severity, []).append(anomaly)

        report = {
            'total_anomalies': len(self.anomalies),
            'by_category': {cat: len(anoms) for cat, anoms in by_category.items()},
            'by_severity': {sev: len(anoms) for sev, anoms in by_severity.items()},
            'anomalies': [asdict(a) for a in self.anomalies],
        }

        return report

    def print_summary(self):
        """Print summary to console"""
        by_category = {}
        by_severity = {}

        for anomaly in self.anomalies:
            by_category.setdefault(anomaly.category, []).append(anomaly)
            by_severity.setdefault(anomaly.severity, []).append(anomaly)

        print("="*60)
        print("ANOMALY SUMMARY")
        print("="*60)

        print("\nBy Category:")
        for category, anomalies in sorted(by_category.items()):
            print(f"  {category}: {len(anomalies)}")

        print("\nBy Severity:")
        for severity, anomalies in sorted(by_severity.items(), reverse=True):
            print(f"  {severity}: {len(anomalies)}")

        print("\nTop 10 Anomalies:")
        for i, anomaly in enumerate(self.anomalies[:10], 1):
            loc = f"{anomaly.file_path}:{anomaly.line_number}" if anomaly.line_number else anomaly.file_path
            print(f"  {i}. [{anomaly.severity}] {anomaly.description}")
            print(f"     {loc}")
            print(f"     → {anomaly.recommendation}")
            print()


def main():
    parser = argparse.ArgumentParser(description='NK Cell anomaly detection')
    parser.add_argument('--baseline', type=str,
                        default='/Users/byron/projects/.immunos/baseline.json',
                        help='Baseline snapshot file')
    parser.add_argument('--output', type=str,
                        default='/Users/byron/projects/.immunos/nk_scan.json',
                        help='Output report file')
    parser.add_argument('--mode', type=str,
                        choices=['pattern', 'negsel', 'hybrid'],
                        default='pattern',
                        help='Detection mode: pattern (regex), negsel (ML), or hybrid (both)')
    parser.add_argument('--negsel', action='store_true',
                        help='Enable negative selection ML detection')
    parser.add_argument('--ollama', action='store_true',
                        help='Enable Ollama semantic analysis (slower)')
    parser.add_argument('--ollama-url', type=str,
                        default='http://localhost:11434',
                        help='Ollama API URL')
    parser.add_argument('--projects-root', type=str,
                        default='/Users/byron/projects',
                        help='Projects root directory')

    args = parser.parse_args()

    # Determine mode
    if args.negsel or args.mode != 'pattern':
        detection_mode = args.mode
        enable_negsel = True
    else:
        detection_mode = 'pattern'
        enable_negsel = False

    # Load baseline
    baseline_path = Path(args.baseline)
    if not baseline_path.exists():
        print(f"Error: Baseline not found: {baseline_path}")
        print("Run immunos_baseline.py first")
        sys.exit(1)

    with open(baseline_path, 'r') as f:
        baseline = json.load(f)

    print(f"Loaded baseline: {baseline['file_count']} files from {baseline['scan_date']}")
    print()

    # Create scanner and run
    scanner = NKCellScanner(
        ollama_url=args.ollama_url,
        mode=detection_mode,
        enable_negsel=enable_negsel
    )
    scanner.scan_directory(Path(args.projects_root), baseline, use_ollama=args.ollama)

    # Generate report
    report = scanner.generate_report()

    # Save report
    output_path = Path(args.output)
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)

    print(f"✓ Report saved to {output_path}")
    print()

    # Print summary
    scanner.print_summary()


if __name__ == '__main__':
    main()
