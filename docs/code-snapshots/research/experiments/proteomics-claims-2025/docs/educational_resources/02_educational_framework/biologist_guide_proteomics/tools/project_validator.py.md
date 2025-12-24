---
source: /Users/byron/projects/research/experiments/proteomics-claims-2025/docs/educational_resources/02_educational_framework/biologist_guide_proteomics/tools/project_validator.py
relative: research/experiments/proteomics-claims-2025/docs/educational_resources/02_educational_framework/biologist_guide_proteomics/tools/project_validator.py
generated_at: 2025-12-23 10:28
---

````python
#!/usr/bin/env python3
"""
Project Validation Script for Biologist Guide to Proteomics

This script validates the integrity and completeness of the proteomics
analysis educational framework.
"""

import os
import sys
import re
import json
from pathlib import Path
from typing import Dict, List, Tuple, Any
import subprocess


class ProjectValidator:
    """Comprehensive project validation and quality assurance"""

    def __init__(self, project_root: str = None):
        """Initialize validator with project root directory"""
        if project_root is None:
            project_root = Path(__file__).parent.parent
        self.project_root = Path(project_root)
        self.issues = []
        self.warnings = []
        self.stats = {}

    def validate_directory_structure(self) -> Dict:
        """Validate expected directory structure exists"""
        print("🔍 Validating directory structure...")

        expected_dirs = [
            '00_getting_started',
            '01_background_knowledge',
            '02_software_setup',
            '03_data_understanding',
            '04_group1_analyses',
            '05_group2_analyses',
            '06_resources_and_support',
            '07_practice_exercises',
            'data',
            'tools',
            'tutorials'
        ]

        missing_dirs = []
        existing_dirs = []

        for dir_name in expected_dirs:
            dir_path = self.project_root / dir_name
            if dir_path.exists() and dir_path.is_dir():
                existing_dirs.append(dir_name)
            else:
                missing_dirs.append(dir_name)
                self.issues.append(f"Missing required directory: {dir_name}")

        return {
            'expected': len(expected_dirs),
            'existing': len(existing_dirs),
            'missing': missing_dirs,
            'found': existing_dirs
        }

    def validate_core_files(self) -> Dict:
        """Validate core project files exist"""
        print("📋 Validating core files...")

        required_files = [
            'README.md',
            'QUICKSTART.md',
            'PROJECT_SUMMARY.md',
            'INDEX.md',
            'requirements.txt',
            'SETUP_GUIDE.md',
            'MULTI_LEVEL_ANALYSIS_GUIDE.md'
        ]

        missing_files = []
        existing_files = []
        file_stats = {}

        for file_name in required_files:
            file_path = self.project_root / file_name
            if file_path.exists():
                existing_files.append(file_name)
                file_stats[file_name] = {
                    'size': file_path.stat().st_size,
                    'lines': self._count_lines(file_path)
                }
            else:
                missing_files.append(file_name)
                self.issues.append(f"Missing required file: {file_name}")

        return {
            'required': len(required_files),
            'existing': len(existing_files),
            'missing': missing_files,
            'stats': file_stats
        }

    def validate_markdown_files(self) -> Dict:
        """Validate markdown files for common issues"""
        print("📝 Validating markdown files...")

        md_files = list(self.project_root.rglob("*.md"))
        md_stats = {
            'total_files': len(md_files),
            'total_lines': 0,
            'broken_links': [],
            'missing_headers': [],
            'todo_items': []
        }

        for md_file in md_files:
            try:
                content = md_file.read_text(encoding='utf-8')
                lines = content.split('\n')
                md_stats['total_lines'] += len(lines)

                # Check for broken internal links
                for line_num, line in enumerate(lines, 1):
                    # Find markdown links
                    link_pattern = r'\[([^\]]+)\]\(([^)]+)\)'
                    for match in re.finditer(link_pattern, line):
                        link_text, link_url = match.groups()

                        # Check internal links (not URLs)
                        if not link_url.startswith(('http', 'mailto:', '#')):
                            # Resolve relative to current file
                            link_path = md_file.parent / link_url
                            if not link_path.exists():
                                md_stats['broken_links'].append({
                                    'file': str(md_file.relative_to(self.project_root)),
                                    'line': line_num,
                                    'link': link_url,
                                    'text': link_text
                                })

                    # Check for TODO/FIXME items
                    if any(keyword in line.upper() for keyword in ['TODO', 'FIXME', 'BUG', 'ISSUE']):
                        md_stats['todo_items'].append({
                            'file': str(md_file.relative_to(self.project_root)),
                            'line': line_num,
                            'content': line.strip()
                        })

                # Check for proper headers
                if not content.startswith('#'):
                    md_stats['missing_headers'].append(str(md_file.relative_to(self.project_root)))

            except UnicodeDecodeError:
                self.warnings.append(f"Could not read file (encoding issue): {md_file}")
            except Exception as e:
                self.warnings.append(f"Error processing {md_file}: {str(e)}")

        return md_stats

    def validate_python_files(self) -> Dict:
        """Validate Python files in tools directory"""
        print("🐍 Validating Python files...")

        py_files = list(self.project_root.rglob("*.py"))
        py_stats = {
            'total_files': len(py_files),
            'total_lines': 0,
            'syntax_errors': [],
            'import_errors': [],
            'missing_docstrings': []
        }

        for py_file in py_files:
            try:
                content = py_file.read_text(encoding='utf-8')
                lines = content.split('\n')
                py_stats['total_lines'] += len(lines)

                # Check syntax
                try:
                    compile(content, str(py_file), 'exec')
                except SyntaxError as e:
                    py_stats['syntax_errors'].append({
                        'file': str(py_file.relative_to(self.project_root)),
                        'line': e.lineno,
                        'error': str(e)
                    })

                # Check for missing module docstrings
                if not content.strip().startswith('"""') and not content.strip().startswith("'''"):
                    if not content.strip().startswith('#'):  # Allow shebang
                        py_stats['missing_docstrings'].append(str(py_file.relative_to(self.project_root)))

            except Exception as e:
                self.warnings.append(f"Error processing {py_file}: {str(e)}")

        return py_stats

    def validate_data_files(self) -> Dict:
        """Validate data directory and sample files"""
        print("📊 Validating data files...")

        data_dir = self.project_root / 'data'
        data_stats = {
            'data_dir_exists': data_dir.exists(),
            'csv_files': [],
            'other_files': [],
            'total_size': 0
        }

        if data_dir.exists():
            for file_path in data_dir.iterdir():
                if file_path.is_file():
                    file_info = {
                        'name': file_path.name,
                        'size': file_path.stat().st_size,
                        'extension': file_path.suffix
                    }
                    data_stats['total_size'] += file_info['size']

                    if file_path.suffix == '.csv':
                        data_stats['csv_files'].append(file_info)
                    else:
                        data_stats['other_files'].append(file_info)
        else:
            self.issues.append("Data directory does not exist")

        return data_stats

    def check_dependencies(self) -> Dict:
        """Check if required Python packages are available"""
        print("📦 Checking dependencies...")

        requirements_file = self.project_root / 'requirements.txt'
        dep_stats = {
            'requirements_exists': requirements_file.exists(),
            'packages_checked': 0,
            'packages_available': 0,
            'missing_packages': [],
            'import_errors': []
        }

        if requirements_file.exists():
            try:
                requirements = requirements_file.read_text().split('\n')
                requirements = [req.strip() for req in requirements if req.strip() and not req.strip().startswith('#')]

                for req in requirements:
                    dep_stats['packages_checked'] += 1
                    # Extract package name (before >= or ==)
                    package_name = re.split('[>=<!]', req)[0].strip()

                    try:
                        __import__(package_name)
                        dep_stats['packages_available'] += 1
                    except ImportError:
                        dep_stats['missing_packages'].append(package_name)
                        dep_stats['import_errors'].append(f"Could not import: {package_name}")

            except Exception as e:
                self.warnings.append(f"Error reading requirements.txt: {str(e)}")

        return dep_stats

    def generate_summary_report(self) -> str:
        """Generate comprehensive validation report"""
        print("\n📋 Generating validation report...")

        # Run all validations
        dir_results = self.validate_directory_structure()
        file_results = self.validate_core_files()
        md_results = self.validate_markdown_files()
        py_results = self.validate_python_files()
        data_results = self.validate_data_files()
        dep_results = self.check_dependencies()

        # Calculate overall health score
        total_checks = 0
        passed_checks = 0

        # Directory structure (weight: 20%)
        total_checks += dir_results['expected']
        passed_checks += dir_results['existing']

        # Core files (weight: 20%)
        total_checks += file_results['required']
        passed_checks += file_results['existing']

        # Dependencies (weight: 20%)
        if dep_results['packages_checked'] > 0:
            total_checks += dep_results['packages_checked']
            passed_checks += dep_results['packages_available']

        # Python syntax (weight: 20%)
        if py_results['total_files'] > 0:
            total_checks += py_results['total_files']
            passed_checks += (py_results['total_files'] - len(py_results['syntax_errors']))

        # Data files (weight: 10%)
        total_checks += 1  # Data directory exists
        passed_checks += (1 if data_results['data_dir_exists'] else 0)

        # Markdown quality (weight: 10%)
        if md_results['total_files'] > 0:
            total_checks += md_results['total_files']
            passed_checks += (md_results['total_files'] - len(md_results['missing_headers']))

        health_score = (passed_checks / total_checks * 100) if total_checks > 0 else 0

        # Generate report
        report = f"""
========================================
🔬 PROJECT VALIDATION REPORT
========================================

📊 OVERALL HEALTH SCORE: {health_score:.1f}%

🗂️  DIRECTORY STRUCTURE
   ✅ Found: {dir_results['existing']}/{dir_results['expected']} directories
   ❌ Missing: {len(dir_results['missing'])} directories

📋 CORE FILES
   ✅ Found: {file_results['existing']}/{file_results['required']} files
   📈 Total documentation: {sum(stats['lines'] for stats in file_results['stats'].values())} lines

📝 MARKDOWN FILES
   📄 Total files: {md_results['total_files']}
   📏 Total lines: {md_results['total_lines']:,}
   🔗 Broken links: {len(md_results['broken_links'])}
   ⚠️  TODO items: {len(md_results['todo_items'])}

🐍 PYTHON FILES
   📄 Total files: {py_results['total_files']}
   📏 Total lines: {py_results['total_lines']:,}
   ❌ Syntax errors: {len(py_results['syntax_errors'])}
   📚 Missing docstrings: {len(py_results['missing_docstrings'])}

📊 DATA FILES
   📁 Data directory: {'✅ Exists' if data_results['data_dir_exists'] else '❌ Missing'}
   📈 CSV files: {len(data_results['csv_files'])}
   💾 Total size: {data_results['total_size']:,} bytes

📦 DEPENDENCIES
   📋 Requirements file: {'✅ Exists' if dep_results['requirements_exists'] else '❌ Missing'}
   ✅ Available packages: {dep_results['packages_available']}/{dep_results['packages_checked']}
   ❌ Missing packages: {len(dep_results['missing_packages'])}

🚨 CRITICAL ISSUES: {len(self.issues)}
⚠️  WARNINGS: {len(self.warnings)}

"""

        # Add detailed issue reporting
        if self.issues:
            report += "\n🚨 CRITICAL ISSUES:\n"
            for i, issue in enumerate(self.issues, 1):
                report += f"   {i}. {issue}\n"

        if self.warnings:
            report += "\n⚠️  WARNINGS:\n"
            for i, warning in enumerate(self.warnings, 1):
                report += f"   {i}. {warning}\n"

        # Add broken links if any
        if md_results['broken_links']:
            report += f"\n🔗 BROKEN LINKS ({len(md_results['broken_links'])}):\n"
            for link in md_results['broken_links'][:10]:  # Show first 10
                report += f"   📄 {link['file']}:{link['line']} → {link['link']}\n"
            if len(md_results['broken_links']) > 10:
                report += f"   ... and {len(md_results['broken_links']) - 10} more\n"

        # Add TODO items summary
        if md_results['todo_items']:
            report += f"\n📝 TODO/FIXME ITEMS ({len(md_results['todo_items'])}):\n"
            for todo in md_results['todo_items'][:5]:  # Show first 5
                report += f"   📄 {todo['file']}:{todo['line']} → {todo['content'][:60]}...\n"
            if len(md_results['todo_items']) > 5:
                report += f"   ... and {len(md_results['todo_items']) - 5} more\n"

        # Add recommendations
        report += "\n💡 RECOMMENDATIONS:\n"

        if health_score >= 90:
            report += "   🎉 Excellent! Project is in great shape.\n"
        elif health_score >= 75:
            report += "   👍 Good project health. Address critical issues.\n"
        elif health_score >= 50:
            report += "   ⚠️  Moderate issues. Focus on critical fixes.\n"
        else:
            report += "   🚨 Significant issues need attention.\n"

        if dir_results['missing']:
            report += "   📁 Create missing directories for complete structure\n"

        if len(md_results['broken_links']) > 0:
            report += "   🔗 Fix broken internal links for better navigation\n"

        if dep_results['missing_packages']:
            report += "   📦 Install missing dependencies with: pip install -r requirements.txt\n"

        if len(md_results['todo_items']) > 10:
            report += "   📝 Address TODO items to improve completeness\n"

        report += "\n" + "="*40 + "\n"

        return report

    def _count_lines(self, file_path: Path) -> int:
        """Count lines in a text file"""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                return sum(1 for _ in f)
        except:
            return 0

    def run_validation(self, save_report: bool = True) -> str:
        """Run complete validation and return report"""
        print("🚀 Starting project validation...\n")

        report = self.generate_summary_report()

        if save_report:
            report_file = self.project_root / 'VALIDATION_REPORT.md'
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write(f"# Project Validation Report\n\n```\n{report}\n```\n")
            print(f"📄 Report saved to: {report_file}")

        return report


def main():
    """Main entry point for validation script"""
    import argparse

    parser = argparse.ArgumentParser(description='Validate proteomics analysis project')
    parser.add_argument('--project-root', help='Project root directory', default=None)
    parser.add_argument('--no-save', action='store_true', help='Do not save report file')
    parser.add_argument('--quiet', action='store_true', help='Minimal output')

    args = parser.parse_args()

    validator = ProjectValidator(args.project_root)

    if not args.quiet:
        print("🔬 Biologist Guide to Proteomics - Project Validator")
        print("="*50)

    report = validator.run_validation(save_report=not args.no_save)

    if not args.quiet:
        print(report)

    # Exit with error code if critical issues found
    if validator.issues:
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()
````
