---
source: /Users/byron/projects/scripts/immunos_journal.py
relative: scripts/immunos_journal.py
generated_at: 2025-12-23 10:28
---

````python
#!/usr/bin/env python3
"""
IMMUNOS Daily Journal Generator

Generate daily journal entry with scan findings from baseline and NK Cell scans.

Usage:
    python immunos_journal.py [--date 2025-12-11]
"""

import json
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, List


class JournalGenerator:
    """Generate daily journal from scan results"""

    def __init__(self, immunos_dir: Path):
        self.immunos_dir = immunos_dir
        self.journal_dir = immunos_dir / "journal"
        self.journal_dir.mkdir(exist_ok=True)

    def load_baseline(self) -> Dict:
        """Load baseline snapshot"""
        baseline_path = self.immunos_dir / "baseline.json"
        with open(baseline_path, 'r') as f:
            return json.load(f)

    def load_nk_scan(self) -> Dict:
        """Load NK Cell scan results"""
        nk_scan_path = self.immunos_dir / "nk_scan.json"
        with open(nk_scan_path, 'r') as f:
            return json.load(f)

    def generate_todo_summary(self, date_str: str) -> str:
        """Generate todo summary for journal"""
        todo_path = Path('/Users/byron/projects/todo/.index/todos.json')

        if not todo_path.exists():
            return "## 📋 Todo Summary\n\nNo todos tracked.\n\n"

        try:
            with open(todo_path, 'r') as f:
                todo_index = json.load(f)

            todos = todo_index.get('todos', {})
            stats = todo_index.get('stats', {})
            today = datetime.strptime(date_str, '%Y-%m-%d').date()

            # Find todos completed today
            completed_today = []
            for todo_id, todo_data in todos.items():
                if todo_data['status'] == 'done' and todo_data.get('completed'):
                    completed_date = datetime.fromisoformat(todo_data['completed']).date()
                    if completed_date == today:
                        completed_today.append(todo_data)

            # Find todos created today
            created_today = []
            for todo_id, todo_data in todos.items():
                created_date = datetime.fromisoformat(todo_data['created']).date()
                if created_date == today:
                    created_today.append(todo_data)

            # Find overdue
            overdue = []
            for todo_id, todo_data in todos.items():
                if todo_data['status'] != 'done' and todo_data.get('due'):
                    due_date = datetime.fromisoformat(todo_data['due']).date()
                    if due_date < today:
                        days_overdue = (today - due_date).days
                        overdue.append({**todo_data, 'days_overdue': days_overdue})

            # Find due soon
            due_this_week = []
            for todo_id, todo_data in todos.items():
                if todo_data['status'] != 'done' and todo_data.get('due'):
                    due_date = datetime.fromisoformat(todo_data['due']).date()
                    days_until = (due_date - today).days
                    if 0 <= days_until <= 7:
                        due_this_week.append({**todo_data, 'days_until': days_until})

            # Generate markdown
            md = "## 📋 Todo Summary\n\n"
            md += "### Statistics\n"
            md += f"- **Total Active**: {stats.get('total', 0) - stats.get('by_status', {}).get('done', 0)}\n"
            md += f"- **Completed Today**: {len(completed_today)}\n"
            md += f"- **Created Today**: {len(created_today)}\n"
            md += f"- **Overdue**: {len(overdue)}\n"
            md += f"- **Due This Week**: {len(due_this_week)}\n"
            md += "\n"

            # Status breakdown
            md += "### By Status\n"
            for status, count in stats.get('by_status', {}).items():
                md += f"- **{status.title()}**: {count}\n"
            md += "\n"

            if completed_today:
                md += "### ✅ Completed Today\n"
                for todo in completed_today:
                    project_str = f" [{todo['project']}]" if todo.get('project') else ""
                    md += f"- **{todo['id']}**: {todo['title']}{project_str}\n"
                md += "\n"

            if overdue:
                md += "### 🔴 Overdue Tasks\n"
                for todo in sorted(overdue, key=lambda x: x['days_overdue'], reverse=True)[:5]:
                    project_str = f" [{todo['project']}]" if todo.get('project') else ""
                    md += f"- **{todo['id']}**: {todo['title']}{project_str} ({todo['days_overdue']} days)\n"
                md += "\n"

            if due_this_week:
                md += "### ⚡ Due This Week\n"
                for todo in sorted(due_this_week, key=lambda x: x['days_until'])[:5]:
                    project_str = f" [{todo['project']}]" if todo.get('project') else ""
                    days_str = "today" if todo['days_until'] == 0 else f"in {todo['days_until']} days"
                    md += f"- **{todo['id']}**: {todo['title']}{project_str} ({days_str})\n"
                md += "\n"

            # Quick commands
            md += "### ⚡ Quick Commands\n"
            md += "```bash\n"
            md += "python scripts/immunos_todo.py list --today\n"
            md += "python scripts/immunos_todo.py next\n"
            md += "python scripts/immunos_todo.py complete TODO-ID\n"
            md += "```\n\n"

            return md

        except Exception as e:
            return f"## 📋 Todo Summary\n\nError loading todos: {e}\n\n"

    def format_anomaly(self, anomaly: Dict, index: int) -> str:
        """Format a single anomaly for display"""
        emoji_map = {
            'security': '🔴',
            'code_quality': '🟡',
            'structural': '🔵',
            'documentation': '📄',
        }

        emoji = emoji_map.get(anomaly['category'], '⚪')
        loc = f"{anomaly['file_path']}:{anomaly['line_number']}" if anomaly.get('line_number') else anomaly['file_path']

        output = f"{index}. **{anomaly['severity']}: {anomaly['description']}**\n"
        output += f"   - File: `{loc}`\n"
        if anomaly.get('pattern'):
            output += f"   - Pattern: `{anomaly['pattern']}`\n"
        output += f"   - Recommendation: {anomaly['recommendation']}\n"

        if anomaly.get('ollama_analysis'):
            output += f"   - Ollama Analysis: {anomaly['ollama_analysis']}\n"

        return output

    def generate_entry(self, date_str: str) -> str:
        """Generate journal entry for date"""
        baseline = self.load_baseline()
        nk_scan = self.load_nk_scan()

        # Sort anomalies by severity
        severity_order = {'HIGH': 0, 'MEDIUM': 1, 'LOW': 2}
        anomalies = sorted(
            nk_scan['anomalies'],
            key=lambda a: (severity_order.get(a['severity'], 3), a['category'])
        )

        # Separate by category
        security = [a for a in anomalies if a['category'] == 'security']
        code_quality = [a for a in anomalies if a['category'] == 'code_quality']
        structural = [a for a in anomalies if a['category'] == 'structural']
        documentation = [a for a in anomalies if a['category'] == 'documentation']

        # Calculate health score (100 - weighted anomalies)
        security_penalty = len([a for a in security if a['severity'] == 'HIGH']) * 8
        security_penalty += len([a for a in security if a['severity'] == 'MEDIUM']) * 3
        quality_penalty = len([a for a in code_quality if a['severity'] == 'MEDIUM']) * 1
        health_score = max(0, 100 - security_penalty - quality_penalty - len(structural) * 0.5)

        # Generate markdown with todo summary
        todo_summary = self.generate_todo_summary(date_str)

        md = f"""# IMMUNOS Daily Scan - {date_str}

{todo_summary}## Summary
- **Scan Duration**: ~30 seconds
- **Files Scanned**: {baseline['file_count']:,}
- **Directories**: {baseline['dir_count']:,}
- **Total Size**: {baseline['total_size_mb']:.1f} MB
- **Anomalies Found**: {nk_scan['total_anomalies']} ({nk_scan['by_severity'].get('HIGH', 0)} high, {nk_scan['by_severity'].get('MEDIUM', 0)} medium, {nk_scan['by_severity'].get('LOW', 0)} low)

## Baseline Information
- **Scan Date**: {baseline['scan_date']}
- **Projects Root**: `{baseline['projects_root']}`
- **Paper Directories**: {len(baseline['known_patterns'].get('paper_directories', []))}
- **Protocol Files**: {len(baseline['known_patterns'].get('protocols', []))}
- **Scripts**: {len(baseline['known_patterns'].get('scripts', []))}

## Anomalies Detected

### 🔴 Security Issues ({len(security)})

"""

        if security:
            for i, anomaly in enumerate(security[:10], 1):  # Top 10
                md += self.format_anomaly(anomaly, i) + "\n"
            if len(security) > 10:
                md += f"\n*... and {len(security) - 10} more security issues*\n"
        else:
            md += "✓ No security issues detected\n"

        md += f"""
### 🟡 Code Quality Issues ({len(code_quality)})

"""

        if code_quality:
            # Show top 5 MEDIUM, then summary of rest
            medium = [a for a in code_quality if a['severity'] == 'MEDIUM']
            low = [a for a in code_quality if a['severity'] == 'LOW']

            for i, anomaly in enumerate(medium[:5], 1):
                md += self.format_anomaly(anomaly, i) + "\n"

            if len(medium) > 5:
                md += f"\n*... and {len(medium) - 5} more MEDIUM code quality issues*\n"

            if low:
                md += f"\n**LOW severity**: {len(low)} issues (mostly TODO markers)\n"
        else:
            md += "✓ No code quality issues detected\n"

        md += f"""
### 🔵 Structural Issues ({len(structural)})

"""

        if structural:
            for i, anomaly in enumerate(structural[:5], 1):  # Top 5
                md += self.format_anomaly(anomaly, i) + "\n"
            if len(structural) > 5:
                md += f"\n*... and {len(structural) - 5} more structural issues*\n"
        else:
            md += "✓ No structural issues detected\n"

        md += f"""
### 📄 Documentation Issues ({len(documentation)})

"""

        if documentation:
            for i, anomaly in enumerate(documentation, 1):
                md += self.format_anomaly(anomaly, i) + "\n"
        else:
            md += "✓ Documentation complete\n"

        md += f"""
## Suggestions

### Security
"""
        if security:
            md += """1. **URGENT**: Review and remove hardcoded API keys and secrets
2. Set up pre-commit hook to prevent secrets in commits
3. Audit file permissions on sensitive files
4. Consider using a secrets manager (e.g., 1Password, AWS Secrets Manager)
"""
        else:
            md += "✓ No immediate security actions needed\n"

        md += """
### Code Quality
"""
        if code_quality:
            md += f"""1. Refactor {len([a for a in code_quality if 'long_function' in a.get('pattern', '')])} long functions (>100 lines)
2. Address {len([a for a in code_quality if 'todo' in a.get('pattern', '')])} TODO/FIXME markers
3. Add type hints for better maintainability
4. Consider extracting duplicate code into utilities
"""
        else:
            md += "✓ Code quality is good\n"

        md += """
### Organization
"""
        if structural:
            md += f"""1. Standardize naming conventions ({len(structural)} inconsistencies found)
2. Review directory structure for clarity
3. Ensure consistent patterns across projects
"""
        else:
            md += "✓ Project structure is consistent\n"

        md += """
### Documentation
"""
        if documentation:
            md += f"""1. Add docstrings to {len(documentation)} files
2. Ensure README files in all project directories
3. Update outdated timestamps
"""
        else:
            md += "✓ Documentation is complete\n"

        md += f"""
## Project Health Score: {health_score:.1f}/100

### Scoring Breakdown
- **Security**: {100 - security_penalty:.0f}/100 ({len(security)} issues, -{security_penalty} points)
- **Code Quality**: {100 - quality_penalty:.0f}/100 ({len([a for a in code_quality if a['severity'] != 'LOW'])} issues, -{quality_penalty} points)
- **Structure**: {100 - len(structural)*0.5:.0f}/100 ({len(structural)} issues, -{len(structural)*0.5:.1f} points)
- **Documentation**: {100 - len(documentation)*2:.0f}/100 ({len(documentation)} issues, -{len(documentation)*2} points)

## Next Actions

1. **Immediate** (Security):
"""

        high_security = [a for a in security if a['severity'] == 'HIGH']
        if high_security:
            for i, anomaly in enumerate(high_security[:3], 1):
                md += f"   - [ ] Fix: {anomaly['file_path']} - {anomaly['description']}\n"
        else:
            md += "   - [x] No immediate security actions\n"

        md += """
2. **Short-term** (Code Quality):
"""
        medium_quality = [a for a in code_quality if a['severity'] == 'MEDIUM']
        if medium_quality:
            for i, anomaly in enumerate(medium_quality[:3], 1):
                md += f"   - [ ] Refactor: {anomaly['file_path']}\n"
        else:
            md += "   - [x] Code quality is satisfactory\n"

        md += """
3. **Long-term** (Maintenance):
   - [ ] Standardize naming conventions across projects
   - [ ] Add comprehensive documentation
   - [ ] Set up automated quality checks
   - [ ] Regular security audits

---

**Generated by**: IMMUNOS Project Scanner
**NK Cell Scan**: Completed
**B Cell Verification**: Pending
**Next Scan**: Tomorrow
"""

        return md

    def save_entry(self, date_str: str, content: str):
        """Save journal entry to file"""
        filename = f"{date_str}.md"
        filepath = self.journal_dir / filename

        with open(filepath, 'w') as f:
            f.write(content)

        print(f"✓ Journal entry saved to {filepath}")
        print(f"  Size: {len(content)} bytes")

        # Update index
        self.update_index()

    def update_index(self):
        """Update journal index"""
        index_path = self.journal_dir / "INDEX.md"

        # List all journal entries
        entries = sorted(self.journal_dir.glob("*.md"), reverse=True)
        entries = [e for e in entries if e.name != "INDEX.md"]

        content = "# IMMUNOS Daily Journal Index\n\n"
        content += "Daily security and quality scans of ~/projects/\n\n"
        content += "## Recent Entries\n\n"

        for entry in entries[:30]:  # Last 30 days
            date = entry.stem
            content += f"- [{date}](./{entry.name})\n"

        with open(index_path, 'w') as f:
            f.write(content)


def main():
    parser = argparse.ArgumentParser(description='Generate daily journal entry')
    parser.add_argument('--date', type=str,
                        default=datetime.now().strftime("%Y-%m-%d"),
                        help='Date for journal entry (YYYY-MM-DD)')
    parser.add_argument('--immunos-dir', type=str,
                        default='/Users/byron/projects/.immunos',
                        help='IMMUNOS data directory')

    args = parser.parse_args()

    immunos_dir = Path(args.immunos_dir)
    if not immunos_dir.exists():
        print(f"Error: IMMUNOS directory not found: {immunos_dir}")
        print("Run immunos_baseline.py and immunos_nk_scan.py first")
        sys.exit(1)

    # Generate journal entry
    generator = JournalGenerator(immunos_dir)
    content = generator.generate_entry(args.date)
    generator.save_entry(args.date, content)

    print()
    print("Journal entry generated successfully!")
    print(f"View at: {immunos_dir}/journal/{args.date}.md")


if __name__ == '__main__':
    import sys
    main()

````
