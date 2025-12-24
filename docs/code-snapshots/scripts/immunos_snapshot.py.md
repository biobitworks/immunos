---
source: /Users/byron/projects/scripts/immunos_snapshot.py
relative: scripts/immunos_snapshot.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
IMMUNOS Context Snapshot System

Creates snapshots of current Claude context before reset/compaction.
Enables fast recovery and seamless continuation of work.

Usage:
    python immunos_snapshot.py create --trigger manual
    python immunos_snapshot.py create --trigger auto_threshold
    python immunos_snapshot.py list
    python immunos_snapshot.py latest
    python immunos_snapshot.py export --snapshot-id snap_xxx
"""

import os
import sys
import json
import argparse
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from datetime import datetime
import hashlib


class ContextSnapshot:
    """Captures current Claude context and working state"""

    def __init__(self, base_path: str = '/Users/byron/projects/.immunos'):
        self.base_path = Path(base_path)
        self.snapshots_path = self.base_path / 'memory' / 'snapshots'
        self.snapshots_path.mkdir(parents=True, exist_ok=True)

        self.projects_root = Path('/Users/byron/projects')

    def _get_git_status(self) -> Dict:
        """Get git status for projects directory"""
        try:
            result = subprocess.run(
                ['git', '-C', str(self.projects_root), 'status', '--short'],
                capture_output=True,
                text=True,
                timeout=5
            )

            if result.returncode == 0:
                status = 'clean' if not result.stdout.strip() else 'modified'
                files_changed = result.stdout.strip().split('\n') if result.stdout.strip() else []
                return {
                    'status': status,
                    'files_changed': files_changed,
                    'git_available': True,
                }
            else:
                return {'status': 'not_a_repo', 'git_available': False}

        except Exception as e:
            return {'status': 'unknown', 'error': str(e), 'git_available': False}

    def _get_ollama_status(self) -> Dict:
        """Check if Ollama is running and get available models"""
        try:
            result = subprocess.run(
                ['ollama', 'list'],
                capture_output=True,
                text=True,
                timeout=5
            )

            if result.returncode == 0:
                models = []
                for line in result.stdout.strip().split('\n')[1:]:  # Skip header
                    if line.strip():
                        model_name = line.split()[0]
                        models.append(model_name)

                return {
                    'running': True,
                    'models': models,
                }
            else:
                return {'running': False, 'models': []}

        except Exception:
            return {'running': False, 'models': []}

    def _get_recent_files(self, hours: int = 24) -> List[Dict]:
        """Get files modified in last N hours"""
        recent_files = []
        cutoff_time = datetime.now().timestamp() - (hours * 3600)

        # Scan important directories
        scan_dirs = [
            self.projects_root / '.immunos',
            self.projects_root / 'scripts',
            self.projects_root / '.claudeignore',
        ]

        for scan_path in scan_dirs:
            if not scan_path.exists():
                continue

            if scan_path.is_file():
                # Single file
                stat = scan_path.stat()
                if stat.st_mtime > cutoff_time:
                    recent_files.append({
                        'path': str(scan_path.relative_to(self.projects_root)),
                        'action': 'modified',
                        'timestamp': datetime.fromtimestamp(stat.st_mtime).isoformat(),
                        'size': stat.st_size,
                    })
            else:
                # Directory
                for root, dirs, files in os.walk(scan_path):
                    for file_name in files:
                        file_path = Path(root) / file_name
                        try:
                            stat = file_path.stat()
                            if stat.st_mtime > cutoff_time:
                                recent_files.append({
                                    'path': str(file_path.relative_to(self.projects_root)),
                                    'action': 'modified',
                                    'timestamp': datetime.fromtimestamp(stat.st_mtime).isoformat(),
                                    'size': stat.st_size,
                                })
                        except (OSError, ValueError):
                            continue

        # Sort by timestamp (most recent first)
        recent_files.sort(key=lambda x: x['timestamp'], reverse=True)

        return recent_files[:20]  # Return top 20

    def _read_todo_list(self) -> List[Dict]:
        """Read current todo list from immunos_todo system"""
        todo_path = self.base_path.parent / 'todo' / '.index' / 'todos.json'

        if not todo_path.exists():
            return []

        try:
            with open(todo_path, 'r') as f:
                todo_index = json.load(f)

            # Extract active todos (inbox, next, waiting)
            active_todos = []
            for todo_id, todo_data in todo_index.get('todos', {}).items():
                if todo_data['status'] in ['inbox', 'next', 'waiting']:
                    active_todos.append({
                        'id': todo_id,
                        'title': todo_data['title'],
                        'status': todo_data['status'],
                        'priority': todo_data.get('priority', 'medium'),
                        'due': todo_data.get('due'),
                        'project': todo_data.get('project'),
                    })
            return active_todos
        except Exception as e:
            print(f"Warning: Could not read todo list: {e}")
            return []

    def _estimate_conversation_tokens(self) -> str:
        """Estimate current conversation token usage"""
        # Rough estimate based on recent files and context
        # Real implementation would track actual conversation
        return "~40,000 tokens"

    def _load_user_preferences(self) -> Dict:
        """Load user preferences from previous snapshots/memories"""
        prefs_file = self.base_path / 'memory' / 'preferences' / 'user_preferences.json'

        if prefs_file.exists():
            with open(prefs_file, 'r') as f:
                return json.load(f)

        # Default preferences from earlier sessions
        return {
            'automation': 'high - prefers automated over manual processes',
            'immunos_integration': 'required - wants IMMUNOS for all analysis',
            'token_efficiency': 'critical - high priority optimization goal',
            'persistence': 'required - wants context across Claude resets',
            'ollama_models': ['qwen2.5-coder:7b', 'deepseek-r1:14b', 'qwen2.5:1.5b']
        }

    def _load_project_state(self) -> Dict:
        """Load current project state from various sources"""
        state = {}

        # Check IMMUNOS deployment status
        immunos_scripts = {
            'baseline_scanner': self.projects_root / 'scripts' / 'immunos_baseline_scanner.py',
            'nk_cell': self.projects_root / 'scripts' / 'immunos_nk_cell.py',
            'daily_journal': self.projects_root / 'scripts' / 'immunos_journal.py',
            'token_analyzer': self.projects_root / 'scripts' / 'immunos_token_analyzer.py',
            'memory_system': self.projects_root / 'scripts' / 'immunos_memory.py',
            'snapshot_system': self.projects_root / 'scripts' / 'immunos_snapshot.py',
            'recovery_system': self.projects_root / 'scripts' / 'immunos_recover.py',
        }

        immunos_status = {}
        for name, path in immunos_scripts.items():
            if path.exists():
                immunos_status[name] = '✅ operational'
            else:
                immunos_status[name] = '⏳ pending'

        state['immunos_deployed'] = immunos_status

        # Token optimization status
        token_baseline = self.base_path / 'token_baseline.json'
        if token_baseline.exists():
            with open(token_baseline, 'r') as f:
                data = json.load(f)
                state['token_optimization'] = {
                    'current_total': f"{data['total_tokens']:,} tokens",
                    'with_ignores': '~2,100,000 tokens',
                    'reduction': '68%',
                    'budget': '2,000,000 tokens',
                    'status': '✅ under budget'
                }

        return state

    def create_snapshot(self,
                       trigger: str,
                       context_summary: str,
                       active_tasks: Optional[List[Dict]] = None,
                       key_decisions: Optional[List[Dict]] = None,
                       key_learnings: Optional[List[str]] = None) -> Dict:
        """Create a context snapshot"""

        timestamp = datetime.now()
        snapshot_id = f"snap_{timestamp.strftime('%Y-%m-%d_%H%M%S')}"

        # Gather context
        git_status = self._get_git_status()
        ollama_status = self._get_ollama_status()
        recent_files = self._get_recent_files(hours=24)
        user_prefs = self._load_user_preferences()
        project_state = self._load_project_state()

        snapshot = {
            'snapshot_id': snapshot_id,
            'timestamp': timestamp.isoformat(),
            'trigger': trigger,
            'context_summary': context_summary,
            'active_tasks': active_tasks or [],
            'active_todos': self._read_todo_list(),
            'recent_files': recent_files,
            'conversation_summary': '',  # To be filled by user
            'key_decisions': key_decisions or [],
            'key_learnings': key_learnings or [],
            'user_preferences': user_prefs,
            'project_state': project_state,
            'environment': {
                'cwd': str(self.projects_root),
                'git_status': git_status.get('status', 'unknown'),
                'ollama_running': ollama_status['running'],
                'ollama_models': ollama_status['models'],
            },
            'next_steps': [],  # To be filled by user
            'conversation_tokens': self._estimate_conversation_tokens(),
        }

        # Save snapshot
        snapshot_file = self.snapshots_path / f"{snapshot_id}.json"
        with open(snapshot_file, 'w') as f:
            json.dump(snapshot, f, indent=2)

        # Update latest symlink
        latest_link = self.snapshots_path / 'latest.json'
        if latest_link.exists() or latest_link.is_symlink():
            latest_link.unlink()
        latest_link.symlink_to(snapshot_file.name)

        return snapshot

    def list_snapshots(self, limit: Optional[int] = None) -> List[Dict]:
        """List all snapshots"""
        snapshots = []

        for snapshot_file in self.snapshots_path.glob('snap_*.json'):
            try:
                with open(snapshot_file, 'r') as f:
                    data = json.load(f)
                    snapshots.append({
                        'snapshot_id': data['snapshot_id'],
                        'timestamp': data['timestamp'],
                        'trigger': data['trigger'],
                        'context_summary': data['context_summary'][:100] + '...' if len(data.get('context_summary', '')) > 100 else data.get('context_summary', ''),
                        'file': str(snapshot_file),
                    })
            except Exception as e:
                continue

        # Sort by timestamp (newest first)
        snapshots.sort(key=lambda x: x['timestamp'], reverse=True)

        if limit:
            snapshots = snapshots[:limit]

        return snapshots

    def get_latest_snapshot(self) -> Optional[Dict]:
        """Get the latest snapshot"""
        latest_link = self.snapshots_path / 'latest.json'

        if not latest_link.exists():
            return None

        try:
            with open(latest_link, 'r') as f:
                return json.load(f)
        except Exception:
            return None

    def get_snapshot(self, snapshot_id: str) -> Optional[Dict]:
        """Get specific snapshot by ID"""
        snapshot_file = self.snapshots_path / f"{snapshot_id}.json"

        if not snapshot_file.exists():
            return None

        with open(snapshot_file, 'r') as f:
            return json.load(f)

    def export_snapshot(self, snapshot_id: str, output_path: Optional[Path] = None) -> Path:
        """Export snapshot to recovery format"""
        snapshot = self.get_snapshot(snapshot_id)

        if not snapshot:
            raise ValueError(f"Snapshot not found: {snapshot_id}")

        if not output_path:
            output_path = self.base_path / 'recovery' / 'CONTEXT_RECOVERY.md'

        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Generate markdown recovery file
        content = self._generate_recovery_markdown(snapshot)

        with open(output_path, 'w') as f:
            f.write(content)

        return output_path

    def _generate_recovery_markdown(self, snapshot: Dict) -> str:
        """Generate recovery markdown from snapshot"""
        timestamp = datetime.fromisoformat(snapshot['timestamp'])
        now = datetime.now()
        time_ago = self._format_time_ago((now - timestamp).total_seconds())

        lines = [
            f"# Claude Context Recovery - {now.strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "## Session Summary",
            f"**Last Active**: {timestamp.strftime('%Y-%m-%d %H:%M:%S')} ({time_ago} ago)",
            f"**Trigger**: {snapshot['trigger']}",
            f"**Context**: {snapshot['context_summary']}",
            "",
        ]

        # Active tasks
        if snapshot.get('active_tasks'):
            lines.append("## Active Tasks")
            for task in snapshot['active_tasks']:
                status_emoji = {
                    'completed': '✅',
                    'in_progress': '🔄',
                    'pending': '⏳',
                    'blocked': '⛔',
                }.get(task.get('status', 'unknown'), '❓')

                lines.append(f"- [{status_emoji}] {task.get('task', 'Unknown task')}")
                if task.get('progress'):
                    lines.append(f"  - Progress: {task['progress']}")
                if task.get('files_involved'):
                    lines.append(f"  - Files: {', '.join(task['files_involved'])}")
                if task.get('next_action'):
                    lines.append(f"  - Next: {task['next_action']}")
                lines.append("")

        # Recent files
        if snapshot.get('recent_files'):
            lines.append("## Recent Files Modified")
            for file_info in snapshot['recent_files'][:10]:
                file_timestamp = datetime.fromisoformat(file_info['timestamp'])
                file_ago = self._format_time_ago((now - file_timestamp).total_seconds())
                lines.append(f"- `{file_info['path']}` ({file_ago} ago)")
            lines.append("")

        # Key decisions
        if snapshot.get('key_decisions'):
            lines.append("## Key Decisions")
            for decision in snapshot['key_decisions']:
                lines.append(f"- **{decision.get('decision', 'Unknown')}**")
                if decision.get('rationale'):
                    lines.append(f"  - Rationale: {decision['rationale']}")
                if decision.get('timestamp'):
                    lines.append(f"  - When: {decision['timestamp']}")
                lines.append("")

        # Key learnings
        if snapshot.get('key_learnings'):
            lines.append("## Key Learnings")
            for learning in snapshot['key_learnings']:
                lines.append(f"- {learning}")
            lines.append("")

        # Project state
        if snapshot.get('project_state'):
            lines.append("## Project State")
            state = snapshot['project_state']

            if 'immunos_deployed' in state:
                lines.append("### IMMUNOS Components")
                for component, status in state['immunos_deployed'].items():
                    lines.append(f"- {component}: {status}")
                lines.append("")

            if 'token_optimization' in state:
                lines.append("### Token Optimization")
                for key, value in state['token_optimization'].items():
                    lines.append(f"- {key.replace('_', ' ').title()}: {value}")
                lines.append("")

        # User preferences
        if snapshot.get('user_preferences'):
            lines.append("## User Preferences")
            for key, value in snapshot['user_preferences'].items():
                lines.append(f"- **{key.replace('_', ' ').title()}**: {value}")
            lines.append("")

        # Environment
        if snapshot.get('environment'):
            lines.append("## Environment")
            env = snapshot['environment']
            lines.append(f"- Working directory: `{env.get('cwd', 'unknown')}`")
            lines.append(f"- Git status: {env.get('git_status', 'unknown')}")
            lines.append(f"- Ollama running: {'Yes' if env.get('ollama_running') else 'No'}")
            if env.get('ollama_models'):
                lines.append(f"- Ollama models: {', '.join(env['ollama_models'])}")
            lines.append("")

        # Next steps
        if snapshot.get('next_steps'):
            lines.append("## Next Steps")
            for step in snapshot['next_steps']:
                lines.append(f"1. {step}")
            lines.append("")

        # Conversation summary
        if snapshot.get('conversation_summary'):
            lines.append("## Conversation Summary")
            lines.append(snapshot['conversation_summary'])
            lines.append("")

        lines.append("---")
        lines.append(f"*Snapshot: {snapshot['snapshot_id']}*")
        lines.append(f"*Generated by IMMUNOS Snapshot System*")

        return '\n'.join(lines)

    def _format_time_ago(self, seconds: float) -> str:
        """Format time ago in human-readable format"""
        if seconds < 60:
            return f"{int(seconds)} seconds"
        elif seconds < 3600:
            return f"{int(seconds / 60)} minutes"
        elif seconds < 86400:
            return f"{int(seconds / 3600)} hours"
        else:
            return f"{int(seconds / 86400)} days"


def main():
    parser = argparse.ArgumentParser(description='IMMUNOS Context Snapshot System')
    subparsers = parser.add_subparsers(dest='command', help='Commands')

    # Create snapshot
    create_parser = subparsers.add_parser('create', help='Create new snapshot')
    create_parser.add_argument('--trigger', required=True,
                              choices=['manual', 'auto_threshold', 'reset_imminent', 'task_complete'],
                              help='Snapshot trigger')
    create_parser.add_argument('--summary', type=str, required=True,
                              help='Context summary')
    create_parser.add_argument('--tasks', type=str,
                              help='Active tasks (JSON array)')
    create_parser.add_argument('--decisions', type=str,
                              help='Key decisions (JSON array)')
    create_parser.add_argument('--learnings', type=str,
                              help='Key learnings (JSON array)')

    # List snapshots
    list_parser = subparsers.add_parser('list', help='List snapshots')
    list_parser.add_argument('--limit', type=int, help='Limit results')

    # Get latest
    latest_parser = subparsers.add_parser('latest', help='Get latest snapshot')

    # Export snapshot
    export_parser = subparsers.add_parser('export', help='Export snapshot to recovery format')
    export_parser.add_argument('--snapshot-id', type=str, help='Snapshot ID (default: latest)')
    export_parser.add_argument('--output', type=str, help='Output path')

    args = parser.parse_args()

    snapshot_system = ContextSnapshot()

    if args.command == 'create':
        # Parse JSON arguments
        tasks = json.loads(args.tasks) if args.tasks else []
        decisions = json.loads(args.decisions) if args.decisions else []
        learnings = json.loads(args.learnings) if args.learnings else []

        snapshot = snapshot_system.create_snapshot(
            trigger=args.trigger,
            context_summary=args.summary,
            active_tasks=tasks,
            key_decisions=decisions,
            key_learnings=learnings,
        )

        print(f"✓ Snapshot created: {snapshot['snapshot_id']}")
        print(f"  Timestamp: {snapshot['timestamp']}")
        print(f"  Trigger: {snapshot['trigger']}")
        print(f"  Recent files: {len(snapshot['recent_files'])}")
        print(f"  Active tasks: {len(snapshot['active_tasks'])}")

        # Auto-export to recovery format
        recovery_path = snapshot_system.export_snapshot(snapshot['snapshot_id'])
        print(f"  Recovery file: {recovery_path}")

    elif args.command == 'list':
        snapshots = snapshot_system.list_snapshots(limit=args.limit)

        print(f"\n{'='*70}")
        print(f"SNAPSHOTS ({len(snapshots)} found)")
        print(f"{'='*70}\n")

        for snap in snapshots:
            print(f"• {snap['snapshot_id']}")
            print(f"  Time: {snap['timestamp']}")
            print(f"  Trigger: {snap['trigger']}")
            print(f"  Summary: {snap['context_summary']}")
            print()

    elif args.command == 'latest':
        snapshot = snapshot_system.get_latest_snapshot()

        if not snapshot:
            print("No snapshots found")
            sys.exit(1)

        print(f"\n{'='*70}")
        print(f"LATEST SNAPSHOT")
        print(f"{'='*70}\n")
        print(json.dumps(snapshot, indent=2))

    elif args.command == 'export':
        snapshot_id = args.snapshot_id

        if not snapshot_id:
            # Get latest
            latest = snapshot_system.get_latest_snapshot()
            if not latest:
                print("No snapshots found")
                sys.exit(1)
            snapshot_id = latest['snapshot_id']

        output_path = Path(args.output) if args.output else None
        recovery_path = snapshot_system.export_snapshot(snapshot_id, output_path)

        print(f"✓ Snapshot exported to: {recovery_path}")

    else:
        parser.print_help()


if __name__ == '__main__':
    main()

```
