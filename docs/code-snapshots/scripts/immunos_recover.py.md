---
source: /Users/byron/projects/scripts/immunos_recover.py
relative: scripts/immunos_recover.py
generated_at: 2025-12-23 10:28
---

````python
#!/usr/bin/env python3
"""
IMMUNOS Context Recovery System

Quickly restore Claude context after reset or compaction.
Generates recovery file from latest snapshot in <5 seconds.

Usage:
    python immunos_recover.py                    # Recover from latest snapshot
    python immunos_recover.py --snapshot snap_xxx # Recover from specific snapshot
    python immunos_recover.py --preview          # Preview without writing
    python immunos_recover.py --auto             # Auto-detect and recover
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
import subprocess


class ContextRecovery:
    """Fast context recovery from snapshots"""

    def __init__(self, base_path: str = '/Users/byron/projects/.immunos'):
        self.base_path = Path(base_path)
        self.snapshots_path = self.base_path / 'memory' / 'snapshots'
        self.recovery_path = self.base_path / 'recovery'
        self.recovery_path.mkdir(parents=True, exist_ok=True)

        self.projects_root = Path('/Users/byron/projects')

    def get_latest_snapshot(self) -> Optional[Dict]:
        """Get the latest snapshot"""
        latest_link = self.snapshots_path / 'latest.json'

        if not latest_link.exists():
            # No symlink, find most recent file
            snapshots = sorted(
                self.snapshots_path.glob('snap_*.json'),
                key=lambda p: p.stat().st_mtime,
                reverse=True
            )

            if not snapshots:
                return None

            with open(snapshots[0], 'r') as f:
                return json.load(f)

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

    def load_relevant_memories(self, snapshot: Dict, max_memories: int = 10) -> List[Dict]:
        """Load relevant memories from memory system"""
        memories = []

        try:
            # Import memory system
            sys.path.insert(0, str(self.projects_root / 'scripts'))
            from immunos_memory import MemoryStore

            store = MemoryStore()

            # Get high-priority memories
            high_priority = store.list_memories(priority='high', limit=5)
            memories.extend([m.to_dict() for m in high_priority])

            # Get recent memories
            recent = store.list_memories(min_relevance=0.5, limit=5)
            memories.extend([m.to_dict() for m in recent if m.to_dict() not in memories])

        except Exception as e:
            # Memory system not available or error
            pass

        return memories[:max_memories]

    def _parse_timestamp(self, timestamp_str: str) -> datetime:
        """Parse timestamp handling both ISO format and 'Z' suffix"""
        # Replace 'Z' with '+00:00' for Python's fromisoformat
        if timestamp_str.endswith('Z'):
            timestamp_str = timestamp_str[:-1] + '+00:00'
        dt = datetime.fromisoformat(timestamp_str)
        # Convert to naive datetime (remove timezone info for comparison)
        if dt.tzinfo is not None:
            dt = dt.replace(tzinfo=None)
        return dt

    def generate_recovery_markdown(self, snapshot: Dict, memories: Optional[List[Dict]] = None) -> str:
        """Generate comprehensive recovery markdown"""
        timestamp = self._parse_timestamp(snapshot['timestamp'])
        now = datetime.now()
        time_ago = self._format_time_ago((now - timestamp).total_seconds())

        lines = [
            f"# Claude Context Recovery",
            f"",
            f"**Recovery Date**: {now.strftime('%Y-%m-%d %H:%M:%S')}",
            f"**Snapshot**: {snapshot['snapshot_id']}",
            f"**Last Active**: {timestamp.strftime('%Y-%m-%d %H:%M:%S')} ({time_ago} ago)",
            f"",
            "---",
            "",
            "## 🎯 Current Focus",
            "",
            f"{snapshot.get('context_summary', 'No summary available')}",
            "",
        ]

        # Active tasks
        if snapshot.get('active_tasks'):
            lines.append("## ✅ Active Tasks")
            lines.append("")

            for task in snapshot['active_tasks']:
                status_emoji = {
                    'completed': '✅',
                    'in_progress': '🔄',
                    'pending': '⏳',
                    'blocked': '⛔',
                }.get(task.get('status', 'unknown'), '❓')

                lines.append(f"### {status_emoji} {task.get('task', 'Unknown task')}")

                if task.get('progress'):
                    lines.append(f"- **Progress**: {task['progress']}")

                if task.get('status'):
                    lines.append(f"- **Status**: {task['status']}")

                if task.get('files_involved'):
                    lines.append(f"- **Files**:")
                    for file in task['files_involved']:
                        lines.append(f"  - `{file}`")

                if task.get('next_action'):
                    lines.append(f"- **Next Action**: {task['next_action']}")

                lines.append("")

        # Active todos
        if snapshot.get('active_todos'):
            lines.append("## 📋 Active Todos")
            lines.append("")

            todos = snapshot['active_todos']

            # Categorize by urgency
            overdue = []
            due_today = []
            due_soon = []

            for todo in todos:
                if todo.get('due'):
                    due_date_str = todo['due']
                    # Handle both full ISO format and date-only format
                    if 'T' in due_date_str:
                        due_date = datetime.fromisoformat(due_date_str.replace('Z', '+00:00'))
                    else:
                        due_date = datetime.strptime(due_date_str, '%Y-%m-%d')

                    days_until = (due_date.date() - now.date()).days

                    if days_until < 0:
                        overdue.append((todo, abs(days_until)))
                    elif days_until == 0:
                        due_today.append(todo)
                    elif days_until <= 7:
                        due_soon.append((todo, days_until))

            # Display sections
            if overdue:
                lines.append("### 🔴 Overdue")
                for todo, days in overdue:
                    project_str = f" [{todo['project']}]" if todo.get('project') else ""
                    lines.append(f"- **{todo['id']}**: {todo['title']}{project_str} ({days} days overdue)")
                lines.append("")

            if due_today:
                lines.append("### 🟡 Due Today")
                for todo in due_today:
                    project_str = f" [{todo['project']}]" if todo.get('project') else ""
                    priority_str = todo.get('priority', 'medium')
                    lines.append(f"- **{todo['id']}**: {todo['title']}{project_str} [{priority_str}]")
                lines.append("")

            if due_soon:
                lines.append("### 🟢 Due This Week")
                for todo, days in due_soon[:5]:
                    project_str = f" [{todo['project']}]" if todo.get('project') else ""
                    lines.append(f"- **{todo['id']}**: {todo['title']}{project_str} (in {days} days)")
                lines.append("")

            # Next actions
            next_actions = [t for t in todos if t.get('status') == 'next']
            if next_actions:
                lines.append("### ⚡ Next Actions")
                for todo in next_actions[:5]:
                    project_str = f" [{todo.get('project', 'general')}]" if todo.get('project') else ""
                    lines.append(f"- **{todo['id']}**: {todo['title']}{project_str}")
                lines.append("")

            # Quick commands
            lines.append("### ⚡ Quick Todo Commands")
            lines.append("```bash")
            lines.append("python scripts/immunos_todo.py list --overdue")
            lines.append("python scripts/immunos_todo.py next")
            lines.append("python scripts/immunos_todo.py complete TODO-ID")
            lines.append("```")
            lines.append("")

        # Conversation summary
        if snapshot.get('conversation_summary'):
            lines.append("## 💬 Conversation Summary")
            lines.append("")
            lines.append(snapshot['conversation_summary'])
            lines.append("")

        # Recent work
        if snapshot.get('recent_files'):
            lines.append("## 📝 Recent Work")
            lines.append("")
            lines.append("Files modified in last 24 hours:")
            lines.append("")

            # Group by importance
            important_files = []
            other_files = []

            for file_info in snapshot['recent_files']:
                if any(keyword in file_info['path'] for keyword in ['.immunos', 'scripts/', '.claudeignore']):
                    important_files.append(file_info)
                else:
                    other_files.append(file_info)

            if important_files:
                lines.append("### Key Files")
                for file_info in important_files[:10]:
                    file_timestamp = self._parse_timestamp(file_info['timestamp'])
                    file_ago = self._format_time_ago((now - file_timestamp).total_seconds())
                    lines.append(f"- `{file_info['path']}` - {file_ago} ago")
                lines.append("")

            if other_files and len(other_files) <= 5:
                lines.append("### Other Changes")
                for file_info in other_files[:5]:
                    file_timestamp = self._parse_timestamp(file_info['timestamp'])
                    file_ago = self._format_time_ago((now - file_timestamp).total_seconds())
                    lines.append(f"- `{file_info['path']}` - {file_ago} ago")
                lines.append("")

        # Key decisions
        if snapshot.get('key_decisions'):
            lines.append("## 🔑 Key Decisions")
            lines.append("")

            for decision in snapshot['key_decisions']:
                lines.append(f"### {decision.get('decision', 'Unknown')}")

                if decision.get('rationale'):
                    lines.append(f"**Rationale**: {decision['rationale']}")
                    lines.append("")

                if decision.get('timestamp'):
                    dec_time = self._parse_timestamp(decision['timestamp'])
                    dec_ago = self._format_time_ago((now - dec_time).total_seconds())
                    lines.append(f"*Made {dec_ago} ago*")
                    lines.append("")

        # Key learnings
        if snapshot.get('key_learnings'):
            lines.append("## 💡 Key Learnings")
            lines.append("")

            for learning in snapshot['key_learnings']:
                lines.append(f"- {learning}")

            lines.append("")

        # Relevant memories
        if memories:
            lines.append("## 🧠 Relevant Memories")
            lines.append("")

            for memory in memories[:5]:
                mem_time = self._parse_timestamp(memory['timestamp'])
                mem_ago = self._format_time_ago((now - mem_time).total_seconds())

                lines.append(f"### {memory.get('type', 'memory').title()} ({mem_ago} ago)")
                lines.append(f"**Priority**: {memory.get('priority', 'unknown')} | **Relevance**: {memory.get('relevance_score', 0):.2f}")

                content = memory.get('content', {})
                if isinstance(content, dict):
                    if 'topic' in content:
                        lines.append(f"**Topic**: {content['topic']}")
                    if 'summary' in content:
                        lines.append(f"{content['summary']}")
                    if 'key_points' in content and content['key_points']:
                        lines.append("**Key Points**:")
                        for point in content['key_points'][:3]:
                            lines.append(f"- {point}")
                else:
                    lines.append(f"{content}")

                lines.append("")

        # Project state
        if snapshot.get('project_state'):
            lines.append("## 📊 Project State")
            lines.append("")

            state = snapshot['project_state']

            if 'immunos_deployed' in state:
                lines.append("### IMMUNOS Components")
                lines.append("")
                for component, status in state['immunos_deployed'].items():
                    component_name = component.replace('_', ' ').title()
                    lines.append(f"- **{component_name}**: {status}")
                lines.append("")

            if 'token_optimization' in state:
                lines.append("### Token Optimization Status")
                lines.append("")
                for key, value in state['token_optimization'].items():
                    key_name = key.replace('_', ' ').title()
                    lines.append(f"- **{key_name}**: {value}")
                lines.append("")

        # User preferences
        if snapshot.get('user_preferences'):
            lines.append("## 👤 User Preferences")
            lines.append("")

            prefs = snapshot['user_preferences']
            for key, value in prefs.items():
                key_name = key.replace('_', ' ').title()
                lines.append(f"- **{key_name}**: {value}")

            lines.append("")

        # Environment
        if snapshot.get('environment'):
            lines.append("## 🖥️ Environment")
            lines.append("")

            env = snapshot['environment']
            lines.append(f"- **Working Directory**: `{env.get('cwd', 'unknown')}`")
            lines.append(f"- **Git Status**: {env.get('git_status', 'unknown')}")
            lines.append(f"- **Ollama Running**: {'✅ Yes' if env.get('ollama_running') else '❌ No'}")

            if env.get('ollama_models'):
                lines.append(f"- **Available Models**: {', '.join(env['ollama_models'])}")

            lines.append("")

        # Next steps
        if snapshot.get('next_steps'):
            lines.append("## ⏭️ Next Steps")
            lines.append("")

            for i, step in enumerate(snapshot['next_steps'], 1):
                lines.append(f"{i}. {step}")

            lines.append("")

        # User context from claude.md
        claude_md_path = self.projects_root / 'claude.md'
        if claude_md_path.exists():
            lines.append("## 👤 User Context")
            lines.append("")
            try:
                with open(claude_md_path, 'r') as f:
                    content = f.read()
                    # Extract personal info section
                    if '## Personal Information' in content:
                        personal_section = content.split('## Personal Information')[1].split('---')[0].strip()
                        lines.append("**Personal Information**: See `claude.md` for full details")
                        lines.append("")
                        # Extract key info
                        for line in personal_section.split('\n'):
                            if line.startswith('**Name**:') or line.startswith('**Email**:') or line.startswith('**LinkedIn**:'):
                                lines.append(line)
                        lines.append("")
                    lines.append("**Full Context File**: `~/projects/claude.md`")
                    lines.append("")
            except Exception:
                lines.append("**Context File**: `~/projects/claude.md` (available)")
                lines.append("")

        # Quick commands
        lines.append("## ⚡ Quick Commands")
        lines.append("")
        lines.append("```bash")
        lines.append("# View user context")
        lines.append("cat ~/projects/claude.md")
        lines.append("")
        lines.append("# Continue current work")
        lines.append("cd ~/projects")
        lines.append("")
        lines.append("# View recent journal")
        lines.append("cat .immunos/journal/$(date +%Y-%m-%d).md")
        lines.append("")
        lines.append("# Check token status")
        lines.append("python scripts/immunos_token_analyzer.py --top 10")
        lines.append("")
        lines.append("# View memory status")
        lines.append("python scripts/immunos_memory.py stats")
        lines.append("")
        lines.append("# Create new snapshot")
        lines.append("python scripts/immunos_snapshot.py create --trigger manual --summary 'Context update'")
        lines.append("```")
        lines.append("")

        # Footer
        lines.append("---")
        lines.append("")
        lines.append(f"*Generated by IMMUNOS Recovery System*")
        lines.append(f"*Snapshot: {snapshot['snapshot_id']}*")
        lines.append(f"*Recovery Time: <5 seconds*")

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

    def recover(self, snapshot_id: Optional[str] = None, preview: bool = False) -> Dict:
        """Recover context from snapshot"""
        # Load snapshot
        if snapshot_id:
            snapshot = self.get_snapshot(snapshot_id)
            if not snapshot:
                raise ValueError(f"Snapshot not found: {snapshot_id}")
        else:
            snapshot = self.get_latest_snapshot()
            if not snapshot:
                raise ValueError("No snapshots available for recovery")

        # Load relevant memories
        memories = self.load_relevant_memories(snapshot)

        # Generate recovery markdown
        recovery_content = self.generate_recovery_markdown(snapshot, memories)

        if not preview:
            # Write recovery file
            recovery_file = self.recovery_path / 'CONTEXT_RECOVERY.md'
            with open(recovery_file, 'w') as f:
                f.write(recovery_content)

            # Also create a quick-start script
            self._create_quickstart_script(snapshot)

        return {
            'snapshot_id': snapshot['snapshot_id'],
            'snapshot_timestamp': snapshot['timestamp'],
            'recovery_content': recovery_content,
            'memories_loaded': len(memories),
            'preview': preview,
        }

    def _create_quickstart_script(self, snapshot: Dict):
        """Create a quick-start bash script"""
        script_path = self.recovery_path / 'quick_start.sh'

        lines = [
            "#!/bin/bash",
            "# IMMUNOS Quick Start - Auto-generated",
            "",
            "echo '🧬 IMMUNOS Context Recovery'",
            "echo '============================'",
            "echo ''",
            "",
            "# Display recovery file",
            "if [ -f ~/.immunos/recovery/CONTEXT_RECOVERY.md ]; then",
            "    cat ~/.immunos/recovery/CONTEXT_RECOVERY.md",
            "else",
            "    echo 'Recovery file not found'",
            "    exit 1",
            "fi",
            "",
            "echo ''",
            "echo 'Ready to continue work!'",
        ]

        with open(script_path, 'w') as f:
            f.write('\n'.join(lines))

        script_path.chmod(0o755)

    def auto_recover(self) -> bool:
        """Auto-detect if recovery needed and recover"""
        recovery_file = self.recovery_path / 'CONTEXT_RECOVERY.md'

        # Check if recovery file is stale (>1 hour old)
        if recovery_file.exists():
            age = datetime.now().timestamp() - recovery_file.stat().st_mtime
            if age < 3600:  # Less than 1 hour old
                print(f"✓ Recovery file is recent ({int(age / 60)} minutes old)")
                return False

        # Check if there are snapshots
        snapshots = list(self.snapshots_path.glob('snap_*.json'))
        if not snapshots:
            print("No snapshots available for recovery")
            return False

        # Perform recovery
        print("Auto-recovering from latest snapshot...")
        result = self.recover()

        print(f"✓ Context recovered from {result['snapshot_id']}")
        print(f"  Snapshot age: {result['snapshot_timestamp']}")
        print(f"  Memories loaded: {result['memories_loaded']}")
        print(f"  Recovery file: {self.recovery_path / 'CONTEXT_RECOVERY.md'}")

        return True


def main():
    parser = argparse.ArgumentParser(description='IMMUNOS Context Recovery')
    parser.add_argument('--snapshot', type=str, help='Specific snapshot ID to recover from')
    parser.add_argument('--preview', action='store_true', help='Preview without writing files')
    parser.add_argument('--auto', action='store_true', help='Auto-detect and recover if needed')

    args = parser.parse_args()

    recovery = ContextRecovery()

    try:
        if args.auto:
            recovery.auto_recover()
        else:
            result = recovery.recover(snapshot_id=args.snapshot, preview=args.preview)

            if args.preview:
                print(result['recovery_content'])
            else:
                print(f"✓ Context recovered successfully!")
                print(f"  Snapshot: {result['snapshot_id']}")
                print(f"  Timestamp: {result['snapshot_timestamp']}")
                print(f"  Memories: {result['memories_loaded']}")
                print(f"  Recovery file: /Users/byron/projects/.immunos/recovery/CONTEXT_RECOVERY.md")
                print("")
                print("To view recovery context:")
                print("  cat ~/.immunos/recovery/CONTEXT_RECOVERY.md")

    except Exception as e:
        print(f"Error during recovery: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()

````
