#!/usr/bin/env python3
"""
IMMUNOS Todo System - GTD Workflow Manager

A Getting Things Done (GTD) todo system with:
- Daily file organization (todo/YYYY-MM-DD.md)
- GTD workflow (inbox → next → waiting → someday → done)
- CLI access without Claude
- Full IMMUNOS integration (snapshots, recovery, journal, NK scan)

Usage:
    python immunos_todo.py add "Task" --due tomorrow --priority high
    python immunos_todo.py list --status next
    python immunos_todo.py next
    python immunos_todo.py complete TODO-ID
    python immunos_todo.py stats
"""

import json
import argparse
import sys
from pathlib import Path
from datetime import datetime, timedelta
from typing import Dict, List, Optional
from dataclasses import dataclass, asdict
from enum import Enum


class Status(Enum):
    """GTD workflow statuses"""
    INBOX = "inbox"
    NEXT = "next"
    WAITING = "waiting"
    SOMEDAY = "someday"
    DONE = "done"


class Priority(Enum):
    """Task priorities"""
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"


@dataclass
class TodoItem:
    """Single todo item with GTD workflow support"""
    id: str
    title: str
    status: str
    created: str
    due: Optional[str] = None
    completed: Optional[str] = None
    priority: str = "medium"
    project: Optional[str] = None
    tags: List[str] = None
    description: Optional[str] = None
    file: Optional[str] = None

    def __post_init__(self):
        if self.tags is None:
            self.tags = []

    def to_dict(self) -> Dict:
        """Convert to dictionary"""
        return asdict(self)

    def is_overdue(self) -> bool:
        """Check if todo is overdue"""
        if not self.due or self.status == Status.DONE.value:
            return False

        due_date = datetime.fromisoformat(self.due)
        return datetime.now() > due_date

    def days_until_due(self) -> Optional[int]:
        """Calculate days until due date"""
        if not self.due:
            return None

        due_date = datetime.fromisoformat(self.due).date()
        today = datetime.now().date()
        return (due_date - today).days

    def to_markdown(self) -> str:
        """Convert to markdown format"""
        lines = []
        lines.append(f"### {self.id}: {self.title}")
        lines.append(f"- **Status**: {self.status}")
        lines.append(f"- **Created**: {datetime.fromisoformat(self.created).strftime('%Y-%m-%d %H:%M')}")

        if self.due:
            lines.append(f"- **Due**: {datetime.fromisoformat(self.due).strftime('%Y-%m-%d')}")

        if self.completed:
            lines.append(f"- **Completed**: {datetime.fromisoformat(self.completed).strftime('%Y-%m-%d %H:%M')}")

        lines.append(f"- **Priority**: {self.priority}")

        if self.project:
            lines.append(f"- **Project**: {self.project}")

        if self.tags:
            tags_str = " ".join([f"#{tag}" for tag in self.tags])
            lines.append(f"- **Tags**: {tags_str}")

        if self.description:
            lines.append(f"\n{self.description}")

        return "\n".join(lines)


class TodoStore:
    """Storage and retrieval with index management"""

    def __init__(self, base_path: Path):
        self.base_path = base_path
        self.todo_dir = base_path / "todo"
        self.daily_dir = self.todo_dir / "daily"
        self.index_dir = self.todo_dir / ".index"
        self.archive_dir = self.todo_dir / "archive"
        self.templates_dir = self.todo_dir / "templates"

        self.index_file = self.index_dir / "todos.json"
        self.stats_file = self.index_dir / "stats.json"

        self._ensure_directories()
        self._load_index()

    def _ensure_directories(self):
        """Create directory structure if needed"""
        self.daily_dir.mkdir(parents=True, exist_ok=True)
        self.index_dir.mkdir(parents=True, exist_ok=True)
        self.archive_dir.mkdir(parents=True, exist_ok=True)
        self.templates_dir.mkdir(parents=True, exist_ok=True)

    def _load_index(self):
        """Load index from JSON file"""
        if self.index_file.exists():
            with open(self.index_file, 'r') as f:
                data = json.load(f)
                self.todos = data.get('todos', {})
                self.last_updated = data.get('last_updated', datetime.now().isoformat())
        else:
            self.todos = {}
            self.last_updated = datetime.now().isoformat()
            self._save_index()

    def _save_index(self):
        """Save index to JSON file"""
        self.last_updated = datetime.now().isoformat()

        data = {
            'last_updated': self.last_updated,
            'todos': self.todos,
            'stats': self._calculate_stats()
        }

        with open(self.index_file, 'w') as f:
            json.dump(data, f, indent=2)

    def _calculate_stats(self) -> Dict:
        """Calculate statistics"""
        stats = {
            'total': len(self.todos),
            'by_status': {},
            'by_priority': {},
            'overdue': 0,
            'due_today': 0,
            'due_this_week': 0
        }

        today = datetime.now().date()
        week_from_now = today + timedelta(days=7)

        for todo_data in self.todos.values():
            # Status counts
            status = todo_data['status']
            stats['by_status'][status] = stats['by_status'].get(status, 0) + 1

            # Priority counts
            priority = todo_data['priority']
            stats['by_priority'][priority] = stats['by_priority'].get(priority, 0) + 1

            # Due date analysis
            if todo_data['status'] != Status.DONE.value and todo_data.get('due'):
                due_date = datetime.fromisoformat(todo_data['due']).date()

                if due_date < today:
                    stats['overdue'] += 1
                elif due_date == today:
                    stats['due_today'] += 1
                elif due_date <= week_from_now:
                    stats['due_this_week'] += 1

        return stats

    def _generate_todo_id(self, date: datetime) -> str:
        """Generate unique todo ID"""
        date_str = date.strftime('%Y%m%d')
        prefix = f"TODO-{date_str}"

        # Find next available number for this date
        existing_ids = [tid for tid in self.todos.keys() if tid.startswith(prefix)]
        if not existing_ids:
            return f"{prefix}-001"

        # Extract numbers and find max
        numbers = [int(tid.split('-')[-1]) for tid in existing_ids]
        next_num = max(numbers) + 1
        return f"{prefix}-{next_num:03d}"

    def _get_daily_file(self, date: datetime) -> Path:
        """Get daily file path for a date"""
        return self.daily_dir / f"{date.strftime('%Y-%m-%d')}.md"

    def _update_daily_file(self, todo: TodoItem):
        """Update daily markdown file with todo"""
        created_date = datetime.fromisoformat(todo.created)
        daily_file = self._get_daily_file(created_date)

        # Load existing file or create new
        if daily_file.exists():
            with open(daily_file, 'r') as f:
                content = f.read()
        else:
            # Create new daily file from template
            content = self._generate_daily_file_header(created_date)

        # Parse frontmatter and sections
        sections = self._parse_daily_file(content)

        # Add or update todo in appropriate section
        status_section = self._get_status_section_name(todo.status)
        if status_section not in sections:
            sections[status_section] = []

        # Remove todo from all sections (in case of status change)
        for section_name in sections:
            sections[section_name] = [t for t in sections[section_name] if not t.startswith(f"### {todo.id}:")]

        # Add to correct section
        sections[status_section].append(todo.to_markdown())

        # Write updated file
        new_content = self._generate_daily_file(created_date, sections)
        with open(daily_file, 'w') as f:
            f.write(new_content)

    def _generate_daily_file_header(self, date: datetime) -> str:
        """Generate header for new daily file"""
        return f"""---
date: {date.strftime('%Y-%m-%d')}
total_todos: 0
completed: 0
---

# Todo List - {date.strftime('%Y-%m-%d')}

"""

    def _parse_daily_file(self, content: str) -> Dict[str, List[str]]:
        """Parse daily file into sections"""
        sections = {
            'Inbox': [],
            'Next Actions': [],
            'Waiting': [],
            'Someday': [],
            'Completed': []
        }

        current_section = None
        current_todo_lines = []

        lines = content.split('\n')
        for line in lines:
            if line.startswith('## '):
                # Save previous todo
                if current_section and current_todo_lines:
                    sections[current_section].append('\n'.join(current_todo_lines))
                    current_todo_lines = []

                # Update section
                section_name = line[3:].strip()
                if section_name in sections:
                    current_section = section_name
            elif line.startswith('### TODO-') and current_section:
                # Save previous todo
                if current_todo_lines:
                    sections[current_section].append('\n'.join(current_todo_lines))
                current_todo_lines = [line]
            elif current_section and current_todo_lines:
                current_todo_lines.append(line)

        # Save last todo
        if current_section and current_todo_lines:
            sections[current_section].append('\n'.join(current_todo_lines))

        return sections

    def _get_status_section_name(self, status: str) -> str:
        """Map status to section name"""
        mapping = {
            Status.INBOX.value: 'Inbox',
            Status.NEXT.value: 'Next Actions',
            Status.WAITING.value: 'Waiting',
            Status.SOMEDAY.value: 'Someday',
            Status.DONE.value: 'Completed'
        }
        return mapping.get(status, 'Inbox')

    def _generate_daily_file(self, date: datetime, sections: Dict[str, List[str]]) -> str:
        """Generate complete daily file content"""
        # Count todos
        total = sum(len(todos) for todos in sections.values())
        completed = len(sections.get('Completed', []))

        lines = [
            "---",
            f"date: {date.strftime('%Y-%m-%d')}",
            f"total_todos: {total}",
            f"completed: {completed}",
            "---",
            "",
            f"# Todo List - {date.strftime('%Y-%m-%d')}",
            ""
        ]

        # Add sections
        for section_name, todos in sections.items():
            if todos:
                lines.append(f"## {section_name}")
                lines.append("")
                for todo_md in todos:
                    lines.append(todo_md)
                    lines.append("")

        return '\n'.join(lines)

    def add_todo(self, title: str, status: str = Status.INBOX.value,
                 due: Optional[str] = None, priority: str = Priority.MEDIUM.value,
                 project: Optional[str] = None, tags: Optional[List[str]] = None,
                 description: Optional[str] = None) -> TodoItem:
        """Add a new todo"""
        now = datetime.now()
        todo_id = self._generate_todo_id(now)

        # Parse due date if provided
        if due:
            due = self._parse_due_date(due)

        todo = TodoItem(
            id=todo_id,
            title=title,
            status=status,
            created=now.isoformat(),
            due=due,
            priority=priority,
            project=project,
            tags=tags or [],
            description=description,
            file=f"daily/{now.strftime('%Y-%m-%d')}.md"
        )

        # Save to index
        self.todos[todo_id] = todo.to_dict()
        self._save_index()

        # Update daily file
        self._update_daily_file(todo)

        return todo

    def _parse_due_date(self, due_str: str) -> str:
        """Parse due date string (tomorrow, YYYY-MM-DD, etc.)"""
        due_str = due_str.lower().strip()

        if due_str == 'today':
            due_date = datetime.now()
        elif due_str == 'tomorrow':
            due_date = datetime.now() + timedelta(days=1)
        elif due_str.startswith('+'):
            # +N days
            days = int(due_str[1:])
            due_date = datetime.now() + timedelta(days=days)
        else:
            # Assume YYYY-MM-DD format
            due_date = datetime.strptime(due_str, '%Y-%m-%d')

        # Set to end of day
        due_date = due_date.replace(hour=23, minute=59, second=59)
        return due_date.isoformat()

    def get_todo(self, todo_id: str) -> Optional[TodoItem]:
        """Get a todo by ID"""
        todo_data = self.todos.get(todo_id)
        if not todo_data:
            return None

        return TodoItem(**todo_data)

    def list_todos(self, status: Optional[str] = None, project: Optional[str] = None,
                   priority: Optional[str] = None, overdue: bool = False,
                   due_today: bool = False) -> List[TodoItem]:
        """List todos with filters"""
        todos = []

        for todo_data in self.todos.values():
            todo = TodoItem(**todo_data)

            # Apply filters
            if status and todo.status != status:
                continue

            if project and todo.project != project:
                continue

            if priority and todo.priority != priority:
                continue

            if overdue and not todo.is_overdue():
                continue

            if due_today:
                days_until = todo.days_until_due()
                if days_until != 0:
                    continue

            todos.append(todo)

        # Sort by due date, then priority
        def sort_key(t):
            priority_order = {Priority.HIGH.value: 0, Priority.MEDIUM.value: 1, Priority.LOW.value: 2}
            return (
                t.due or '9999-12-31',
                priority_order.get(t.priority, 3)
            )

        return sorted(todos, key=sort_key)

    def update_todo(self, todo_id: str, **kwargs) -> Optional[TodoItem]:
        """Update a todo"""
        if todo_id not in self.todos:
            return None

        # Update fields
        for key, value in kwargs.items():
            if value is not None:
                if key == 'due' and isinstance(value, str):
                    value = self._parse_due_date(value)
                self.todos[todo_id][key] = value

        self._save_index()

        # Update daily file
        todo = TodoItem(**self.todos[todo_id])
        self._update_daily_file(todo)

        return todo

    def move_todo(self, todo_id: str, new_status: str) -> Optional[TodoItem]:
        """Move todo to new status"""
        return self.update_todo(todo_id, status=new_status)

    def complete_todo(self, todo_id: str) -> Optional[TodoItem]:
        """Mark todo as completed"""
        return self.update_todo(
            todo_id,
            status=Status.DONE.value,
            completed=datetime.now().isoformat()
        )


class TodoCLI:
    """CLI interface for todo management"""

    def __init__(self):
        self.base_path = Path('/Users/byron/projects')
        self.store = TodoStore(self.base_path)

    def run(self, args: List[str]):
        """Run CLI"""
        parser = argparse.ArgumentParser(
            description='IMMUNOS GTD Todo System',
            formatter_class=argparse.RawDescriptionHelpFormatter
        )

        subparsers = parser.add_subparsers(dest='command', help='Commands')

        # Add command
        add_parser = subparsers.add_parser('add', help='Add a new todo')
        add_parser.add_argument('title', help='Todo title')
        add_parser.add_argument('--due', help='Due date (today, tomorrow, +N, YYYY-MM-DD)')
        add_parser.add_argument('--priority', choices=['high', 'medium', 'low'], default='medium')
        add_parser.add_argument('--project', help='Project name')
        add_parser.add_argument('--tags', help='Comma-separated tags')
        add_parser.add_argument('--description', help='Todo description')

        # List command
        list_parser = subparsers.add_parser('list', help='List todos')
        list_parser.add_argument('--status', choices=['inbox', 'next', 'waiting', 'someday', 'done'])
        list_parser.add_argument('--project', help='Filter by project')
        list_parser.add_argument('--priority', choices=['high', 'medium', 'low'])
        list_parser.add_argument('--overdue', action='store_true', help='Show only overdue')
        list_parser.add_argument('--today', action='store_true', help='Show due today')

        # Shortcuts
        subparsers.add_parser('next', help='Show next actions')
        subparsers.add_parser('overdue', help='Show overdue todos')
        subparsers.add_parser('inbox', help='Show inbox')

        # Move command
        move_parser = subparsers.add_parser('move', help='Move todo to new status')
        move_parser.add_argument('todo_id', help='Todo ID')
        move_parser.add_argument('status', choices=['inbox', 'next', 'waiting', 'someday'])

        # Complete command
        complete_parser = subparsers.add_parser('complete', help='Mark todo as done')
        complete_parser.add_argument('todo_id', help='Todo ID')

        # Stats command
        subparsers.add_parser('stats', help='Show statistics')

        # Parse args
        parsed = parser.parse_args(args)

        if not parsed.command:
            parser.print_help()
            return

        # Execute command
        if parsed.command == 'add':
            self.cmd_add(parsed)
        elif parsed.command == 'list':
            self.cmd_list(parsed)
        elif parsed.command == 'next':
            self.cmd_next()
        elif parsed.command == 'overdue':
            self.cmd_overdue()
        elif parsed.command == 'inbox':
            self.cmd_inbox()
        elif parsed.command == 'move':
            self.cmd_move(parsed)
        elif parsed.command == 'complete':
            self.cmd_complete(parsed)
        elif parsed.command == 'stats':
            self.cmd_stats()

    def cmd_add(self, args):
        """Add todo command"""
        tags = args.tags.split(',') if args.tags else []

        todo = self.store.add_todo(
            title=args.title,
            due=args.due,
            priority=args.priority,
            project=args.project,
            tags=tags,
            description=args.description
        )

        print(f"✅ Created: {todo.id}")
        print(f"   {todo.title}")
        if todo.due:
            due_str = datetime.fromisoformat(todo.due).strftime('%Y-%m-%d')
            print(f"   Due: {due_str}")

    def cmd_list(self, args):
        """List todos command"""
        todos = self.store.list_todos(
            status=args.status,
            project=args.project,
            priority=args.priority,
            overdue=args.overdue,
            due_today=args.today
        )

        if not todos:
            print("No todos found.")
            return

        self._print_todo_list(todos)

    def cmd_next(self):
        """Show next actions"""
        todos = self.store.list_todos(status=Status.NEXT.value)
        print("⚡ Next Actions:")
        print()
        self._print_todo_list(todos)

    def cmd_overdue(self):
        """Show overdue todos"""
        todos = self.store.list_todos(overdue=True)
        print("🔴 Overdue:")
        print()
        self._print_todo_list(todos)

    def cmd_inbox(self):
        """Show inbox"""
        todos = self.store.list_todos(status=Status.INBOX.value)
        print("📥 Inbox:")
        print()
        self._print_todo_list(todos)

    def cmd_move(self, args):
        """Move todo command"""
        todo = self.store.move_todo(args.todo_id, args.status)

        if not todo:
            print(f"❌ Todo not found: {args.todo_id}")
            return

        print(f"✅ Moved to {args.status}: {todo.id}")
        print(f"   {todo.title}")

    def cmd_complete(self, args):
        """Complete todo command"""
        todo = self.store.complete_todo(args.todo_id)

        if not todo:
            print(f"❌ Todo not found: {args.todo_id}")
            return

        print(f"✅ Completed: {todo.id}")
        print(f"   {todo.title}")

    def cmd_stats(self):
        """Show statistics"""
        stats = self.store._calculate_stats()

        print("📊 Todo Statistics")
        print()
        print(f"Total: {stats['total']}")
        print()

        print("By Status:")
        for status, count in stats['by_status'].items():
            print(f"  {status}: {count}")
        print()

        print("By Priority:")
        for priority, count in stats['by_priority'].items():
            print(f"  {priority}: {count}")
        print()

        print("Due Dates:")
        print(f"  Overdue: {stats['overdue']}")
        print(f"  Due today: {stats['due_today']}")
        print(f"  Due this week: {stats['due_this_week']}")

    def _print_todo_list(self, todos: List[TodoItem]):
        """Print formatted todo list"""
        for todo in todos:
            # Priority indicator
            priority_icon = {
                Priority.HIGH.value: '🔴',
                Priority.MEDIUM.value: '🟡',
                Priority.LOW.value: '🟢'
            }.get(todo.priority, '⚪')

            # Status
            status_short = todo.status[0].upper()

            # Due date
            due_str = ''
            if todo.due:
                days = todo.days_until_due()
                if days is not None:
                    if days < 0:
                        due_str = f" (overdue {abs(days)}d)"
                    elif days == 0:
                        due_str = " (due today)"
                    else:
                        due_str = f" (due in {days}d)"

            # Project
            project_str = f" [{todo.project}]" if todo.project else ""

            print(f"{priority_icon} {todo.id} ({status_short}): {todo.title}{project_str}{due_str}")


def main():
    """Main entry point"""
    cli = TodoCLI()
    cli.run(sys.argv[1:])


if __name__ == '__main__':
    main()
