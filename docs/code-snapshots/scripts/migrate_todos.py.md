---
source: /Users/byron/projects/scripts/migrate_todos.py
relative: scripts/migrate_todos.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
IMMUNOS Todo Migration Script

Migrate existing TODO.md files to the new GTD todo system.

Usage:
    python migrate_todos.py <TODO_FILE> --project <PROJECT_NAME> [--dry-run]

Example:
    python migrate_todos.py prion-clock/TODO.md --project prion-clock --dry-run
    python migrate_todos.py immunos81/TODO.md --project immunos
"""

import re
import sys
import argparse
from pathlib import Path
from datetime import datetime, timedelta
from typing import List, Dict, Tuple
import json

# Import the todo system
sys.path.insert(0, str(Path(__file__).parent))
from immunos_todo import TodoStore, Priority


class TodoMigrator:
    """Migrate TODO.md files to GTD system"""

    # Priority emoji to priority mapping
    PRIORITY_MAP = {
        '🔴': Priority.HIGH.value,      # Critical
        '🟡': Priority.HIGH.value,      # High
        '🟢': Priority.MEDIUM.value,    # Medium
        '🔵': Priority.LOW.value,       # Optional/Low
    }

    # Section to GTD status mapping
    STATUS_MAP = {
        'critical': 'next',
        'high': 'next',
        'medium': 'next',
        'optional': 'someday',
        'future': 'someday',
        'completed': 'done',
    }

    def __init__(self, todo_file: Path, project: str, dry_run: bool = False):
        self.todo_file = todo_file
        self.project = project
        self.dry_run = dry_run
        self.store = TodoStore(Path('/Users/byron/projects'))

        self.todos_to_create: List[Dict] = []
        self.migration_log: List[str] = []

    def parse_section_header(self, line: str) -> Tuple[str, str]:
        """Parse section header to get priority and name"""
        # Match: ## 🔴 Critical Priority
        match = re.match(r'##\s*([🔴🟡🟢🔵])?\s*(.+)', line)
        if match:
            emoji, name = match.groups()
            priority_emoji = emoji or '🟢'
            section_name = name.strip().lower()

            # Extract priority keyword
            for keyword in ['critical', 'high', 'medium', 'optional', 'low', 'future', 'completed']:
                if keyword in section_name:
                    return keyword, self.PRIORITY_MAP.get(priority_emoji, Priority.MEDIUM.value)

            return 'medium', self.PRIORITY_MAP.get(priority_emoji, Priority.MEDIUM.value)

        return 'medium', Priority.MEDIUM.value

    def parse_checkbox_item(self, line: str) -> Dict:
        """Parse a checkbox line to extract todo information"""
        # Match: - [ ] **Task title**
        # Or:    - [x] **Task title**
        # Or:    - [ ] Task title (no bold)

        checkbox_match = re.match(r'-\s*\[([ xX])\]\s*(.+)', line)
        if not checkbox_match:
            return None

        checked, content = checkbox_match.groups()
        completed = checked.lower() == 'x'

        # Extract title (remove bold markers if present)
        title_match = re.match(r'\*\*(.+?)\*\*', content)
        if title_match:
            title = title_match.group(1)
            # Get description after title
            description = content[len(title) + 4:].strip()
        else:
            # No bold, use first line as title
            title = content.split('\n')[0].strip()
            description = None

        # Clean up title
        title = title.strip()
        if title.startswith('**') and title.endswith('**'):
            title = title[2:-2]

        return {
            'title': title,
            'description': description,
            'completed': completed,
        }

    def parse_todo_file(self) -> List[Dict]:
        """Parse TODO.md file and extract todos"""
        if not self.todo_file.exists():
            print(f"❌ File not found: {self.todo_file}")
            sys.exit(1)

        with open(self.todo_file, 'r') as f:
            lines = f.readlines()

        current_section = 'medium'
        current_priority = Priority.MEDIUM.value
        todos = []

        i = 0
        while i < len(lines):
            line = lines[i].rstrip()

            # Check for section header
            if line.startswith('##'):
                current_section, current_priority = self.parse_section_header(line)
                self.migration_log.append(f"Found section: {line} → status={current_section}, priority={current_priority}")

            # Check for checkbox item
            elif line.strip().startswith('- ['):
                todo_info = self.parse_checkbox_item(line)

                if todo_info:
                    # Collect multi-line description
                    description_lines = []
                    if todo_info.get('description'):
                        description_lines.append(todo_info['description'])

                    # Look ahead for continuation lines
                    j = i + 1
                    while j < len(lines):
                        next_line = lines[j].rstrip()
                        # Stop at next todo or section
                        if next_line.startswith('##') or next_line.strip().startswith('- ['):
                            break
                        # Add continuation line
                        if next_line.strip() and next_line.startswith('  '):
                            description_lines.append(next_line.strip())
                        elif not next_line.strip():
                            break
                        j += 1

                    description = '\n'.join(description_lines) if description_lines else None

                    # Determine status
                    if todo_info['completed']:
                        status = 'done'
                    else:
                        status = self.STATUS_MAP.get(current_section, 'next')

                    # Determine due date (heuristic: high priority = 1 week, medium = 2 weeks, low = 1 month)
                    if status == 'next':
                        if current_priority == Priority.HIGH.value:
                            due_days = 7
                        elif current_priority == Priority.MEDIUM.value:
                            due_days = 14
                        else:
                            due_days = 30
                        due_date = (datetime.now() + timedelta(days=due_days)).isoformat()
                    else:
                        due_date = None

                    todos.append({
                        'title': todo_info['title'],
                        'description': description,
                        'status': status,
                        'priority': current_priority,
                        'due': due_date,
                        'project': self.project,
                        'completed': todo_info['completed'],
                    })

                    self.migration_log.append(f"  Parsed: [{status}] {todo_info['title'][:50]}...")

            i += 1

        return todos

    def create_todos(self, todos: List[Dict]):
        """Create todos in the GTD system"""
        created_count = 0
        skipped_count = 0

        for todo_data in todos:
            # Extract tags from description
            tags = []
            if todo_data.get('description'):
                # Look for hashtags
                hashtags = re.findall(r'#(\w+)', todo_data['description'])
                tags.extend(hashtags)

            # Add project as tag
            tags.append(todo_data['project'])
            tags.append('migrated')

            if self.dry_run:
                print(f"  [DRY RUN] Would create: {todo_data['title']}")
                print(f"    Status: {todo_data['status']}")
                print(f"    Priority: {todo_data['priority']}")
                if todo_data.get('due'):
                    due_str = datetime.fromisoformat(todo_data['due']).strftime('%Y-%m-%d')
                    print(f"    Due: {due_str}")
                created_count += 1
            else:
                try:
                    # Create the todo
                    if todo_data['completed']:
                        # Create as completed
                        new_todo = self.store.add_todo(
                            title=todo_data['title'],
                            status='done',
                            due=todo_data.get('due'),
                            priority=todo_data['priority'],
                            project=todo_data['project'],
                            tags=tags,
                            description=todo_data.get('description')
                        )
                        # Mark as completed
                        self.store.update_todo(
                            new_todo.id,
                            completed=datetime.now().isoformat()
                        )
                    else:
                        new_todo = self.store.add_todo(
                            title=todo_data['title'],
                            status=todo_data['status'],
                            due=todo_data.get('due'),
                            priority=todo_data['priority'],
                            project=todo_data['project'],
                            tags=tags,
                            description=todo_data.get('description')
                        )

                    print(f"  ✅ Created: {new_todo.id} - {todo_data['title'][:60]}")
                    created_count += 1
                except Exception as e:
                    print(f"  ❌ Failed: {todo_data['title'][:60]}")
                    print(f"     Error: {e}")
                    skipped_count += 1

        return created_count, skipped_count

    def archive_original(self):
        """Archive the original TODO.md file"""
        if self.dry_run:
            print(f"\n[DRY RUN] Would archive: {self.todo_file} → {self.todo_file}.migrated")
            return

        archive_path = self.todo_file.with_suffix('.md.migrated')
        self.todo_file.rename(archive_path)
        print(f"\n✅ Original archived to: {archive_path}")

    def save_migration_log(self):
        """Save migration log"""
        log_file = self.todo_file.parent / f"migration_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"

        if self.dry_run:
            print(f"\n[DRY RUN] Would save log to: {log_file}")
            return

        with open(log_file, 'w') as f:
            f.write(f"Migration Log - {datetime.now().isoformat()}\n")
            f.write(f"Source: {self.todo_file}\n")
            f.write(f"Project: {self.project}\n")
            f.write("\n" + "="*60 + "\n\n")
            f.write('\n'.join(self.migration_log))

        print(f"📄 Log saved to: {log_file}")

    def migrate(self):
        """Run the migration"""
        print(f"{'='*60}")
        print(f"TODO Migration: {self.todo_file}")
        print(f"Project: {self.project}")
        print(f"Dry Run: {self.dry_run}")
        print(f"{'='*60}\n")

        # Parse the file
        print("📖 Parsing TODO.md file...")
        todos = self.parse_todo_file()
        print(f"   Found {len(todos)} todos")

        # Show summary
        by_status = {}
        by_priority = {}
        for todo in todos:
            by_status[todo['status']] = by_status.get(todo['status'], 0) + 1
            by_priority[todo['priority']] = by_priority.get(todo['priority'], 0) + 1

        print(f"\n📊 Summary:")
        print(f"   By Status: {dict(by_status)}")
        print(f"   By Priority: {dict(by_priority)}")
        print()

        # Create todos
        print("📝 Creating todos in GTD system...\n")
        created, skipped = self.create_todos(todos)

        print(f"\n{'='*60}")
        print(f"Migration {'Preview' if self.dry_run else 'Complete'}!")
        print(f"{'='*60}")
        print(f"  Created: {created}")
        print(f"  Skipped: {skipped}")
        print(f"  Total: {len(todos)}")

        if not self.dry_run:
            # Archive original
            self.archive_original()

            # Save log
            self.save_migration_log()

            print(f"\n✅ Migration complete!")
            print(f"\nNext steps:")
            print(f"  python scripts/immunos_todo.py list --project {self.project}")
            print(f"  python scripts/immunos_todo.py next")
        else:
            print(f"\n💡 This was a dry run. Run without --dry-run to actually migrate.")


def main():
    parser = argparse.ArgumentParser(
        description='Migrate TODO.md files to GTD system',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Preview migration
  python migrate_todos.py prion-clock/TODO.md --project prion-clock --dry-run

  # Execute migration
  python migrate_todos.py prion-clock/TODO.md --project prion-clock
  python migrate_todos.py immunos81/TODO.md --project immunos
"""
    )

    parser.add_argument('todo_file', type=str, help='Path to TODO.md file')
    parser.add_argument('--project', type=str, required=True, help='Project name')
    parser.add_argument('--dry-run', action='store_true', help='Preview without making changes')

    args = parser.parse_args()

    # Resolve path
    todo_file = Path(args.todo_file)
    if not todo_file.is_absolute():
        todo_file = Path('/Users/byron/projects') / todo_file

    # Run migration
    migrator = TodoMigrator(todo_file, args.project, args.dry_run)
    migrator.migrate()


if __name__ == '__main__':
    main()

```
