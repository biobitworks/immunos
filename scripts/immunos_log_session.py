#!/usr/bin/env python3
"""
IMMUNOS Universal Session Logger

Logs AI model sessions to universal conversation database.
Can be called by any model (Claude, Ollama models, etc.)

Usage:
    python immunos_log_session.py --model "claude-sonnet-4.5" --summary "Fixed Prion Clock citations"
    python immunos_log_session.py --model "qwen-coder" --summary "Implemented NK scanner" --files "scripts/immunos_nk_scan.py"
"""

import json
import argparse
from datetime import datetime
from pathlib import Path
from typing import List, Optional


class SessionLogger:
    """Log AI model sessions to universal database"""

    def __init__(self, base_path: str = '/Users/byron/projects'):
        self.base_path = Path(base_path)
        self.log_file = self.base_path / '.immunos' / 'conversations' / 'universal-log.jsonl'
        self.log_file.parent.mkdir(parents=True, exist_ok=True)

    def log_session(
        self,
        model: str,
        summary: str,
        tags: Optional[List[str]] = None,
        files_modified: Optional[List[str]] = None,
        priority: str = 'medium',
        session_id: Optional[str] = None
    ):
        """Log a session to the universal conversation database"""

        # Generate session ID if not provided
        if not session_id:
            model_short = model.split('-')[0].lower()  # e.g., "claude" from "claude-sonnet-4.5"
            session_id = f"{model_short}-{datetime.now().strftime('%Y%m%d-%H%M%S')}"

        # Create log entry
        entry = {
            'timestamp': datetime.utcnow().isoformat() + 'Z',
            'model': model,
            'session_id': session_id,
            'user': 'Byron P. Lee',
            'summary': summary,
            'priority': priority,
            'tags': tags or [],
            'files_modified': files_modified or []
        }

        # Append to log file (JSON Lines format)
        with open(self.log_file, 'a') as f:
            f.write(json.dumps(entry) + '\n')

        return session_id

    def get_recent_sessions(self, model: Optional[str] = None, limit: int = 10):
        """Get recent sessions, optionally filtered by model"""
        if not self.log_file.exists():
            return []

        sessions = []
        with open(self.log_file, 'r') as f:
            for line in f:
                if line.strip():
                    entry = json.loads(line)
                    if model is None or entry['model'] == model:
                        sessions.append(entry)

        # Return most recent first
        return sorted(sessions, key=lambda x: x['timestamp'], reverse=True)[:limit]

    def search_sessions(self, query: str):
        """Search sessions by summary or tags"""
        if not self.log_file.exists():
            return []

        matches = []
        query_lower = query.lower()

        with open(self.log_file, 'r') as f:
            for line in f:
                if line.strip():
                    entry = json.loads(line)
                    # Search in summary and tags
                    if (query_lower in entry['summary'].lower() or
                        any(query_lower in tag.lower() for tag in entry['tags'])):
                        matches.append(entry)

        return sorted(matches, key=lambda x: x['timestamp'], reverse=True)


def main():
    parser = argparse.ArgumentParser(description='IMMUNOS Universal Session Logger')
    parser.add_argument('--model', required=True, help='Model name (e.g., "claude-sonnet-4.5", "qwen-coder")')
    parser.add_argument('--summary', required=True, help='One-line summary of session')
    parser.add_argument('--tags', help='Comma-separated tags (e.g., "prion-clock,citation,fix")')
    parser.add_argument('--files', help='Comma-separated list of files modified')
    parser.add_argument('--priority', choices=['high', 'medium', 'low'], default='medium')
    parser.add_argument('--session-id', help='Optional custom session ID')

    # Query options
    parser.add_argument('--recent', type=int, metavar='N', help='Show N recent sessions')
    parser.add_argument('--search', help='Search sessions by keyword')
    parser.add_argument('--filter-model', help='Filter by model name')

    args = parser.parse_args()

    logger = SessionLogger()

    # Query mode
    if args.recent or args.search:
        if args.search:
            sessions = logger.search_sessions(args.search)
        else:
            sessions = logger.get_recent_sessions(
                model=args.filter_model,
                limit=args.recent or 10
            )

        if not sessions:
            print("No sessions found")
            return

        print(f"Found {len(sessions)} session(s):\n")
        for session in sessions:
            timestamp = datetime.fromisoformat(session['timestamp'].replace('Z', '+00:00'))
            print(f"[{timestamp.strftime('%Y-%m-%d %H:%M')}] {session['model']}")
            print(f"  ID: {session['session_id']}")
            print(f"  Summary: {session['summary']}")
            if session['tags']:
                print(f"  Tags: {', '.join(session['tags'])}")
            if session['files_modified']:
                print(f"  Files: {', '.join(session['files_modified'])}")
            print()

        return

    # Log mode
    tags = args.tags.split(',') if args.tags else []
    files = args.files.split(',') if args.files else []

    session_id = logger.log_session(
        model=args.model,
        summary=args.summary,
        tags=tags,
        files_modified=files,
        priority=args.priority,
        session_id=args.session_id
    )

    print(f"✓ Session logged: {session_id}")
    print(f"  Model: {args.model}")
    print(f"  Summary: {args.summary}")
    if tags:
        print(f"  Tags: {', '.join(tags)}")
    if files:
        print(f"  Files modified: {', '.join(files)}")


if __name__ == '__main__':
    main()
