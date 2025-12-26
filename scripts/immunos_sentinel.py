#!/usr/bin/env python3
"""
IMMUNOS Sentinel Agent - Log Watcher & Ticket Generator

Continuously monitors log files, stdout/stderr streams, and server output.
Detects anomalies, creates tickets, and routes to appropriate IMMUNOS agents.

Usage:
    # Watch a log file
    python immunos_sentinel.py watch /var/log/app.log

    # Watch dashboard logs
    python immunos_sentinel.py watch-dashboard

    # Start as background service
    python immunos_sentinel.py service --port 5002
"""

import os
import sys
import re
import json
import time
import hashlib
import argparse
import threading
import sqlite3
from pathlib import Path
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple, Callable
from dataclasses import dataclass, field, asdict
from enum import Enum
from collections import deque

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))


class TicketSeverity(Enum):
    CRITICAL = "critical"
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    INFO = "info"


class TicketStatus(Enum):
    OPEN = "open"
    IN_PROGRESS = "in_progress"
    RESOLVED = "resolved"
    CLOSED = "closed"
    IGNORED = "ignored"


class TicketCategory(Enum):
    ERROR = "error"
    SECURITY = "security"
    PERFORMANCE = "performance"
    WARNING = "warning"
    ANOMALY = "anomaly"
    INFO = "info"


@dataclass
class LogPattern:
    """Pattern to match in logs"""
    name: str
    pattern: str
    regex: re.Pattern
    severity: TicketSeverity
    category: TicketCategory
    description: str
    agent: str = "nk_cell"  # Which agent should handle

    def __post_init__(self):
        if isinstance(self.pattern, str):
            self.regex = re.compile(self.pattern, re.IGNORECASE)


@dataclass
class LogEntry:
    """Parsed log entry"""
    timestamp: datetime
    level: str
    message: str
    source: str
    line_number: int
    raw: str
    metadata: Dict = field(default_factory=dict)


# Token cost estimates per model (per 1K tokens)
MODEL_TOKEN_COSTS = {
    'claude-opus': {'input': 0.015, 'output': 0.075},
    'claude-sonnet': {'input': 0.003, 'output': 0.015},
    'claude-haiku': {'input': 0.00025, 'output': 0.00125},
    'gpt-4': {'input': 0.03, 'output': 0.06},
    'gpt-4-turbo': {'input': 0.01, 'output': 0.03},
    'gpt-3.5-turbo': {'input': 0.0005, 'output': 0.0015},
    'ollama-qwen2.5-coder': {'input': 0.0, 'output': 0.0},  # Free local
    'ollama-deepseek-r1': {'input': 0.0, 'output': 0.0},
    'ollama-qwen2.5': {'input': 0.0, 'output': 0.0},
}

# Average tokens per task type
TASK_TOKEN_ESTIMATES = {
    'error_analysis': {'input': 500, 'output': 300},
    'security_scan': {'input': 800, 'output': 500},
    'code_review': {'input': 1000, 'output': 600},
    'log_analysis': {'input': 400, 'output': 250},
    'ticket_triage': {'input': 200, 'output': 150},
    'resolution': {'input': 600, 'output': 400},
}


@dataclass
class TokenEstimate:
    """Token usage estimate for a task"""
    input_tokens: int = 0
    output_tokens: int = 0
    model: str = "ollama-qwen2.5-coder"
    estimated_cost: float = 0.0
    task_type: str = "ticket_triage"

    def calculate_cost(self) -> float:
        """Calculate cost based on model"""
        costs = MODEL_TOKEN_COSTS.get(self.model, {'input': 0, 'output': 0})
        self.estimated_cost = (
            (self.input_tokens / 1000) * costs['input'] +
            (self.output_tokens / 1000) * costs['output']
        )
        return self.estimated_cost

    def to_dict(self) -> Dict:
        return {
            'input_tokens': self.input_tokens,
            'output_tokens': self.output_tokens,
            'total_tokens': self.input_tokens + self.output_tokens,
            'model': self.model,
            'estimated_cost': self.estimated_cost,
            'task_type': self.task_type
        }

    @classmethod
    def for_task(cls, task_type: str, model: str = "ollama-qwen2.5-coder") -> 'TokenEstimate':
        """Create estimate for a task type"""
        tokens = TASK_TOKEN_ESTIMATES.get(task_type, {'input': 300, 'output': 200})
        estimate = cls(
            input_tokens=tokens['input'],
            output_tokens=tokens['output'],
            model=model,
            task_type=task_type
        )
        estimate.calculate_cost()
        return estimate


@dataclass
class Ticket:
    """Issue ticket created from log detection"""
    id: Optional[int] = None
    title: str = ""
    description: str = ""
    severity: TicketSeverity = TicketSeverity.MEDIUM
    category: TicketCategory = TicketCategory.ERROR
    status: TicketStatus = TicketStatus.OPEN
    source_file: str = ""
    source_line: int = 0
    pattern_matched: str = ""
    log_entries: List[str] = field(default_factory=list)
    assigned_agent: str = ""
    created_at: datetime = field(default_factory=datetime.now)
    updated_at: datetime = field(default_factory=datetime.now)
    resolved_at: Optional[datetime] = None
    resolution_note: str = ""
    context: Dict = field(default_factory=dict)
    token_estimate: Optional[TokenEstimate] = None
    actual_tokens_used: int = 0

    def estimate_tokens(self, model: str = "ollama-qwen2.5-coder"):
        """Estimate tokens needed to resolve this ticket"""
        # Map category to task type
        task_map = {
            TicketCategory.ERROR: 'error_analysis',
            TicketCategory.SECURITY: 'security_scan',
            TicketCategory.PERFORMANCE: 'log_analysis',
            TicketCategory.WARNING: 'ticket_triage',
            TicketCategory.ANOMALY: 'error_analysis',
            TicketCategory.INFO: 'ticket_triage',
        }
        task_type = task_map.get(self.category, 'ticket_triage')

        # Adjust based on severity
        multiplier = {
            TicketSeverity.CRITICAL: 2.0,
            TicketSeverity.HIGH: 1.5,
            TicketSeverity.MEDIUM: 1.0,
            TicketSeverity.LOW: 0.7,
            TicketSeverity.INFO: 0.5,
        }.get(self.severity, 1.0)

        estimate = TokenEstimate.for_task(task_type, model)
        estimate.input_tokens = int(estimate.input_tokens * multiplier)
        estimate.output_tokens = int(estimate.output_tokens * multiplier)
        estimate.calculate_cost()
        self.token_estimate = estimate
        return estimate

    def to_dict(self) -> Dict:
        return {
            'id': self.id,
            'title': self.title,
            'description': self.description,
            'severity': self.severity.value,
            'category': self.category.value,
            'status': self.status.value,
            'source_file': self.source_file,
            'source_line': self.source_line,
            'pattern_matched': self.pattern_matched,
            'log_entries': self.log_entries,
            'assigned_agent': self.assigned_agent,
            'created_at': self.created_at.isoformat(),
            'updated_at': self.updated_at.isoformat(),
            'resolved_at': self.resolved_at.isoformat() if self.resolved_at else None,
            'resolution_note': self.resolution_note,
            'context': self.context,
            'token_estimate': self.token_estimate.to_dict() if self.token_estimate else None,
            'actual_tokens_used': self.actual_tokens_used
        }


class SentinelPatterns:
    """Built-in patterns for common log issues"""

    # Error patterns
    ERROR_PATTERNS = [
        LogPattern(
            name="python_exception",
            pattern=r"(Traceback|Exception|Error):\s*(.+)",
            regex=None,
            severity=TicketSeverity.HIGH,
            category=TicketCategory.ERROR,
            description="Python exception detected",
            agent="nk_cell"
        ),
        LogPattern(
            name="http_500",
            pattern=r"HTTP/\d\.\d\"\s+5\d{2}",
            regex=None,
            severity=TicketSeverity.HIGH,
            category=TicketCategory.ERROR,
            description="HTTP 500 server error",
            agent="nk_cell"
        ),
        LogPattern(
            name="http_4xx",
            pattern=r"HTTP/\d\.\d\"\s+4\d{2}",
            regex=None,
            severity=TicketSeverity.LOW,
            category=TicketCategory.WARNING,
            description="HTTP 4xx client error",
            agent="b_cell"
        ),
        LogPattern(
            name="connection_refused",
            pattern=r"(Connection refused|ECONNREFUSED|ConnectionError)",
            regex=None,
            severity=TicketSeverity.MEDIUM,
            category=TicketCategory.ERROR,
            description="Connection refused error",
            agent="nk_cell"
        ),
        LogPattern(
            name="timeout",
            pattern=r"(timeout|timed out|TimeoutError)",
            regex=None,
            severity=TicketSeverity.MEDIUM,
            category=TicketCategory.PERFORMANCE,
            description="Timeout detected",
            agent="nk_cell"
        ),
    ]

    # Security patterns
    SECURITY_PATTERNS = [
        LogPattern(
            name="auth_failure",
            pattern=r"(authentication failed|invalid password|login failed|unauthorized)",
            regex=None,
            severity=TicketSeverity.HIGH,
            category=TicketCategory.SECURITY,
            description="Authentication failure",
            agent="b_cell"
        ),
        LogPattern(
            name="sql_injection",
            pattern=r"(SELECT.*FROM.*WHERE|DROP TABLE|UNION SELECT|OR 1=1)",
            regex=None,
            severity=TicketSeverity.CRITICAL,
            category=TicketCategory.SECURITY,
            description="Possible SQL injection attempt",
            agent="b_cell"
        ),
        LogPattern(
            name="xss_attempt",
            pattern=r"(<script|javascript:|onerror=|onload=)",
            regex=None,
            severity=TicketSeverity.CRITICAL,
            category=TicketCategory.SECURITY,
            description="Possible XSS attempt",
            agent="b_cell"
        ),
        LogPattern(
            name="path_traversal",
            pattern=r"(\.\./|\.\.\\|%2e%2e)",
            regex=None,
            severity=TicketSeverity.CRITICAL,
            category=TicketCategory.SECURITY,
            description="Possible path traversal attempt",
            agent="b_cell"
        ),
        LogPattern(
            name="sensitive_data",
            pattern=r"(password|api_key|secret|token|credential).*=.*['\"][^'\"]+['\"]",
            regex=None,
            severity=TicketSeverity.HIGH,
            category=TicketCategory.SECURITY,
            description="Possible sensitive data in logs",
            agent="b_cell"
        ),
    ]

    # Performance patterns
    PERFORMANCE_PATTERNS = [
        LogPattern(
            name="slow_query",
            pattern=r"(slow query|query took|execution time:?\s*\d+\s*(ms|s|seconds))",
            regex=None,
            severity=TicketSeverity.MEDIUM,
            category=TicketCategory.PERFORMANCE,
            description="Slow query detected",
            agent="dendritic"
        ),
        LogPattern(
            name="memory_warning",
            pattern=r"(memory warning|out of memory|MemoryError|heap)",
            regex=None,
            severity=TicketSeverity.HIGH,
            category=TicketCategory.PERFORMANCE,
            description="Memory issue detected",
            agent="nk_cell"
        ),
        LogPattern(
            name="high_cpu",
            pattern=r"(cpu.*high|cpu usage|load average)",
            regex=None,
            severity=TicketSeverity.MEDIUM,
            category=TicketCategory.PERFORMANCE,
            description="High CPU usage",
            agent="dendritic"
        ),
    ]

    # Warning patterns
    WARNING_PATTERNS = [
        LogPattern(
            name="deprecation",
            pattern=r"(deprecated|deprecation warning|will be removed)",
            regex=None,
            severity=TicketSeverity.LOW,
            category=TicketCategory.WARNING,
            description="Deprecation warning",
            agent="dendritic"
        ),
        LogPattern(
            name="werkzeug_warning",
            pattern=r"WARNING.*production deployment",
            regex=None,
            severity=TicketSeverity.INFO,
            category=TicketCategory.WARNING,
            description="Werkzeug production warning",
            agent="dendritic"
        ),
    ]

    @classmethod
    def get_all_patterns(cls) -> List[LogPattern]:
        """Get all built-in patterns"""
        patterns = []
        patterns.extend(cls.ERROR_PATTERNS)
        patterns.extend(cls.SECURITY_PATTERNS)
        patterns.extend(cls.PERFORMANCE_PATTERNS)
        patterns.extend(cls.WARNING_PATTERNS)
        return patterns


class TicketStore:
    """SQLite-backed ticket storage"""

    def __init__(self, db_path: str):
        self.db_path = db_path
        self._init_db()

    def _init_db(self):
        """Initialize ticket tables"""
        conn = sqlite3.connect(self.db_path)
        conn.executescript('''
            -- Tickets Table
            CREATE TABLE IF NOT EXISTS tickets (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                title TEXT NOT NULL,
                description TEXT,
                severity TEXT NOT NULL,
                category TEXT NOT NULL,
                status TEXT NOT NULL DEFAULT 'open',
                source_file TEXT,
                source_line INTEGER,
                pattern_matched TEXT,
                log_entries_json TEXT,
                assigned_agent TEXT,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                updated_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                resolved_at DATETIME,
                resolution_note TEXT,
                context_json TEXT,
                estimated_input_tokens INTEGER DEFAULT 0,
                estimated_output_tokens INTEGER DEFAULT 0,
                estimated_model TEXT,
                estimated_cost REAL DEFAULT 0.0,
                actual_tokens_used INTEGER DEFAULT 0
            );

            CREATE INDEX IF NOT EXISTS idx_tickets_status ON tickets(status);
            CREATE INDEX IF NOT EXISTS idx_tickets_severity ON tickets(severity);
            CREATE INDEX IF NOT EXISTS idx_tickets_created ON tickets(created_at DESC);

            -- Sentinel Log Stream
            CREATE TABLE IF NOT EXISTS sentinel_log_stream (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                source TEXT NOT NULL,
                level TEXT,
                message TEXT NOT NULL,
                raw_line TEXT,
                line_number INTEGER,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP
            );

            CREATE INDEX IF NOT EXISTS idx_log_stream_source ON sentinel_log_stream(source, created_at DESC);

            -- Sentinel Watchers
            CREATE TABLE IF NOT EXISTS sentinel_watchers (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT NOT NULL UNIQUE,
                source_type TEXT NOT NULL,
                source_path TEXT,
                enabled BOOLEAN DEFAULT 1,
                patterns_json TEXT,
                last_position INTEGER DEFAULT 0,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
                updated_at DATETIME DEFAULT CURRENT_TIMESTAMP
            );
        ''')
        conn.commit()
        conn.close()

    def create_ticket(self, ticket: Ticket) -> int:
        """Create a new ticket"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Get token estimate values
        est_input = ticket.token_estimate.input_tokens if ticket.token_estimate else 0
        est_output = ticket.token_estimate.output_tokens if ticket.token_estimate else 0
        est_model = ticket.token_estimate.model if ticket.token_estimate else None
        est_cost = ticket.token_estimate.estimated_cost if ticket.token_estimate else 0.0

        cursor.execute('''
            INSERT INTO tickets (
                title, description, severity, category, status,
                source_file, source_line, pattern_matched,
                log_entries_json, assigned_agent, context_json,
                estimated_input_tokens, estimated_output_tokens,
                estimated_model, estimated_cost, actual_tokens_used
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (
            ticket.title,
            ticket.description,
            ticket.severity.value,
            ticket.category.value,
            ticket.status.value,
            ticket.source_file,
            ticket.source_line,
            ticket.pattern_matched,
            json.dumps(ticket.log_entries),
            ticket.assigned_agent,
            json.dumps(ticket.context),
            est_input,
            est_output,
            est_model,
            est_cost,
            ticket.actual_tokens_used
        ))
        ticket_id = cursor.lastrowid
        conn.commit()
        conn.close()
        return ticket_id

    def get_ticket(self, ticket_id: int) -> Optional[Ticket]:
        """Get a ticket by ID"""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute('SELECT * FROM tickets WHERE id = ?', (ticket_id,))
        row = cursor.fetchone()
        conn.close()

        if row:
            return self._row_to_ticket(row)
        return None

    def get_tickets(self, status: Optional[str] = None,
                    severity: Optional[str] = None,
                    limit: int = 100) -> List[Ticket]:
        """Get tickets with optional filters"""
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        query = 'SELECT * FROM tickets WHERE 1=1'
        params = []

        if status:
            query += ' AND status = ?'
            params.append(status)
        if severity:
            query += ' AND severity = ?'
            params.append(severity)

        query += ' ORDER BY created_at DESC LIMIT ?'
        params.append(limit)

        cursor.execute(query, params)
        rows = cursor.fetchall()
        conn.close()

        return [self._row_to_ticket(row) for row in rows]

    def update_ticket(self, ticket_id: int, updates: Dict) -> bool:
        """Update a ticket"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        set_clauses = []
        params = []
        for key, value in updates.items():
            if key in ['log_entries', 'context']:
                set_clauses.append(f'{key}_json = ?')
                params.append(json.dumps(value))
            elif key == 'severity':
                set_clauses.append('severity = ?')
                params.append(value.value if isinstance(value, TicketSeverity) else value)
            elif key == 'status':
                set_clauses.append('status = ?')
                params.append(value.value if isinstance(value, TicketStatus) else value)
            elif key == 'category':
                set_clauses.append('category = ?')
                params.append(value.value if isinstance(value, TicketCategory) else value)
            else:
                set_clauses.append(f'{key} = ?')
                params.append(value)

        set_clauses.append('updated_at = CURRENT_TIMESTAMP')
        params.append(ticket_id)

        query = f"UPDATE tickets SET {', '.join(set_clauses)} WHERE id = ?"
        cursor.execute(query, params)
        conn.commit()
        conn.close()
        return cursor.rowcount > 0

    def log_entry(self, source: str, level: str, message: str,
                  raw_line: str = "", line_number: int = 0):
        """Log a stream entry for training data"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute('''
            INSERT INTO sentinel_log_stream (source, level, message, raw_line, line_number)
            VALUES (?, ?, ?, ?, ?)
        ''', (source, level, message, raw_line, line_number))
        conn.commit()
        conn.close()

    def _row_to_ticket(self, row) -> Ticket:
        """Convert database row to Ticket object"""
        # Reconstruct token estimate if present
        token_estimate = None
        try:
            if row['estimated_input_tokens'] or row['estimated_output_tokens']:
                token_estimate = TokenEstimate(
                    input_tokens=row['estimated_input_tokens'] or 0,
                    output_tokens=row['estimated_output_tokens'] or 0,
                    model=row['estimated_model'] or 'ollama-qwen2.5-coder',
                    estimated_cost=row['estimated_cost'] or 0.0
                )
        except (KeyError, TypeError):
            pass

        return Ticket(
            id=row['id'],
            title=row['title'],
            description=row['description'] or '',
            severity=TicketSeverity(row['severity']),
            category=TicketCategory(row['category']),
            status=TicketStatus(row['status']),
            source_file=row['source_file'] or '',
            source_line=row['source_line'] or 0,
            pattern_matched=row['pattern_matched'] or '',
            log_entries=json.loads(row['log_entries_json']) if row['log_entries_json'] else [],
            assigned_agent=row['assigned_agent'] or '',
            created_at=datetime.fromisoformat(row['created_at']) if row['created_at'] else datetime.now(),
            updated_at=datetime.fromisoformat(row['updated_at']) if row['updated_at'] else datetime.now(),
            resolved_at=datetime.fromisoformat(row['resolved_at']) if row['resolved_at'] else None,
            resolution_note=row['resolution_note'] or '',
            context=json.loads(row['context_json']) if row['context_json'] else {},
            token_estimate=token_estimate,
            actual_tokens_used=row['actual_tokens_used'] if 'actual_tokens_used' in row.keys() else 0
        )


class Sentinel:
    """Log watcher and ticket generator"""

    def __init__(self, db_path: str, patterns: Optional[List[LogPattern]] = None):
        self.db_path = db_path
        self.store = TicketStore(db_path)
        self.patterns = patterns or SentinelPatterns.get_all_patterns()
        self.running = False
        self.watchers: Dict[str, threading.Thread] = {}
        self.callbacks: List[Callable[[Ticket], None]] = []

        # Deduplication window (don't create duplicate tickets within 5 minutes)
        self.recent_tickets: deque = deque(maxlen=100)
        self.dedup_window = timedelta(minutes=5)

    def add_callback(self, callback: Callable[[Ticket], None]):
        """Add callback for new tickets"""
        self.callbacks.append(callback)

    def _notify_callbacks(self, ticket: Ticket):
        """Notify all callbacks of new ticket"""
        for callback in self.callbacks:
            try:
                callback(ticket)
            except Exception as e:
                print(f"Callback error: {e}")

    def _should_create_ticket(self, pattern: LogPattern, message: str) -> bool:
        """Check if we should create a ticket (deduplication)"""
        ticket_hash = hashlib.md5(f"{pattern.name}:{message[:100]}".encode()).hexdigest()
        now = datetime.now()

        # Check recent tickets
        for recent_hash, timestamp in self.recent_tickets:
            if recent_hash == ticket_hash and now - timestamp < self.dedup_window:
                return False

        self.recent_tickets.append((ticket_hash, now))
        return True

    def process_line(self, line: str, source: str, line_number: int = 0) -> Optional[Ticket]:
        """Process a log line and create ticket if pattern matches"""
        for pattern in self.patterns:
            match = pattern.regex.search(line)
            if match:
                if not self._should_create_ticket(pattern, line):
                    return None

                # Create ticket
                ticket = Ticket(
                    title=f"{pattern.name}: {pattern.description}",
                    description=f"Pattern '{pattern.name}' matched in {source}\n\nMatched: {match.group(0)}",
                    severity=pattern.severity,
                    category=pattern.category,
                    source_file=source,
                    source_line=line_number,
                    pattern_matched=pattern.name,
                    log_entries=[line],
                    assigned_agent=pattern.agent,
                    context={
                        'match': match.group(0),
                        'pattern': pattern.pattern,
                        'groups': match.groups() if match.groups() else []
                    }
                )

                # Estimate tokens for this ticket
                ticket.estimate_tokens()

                # Store ticket
                ticket.id = self.store.create_ticket(ticket)

                # Log entry for training
                self.store.log_entry(
                    source=source,
                    level=pattern.severity.value,
                    message=line[:500],
                    raw_line=line,
                    line_number=line_number
                )

                # Notify callbacks
                self._notify_callbacks(ticket)

                print(f"[SENTINEL] Ticket #{ticket.id}: {ticket.title} ({ticket.severity.value})")
                return ticket

        return None

    def watch_file(self, file_path: str, follow: bool = True):
        """Watch a log file for patterns"""
        path = Path(file_path)
        source = str(path)

        if not path.exists():
            print(f"[SENTINEL] Warning: {file_path} does not exist, waiting...")

        line_number = 0

        while self.running:
            try:
                if not path.exists():
                    time.sleep(1)
                    continue

                with open(path, 'r') as f:
                    # Seek to end if following
                    if follow and line_number == 0:
                        f.seek(0, 2)  # Seek to end

                    while self.running:
                        line = f.readline()
                        if line:
                            line_number += 1
                            self.process_line(line.strip(), source, line_number)
                        else:
                            if not follow:
                                break
                            time.sleep(0.1)

                if not follow:
                    break

            except Exception as e:
                print(f"[SENTINEL] Error watching {file_path}: {e}")
                time.sleep(1)

    def watch_stream(self, stream, source: str = "stream"):
        """Watch a stream (stdout/stderr) for patterns"""
        line_number = 0

        while self.running:
            try:
                line = stream.readline()
                if line:
                    line_number += 1
                    self.process_line(line.strip(), source, line_number)
                else:
                    time.sleep(0.1)
            except Exception as e:
                print(f"[SENTINEL] Error reading stream: {e}")
                break

    def start_watcher(self, name: str, file_path: str, follow: bool = True):
        """Start a file watcher in a separate thread"""
        if name in self.watchers:
            print(f"[SENTINEL] Watcher '{name}' already running")
            return

        self.running = True
        thread = threading.Thread(
            target=self.watch_file,
            args=(file_path, follow),
            daemon=True,
            name=f"sentinel-{name}"
        )
        thread.start()
        self.watchers[name] = thread
        print(f"[SENTINEL] Started watcher '{name}' on {file_path}")

    def stop_watcher(self, name: str):
        """Stop a specific watcher"""
        if name in self.watchers:
            # Thread will stop on next iteration
            del self.watchers[name]
            print(f"[SENTINEL] Stopped watcher '{name}'")

    def stop_all(self):
        """Stop all watchers"""
        self.running = False
        self.watchers.clear()
        print("[SENTINEL] Stopped all watchers")

    def get_stats(self) -> Dict:
        """Get sentinel statistics"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Ticket counts by status
        cursor.execute('''
            SELECT status, COUNT(*) FROM tickets
            GROUP BY status
        ''')
        status_counts = dict(cursor.fetchall())

        # Ticket counts by severity
        cursor.execute('''
            SELECT severity, COUNT(*) FROM tickets
            GROUP BY severity
        ''')
        severity_counts = dict(cursor.fetchall())

        # Recent tickets (last 24h)
        cursor.execute('''
            SELECT COUNT(*) FROM tickets
            WHERE created_at > datetime('now', '-24 hours')
        ''')
        recent_count = cursor.fetchone()[0]

        # Log entries
        cursor.execute('SELECT COUNT(*) FROM sentinel_log_stream')
        log_count = cursor.fetchone()[0]

        conn.close()

        return {
            'watchers': list(self.watchers.keys()),
            'patterns_count': len(self.patterns),
            'status_counts': status_counts,
            'severity_counts': severity_counts,
            'recent_tickets_24h': recent_count,
            'total_log_entries': log_count
        }


def main():
    parser = argparse.ArgumentParser(description='IMMUNOS Sentinel - Log Watcher')
    subparsers = parser.add_subparsers(dest='command', help='Commands')

    # Watch command
    watch_parser = subparsers.add_parser('watch', help='Watch a log file')
    watch_parser.add_argument('file', help='Log file to watch')
    watch_parser.add_argument('--no-follow', action='store_true', help='Process once and exit')

    # Watch dashboard command
    dashboard_parser = subparsers.add_parser('watch-dashboard', help='Watch IMMUNOS dashboard logs')

    # Stats command
    stats_parser = subparsers.add_parser('stats', help='Show sentinel statistics')

    # Tickets command
    tickets_parser = subparsers.add_parser('tickets', help='List tickets')
    tickets_parser.add_argument('--status', help='Filter by status')
    tickets_parser.add_argument('--severity', help='Filter by severity')
    tickets_parser.add_argument('--limit', type=int, default=20, help='Max tickets to show')

    args = parser.parse_args()

    # Database path
    db_path = str(Path(__file__).parent.parent / '.immunos' / 'db' / 'dashboard.db')

    if args.command == 'watch':
        sentinel = Sentinel(db_path)
        sentinel.running = True
        print(f"[SENTINEL] Watching {args.file}...")
        try:
            sentinel.watch_file(args.file, follow=not args.no_follow)
        except KeyboardInterrupt:
            print("\n[SENTINEL] Stopped")

    elif args.command == 'watch-dashboard':
        sentinel = Sentinel(db_path)
        sentinel.running = True

        # Watch common dashboard log locations
        log_dir = Path(__file__).parent / 'logs'
        log_dir.mkdir(exist_ok=True)

        print("[SENTINEL] Watching dashboard logs...")
        print("  Tip: Redirect dashboard output to a file to watch it:")
        print("  python immunos_dashboard.py 2>&1 | tee logs/dashboard.log")

        try:
            sentinel.start_watcher('dashboard', str(log_dir / 'dashboard.log'))
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            sentinel.stop_all()
            print("\n[SENTINEL] Stopped")

    elif args.command == 'stats':
        sentinel = Sentinel(db_path)
        stats = sentinel.get_stats()
        print("\n=== SENTINEL STATISTICS ===")
        print(f"Patterns loaded: {stats['patterns_count']}")
        print(f"Active watchers: {stats['watchers']}")
        print(f"Tickets (24h): {stats['recent_tickets_24h']}")
        print(f"Log entries: {stats['total_log_entries']}")
        print("\nBy Status:", stats['status_counts'])
        print("By Severity:", stats['severity_counts'])

    elif args.command == 'tickets':
        sentinel = Sentinel(db_path)
        tickets = sentinel.store.get_tickets(
            status=args.status,
            severity=args.severity,
            limit=args.limit
        )

        print(f"\n=== TICKETS ({len(tickets)}) ===\n")
        for ticket in tickets:
            status_icon = {
                'open': '🔴',
                'in_progress': '🟡',
                'resolved': '🟢',
                'closed': '⚪',
                'ignored': '⚫'
            }.get(ticket.status.value, '⚪')

            print(f"{status_icon} #{ticket.id} [{ticket.severity.value.upper()}] {ticket.title}")
            print(f"   Category: {ticket.category.value} | Agent: {ticket.assigned_agent}")
            print(f"   Created: {ticket.created_at.strftime('%Y-%m-%d %H:%M')}")
            print()

    else:
        parser.print_help()


if __name__ == '__main__':
    main()
