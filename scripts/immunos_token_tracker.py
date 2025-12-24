#!/usr/bin/env python3
"""
IMMUNOS Token Tracker
=====================

Tracks token usage and calculates cost savings between Claude and Ollama.

Features:
- Record token usage per request
- Track Claude session tokens
- Compare Claude vs Ollama usage
- Calculate cost savings
- Per-model breakdowns
- Routing statistics
"""

import os
import sys
import sqlite3
import uuid
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Any, Optional


class TokenTracker:
    """
    Track and compare Claude vs Ollama token usage

    Manages token budgets, calculates costs, and provides analytics
    for the monitoring dashboard.
    """

    # Claude pricing (per 1K tokens)
    # Source: https://www.anthropic.com/pricing
    CLAUDE_COST_INPUT = 0.003  # $0.003 per 1K input tokens
    CLAUDE_COST_OUTPUT = 0.015  # $0.015 per 1K output tokens

    # Token thresholds for session management
    WARNING_THRESHOLD = 150_000  # Warn user
    AUTO_SWITCH_THRESHOLD = 180_000  # Auto-switch to Ollama

    def __init__(self, db_path: str = None):
        self.db_path = db_path or str(Path.home() / ".immunos" / "db" / "immunos.db")
        self._init_db()

    # ========================================================================
    # DATABASE INITIALIZATION
    # ========================================================================

    def _init_db(self):
        """Initialize database tables for token tracking"""
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)

        with sqlite3.connect(self.db_path) as conn:
            # Claude sessions table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS claude_sessions (
                    session_id TEXT PRIMARY KEY,
                    created_at TEXT NOT NULL,
                    total_tokens INTEGER DEFAULT 0,
                    warning_triggered INTEGER DEFAULT 0,
                    auto_switched INTEGER DEFAULT 0,
                    last_updated TEXT
                )
            """)

            # LLM model usage table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS llm_model_usage (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    timestamp TEXT NOT NULL,
                    model_name TEXT NOT NULL,
                    component TEXT NOT NULL,
                    prompt_tokens INTEGER DEFAULT 0,
                    response_tokens INTEGER DEFAULT 0,
                    response_time_ms INTEGER,
                    provider TEXT NOT NULL,
                    routing_reason TEXT,
                    session_id TEXT
                )
            """)

            # Model routing log table
            conn.execute("""
                CREATE TABLE IF NOT EXISTS model_routing_log (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    timestamp TEXT NOT NULL,
                    task_type TEXT NOT NULL,
                    task_classification TEXT NOT NULL,
                    model_selected TEXT NOT NULL,
                    routing_reason TEXT NOT NULL,
                    token_count INTEGER,
                    duration_ms INTEGER,
                    success INTEGER NOT NULL
                )
            """)

            conn.commit()

    # ========================================================================
    # SESSION MANAGEMENT
    # ========================================================================

    def create_session(self) -> str:
        """
        Create a new Claude session

        Returns:
            Session ID
        """
        session_id = str(uuid.uuid4())

        with sqlite3.connect(self.db_path) as conn:
            conn.execute("""
                INSERT INTO claude_sessions (session_id, created_at, last_updated)
                VALUES (?, ?, ?)
            """, (session_id, datetime.now().isoformat(), datetime.now().isoformat()))
            conn.commit()

        return session_id

    def get_or_create_session(self) -> str:
        """
        Get current active session or create new one

        Returns:
            Session ID
        """
        # Get most recent session (within last 24 hours)
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute("""
                SELECT session_id, created_at, total_tokens
                FROM claude_sessions
                WHERE datetime(created_at) > datetime('now', '-24 hours')
                ORDER BY created_at DESC
                LIMIT 1
            """)

            row = cursor.fetchone()

            if row:
                return row[0]

            # Create new session if none exists
            return self.create_session()

    def get_session_usage(self, session_id: str) -> Dict[str, Any]:
        """
        Get token usage for a specific Claude session

        Args:
            session_id: Session ID

        Returns:
            Dict with total_tokens, warning_triggered, auto_switched, created_at
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute("""
                SELECT total_tokens, warning_triggered, auto_switched, created_at, last_updated
                FROM claude_sessions
                WHERE session_id = ?
            """, (session_id,))

            row = cursor.fetchone()

            if row:
                return {
                    'session_id': session_id,
                    'total_tokens': row[0],
                    'warning_triggered': bool(row[1]),
                    'auto_switched': bool(row[2]),
                    'created_at': row[3],
                    'last_updated': row[4]
                }

            # Session not found
            return {
                'session_id': session_id,
                'total_tokens': 0,
                'warning_triggered': False,
                'auto_switched': False,
                'created_at': None,
                'last_updated': None
            }

    def update_session_tokens(self, session_id: str, tokens_used: int):
        """
        Update session token count

        Args:
            session_id: Session ID
            tokens_used: Number of tokens used in this request
        """
        with sqlite3.connect(self.db_path) as conn:
            # Get current total
            cursor = conn.execute("""
                SELECT total_tokens, warning_triggered, auto_switched
                FROM claude_sessions
                WHERE session_id = ?
            """, (session_id,))

            row = cursor.fetchone()

            if not row:
                # Create session if doesn't exist
                conn.execute("""
                    INSERT INTO claude_sessions (session_id, created_at, total_tokens, last_updated)
                    VALUES (?, ?, ?, ?)
                """, (session_id, datetime.now().isoformat(), tokens_used, datetime.now().isoformat()))
                conn.commit()
                return

            current_total, warning_triggered, auto_switched = row
            new_total = current_total + tokens_used

            # Check thresholds
            new_warning = warning_triggered or (new_total >= self.WARNING_THRESHOLD)
            new_auto_switch = auto_switched or (new_total >= self.AUTO_SWITCH_THRESHOLD)

            # Update
            conn.execute("""
                UPDATE claude_sessions
                SET total_tokens = ?,
                    warning_triggered = ?,
                    auto_switched = ?,
                    last_updated = ?
                WHERE session_id = ?
            """, (new_total, int(new_warning), int(new_auto_switch), datetime.now().isoformat(), session_id))

            conn.commit()

    # ========================================================================
    # USAGE RECORDING
    # ========================================================================

    def record_usage(self,
                    model_name: str,
                    component: str,
                    prompt_tokens: int,
                    response_tokens: int,
                    response_time_ms: int,
                    provider: str,
                    routing_reason: str = None,
                    session_id: str = None):
        """
        Record token usage for a request

        Args:
            model_name: Name of the model used
            component: Component that made the request (e.g., 'nk_scanner', 'training')
            prompt_tokens: Number of prompt tokens
            response_tokens: Number of response tokens
            response_time_ms: Response time in milliseconds
            provider: 'claude' or 'ollama'
            routing_reason: Reason for model selection (optional)
            session_id: Claude session ID (optional)
        """
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("""
                INSERT INTO llm_model_usage
                (timestamp, model_name, component, prompt_tokens, response_tokens,
                 response_time_ms, provider, routing_reason, session_id)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                datetime.now().isoformat(),
                model_name,
                component,
                prompt_tokens,
                response_tokens,
                response_time_ms,
                provider,
                routing_reason,
                session_id
            ))
            conn.commit()

        # If Claude, update session token count
        if provider == 'claude' and session_id:
            self.update_session_tokens(session_id, prompt_tokens + response_tokens)

    # ========================================================================
    # ANALYTICS & COMPARISON
    # ========================================================================

    def get_provider_comparison(self, hours: int = 24) -> Dict[str, Any]:
        """
        Compare Claude vs Ollama usage over a time period

        Args:
            hours: Number of hours to look back (default: 24)

        Returns:
            Dict with claude, ollama, savings_usd, savings_percent, total_requests
        """
        cutoff = (datetime.now() - timedelta(hours=hours)).isoformat()

        with sqlite3.connect(self.db_path) as conn:
            # Claude stats
            cursor = conn.execute("""
                SELECT
                    COUNT(*) as requests,
                    SUM(prompt_tokens) as prompt_tokens,
                    SUM(response_tokens) as response_tokens,
                    AVG(response_time_ms) as avg_latency_ms
                FROM llm_model_usage
                WHERE provider = 'claude' AND timestamp >= ?
            """, (cutoff,))

            row = cursor.fetchone()
            claude_requests = row[0] or 0
            claude_prompt = row[1] or 0
            claude_response = row[2] or 0
            claude_latency = row[3] or 0

            claude_tokens = claude_prompt + claude_response
            claude_cost = (
                (claude_prompt / 1000) * self.CLAUDE_COST_INPUT +
                (claude_response / 1000) * self.CLAUDE_COST_OUTPUT
            )

            # Ollama stats
            cursor = conn.execute("""
                SELECT
                    COUNT(*) as requests,
                    SUM(prompt_tokens + response_tokens) as tokens,
                    AVG(response_time_ms) as avg_latency_ms
                FROM llm_model_usage
                WHERE provider = 'ollama' AND timestamp >= ?
            """, (cutoff,))

            row = cursor.fetchone()
            ollama_requests = row[0] or 0
            ollama_tokens = row[1] or 0
            ollama_latency = row[2] or 0

            # Calculate what it would have cost if all Ollama requests used Claude
            # Assume 50/50 split between input and output tokens for estimation
            estimated_ollama_cost = (
                ((ollama_tokens / 2) / 1000) * self.CLAUDE_COST_INPUT +
                ((ollama_tokens / 2) / 1000) * self.CLAUDE_COST_OUTPUT
            )

            savings_usd = estimated_ollama_cost
            total_cost = claude_cost
            savings_percent = (savings_usd / (savings_usd + total_cost) * 100) if (savings_usd + total_cost) > 0 else 0

            return {
                'claude': {
                    'tokens': claude_tokens,
                    'requests': claude_requests,
                    'cost_usd': round(claude_cost, 4),
                    'avg_latency_ms': round(claude_latency, 2)
                },
                'ollama': {
                    'tokens': ollama_tokens,
                    'requests': ollama_requests,
                    'cost_usd': 0.0,  # Ollama is free
                    'avg_latency_ms': round(ollama_latency, 2)
                },
                'savings_usd': round(savings_usd, 4),
                'savings_percent': round(savings_percent, 2),
                'total_requests': claude_requests + ollama_requests
            }

    def get_model_breakdown(self) -> List[Dict[str, Any]]:
        """
        Get per-model token usage breakdown

        Returns:
            List of dicts with model_name, provider, tokens, requests, avg_latency, cost
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute("""
                SELECT
                    model_name,
                    provider,
                    SUM(prompt_tokens) as prompt_tokens,
                    SUM(response_tokens) as response_tokens,
                    COUNT(*) as requests,
                    AVG(response_time_ms) as avg_latency_ms
                FROM llm_model_usage
                GROUP BY model_name, provider
                ORDER BY requests DESC
            """)

            results = []
            for row in cursor.fetchall():
                model_name, provider, prompt, response, requests, latency = row

                total_tokens = (prompt or 0) + (response or 0)

                # Calculate cost
                if provider == 'claude':
                    cost = (
                        ((prompt or 0) / 1000) * self.CLAUDE_COST_INPUT +
                        ((response or 0) / 1000) * self.CLAUDE_COST_OUTPUT
                    )
                else:
                    cost = 0.0

                results.append({
                    'model_name': model_name,
                    'provider': provider,
                    'tokens': total_tokens,
                    'requests': requests,
                    'avg_latency_ms': round(latency or 0, 2),
                    'cost_usd': round(cost, 4)
                })

            return results

    def get_routing_stats(self) -> Dict[str, Any]:
        """
        Get model routing decision statistics

        Returns:
            Dict with total_routes, by_reason, by_task_type, success_rate
        """
        with sqlite3.connect(self.db_path) as conn:
            # Total routes
            cursor = conn.execute("SELECT COUNT(*) FROM model_routing_log")
            total_routes = cursor.fetchone()[0]

            # By routing reason
            cursor = conn.execute("""
                SELECT routing_reason, COUNT(*) as count
                FROM model_routing_log
                GROUP BY routing_reason
                ORDER BY count DESC
            """)
            by_reason = {row[0]: row[1] for row in cursor.fetchall()}

            # By task type
            cursor = conn.execute("""
                SELECT task_type, COUNT(*) as count
                FROM model_routing_log
                GROUP BY task_type
                ORDER BY count DESC
            """)
            by_task_type = {row[0]: row[1] for row in cursor.fetchall()}

            # Success rate
            cursor = conn.execute("""
                SELECT
                    SUM(success) as successes,
                    COUNT(*) as total
                FROM model_routing_log
            """)
            row = cursor.fetchone()
            successes = row[0] or 0
            total = row[1] or 1
            success_rate = (successes / total) * 100

            return {
                'total_routes': total_routes,
                'by_reason': by_reason,
                'by_task_type': by_task_type,
                'success_rate': round(success_rate, 2)
            }

    def log_routing_decision(self,
                            task_type: str,
                            task_classification: str,
                            model_selected: str,
                            routing_reason: str,
                            token_count: int = None,
                            duration_ms: int = None,
                            success: bool = True):
        """
        Log a model routing decision

        Args:
            task_type: Type of task (e.g., 'code_scan', 'reasoning')
            task_classification: Classification (e.g., 'routine', 'creative', 'critical')
            model_selected: Model that was selected
            routing_reason: Reason for selection
            token_count: Estimated token count (optional)
            duration_ms: Duration in milliseconds (optional)
            success: Whether the routing was successful
        """
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("""
                INSERT INTO model_routing_log
                (timestamp, task_type, task_classification, model_selected, routing_reason,
                 token_count, duration_ms, success)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """, (
                datetime.now().isoformat(),
                task_type,
                task_classification,
                model_selected,
                routing_reason,
                token_count,
                duration_ms,
                int(success)
            ))
            conn.commit()


# ============================================================================
# CLI INTERFACE
# ============================================================================

def main():
    """Command-line interface for token tracker"""
    import argparse

    parser = argparse.ArgumentParser(description='IMMUNOS Token Tracker')
    parser.add_argument('action', choices=['session', 'comparison', 'models', 'routing'],
                       help='Action to perform')
    parser.add_argument('--session-id', help='Session ID (for session action)')
    parser.add_argument('--hours', type=int, default=24, help='Hours to look back (for comparison)')
    parser.add_argument('--db', help='Database path')

    args = parser.parse_args()

    tracker = TokenTracker(db_path=args.db)

    if args.action == 'session':
        if args.session_id:
            usage = tracker.get_session_usage(args.session_id)
        else:
            session_id = tracker.get_or_create_session()
            usage = tracker.get_session_usage(session_id)

        print(f"\nSession: {usage['session_id']}")
        print(f"Total Tokens: {usage['total_tokens']:,}")
        print(f"Warning Triggered: {usage['warning_triggered']}")
        print(f"Auto-Switched: {usage['auto_switched']}")
        print(f"Created: {usage['created_at']}")
        print(f"Last Updated: {usage['last_updated']}")

    elif args.action == 'comparison':
        data = tracker.get_provider_comparison(hours=args.hours)

        print(f"\nClaude vs Ollama Comparison (last {args.hours} hours):\n")

        print("Claude:")
        print(f"  Tokens: {data['claude']['tokens']:,}")
        print(f"  Requests: {data['claude']['requests']}")
        print(f"  Cost: ${data['claude']['cost_usd']:.4f}")
        print(f"  Avg Latency: {data['claude']['avg_latency_ms']:.0f} ms")

        print("\nOllama:")
        print(f"  Tokens: {data['ollama']['tokens']:,}")
        print(f"  Requests: {data['ollama']['requests']}")
        print(f"  Cost: $0.00 (free)")
        print(f"  Avg Latency: {data['ollama']['avg_latency_ms']:.0f} ms")

        print(f"\nSavings: ${data['savings_usd']:.4f} ({data['savings_percent']:.1f}%)")
        print(f"Total Requests: {data['total_requests']}")

    elif args.action == 'models':
        models = tracker.get_model_breakdown()

        print(f"\nPer-Model Breakdown ({len(models)} models):\n")

        for model in models:
            print(f"{model['model_name']} ({model['provider']}):")
            print(f"  Tokens: {model['tokens']:,}")
            print(f"  Requests: {model['requests']}")
            print(f"  Cost: ${model['cost_usd']:.4f}")
            print(f"  Avg Latency: {model['avg_latency_ms']:.0f} ms")
            print()

    elif args.action == 'routing':
        stats = tracker.get_routing_stats()

        print(f"\nRouting Statistics:\n")
        print(f"Total Routes: {stats['total_routes']}")
        print(f"Success Rate: {stats['success_rate']:.1f}%")

        print("\nBy Routing Reason:")
        for reason, count in stats['by_reason'].items():
            print(f"  {reason}: {count}")

        print("\nBy Task Type:")
        for task_type, count in stats['by_task_type'].items():
            print(f"  {task_type}: {count}")


if __name__ == '__main__':
    main()
