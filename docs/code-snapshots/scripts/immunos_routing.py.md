---
source: /Users/byron/projects/scripts/immunos_routing.py
relative: scripts/immunos_routing.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
IMMUNOS Model Router
====================

Intelligent model routing based on task type, token budget, and performance.

Features:
- Token budget awareness (150k warning, 180k auto-switch)
- Task-based routing (code → qwen2.5-coder, reasoning → deepseek-r1)
- Classification-based routing (critical → Claude, routine → Ollama)
- Configurable thresholds
- Routing decision logging
"""

import os
import sys
import sqlite3
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass

# Import token tracker
try:
    from immunos_token_tracker import TokenTracker
except ImportError:
    # Try relative import
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from immunos_token_tracker import TokenTracker


@dataclass
class RoutingDecision:
    """Result of a routing decision"""
    model_name: str
    routing_reason: str
    session_tokens: int
    warning_level: str  # 'normal', 'warning', 'critical'


class ModelRouter:
    """
    Intelligent model routing with token budget awareness

    Routing Logic Priority:
    1. Token budget exceeded (>= 180k) → Always Ollama
    2. Token budget warning (>= 150k) + routine task → Ollama
    3. Critical task → Always Claude
    4. Task-optimized routing → Specific Ollama model
    5. Default → Claude for complex tasks
    """

    # Token thresholds
    WARNING_THRESHOLD = 150_000
    AUTO_SWITCH_THRESHOLD = 180_000

    # Task type to model mappings
    TASK_MAPPINGS = {
        # Code tasks -> qwen2.5-coder
        'code_scan': 'qwen2.5-coder:7b',
        'vulnerability_detection': 'qwen2.5-coder:7b',
        'code_review': 'qwen2.5-coder:7b',
        'code_generation': 'qwen2.5-coder:7b',
        'debug': 'qwen2.5-coder:7b',

        # Reasoning tasks -> deepseek-r1
        'reasoning': 'deepseek-r1:14b',
        'planning': 'deepseek-r1:14b',
        'analysis': 'deepseek-r1:14b',
        'problem_solving': 'deepseek-r1:14b',

        # Simple tasks -> lightweight model
        'file_analysis': 'qwen2.5:1.5b',
        'summarize': 'qwen2.5:1.5b',
        'simple_question': 'qwen2.5:1.5b',
    }

    # Default models
    CLAUDE_MODEL = 'claude-sonnet-4.5'
    DEFAULT_OLLAMA_MODEL = 'qwen2.5-coder:7b'

    def __init__(self, db_path: str = None):
        self.db_path = db_path or str(Path.home() / ".immunos" / "db" / "immunos.db")
        self.token_tracker = TokenTracker(db_path=self.db_path)
        self.session_id = self.token_tracker.get_or_create_session()

        # Load custom configuration if exists
        self._load_config()

    # ========================================================================
    # CONFIGURATION
    # ========================================================================

    def _load_config(self):
        """Load routing configuration from database settings table"""
        try:
            with sqlite3.connect(self.db_path) as conn:
                # Create settings table if doesn't exist
                conn.execute("""
                    CREATE TABLE IF NOT EXISTS settings (
                        key TEXT PRIMARY KEY,
                        value TEXT NOT NULL,
                        updated_at TEXT NOT NULL
                    )
                """)
                conn.commit()

                # Load warning threshold
                cursor = conn.execute("""
                    SELECT value FROM settings WHERE key = 'token_warning_threshold'
                """)
                row = cursor.fetchone()
                if row:
                    self.WARNING_THRESHOLD = int(row[0])

                # Load auto-switch threshold
                cursor = conn.execute("""
                    SELECT value FROM settings WHERE key = 'token_auto_switch_threshold'
                """)
                row = cursor.fetchone()
                if row:
                    self.AUTO_SWITCH_THRESHOLD = int(row[0])

        except Exception as e:
            print(f"[WARNING] Failed to load routing config: {e}")

    def get_routing_config(self) -> Dict[str, Any]:
        """
        Get current routing configuration

        Returns:
            Dict with warning_threshold, auto_switch_threshold, task_mappings
        """
        return {
            'warning_threshold': self.WARNING_THRESHOLD,
            'auto_switch_threshold': self.AUTO_SWITCH_THRESHOLD,
            'task_mappings': self.TASK_MAPPINGS.copy(),
            'claude_model': self.CLAUDE_MODEL,
            'default_ollama_model': self.DEFAULT_OLLAMA_MODEL
        }

    def update_routing_config(self, config: Dict[str, Any]) -> bool:
        """
        Update routing configuration

        Args:
            config: Dict with warning_threshold, auto_switch_threshold, task_mappings

        Returns:
            True if successful
        """
        try:
            with sqlite3.connect(self.db_path) as conn:
                from datetime import datetime

                # Update warning threshold
                if 'warning_threshold' in config:
                    conn.execute("""
                        INSERT OR REPLACE INTO settings (key, value, updated_at)
                        VALUES ('token_warning_threshold', ?, ?)
                    """, (str(config['warning_threshold']), datetime.now().isoformat()))

                    self.WARNING_THRESHOLD = config['warning_threshold']

                # Update auto-switch threshold
                if 'auto_switch_threshold' in config:
                    conn.execute("""
                        INSERT OR REPLACE INTO settings (key, value, updated_at)
                        VALUES ('token_auto_switch_threshold', ?, ?)
                    """, (str(config['auto_switch_threshold']), datetime.now().isoformat()))

                    self.AUTO_SWITCH_THRESHOLD = config['auto_switch_threshold']

                # Update task mappings
                if 'task_mappings' in config:
                    import json
                    conn.execute("""
                        INSERT OR REPLACE INTO settings (key, value, updated_at)
                        VALUES ('task_mappings', ?, ?)
                    """, (json.dumps(config['task_mappings']), datetime.now().isoformat()))

                    self.TASK_MAPPINGS.update(config['task_mappings'])

                conn.commit()
                return True

        except Exception as e:
            print(f"[ERROR] Failed to update routing config: {e}")
            return False

    # ========================================================================
    # ROUTING LOGIC
    # ========================================================================

    def route_task(self,
                   task_type: str,
                   task_classification: str,
                   estimated_tokens: int = None) -> Tuple[str, str]:
        """
        Determine which model to use for a task

        Args:
            task_type: Type of task (e.g., 'code_scan', 'reasoning', 'creative')
            task_classification: Classification ('routine', 'creative', 'critical')
            estimated_tokens: Estimated token count for this task (optional)

        Returns:
            Tuple of (model_name, routing_reason)
        """
        # Get current session status
        session = self.token_tracker.get_session_usage(self.session_id)
        session_tokens = session['total_tokens']

        # PRIORITY 1: Token budget exceeded (>= 180k)
        if session_tokens >= self.AUTO_SWITCH_THRESHOLD:
            model = self._select_ollama_model(task_type)

            self.token_tracker.log_routing_decision(
                task_type=task_type,
                task_classification=task_classification,
                model_selected=model,
                routing_reason='token_budget_exceeded',
                token_count=estimated_tokens
            )

            return (model, 'token_budget_exceeded')

        # PRIORITY 2: Token budget warning (>= 150k) + routine task
        if session_tokens >= self.WARNING_THRESHOLD and task_classification == 'routine':
            model = self._select_ollama_model(task_type)

            self.token_tracker.log_routing_decision(
                task_type=task_type,
                task_classification=task_classification,
                model_selected=model,
                routing_reason='token_budget_warning',
                token_count=estimated_tokens
            )

            return (model, 'token_budget_warning')

        # PRIORITY 3: Critical task → Always Claude
        if task_classification == 'critical':
            self.token_tracker.log_routing_decision(
                task_type=task_type,
                task_classification=task_classification,
                model_selected=self.CLAUDE_MODEL,
                routing_reason='task_critical',
                token_count=estimated_tokens
            )

            return (self.CLAUDE_MODEL, 'task_critical')

        # PRIORITY 4: Task-optimized routing
        if task_type in self.TASK_MAPPINGS:
            model = self.TASK_MAPPINGS[task_type]

            self.token_tracker.log_routing_decision(
                task_type=task_type,
                task_classification=task_classification,
                model_selected=model,
                routing_reason='task_optimized',
                token_count=estimated_tokens
            )

            return (model, 'task_optimized')

        # PRIORITY 5: Default - Claude for creative/complex tasks
        self.token_tracker.log_routing_decision(
            task_type=task_type,
            task_classification=task_classification,
            model_selected=self.CLAUDE_MODEL,
            routing_reason='default_complex',
            token_count=estimated_tokens
        )

        return (self.CLAUDE_MODEL, 'default_complex')

    def route_task_detailed(self,
                          task_type: str,
                          task_classification: str,
                          estimated_tokens: int = None) -> RoutingDecision:
        """
        Get detailed routing decision with warning level

        Args:
            task_type: Type of task
            task_classification: Classification
            estimated_tokens: Estimated token count

        Returns:
            RoutingDecision object with model, reason, session_tokens, warning_level
        """
        model, reason = self.route_task(task_type, task_classification, estimated_tokens)

        session = self.token_tracker.get_session_usage(self.session_id)
        session_tokens = session['total_tokens']

        # Determine warning level
        if session_tokens >= self.AUTO_SWITCH_THRESHOLD:
            warning_level = 'critical'
        elif session_tokens >= self.WARNING_THRESHOLD:
            warning_level = 'warning'
        else:
            warning_level = 'normal'

        return RoutingDecision(
            model_name=model,
            routing_reason=reason,
            session_tokens=session_tokens,
            warning_level=warning_level
        )

    def _select_ollama_model(self, task_type: str) -> str:
        """
        Select appropriate Ollama model for task type

        Args:
            task_type: Type of task

        Returns:
            Ollama model name
        """
        # Check if task type has specific mapping
        if task_type in self.TASK_MAPPINGS:
            return self.TASK_MAPPINGS[task_type]

        # Default to general code model
        return self.DEFAULT_OLLAMA_MODEL

    # ========================================================================
    # SESSION MANAGEMENT
    # ========================================================================

    def get_session_status(self) -> Dict[str, Any]:
        """
        Get current session status with routing recommendations

        Returns:
            Dict with session info, warning level, recommended model
        """
        session = self.token_tracker.get_session_usage(self.session_id)

        # Determine warning level
        if session['total_tokens'] >= self.AUTO_SWITCH_THRESHOLD:
            warning_level = 'critical'
            recommended_model = self.DEFAULT_OLLAMA_MODEL
            message = 'Token budget exceeded. All tasks routed to Ollama.'
        elif session['total_tokens'] >= self.WARNING_THRESHOLD:
            warning_level = 'warning'
            recommended_model = 'mixed'
            message = 'Approaching token limit. Routine tasks routed to Ollama.'
        else:
            warning_level = 'normal'
            recommended_model = self.CLAUDE_MODEL
            message = 'Token budget normal. Using optimal routing.'

        return {
            'session_id': self.session_id,
            'total_tokens': session['total_tokens'],
            'warning_threshold': self.WARNING_THRESHOLD,
            'auto_switch_threshold': self.AUTO_SWITCH_THRESHOLD,
            'warning_level': warning_level,
            'recommended_model': recommended_model,
            'message': message,
            'percent_used': (session['total_tokens'] / 200000) * 100  # Assuming 200k limit
        }

    def reset_session(self) -> str:
        """
        Create a new session (effectively resets token budget)

        Returns:
            New session ID
        """
        self.session_id = self.token_tracker.create_session()
        return self.session_id


# ============================================================================
# CLI INTERFACE
# ============================================================================

def main():
    """Command-line interface for model router"""
    import argparse

    parser = argparse.ArgumentParser(description='IMMUNOS Model Router')
    parser.add_argument('action', choices=['route', 'status', 'config', 'reset'],
                       help='Action to perform')
    parser.add_argument('--task-type', help='Task type (for route action)')
    parser.add_argument('--classification', choices=['routine', 'creative', 'critical'],
                       help='Task classification (for route action)')
    parser.add_argument('--tokens', type=int, help='Estimated token count')
    parser.add_argument('--db', help='Database path')

    args = parser.parse_args()

    router = ModelRouter(db_path=args.db)

    if args.action == 'route':
        if not args.task_type or not args.classification:
            print("[ERROR] --task-type and --classification required for route action")
            sys.exit(1)

        decision = router.route_task_detailed(
            task_type=args.task_type,
            task_classification=args.classification,
            estimated_tokens=args.tokens
        )

        print(f"\nRouting Decision:")
        print(f"Model: {decision.model_name}")
        print(f"Reason: {decision.routing_reason}")
        print(f"Session Tokens: {decision.session_tokens:,}")
        print(f"Warning Level: {decision.warning_level}")

    elif args.action == 'status':
        status = router.get_session_status()

        print(f"\nSession Status:")
        print(f"Session ID: {status['session_id']}")
        print(f"Total Tokens: {status['total_tokens']:,} ({status['percent_used']:.1f}%)")
        print(f"Warning Level: {status['warning_level']}")
        print(f"Recommended Model: {status['recommended_model']}")
        print(f"\n{status['message']}")

        print(f"\nThresholds:")
        print(f"  Warning: {status['warning_threshold']:,}")
        print(f"  Auto-Switch: {status['auto_switch_threshold']:,}")

    elif args.action == 'config':
        config = router.get_routing_config()

        print(f"\nRouting Configuration:")
        print(f"\nThresholds:")
        print(f"  Warning: {config['warning_threshold']:,}")
        print(f"  Auto-Switch: {config['auto_switch_threshold']:,}")

        print(f"\nModels:")
        print(f"  Claude: {config['claude_model']}")
        print(f"  Default Ollama: {config['default_ollama_model']}")

        print(f"\nTask Mappings:")
        for task, model in config['task_mappings'].items():
            print(f"  {task}: {model}")

    elif args.action == 'reset':
        session_id = router.reset_session()
        print(f"\nSession reset. New session ID: {session_id}")


if __name__ == '__main__':
    main()

```
