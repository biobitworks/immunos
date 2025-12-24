---
source: /Users/byron/projects/scripts/immunos_handoff.py
relative: scripts/immunos_handoff.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
IMMUNOS Context Handoff System
===============================

Automatically saves conversation state when Claude approaches token limits.
Allows seamless continuation with Ollama models and resumption when Claude returns.

Purpose: IMMUNOS acts as negative selection to Claude - validating outputs,
reducing hallucinations, maintaining scientific rigor (RAIT principles).
"""
import json
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Optional

# Add scripts to path
sys.path.append(str(Path(__file__).parent))
from immunos_token_tracker import TokenTracker


class ContextHandoff:
    """Manages context handoff between Claude and IMMUNOS (Ollama models)"""

    def __init__(self, handoff_dir: str = None):
        self.handoff_dir = Path(handoff_dir or Path.home() / ".immunos" / "handoffs")
        self.handoff_dir.mkdir(parents=True, exist_ok=True)
        self.tracker = TokenTracker()
        self.current_session_id = self.tracker.get_or_create_session()

    def save_handoff(self,
                    conversation_history: List[Dict[str, str]],
                    current_task: str,
                    files_being_worked_on: List[str],
                    next_steps: List[str],
                    context: Dict[str, Any] = None,
                    reason: str = "token_limit_approaching") -> str:
        """
        Save current state for handoff to IMMUNOS

        Args:
            conversation_history: List of {role, content} messages
            current_task: Description of what's being worked on
            files_being_worked_on: Paths to files being edited
            next_steps: What should happen next
            context: Additional context (variables, state, etc.)
            reason: Why handoff is happening

        Returns:
            Path to handoff file
        """
        timestamp = datetime.now().isoformat()
        handoff_id = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Get current token usage
        session = self.tracker.get_session_usage(self.current_session_id)

        handoff_data = {
            "handoff_id": handoff_id,
            "timestamp": timestamp,
            "reason": reason,
            "session_id": self.current_session_id,
            "token_usage": {
                "total_tokens": session.get('total_tokens', 0),
                "warning_threshold": 150000,
                "auto_switch_threshold": 180000,
                "percentage_used": (session.get('total_tokens', 0) / 200000) * 100
            },
            "conversation_history": conversation_history[-50:],  # Last 50 messages
            "current_task": current_task,
            "files_being_worked_on": files_being_worked_on,
            "next_steps": next_steps,
            "context": context or {},
            "resumption_instructions": self._generate_resumption_instructions(
                current_task, next_steps
            )
        }

        # Save to file
        handoff_file = self.handoff_dir / f"handoff_{handoff_id}.json"
        with open(handoff_file, 'w') as f:
            json.dump(handoff_data, f, indent=2)

        # Save latest link
        latest_link = self.handoff_dir / "latest_handoff.json"
        latest_link.unlink(missing_ok=True)
        latest_link.symlink_to(handoff_file.name)

        return str(handoff_file)

    def load_latest_handoff(self) -> Optional[Dict[str, Any]]:
        """Load the most recent handoff"""
        latest_link = self.handoff_dir / "latest_handoff.json"

        if not latest_link.exists():
            return None

        with open(latest_link, 'r') as f:
            return json.load(f)

    def load_handoff(self, handoff_id: str) -> Optional[Dict[str, Any]]:
        """Load a specific handoff by ID"""
        handoff_file = self.handoff_dir / f"handoff_{handoff_id}.json"

        if not handoff_file.exists():
            return None

        with open(handoff_file, 'r') as f:
            return json.load(f)

    def list_handoffs(self) -> List[Dict[str, Any]]:
        """List all available handoffs"""
        handoffs = []
        for file in sorted(self.handoff_dir.glob("handoff_*.json"), reverse=True):
            with open(file, 'r') as f:
                data = json.load(f)
                handoffs.append({
                    'handoff_id': data['handoff_id'],
                    'timestamp': data['timestamp'],
                    'current_task': data['current_task'],
                    'reason': data['reason'],
                    'token_usage': data['token_usage']['total_tokens']
                })
        return handoffs

    def _generate_resumption_instructions(self, current_task: str,
                                         next_steps: List[str]) -> str:
        """Generate instructions for resuming the task"""
        return f"""
# RESUMPTION INSTRUCTIONS FOR CLAUDE

## Current Task
{current_task}

## Next Steps
{chr(10).join(f'{i+1}. {step}' for i, step in enumerate(next_steps))}

## Context
This handoff was created because Claude approached token limits. The task should
be continued by IMMUNOS (Ollama models) and can be resumed by Claude when available.

IMMUNOS operates as a negative selection system to Claude, providing:
- Validation of outputs
- Hallucination detection
- Scientific rigor (RAIT: Rigorous, Accurate, Interpretable, Transparent)
- Multi-modal analysis (text, images, future: video)

When resuming, verify IMMUNOS's work and continue from the last completed step.
"""

    def should_handoff(self, warning_only: bool = False) -> tuple[bool, str]:
        """
        Check if handoff should occur based on token usage

        Args:
            warning_only: If True, only return warning without forcing handoff

        Returns:
            (should_handoff, reason)
        """
        session = self.tracker.get_session_usage(self.current_session_id)
        tokens = session.get('total_tokens', 0)

        if tokens >= 180000:
            return (True, "auto_switch_threshold_reached")
        elif tokens >= 150000:
            if warning_only:
                return (False, "warning_threshold_reached")
            else:
                return (True, "approaching_limit")
        else:
            return (False, "normal")


def auto_handoff_decorator(func):
    """
    Decorator to automatically trigger handoff when approaching limits
    """
    def wrapper(*args, **kwargs):
        handoff = ContextHandoff()
        should_handoff, reason = handoff.should_handoff()

        if should_handoff:
            print(f"\n⚠️  Token limit approaching ({reason})")
            print("Consider creating handoff before continuing...")

        return func(*args, **kwargs)

    return wrapper


if __name__ == '__main__':
    # CLI interface
    import argparse

    parser = argparse.ArgumentParser(description='IMMUNOS Context Handoff Manager')
    parser.add_argument('command', choices=['save', 'load', 'list', 'check'],
                       help='Command to execute')
    parser.add_argument('--task', help='Current task description')
    parser.add_argument('--files', nargs='+', help='Files being worked on')
    parser.add_argument('--steps', nargs='+', help='Next steps')
    parser.add_argument('--handoff-id', help='Handoff ID to load')

    args = parser.parse_args()
    handoff = ContextHandoff()

    if args.command == 'save':
        if not args.task:
            print("ERROR: --task required for save")
            sys.exit(1)

        filepath = handoff.save_handoff(
            conversation_history=[],
            current_task=args.task,
            files_being_worked_on=args.files or [],
            next_steps=args.steps or [],
            reason="manual_save"
        )
        print(f"✅ Handoff saved: {filepath}")

    elif args.command == 'load':
        if args.handoff_id:
            data = handoff.load_handoff(args.handoff_id)
        else:
            data = handoff.load_latest_handoff()

        if not data:
            print("No handoff found")
            sys.exit(1)

        print(json.dumps(data, indent=2))

    elif args.command == 'list':
        handoffs = handoff.list_handoffs()
        print(f"\n{'='*70}")
        print("Available Handoffs")
        print(f"{'='*70}\n")

        for h in handoffs:
            print(f"ID: {h['handoff_id']}")
            print(f"  Time: {h['timestamp']}")
            print(f"  Task: {h['current_task']}")
            print(f"  Reason: {h['reason']}")
            print(f"  Tokens: {h['token_usage']:,}")
            print()

    elif args.command == 'check':
        should_handoff, reason = handoff.should_handoff()
        session = handoff.tracker.get_session_usage(handoff.current_session_id)
        tokens = session.get('total_tokens', 0)

        print(f"\n{'='*70}")
        print("IMMUNOS Handoff Status")
        print(f"{'='*70}\n")
        print(f"Session ID: {handoff.current_session_id}")
        print(f"Token Usage: {tokens:,} / 200,000 ({(tokens/200000)*100:.1f}%)")
        print(f"Should Handoff: {'YES' if should_handoff else 'NO'}")
        print(f"Reason: {reason}")
        print()

```
