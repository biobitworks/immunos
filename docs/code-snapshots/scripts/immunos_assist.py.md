---
source: /Users/byron/projects/scripts/immunos_assist.py
relative: scripts/immunos_assist.py
generated_at: 2025-12-23 10:28
---

````python
#!/usr/bin/env python3
"""
IMMUNOS Assistant - Route tasks to appropriate models
Saves Claude tokens by using local Ollama models for routine work
"""
import sys
import requests
import json
from pathlib import Path
from datetime import datetime

# Add scripts to path
sys.path.append(str(Path(__file__).parent))
from immunos_routing import ModelRouter
from immunos_token_tracker import TokenTracker

class ImmunosAssist:
    """Helper to route work to appropriate models"""

    def __init__(self):
        self.router = ModelRouter()
        self.tracker = TokenTracker()
        self.ollama_url = "http://localhost:11434"

    def estimate_tokens(self, text: str) -> int:
        """Rough token estimate (1 token ≈ 4 chars)"""
        return len(text) // 4

    def ask_ollama(self, model: str, prompt: str, system: str = None) -> str:
        """Ask Ollama model a question"""
        messages = []
        if system:
            messages.append({"role": "system", "content": system})
        messages.append({"role": "user", "content": prompt})

        try:
            response = requests.post(
                f"{self.ollama_url}/api/chat",
                json={
                    "model": model,
                    "messages": messages,
                    "stream": False
                },
                timeout=60
            )

            if response.status_code == 200:
                result = response.json()
                return result['message']['content']
            else:
                return f"ERROR: Ollama returned {response.status_code}"

        except Exception as e:
            return f"ERROR: {str(e)}"

    def route_and_execute(self, task_type: str, task_classification: str,
                         prompt: str, system: str = None) -> dict:
        """Route task and execute on appropriate model"""

        # Get routing decision
        model, reason = self.router.route_task(task_type, task_classification)

        start_time = datetime.now()

        if model.startswith('claude'):
            # Return to Claude - can't execute Claude from here
            return {
                'model': model,
                'reason': reason,
                'result': None,
                'message': 'Route to Claude - creative/complex task'
            }
        else:
            # Execute on Ollama
            result = self.ask_ollama(model, prompt, system)
            duration_ms = (datetime.now() - start_time).total_seconds() * 1000

            # Track usage
            prompt_tokens = self.estimate_tokens(prompt + (system or ''))
            response_tokens = self.estimate_tokens(result)

            self.tracker.record_usage(
                model_name=model,
                component='immunos_assist',
                prompt_tokens=prompt_tokens,
                response_tokens=response_tokens,
                response_time_ms=int(duration_ms),
                provider='ollama',
                routing_reason=reason
            )

            return {
                'model': model,
                'reason': reason,
                'result': result,
                'duration_ms': duration_ms,
                'tokens': prompt_tokens + response_tokens
            }


# Common helper functions for Claude Code to use

def code_review(code: str, language: str = "python") -> str:
    """Review code using qwen2.5-coder"""
    assist = ImmunosAssist()
    response = assist.route_and_execute(
        task_type='code_review',
        task_classification='routine',
        prompt=f"Review this {language} code for bugs, style issues, and improvements:\n\n```{language}\n{code}\n```",
        system="You are an expert code reviewer. Provide concise, actionable feedback."
    )
    return response['result'] if response['result'] else "Route to Claude"


def analyze_file(filepath: str, question: str = "Summarize this file") -> str:
    """Analyze a file using appropriate model"""
    assist = ImmunosAssist()

    try:
        with open(filepath, 'r') as f:
            content = f.read()
    except Exception as e:
        return f"ERROR reading file: {e}"

    response = assist.route_and_execute(
        task_type='file_analysis',
        task_classification='routine',
        prompt=f"{question}\n\nFile: {filepath}\n\n```\n{content[:4000]}\n```",
        system="You are a helpful file analyzer."
    )
    return response['result'] if response['result'] else "Route to Claude"


def simple_code_task(description: str, language: str = "python") -> str:
    """Generate simple code using qwen2.5-coder"""
    assist = ImmunosAssist()
    response = assist.route_and_execute(
        task_type='code_generation',
        task_classification='routine',
        prompt=f"Write {language} code for: {description}",
        system=f"You are an expert {language} programmer. Write clean, concise code."
    )
    return response['result'] if response['result'] else "Route to Claude"


def reason_about(problem: str) -> str:
    """Use deepseek-r1 for reasoning tasks"""
    assist = ImmunosAssist()
    response = assist.route_and_execute(
        task_type='reasoning',
        task_classification='routine',
        prompt=problem,
        system="You are an expert problem solver. Think step by step."
    )
    return response['result'] if response['result'] else "Route to Claude"


if __name__ == '__main__':
    # CLI interface
    if len(sys.argv) < 3:
        print("Usage: python3 immunos_assist.py <command> <args...>")
        print("\nCommands:")
        print("  review <file>              - Review code")
        print("  analyze <file> [question]  - Analyze file")
        print("  code <description>         - Generate code")
        print("  reason <problem>           - Reason about problem")
        sys.exit(1)

    command = sys.argv[1]

    if command == 'review' and len(sys.argv) >= 3:
        filepath = sys.argv[2]
        with open(filepath, 'r') as f:
            code = f.read()
        print(code_review(code))

    elif command == 'analyze' and len(sys.argv) >= 3:
        filepath = sys.argv[2]
        question = ' '.join(sys.argv[3:]) if len(sys.argv) > 3 else "Summarize this file"
        print(analyze_file(filepath, question))

    elif command == 'code' and len(sys.argv) >= 3:
        description = ' '.join(sys.argv[2:])
        print(simple_code_task(description))

    elif command == 'reason' and len(sys.argv) >= 3:
        problem = ' '.join(sys.argv[2:])
        print(reason_about(problem))

    else:
        print(f"Unknown command: {command}")
        sys.exit(1)

````
