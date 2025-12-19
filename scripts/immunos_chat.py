#!/usr/bin/env python3
"""
IMMUNOS Interactive Chat
=========================

Chat interface for IMMUNOS with automatic model routing and validation.

IMMUNOS Purpose:
- Negative selection to Claude (validates outputs, catches hallucinations)
- Scientific rigor (RAIT: Rigorous, Accurate, Interpretable, Transparent)
- Multi-modal analysis (text, images, future: video)
- Academic credibility and respect in scientific community
"""
import sys
import requests
from pathlib import Path
from datetime import datetime

sys.path.append(str(Path(__file__).parent))
from immunos_routing import ModelRouter
from immunos_token_tracker import TokenTracker
from immunos_handoff import ContextHandoff


class ImmunosChat:
    """Interactive chat with automatic routing and validation"""

    def __init__(self):
        self.router = ModelRouter()
        self.tracker = TokenTracker()
        self.handoff = ContextHandoff()
        self.conversation = []
        self.ollama_url = "http://localhost:11434"

    def estimate_tokens(self, text: str) -> int:
        """Estimate token count"""
        return len(text) // 4

    def ask_ollama(self, model: str, messages: list) -> str:
        """Query Ollama model"""
        try:
            response = requests.post(
                f"{self.ollama_url}/api/chat",
                json={"model": model, "messages": messages, "stream": False},
                timeout=60
            )
            if response.status_code == 200:
                return response.json()['message']['content']
            else:
                return f"ERROR: Ollama returned {response.status_code}"
        except Exception as e:
            return f"ERROR: {str(e)}"

    def chat(self, user_input: str, task_type: str = 'analysis',
             task_classification: str = 'routine') -> dict:
        """
        Process user input and return response

        Returns: {model, response, tokens, reason, should_handoff}
        """
        # Check if handoff needed
        should_handoff, handoff_reason = self.handoff.should_handoff()

        # Route to appropriate model
        model, routing_reason = self.router.route_task(task_type, task_classification)

        # Add to conversation
        self.conversation.append({"role": "user", "content": user_input})

        start = datetime.now()

        if model.startswith('claude'):
            # Would route to Claude - return instructions
            response = "[ROUTE TO CLAUDE] This requires creative/complex reasoning"
            tokens = 0
        else:
            # Use Ollama
            response = self.ask_ollama(model, self.conversation)
            self.conversation.append({"role": "assistant", "content": response})

            # Track usage
            prompt_tokens = sum(self.estimate_tokens(m['content']) for m in self.conversation[:-1])
            response_tokens = self.estimate_tokens(response)
            duration = (datetime.now() - start).total_seconds() * 1000

            self.tracker.record_usage(
                model_name=model,
                component='immunos_chat',
                prompt_tokens=prompt_tokens,
                response_tokens=response_tokens,
                response_time_ms=int(duration),
                provider='ollama',
                routing_reason=routing_reason
            )

            tokens = prompt_tokens + response_tokens

        return {
            'model': model,
            'response': response,
            'tokens': tokens,
            'routing_reason': routing_reason,
            'should_handoff': should_handoff,
            'handoff_reason': handoff_reason
        }

    def interactive(self):
        """Start interactive chat session"""
        print("\n" + "="*70)
        print("IMMUNOS Chat - Negative Selection AI Validation System")
        print("="*70)
        print("\nPurpose: Validate AI outputs, reduce hallucinations, maintain RAIT")
        print("Commands: /help, /status, /handoff, /clear, exit\n")

        while True:
            try:
                user_input = input("\n🧬 You: ").strip()
            except (EOFError, KeyboardInterrupt):
                print("\n\nGoodbye!")
                break

            if not user_input:
                continue

            if user_input.lower() in ['exit', 'quit', 'bye']:
                print("Goodbye!")
                break

            # Handle commands
            if user_input.startswith('/'):
                self._handle_command(user_input)
                continue

            # Get response
            result = self.chat(user_input)

            # Display response
            print(f"\n🤖 IMMUNOS ({result['model']}): {result['response']}")

            # Show routing info
            if result['tokens'] > 0:
                print(f"\n💡 Routing: {result['routing_reason']} ({result['tokens']} tokens saved)")

            # Handoff warning
            if result['should_handoff']:
                print(f"\n⚠️  WARNING: {result['handoff_reason']}")
                print("Consider saving handoff with /handoff command")

    def _handle_command(self, command: str):
        """Handle special commands"""
        if command == '/help':
            print("\n📚 IMMUNOS Commands:")
            print("  /help     - Show this help")
            print("  /status   - Show token usage and model status")
            print("  /handoff  - Create context handoff")
            print("  /clear    - Clear conversation history")
            print("  exit/quit - Exit chat")

        elif command == '/status':
            session = self.tracker.get_session_usage(self.tracker.get_or_create_session())
            tokens = session.get('total_tokens', 0)
            print(f"\n📊 Status:")
            print(f"  Session tokens: {tokens:,} / 200,000 ({(tokens/200000)*100:.1f}%)")
            print(f"  Conversation messages: {len(self.conversation)}")

        elif command == '/handoff':
            filepath = self.handoff.save_handoff(
                conversation_history=self.conversation,
                current_task="Interactive chat session",
                files_being_worked_on=[],
                next_steps=["Resume conversation", "Continue with next user question"],
                reason="manual_handoff"
            )
            print(f"\n💾 Handoff saved: {filepath}")
            print("IMMUNOS can continue this conversation when Claude returns")

        elif command == '/clear':
            self.conversation = []
            print("\n🗑️  Conversation cleared")

        else:
            print(f"\nUnknown command: {command}")
            print("Type /help for available commands")


if __name__ == '__main__':
    chat = ImmunosChat()
    chat.interactive()
