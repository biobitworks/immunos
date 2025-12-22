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
import os
import sys
import json
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
        self.base_path = Path(__file__).parent.parent

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

    def ask_anthropic_messages(self, base_url: str, api_key: str, model: str, messages: list) -> str:
        """Query Anthropic Messages API."""
        if not api_key:
            return "ERROR: Missing ANTHROPIC_AUTH_TOKEN"
        if not model:
            return "ERROR: Missing Anthropic model"

        url = base_url.rstrip("/") + "/v1/messages"
        headers = {
            "Content-Type": "application/json",
            "x-api-key": api_key,
            "anthropic-version": "2023-06-01"
        }

        payload = {
            "model": model,
            "max_tokens": 1024,
            "messages": [{"role": m["role"], "content": m["content"]} for m in messages]
        }

        try:
            response = requests.post(url, headers=headers, json=payload, timeout=60)
            if response.status_code != 200:
                return f"ERROR: Anthropic returned {response.status_code}"

            data = response.json()
            content = data.get("content", [])
            if isinstance(content, list):
                return "".join(block.get("text", "") for block in content if isinstance(block, dict))
            return str(content)
        except Exception as e:
            return f"ERROR: {str(e)}"

    def ask_openai_compatible(self, base_url: str, api_key: str, model: str, messages: list) -> str:
        """Query an OpenAI-compatible endpoint."""
        if not base_url or not model:
            return "ERROR: Missing orchestrator base_url or model"

        url = base_url.rstrip("/") + "/chat/completions"
        headers = {"Content-Type": "application/json"}
        if api_key:
            headers["Authorization"] = f"Bearer {api_key}"

        try:
            response = requests.post(
                url,
                headers=headers,
                json={"model": model, "messages": messages, "stream": False},
                timeout=60
            )
            if response.status_code != 200:
                return f"ERROR: Orchestrator returned {response.status_code}"

            payload = response.json()
            return payload["choices"][0]["message"]["content"]
        except Exception as e:
            return f"ERROR: {str(e)}"

    def _load_orchestrator_config(self) -> dict:
        """Load orchestrator configuration from disk."""
        config_path = self.base_path / ".immunos" / "config" / "orchestrator.json"
        default_config = {
            "connectivity": "online",
            "online_backend": {
                "provider": "claude_code",
                "model": "",
                "base_url": "",
                "api_key_env": "ANTHROPIC_AUTH_TOKEN"
            },
            "offline_backend": {
                "provider": "ollama",
                "model": "qwen2.5-coder:7b",
                "base_url": "http://localhost:11434",
                "api_key_env": ""
            }
        }

        if not config_path.exists():
            config_path.parent.mkdir(parents=True, exist_ok=True)
            with config_path.open("w", encoding="utf-8") as handle:
                json.dump(default_config, handle, indent=2)
            return default_config

        try:
            with config_path.open(encoding="utf-8") as handle:
                data = json.load(handle)
        except json.JSONDecodeError:
            data = {}

        merged = default_config.copy()
        merged.update(data)
        merged["online_backend"].update(data.get("online_backend", {}))
        merged["offline_backend"].update(data.get("offline_backend", {}))
        return merged

    def _resolve_orchestrator_backend(self, connectivity: str = None) -> dict:
        """Resolve active orchestrator backend based on connectivity."""
        config = self._load_orchestrator_config()
        if connectivity in {"online", "offline"}:
            config["connectivity"] = connectivity

        online_backend = config.get("online_backend", {})
        offline_backend = config.get("offline_backend", {})

        def backend_available(backend: dict) -> tuple:
            provider = backend.get("provider")
            if provider == "ollama":
                return True, ""

            if provider == "claude_code":
                api_key_env = backend.get("api_key_env") or "ANTHROPIC_AUTH_TOKEN"
                if not os.getenv(api_key_env):
                    return False, "missing_api_key"
                model = backend.get("model") or os.getenv("ANTHROPIC_MODEL")
                if not model:
                    return False, "missing_model"
                return True, ""

            if provider == "chatgpt":
                api_key_env = backend.get("api_key_env") or "OPENAI_API_KEY"
                if not os.getenv(api_key_env):
                    return False, "missing_api_key"
                model = backend.get("model") or os.getenv("OPENAI_MODEL")
                if not model:
                    return False, "missing_model"
                return True, ""

            if provider == "openrouter":
                api_key_env = backend.get("api_key_env") or "OPENROUTER_API_KEY"
                if not os.getenv(api_key_env):
                    return False, "missing_api_key"
                model = backend.get("model") or os.getenv("OPENROUTER_MODEL")
                if not model:
                    return False, "missing_model"
                return True, ""

            if provider == "local_server":
                if not backend.get("base_url"):
                    return False, "missing_base_url"
                if not backend.get("model"):
                    return False, "missing_model"
                return True, ""

            return False, "unknown_provider"

        backend = online_backend if config.get("connectivity") == "online" else offline_backend
        fallback_reason = ""

        if config.get("connectivity") == "online":
            ok, reason = backend_available(backend)
            if not ok:
                backend = offline_backend
                config["connectivity"] = "offline"
                fallback_reason = reason

        return {
            "connectivity": config.get("connectivity", "online"),
            "provider": backend.get("provider"),
            "model": backend.get("model"),
            "base_url": backend.get("base_url"),
            "api_key_env": backend.get("api_key_env"),
            "fallback_reason": fallback_reason
        }

    def chat(self, user_input: str, task_type: str = 'analysis',
             task_classification: str = 'routine', mode: str = 'orchestrator',
             domain: str = None, connectivity: str = None) -> dict:
        """
        Process user input and return response

        Returns: {model, response, tokens, reason, should_handoff}
        """
        # Check if handoff needed
        should_handoff, handoff_reason = self.handoff.should_handoff()

        if domain:
            task_type = domain

        if mode == "orchestrator":
            backend = self._resolve_orchestrator_backend(connectivity)
            model = backend.get("model") or ""
            routing_reason = f"orchestrator:{backend.get('provider')}:{backend.get('connectivity')}"
            if backend.get("fallback_reason"):
                routing_reason += f":fallback_{backend['fallback_reason']}"
        else:
            model, routing_reason = self.router.route_task(task_type, task_classification)

        # Add to conversation
        self.conversation.append({"role": "user", "content": user_input})

        start = datetime.now()

        if mode == "orchestrator":
            provider = backend.get("provider")
            if provider == "ollama":
                self.ollama_url = backend.get("base_url") or self.ollama_url
                model_name = model or os.getenv("IMMUNOS_OFFLINE_MODEL") or "qwen2.5-coder:7b"
                response = self.ask_ollama(model_name, self.conversation)
                model = model_name
            elif provider == "claude_code":
                base_url = backend.get("base_url") or os.getenv("ANTHROPIC_BASE_URL") or "https://api.anthropic.com"
                api_key_env = backend.get("api_key_env") or "ANTHROPIC_AUTH_TOKEN"
                api_key = os.getenv(api_key_env, "")
                model_name = model or os.getenv("ANTHROPIC_MODEL", "")
                response = self.ask_anthropic_messages(base_url, api_key, model_name, self.conversation)
                model = model_name or model
            else:
                base_url = backend.get("base_url") or os.getenv("IMMUNOS_ORCHESTRATOR_BASE_URL") or os.getenv("ANTHROPIC_BASE_URL") or ""
                api_key_env = backend.get("api_key_env") or "ANTHROPIC_AUTH_TOKEN"
                api_key = os.getenv(api_key_env, "")
                model_name = model or os.getenv("ANTHROPIC_MODEL", "")
                response = self.ask_openai_compatible(base_url, api_key, model_name, self.conversation)
                model = model_name or model
            self.conversation.append({"role": "assistant", "content": response})
            provider_name = provider or "orchestrator"
        elif model.startswith('claude'):
            # Would route to Claude - return instructions
            response = "[ROUTE TO CLAUDE] This requires creative/complex reasoning"
            tokens = 0
            provider_name = "claude"
        else:
            # Use Ollama
            response = self.ask_ollama(model, self.conversation)
            self.conversation.append({"role": "assistant", "content": response})
            provider_name = "ollama"

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
            provider=provider_name,
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
