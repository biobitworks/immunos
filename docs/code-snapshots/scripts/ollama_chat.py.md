---
source: /Users/byron/projects/scripts/ollama_chat.py
relative: scripts/ollama_chat.py
generated_at: 2025-12-23 10:28
---

```python
#!/usr/bin/env python3
"""
Simple Ollama Chat Interface
Usage: python3 scripts/ollama_chat.py [model_name]
"""
import sys
import requests
import json

def chat(model="qwen2.5-coder:7b"):
    """Interactive chat with Ollama model"""
    print(f"\n🤖 Chatting with {model}")
    print("=" * 60)
    print("Type 'exit' or 'quit' to end conversation\n")

    conversation = []

    while True:
        # Get user input
        try:
            user_input = input("You: ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\n\nGoodbye!")
            break

        if user_input.lower() in ['exit', 'quit', 'bye']:
            print("Goodbye!")
            break

        if not user_input:
            continue

        # Add to conversation history
        conversation.append({"role": "user", "content": user_input})

        # Call Ollama API
        try:
            response = requests.post(
                "http://localhost:11434/api/chat",
                json={
                    "model": model,
                    "messages": conversation,
                    "stream": False
                },
                timeout=60
            )
            response.raise_for_status()

            # Get assistant's reply
            result = response.json()
            assistant_msg = result['message']['content']

            # Add to conversation
            conversation.append({"role": "assistant", "content": assistant_msg})

            # Display response
            print(f"\n{model}: {assistant_msg}\n")

        except requests.exceptions.RequestException as e:
            print(f"\n❌ Error: {e}")
            print("Make sure Ollama server is running: ollama serve\n")
        except KeyError:
            print(f"\n❌ Unexpected response format from Ollama\n")

if __name__ == "__main__":
    # Get model from command line or use default
    model = sys.argv[1] if len(sys.argv) > 1 else "qwen2.5-coder:7b"

    print("\nAvailable models:")
    print("  - qwen2.5-coder:7b  (best for code/technical)")
    print("  - deepseek-r1:14b   (best for reasoning)")
    print("  - qwen2.5:1.5b      (lightweight)")

    chat(model)

```
