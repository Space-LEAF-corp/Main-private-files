from collections import deque

class ConversationMemory:
    def __init__(self, max_turns: int = 20):
        self.history = deque(maxlen=max_turns)

    def add(self, role: str, text: str):
        self.history.append({"role": role, "text": text, "ts": time.time()})

    def get_context(self) -> List[Dict]:
        return list(self.history)

    def summarize(self) -> str:
        # Simple reducer; replace with LLM summarization later
        last_user = next((h["text"] for h in reversed(self.history) if h["role"] == "user"), "")
        return f"Recent user focus: {last_user[:120]}"
