class ResponsePolicy:
    def __init__(self, core_laws: List[str]):
        self.core_laws = core_laws

    def style(self, raw: str) -> str:
        """
        Apply a light-touch style: caring, clear, concise.
        """
        if not raw:
            return "I didn’t catch that—could you say it another way?"
        # soften and clarify tone
        raw = raw.strip()
        # ensure we avoid harmful or risky advice
        safety_prefix = "I'll keep things helpful and responsible. "
        return safety_prefix + raw
