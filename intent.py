class IntentDetector:
    def detect(self, user_input: str) -> str:
        text = user_input.lower().strip()
        if any(k in text for k in ["time", "clock", "date"]):
            return "get_time"
        if text.startswith("echo "):
            return "echo"
        if any(k in text for k in ["random", "rand", "dice"]):
            return "random"
        if any(k in text for k in ["help", "commands", "tools"]):
            return "help"
        if any(k in text for k in ["version", "info", "uptime"]):
            return "system_info"
        # default chat intent
        return "chat"
