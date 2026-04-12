import json
import random

class ToolRegistry:
    def __init__(self):
        self._tools: Dict[str, Callable[..., str]] = {}

    def register(self, name: str, func: Callable[..., str], description: str = ""):
        self._tools[name] = func
        setattr(self, name, func)  # optional direct attr access

    def list(self) -> Dict[str, str]:
        return {name: getattr(func, "__doc__", "") or "" for name, func in self._tools.items()}

    def call(self, name: str, *args, **kwargs) -> str:
        if name not in self._tools:
            return f"[Tool Error] No such tool: {name}"
        try:
            return self._tools[name](*args, **kwargs)
        except Exception as e:
            return f"[Tool Error] {name} failed: {e}"


# Example tools

def tool_time():
    """Return the current system time."""
    return time.strftime("%Y-%m-%d %H:%M:%S")

def tool_randint(low: int = 0, high: int = 100):
    """Return a random integer in range [low, high]."""
    return str(random.randint(low, high))

def tool_echo(text: str):
    """Echo the provided text."""
    return text
