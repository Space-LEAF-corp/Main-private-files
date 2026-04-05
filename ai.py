# ai.py

from jarvondisbase import JarvondisBase
from memory import ConversationMemory
from tools import ToolRegistry
from policy import ResponsePolicy
from intent import IntentDetector


class JarvondisAI:
    def __init__(self, base, memory, tools, policy):
        self.base = base
        self.memory = memory
        self.tools = tools
        self.policy = policy
        self.intent = IntentDetector()

    def _handle_get_time(self):
        result = self.tools.call("time")
        return self.policy.style(f"The current time is {result}.")

    def _handle_random(self):
        result = self.tools.call("randint", 1, 100)
        return self.policy.style(f"Your random number is {result}.")

    def _handle_echo(self, user_input):
        payload = user_input[5:].strip()
        result = self.tools.call("echo", payload)
        return self.policy.style(f"Echoing back: {result}")

    def _handle_help(self):
        tools_list = self.tools.list()
        lines = [f"- {name}: {desc}" for name, desc in tools_list.items()]
        help_text = (
            "Here are available tools and commands:\n"
            + "\n".join(lines)
            + "\nYou can try: 'time', 'random', 'echo your text', 'version'."
        )
        return self.policy.style(help_text)

    def _handle_system_info(self):
        info = self.base.system_info()
        text = (
            f"System info: version {info['version']}, "
            f"architecture {info['architecture']}, "
            f"uptime {info['uptime_seconds']}s, "
            f"firewall={info['firewall_name']}."
        )
        return self.policy.style(text)

    def _handle_firewall_status(self):
        fw = self.base.firewall
        text = (
            f"Diamond Firewall Status:\n"
            f" - Name: {fw.get('name')}\n"
            f" - Architecture: {fw.get('architecture')}\n"
            f" - Trust Chain: {fw.get('trust_chain', {}).get('hash_algorithm')}\n"
            f" - Facets: {fw.get('facets', {}).get('count')}"
        )
        return self.policy.style(text)

    def _handle_chat(self):
        summary = self.memory.summarize()
        laws_hint = "; ".join(self.base.core_laws[:3])
        reply = (
            f"I’m hearing you. {summary}. "
            f"Guided by: {laws_hint}. How can I help further?"
        )
        self.memory.add("assistant", reply)
        return self.policy.style(reply)

    def process_input(self, user_input: str) -> str:
        self.memory.add("user", user_input)
        intent = self.intent.detect(user_input)

        if intent == "get_time":
            return self._handle_get_time()
        if intent == "random":
            return self._handle_random()
        if intent == "echo":
            return self._handle_echo(user_input)
        if intent == "help":
            return self._handle_help()
        if intent == "system_info":
            return self._handle_system_info()
        if intent == "firewall_status":
            return self._handle_firewall_status()

        return self._handle_chat()
