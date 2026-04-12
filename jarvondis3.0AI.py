class jarvondis3.0AI:
    def __init__(self, base: jarvondis3.0Base, memory: ConversationMemory, tools: ToolRegistry, policy: ResponsePolicy):
        self.base = base
        self.memory = memory
        self.tools = tools
        self.policy = policy
        self.intent = IntentDetector()

    def process_input(self, user_input: str) -> str:
        """
        Route input: detect intent, optionally call tools, generate response, and apply policy.
        """
        self.memory.add("user", user_input)
        intent = self.intent.detect(user_input)

        if intent == "get_time":
            result = self.tools.call("time")
            return self.policy.style(f"The current time is {result}.")

        if intent == "random":
            result = self.tools.call("randint", 1, 100)
            return self.policy.style(f"Your random number is {result}.")

        if intent == "echo":
            payload = user_input[5:].strip()  # after "echo "
            result = self.tools.call("echo", payload)
            return self.policy.style(f"Echoing back: {result}")

        if intent == "help":
            tools_list = self.tools.list()
            lines = [f"- {name}: {desc}" for name, desc in tools_list.items()]
            help_text = "Here are available tools and commands:\n" + "\n".join(lines) + \
                        "\nYou can try: 'time', 'random', 'echo your text', 'version'."
            return self.policy.style(help_text)

        if intent == "system_info":
            info = self.base.system_info()
            return self.policy.style(f"System info: version {info['version']}, architecture {self.base.architecture}, "
                                     f"uptime {info['uptime_seconds']}s.")

        # default chat: lightweight rule-based reply using memory summary and core laws
        summary = self.memory.summarize()
        laws_hint = "; ".join(self.base.core_laws[:3])  # highlight a subset
        reply = f"I’m hearing you. {summary}. Guided by: {laws_hint}. How can I help further?"
        self.memory.add("assistant", reply)
        return self.policy.style(reply)
