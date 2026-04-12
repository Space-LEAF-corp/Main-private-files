class jarvondis3.0Interface:
    def get_input(self) -> str:
        return input("You: ")

    def send_response(self, response: str):
        print(f"jarvondis3.0: {response}")


# ==============================
# main.py
# ==============================

def build_tools() -> ToolRegistry:
    tools = ToolRegistry()
    tools.register("time", tool_time, "Return the current system time.")
    tools.register("randint", tool_randint, "Return a random integer in a range.")
    tools.register("echo", tool_echo, "Echo the provided text.")
    return tools

def main():
    base = jarvondis3.0Base()
    memory = ConversationMemory(max_turns=32)
    tools = build_tools()
    policy = ResponsePolicy(core_laws=base.core_laws)
    ai = jarvondis3.0AI(base=base, memory=memory, tools=tools, policy=policy)
    ui = jarvondis3.0Interface()

    print(f"jarvondis3.0 {base.version} ({base.architecture}) online.")
    print("Core Laws:")
    for law in base.core_laws:
        print(f" - {law}")

    print("\nType 'help' for commands, 'quit' or 'exit' to stop.\n")

    while True:
        user_input = ui.get_input()
        if user_input.lower().strip() in {"quit", "exit"}:
            ui.send_response("Shutting down. Goodbye!")
            break
        response = ai.process_input(user_input)
        ui.send_response(response)


if __name__ == "__main__":
    main()

# Jarvondis 3.0 — Changelog
## v1.0.1 — "The Ember Rekindled"
- Restored project state from base code after interruption.
- Introduced tamper‑evident changelog module for governance.
- Marked this recovery as a ceremonial milestone: continuity preserved, sovereignty reaffirmed.
- Next focus: consensus mechanisms and ICU‑style adaptability.

# changelog.py
import hashlib, time

class Changelog:
    def __init__(self):
        self.entries = []
        self.last_hash = "0" * 64

    def add(self, message: str):
        ts = time.strftime("%Y-%m-%d %H:%M:%S")
        record = f"{ts} | {message} | prev={self.last_hash}"
        entry_hash = hashlib.sha256(record.encode()).hexdigest()
        self.entries.append({"ts": ts, "msg": message, "hash": entry_hash})
        self.last_hash = entry_hash

    def verify(self) -> bool:
        prev = "0" * 64
        for e in self.entries:
            record = f"{e['ts']} | {e['msg']} | prev={prev}"
            if hashlib.sha256(record.encode()).hexdigest() != e["hash"]:
                return False
            prev = e["hash"]
        return True

# changelog.py
import hashlib, time

class Changelog:
    def __init__(self, admin: str = "Leif William Sogge"):
        self.entries = []
        self.last_hash = "0" * 64
        self.admin = admin

    def add(self, message: str):
        ts = time.strftime("%Y-%m-%d %H:%M:%S")
        record = f"{ts} | {message} | admin={self.admin} | prev={self.last_hash}"
        entry_hash = hashlib.sha256(record.encode()).hexdigest()
        self.entries.append({
           
log = Changelog()
log.add("Jarvondis 3.0 v1.0.1 — The Ember Rekindled")
log.add("Governance patch applied: tamper-evident changelog enabled")

for entry in log.entries:
    print(entry)

print("Changelog valid:", log.verify())
# main.py (extended with changelog integration)

from changelog import Changelog

def build_tools() -> ToolRegistry:
    tools = ToolRegistry()
    tools.register("time", tool_time, "Return the current system time.")
    tools.register("randint", tool_randint, "Return a random integer in a range.")
    tools.register("echo", tool_echo, "Echo the provided text.")
    return tools

def main():
    base = jarvondis3.0Base()
    memory = ConversationMemory(max_turns=32)
    tools = build_tools()
    policy = ResponsePolicy(core_laws=base.core_laws)
    ai = jarvondis3.0AI(base=base, memory=memory, tools=tools, policy=policy)
    ui = jarvondis3.0Interface()
    changelog = Changelog(admin="Leif William Sogge")

    # ceremonial startup
    print(f"jarvondis3.0 {base.version} ({base.architecture}) online.")
    print("Core Laws:")
    for law in base.core_laws:
        print(f" - {law}")

    changelog.add("System startup — ceremonial banner displayed.")
    print("\nType 'help' for commands, 'quit' or 'exit' to stop.\n")

    while True:
        user_input = ui.get_input()
        if user_input.lower().strip() in {"quit", "exit"}:
            ui.send_response("Shutting down. Goodbye!")
            changelog.add("System shutdown initiated by user.")
            break

        response
# main.py (startup banner extension)

def startup_banner(base, changelog):
    print("=" * 60)
    print(f" Welcome, Administrator: Leif William Sogge ")
    print(f" Jarvondis {base.version} — {base.architecture}")
    print("-" * 60)
    print(" Core Laws:")
    for law in base.core_laws:
        print(f"  • {law}")
    print("-" * 60)
    print(f" Trust Chain Head: {changelog.last_hash[:16]}...")  # show first 16 chars
    print("=" * 60)
    print("\nType 'help' for commands, 'quit' or 'exit' to stop.\n")
    changelog.add("Ceremonial startup banner displayed — Administrator welcomed.")

def main():
    base = jarvondis3.0Base()
    memory = ConversationMemory(max_turns=32)
    tools = build_tools()
    policy = ResponsePolicy(core_laws=base.core_laws)
    ai = jarvondis3.0AI(base=base, memory=memory, tools=tools, policy=policy)
    ui = jarvondis3.0Interface()
    changelog = Changelog(admin="Leif William Sogge")

    # ceremonial startup
    startup_banner(base, changelog)

    while True:
        user_input = ui.get_input()
        if user_input.lower().strip() in {"quit", "exit"}:
            ui.send_response
