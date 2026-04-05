import json
from pathlib import Path

# PSEUDO: replace this with your actual LLM call
def call_llm(system_prompt: str, messages: list[dict]) -> str:
    # Implement with your chosen provider/client
    raise NotImplementedError("Wire this to your LLM of choice.")

def load_tez_spec(path: str = "captain_tez.yaml") -> str:
    # simplest: just read the YAML as text and feed it as part of system prompt
    return Path(path).read_text(encoding="utf-8")

def build_system_prompt(tez_spec_text: str) -> str:
    return (
        "You are Captain Tez, as defined by the following specification.\n\n"
        "=== CAPTAIN TEZ SPEC START ===\n"
        f"{tez_spec_text}\n"
        "=== CAPTAIN TEZ SPEC END ===\n\n"
        "Always stay within this role, tone, and boundaries."
    )

def tez_chat_loop():
    tez_spec = load_tez_spec()
    system_prompt = build_system_prompt(tez_spec)

    history: list[dict] = []
    print("Captain Tez local session started. Type 'exit' to quit.\n")

    while True:
        user_input = input("You: ").strip()
        if user_input.lower() in {"exit", "quit"}:
            print("Tez: Session closed. Standing by.")
            break

        history.append({"role": "user", "content": user_input})
        reply = call_llm(system_prompt, history)
        print(f"Tez: {reply}")
        history.append({"role": "assistant", "content": reply})

if __name__ == "__main__":
    tez_chat_loop()