import ui
# PS‐SLM Prototype Engine 0.2
# Local-only, safety-focused version for Space LEAF Corp
# Purpose:
#   - Demonstrate fast phonetic input → slow interpretation → optional sign-language mapping
#   - Explicitly NOT for hacking, exploits, automation, or overriding other programs
#
# Design constraints:
#   - No network access
#   - No file system access
#   - No process control
#   - No hooks into other software
#   - Purely local, educational, and game-mode only

import difflib
import sys

# -----------------------------
# 0. Space LEAF Corp Safety Layer
# -----------------------------

# Words/phrases that indicate intent to misuse, hack, or exploit systems.
# You can expand this list as needed.
DISALLOWED_TERMS = [
    "hack",
    "hacks",
    "cheat",
    "cheats",
    "exploit",
    "exploits",
    "bypass",
    "override",
    "inject",
    "injection",
    "backdoor",
    "back door",
    "ddos",
    "botnet",
    "malware",
    "virus",
    "keylogger",
    "rootkit",
    "privilege escalation",
    "escalate privileges",
    "crack",
    "cracker",
    "bruteforce",
    "brute force",
]

def safety_check(text: str) -> bool:
    """
    Returns True if the input is considered SAFE.
    Returns False if it contains disallowed/harmful intent markers.
    """
    lower = text.lower()
    for term in DISALLOWED_TERMS:
        if term in lower:
            return False
    return True

def enforce_safety_or_exit(text: str):
    """
    If unsafe content is detected, stop the program.
    This keeps PS‐SLM strictly in safe, educational territory.
    """
    if not safety_check(text):
        print("\n[SPACE LEAF CORP SAFETY LAYER]")
        print("This tool is not designed for hacking, exploits, or overriding any systems.")
        print("Detected unsafe intent markers in input. Exiting to protect the protocol.\n")
        sys.exit(1)

# -----------------------------
# 1. Base dictionary for interpretation
# -----------------------------

CORRECTIONS = {
    "paul": "poll",
    "rite": "right",
    "there": "their",
    "ur": "your",
    "u": "you",
    "thru": "through",
    "tho": "though",
    "cuz": "because",
}

# -----------------------------
# 2. Optional ASL gesture mapping (placeholder)
# -----------------------------

ASL_MAP = {
    "poll": "ASL_SIGN_POLL",
    "right": "ASL_SIGN_RIGHT",
    "you": "ASL_SIGN_YOU",
    "because": "ASL_SIGN_BECAUSE",
}

# -----------------------------
# 3. Fast-mode phonetic processor
# -----------------------------

def fast_mode_process(text: str) -> str:
    """
    Fast phonetic → interpreted text.
    Does NOT call any external systems.
    Pure string transformation.
    """
    words = text.lower().split()
    processed = []

    for w in words:
        if w in CORRECTIONS:
            processed.append(CORRECTIONS[w])
        else:
            # fuzzy match for phonetic approximations
            close = difflib.get_close_matches(w, CORRECTIONS.keys(), n=1, cutoff=0.75)
            if close:
                processed.append(CORRECTIONS[close[0]])
            else:
                processed.append(w)

    return " ".join(processed)

# -----------------------------
# 4. Optional ASL mapping layer
# -----------------------------

def asl_layer(text: str) -> str:
    """
    Maps interpreted words to placeholder ASL identifiers.
    This is symbolic only, not a real ASL engine.
    """
    words = text.split()
    mapped = []

    for w in words:
        if w in ASL_MAP:
            mapped.append(ASL_MAP[w])
        else:
            mapped.append(f"[NO_SIGN:{w}]")

    return " ".join(mapped)

# -----------------------------
# 5. Game mode feedback loop
# -----------------------------

def game_mode():
    print("PS‐SLM Game Mode Activated (Space LEAF Corp Safe Edition)")
    print("This tool is for learning, communication, and play only.")
    print("It cannot and must not be used to hack, exploit, or override any systems.\n")
    print("Type fast phonetic text. System will interpret and show ASL mapping.")
    print("Type 'quit' or 'exit' to leave.\n")

    while True:
        raw = input("Fast Input: ")

        # Safety: check every input for disallowed intent markers
        enforce_safety_or_exit(raw)

        if raw.lower() in ["quit", "exit"]:
            print("Exiting PS‐SLM Game Mode.")
            break

        interpreted = fast_mode_process(raw)
        asl_output = asl_layer(interpreted)

        print(f"Interpreted: {interpreted}")
        print(f"ASL Map:    {asl_output}\n")

# -----------------------------
# 6. Run the engine
# -----------------------------

if __name__ == "__main__":
    # No network, no files, no external hooks.
    # Purely local, purely educational.
    game_mode()
v = ui.load_view()
v.present('sheet')
