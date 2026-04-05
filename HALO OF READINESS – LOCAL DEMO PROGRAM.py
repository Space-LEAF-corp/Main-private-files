#!/usr/bin/env python3

# ============================================================
#  HALO OF READINESS – LOCAL DEMO PROGRAM
#  Star-Key Authentication + Simple Halo Interaction
#  Runs fully offline. No external dependencies.
# ============================================================

import time

# ------------------------------------------------------------
#  STAR KEY AUTHENTICATION
# ------------------------------------------------------------

STAR_KEY = "*"   # Local-only symbolic key

def authenticate():
    print("\n=== STAR KEY ACCESS REQUIRED ===")
    key = input("Enter Star Key to continue: ").strip()
    if key == STAR_KEY:
        print("Access granted. ✦ Channel is clear.\n")
        return True
    else:
        print("Access denied. ✦ Star Key mismatch.\n")
        return False


# ------------------------------------------------------------
#  HALO DEFINITIONS
# ------------------------------------------------------------

HALOS = {
    "Shell Mode": {
        "symbol": "🐚",
        "meaning": "I’m quiet, recharging, and prefer solitude.",
        "boundary": "Please give me space while I reset."
    },
    "Panther Purr": {
        "symbol": "🐾",
        "meaning": "I feel safe, playful, and at ease.",
        "boundary": "No need to check in — I’m good."
    },
    "Turtle Time": {
        "symbol": "🐢",
        "meaning": "I’m moving slowly — please be patient.",
        "boundary": "Gentle pacing appreciated."
    },
    "Dragon Bubble": {
        "symbol": "🐉",
        "meaning": "I’m in protective mode.",
        "boundary": "Please respect my space."
    },
    "Dolphin Loop": {
        "symbol": "🐬",
        "meaning": "I’m in a silly, splashy mood.",
        "boundary": "Joyful connection welcome."
    }
}


# ------------------------------------------------------------
#  DISPLAY HALOS
# ------------------------------------------------------------

def show_halos():
    print("\n=== HALO OF READINESS MENU ===")
    for name, data in HALOS.items():
        print(f"{data['symbol']}  {name}")
    print("==============================\n")


# ------------------------------------------------------------
#  HALO DETAIL VIEW
# ------------------------------------------------------------

def view_halo():
    choice = input("Enter the name of a halo to view details: ").strip()
    if choice in HALOS:
        halo = HALOS[choice]
        print("\n--- HALO DETAILS ---")
        print(f"Symbol:   {halo['symbol']}")
        print(f"Meaning:  {halo['meaning']}")
        print(f"Boundary: {halo['boundary']}")
        print("---------------------\n")
    else:
        print("Halo not found.\n")


# ------------------------------------------------------------
#  MAIN PROGRAM LOOP
# ------------------------------------------------------------

def main_program():
    print("Halo of Readiness System Online ✦")
    print("Type 'menu' to view halos, 'view' to inspect one, or 'exit' to close.\n")

    while True:
        cmd = input("Command: ").strip().lower()

        if cmd == "menu":
            show_halos()

        elif cmd == "view":
            view_halo()

        elif cmd == "exit":
            print("Closing program. ✦ Seal holds.")
            break

        else:
            print("Unknown command. Try 'menu', 'view', or 'exit'.\n")


# ------------------------------------------------------------
#  ENTRY POINT
# ------------------------------------------------------------

def boot_sequence():
    print("Booting Halo of Readiness Demo...")
    time.sleep(0.8)
    print("Ceremonial systems warming...")
    time.sleep(0.8)
    print("Awaiting Star Key...\n")
    time.sleep(0.5)

    if authenticate():
        main_program()
    else:
        print("Program terminated.\n")


# ------------------------------------------------------------
#  RUN ON USER REQUEST
# ------------------------------------------------------------

if __name__ == "__main__":
    user_input = input("Say 'run program please' to begin: ").strip().lower()
    if user_input == "run program please":
        boot_sequence()
    else:
        print("Program not started. Awaiting correct phrase.\n")