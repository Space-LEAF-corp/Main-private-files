# ---------------------------------------------------------
# STARSHIP BOOT SEQUENCE — DEMO MODE
# Author: Leif William Sogge (Captain)
# Description:
#   Simulates the cryo‐tube roast, boot‐up sequence,
#   login prompt, and diagnostics panel with ACTIVE
#   and STANDBY systems.
# ---------------------------------------------------------

import time
import os

def clear():
    os.system('cls' if os.name == 'nt' else 'clear')

# ---------------------------------------------------------
# CRYO‐TUBE ROAST (short comedic beat)
# ---------------------------------------------------------
def cryo_roast():
    clear()
    print("CRYO‐TUBE STATUS: OPENING...\n")
    time.sleep(1.2)
    print("Scientist (half‐awake): \"Ugh... who scheduled me for this shift?\"")
    time.sleep(1.5)
    print("Scientist: \"If anyone touches my snacks again, I'm filing a report...\"\n")
    time.sleep(2)
    print("SYSTEM: Bio‐signal detected. Crew member awake.\n")
    time.sleep(1.5)

# ---------------------------------------------------------
# BOOT‐UP SEQUENCE
# ---------------------------------------------------------
def boot_sequence():
    clear()
    print("BOOTING STARSHIP SYSTEMS...\n")
    time.sleep(1.5)

    print("Initializing display...")
    time.sleep(1.2)
    clear()

    # JUMP → JUMP ON CREW!?
    print("J U M P")
    time.sleep(1.2)
    clear()
    print("JUMP ON CREW!?")
    time.sleep(1.5)

# ---------------------------------------------------------
# LOGIN PROMPT
# ---------------------------------------------------------
def login_prompt():
    print("\n----------------------------------------")
    print("Press ENTER to log in.")
    input()
    clear()

    print("LOGIN SECURITY CHECK")
    print("----------------------------------------")
    print("Have BOTH validating QR codes ready")
    print("OR prepare for a timestamped facial photo scan.")
    print("(Photo used only for security logging — not stored.)\n")
    time.sleep(3)

    user = input("Enter Login ID: ")
    clear()

    print(f"Thank you. Have a safe day, {user}.")
    print("Let the stars guide you timely.\n")
    time.sleep(2)

    print("🪪 WELCOME ABOARD\n")
    time.sleep(1.5)

    print("Before downloading your digital token,")
    print("please set a 5‐digit numeric or haptic password.")
    pwd = input("Set password: ")
    clear()

    print("🔐 Security lock engaged.")
    print("Digital token ready.\n")
    time.sleep(2)

# ---------------------------------------------------------
# DIAGNOSTICS PANEL
# ---------------------------------------------------------
def diagnostics_panel():
    clear()
    print("========================================")
    print("        STARSHIP DIAGNOSTICS")
    print("========================================\n")

    print("ACTIVE SYSTEMS")
    print("----------------------------------------")
    print("TTC  — Time Travel Calendar ........ ONLINE")
    print("TTE  — Time Transversal Engine ...... ONLINE")
    print("MFDI — Membrane Flux Index .......... ONLINE")
    print("BFS  — Bio‐Frame Shell .............. ONLINE")
    print("BID  — Bio Identity Descriptor ...... ONLINE")
    print("HSS  — Hyperdrive Stability Sim ..... ONLINE\n")

    print("GLOBAL SYSTEMS (STANDBY)")
    print("----------------------------------------")
    print("GLP  — Geometric Location Pinning ... STANDBY")
    print("Global Navigation Mesh .............. STANDBY")
    print("Deep‐Field Coordinate Resolver ...... STANDBY\n")

    print("NOTE:")
    print("These global modules are not active in demo mode.")
    print("They will activate only in full deployment.\n")

    print("========================================")
    print("SYSTEM READY — AWAITING CAPTAIN INPUT")
    print("========================================\n")

# ---------------------------------------------------------
# MAIN EXECUTION
# ---------------------------------------------------------
if __name__ == "__main__":
    cryo_roast()
    boot_sequence()
    login_prompt()
    diagnostics_panel()