# ---------------------------------------------------------
# STARSHIP BOOT SEQUENCE — DEMO MODE (EXTENDED)
# Author: Leif William Sogge (Captain)
# Description:
#   Simulates:
#     - cryo-tube roast
#     - boot-up sequence
#     - login prompt
#     - diagnostics panel
#     - lost-QR recovery flow
#     - blocked core safety modification attempt
# ---------------------------------------------------------

import time
import os

def clear():
    os.system('cls' if os.name == 'nt' else 'clear')

# ---------------------------------------------------------
# CRYO-TUBE ROAST (short comedic beat)
# ---------------------------------------------------------
def cryo_roast():
    clear()
    print("CRYO-TUBE STATUS: OPENING...\n")
    time.sleep(1.2)
    print("Scientist (half-awake): \"Ugh... who scheduled me for this shift?\"")
    time.sleep(1.5)
    print("Scientist: \"If anyone touches my snacks again, I'm filing a report...\"\n")
    time.sleep(2)
    print("SYSTEM: Bio-signal detected. Crew member awake.\n")
    time.sleep(1.5)

# ---------------------------------------------------------
# BOOT-UP SEQUENCE
# ---------------------------------------------------------
def boot_sequence():
    clear()
    print("BOOTING STARSHIP SYSTEMS...\n")
    time.sleep(1.5)

    print("Initializing display...")
    time.sleep(1.2)
    clear()

    # JUMP -> JUMP ON CREW!?
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
    print("please set a 5-digit numeric or haptic password.")
    pwd = input("Set password: ")
    clear()

    print("🔐 Security lock engaged.")
    print("Digital token ready.\n")
    time.sleep(2)

    return user  # so we can reuse in recovery if needed

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
    print("BFS  — Bio-Frame Shell .............. ONLINE")
    print("BID  — Bio Identity Descriptor ...... ONLINE")
    print("HSS  — Hyperdrive Stability Sim ..... ONLINE\n")

    print("GLOBAL SYSTEMS (STANDBY)")
    print("----------------------------------------")
    print("GLP  — Geometric Location Pinning ... STANDBY")
    print("Global Navigation Mesh .............. STANDBY")
    print("Deep-Field Coordinate Resolver ...... STANDBY\n")

    print("NOTE:")
    print("These global modules are not active in demo mode.")
    print("They will activate only in full deployment.\n")

    print("========================================")
    print("SYSTEM READY — AWAITING CAPTAIN INPUT")
    print("========================================\n")

# ---------------------------------------------------------
# LOST-QR RECOVERY FLOW (DEMO)
# ---------------------------------------------------------
def recover_qr_token():
    clear()
    print("========================================")
    print("      CREDENTIAL RECOVERY TERMINAL")
    print("========================================\n")
    time.sleep(1.5)

    print("SYSTEM: QR token reported lost.")
    print("Entering RECOVERY MODE for this user only.\n")
    time.sleep(2)

    print("NOTE:")
    print("- No jump commands allowed.")
    print("- No system-level changes allowed.")
    print("- Only identity and access tools are available.\n")
    time.sleep(3)

    # Step 1 — Identity re-anchoring via birth data
    print("STEP 1 — IDENTITY RE-ANCHORING\n")
    legal_name = input("Enter legal name (as on record): ")
    birth_time = input("Enter time of birth (HH:MM, 24h): ")
    birth_place = input("Enter place of birth (City, Country): ")
    clear()

    # In a real system, this would check BFS/BID + stored records.
    # Here we just simulate a successful match.
    print("Verifying identity against BID profile...\n")
    time.sleep(2)
    print(f"Identity match confirmed for: {legal_name}")
    print(f"Birth time: {birth_time}")
    print(f"Birth place: {birth_place}\n")
    time.sleep(2)

    # Step 2 — One-time recovery session
    print("STEP 2 — ONE-TIME RECOVERY SESSION\n")
    print("Recovery Mode active for this user only.")
    print("No navigation or core system changes permitted.\n")
    time.sleep(2.5)

    # Step 3 — New QR token issuance
    print("STEP 3 — NEW QR TOKEN ISSUANCE\n")
    time.sleep(1.5)
    print("Generating new digital token...")
    time.sleep(2)
    print("Revoking old token...")
    time.sleep(1.5)
    print("Old token permanently invalidated.\n")
    time.sleep(2)

    # Step 4 — New password/haptic lock
    print("STEP 4 — ACCESS LOCK RESET\n")
    print("Please set a NEW 5-digit numeric or haptic password.")
    new_pwd = input("Set new password: ")
    clear()

    print("🔐 Security lock engaged with new credentials.")
    print("New digital token bound to your BID profile.\n")
    time.sleep(2)

    # Step 5 — Captain notification (simulated)
    print("STEP 5 — CAPTAIN NOTIFICATION\n")
    time.sleep(1.5)
    print("Drafting notification to: spaceleafcorp@outlook.com")
    time.sleep(2)
    print("SUBJECT: CREDENTIAL RECOVERY EVENT")
    print("BODY:")
    print(f"- User: {legal_name}")
    print(f"- Birth time: {birth_time}")
    print(f"- Birth place: {birth_place}")
    print("- Old token: REVOKED")
    print("- New token: ISSUED")
    print("- Mode: DEMO SIMULATION (no real email sent)\n")
    time.sleep(4)

    print("Recovery sequence complete.")
    print("User access restored at previous permission level.\n")
    time.sleep(2)

# ---------------------------------------------------------
# BLOCKED CORE SAFETY LOGIC MODIFICATION ATTEMPT
# ---------------------------------------------------------
def attempt_core_safety_modification():
    clear()
    print("========================================")
    print("   CORE SAFETY MODULE ACCESS TERMINAL")
    print("========================================\n")
    time.sleep(1.5)

    print("WARNING:")
    print("You are attempting to access a protected system area.\n")
    time.sleep(2)

    print("Modules requested:")
    print("- BFS  (Bio-Frame Shell)")
    print("- BID  (Bio Identity Descriptor)")
    print("- HSS  (Hyperdrive Stability Simulation)")
    print("- Physics Constraint Layer\n")
    time.sleep(3)

    print("Verifying authorization...\n")
    time.sleep(2)

    # This is the hard wall
    print("🚫 ACCESS DENIED — CAPTAIN CONSENT REQUIRED 🚫\n")
    time.sleep(2)

    print("These modules are protected by immutable safety laws.")
    print("They cannot be altered, overridden, or bypassed by any user level.\n")
    time.sleep(3)

    print("A notification has been drafted to the Captain:")
    print("   spaceleafcorp@outlook.com\n")
    time.sleep(2)

    print("SUBJECT: BLOCKED CORE SAFETY ACCESS ATTEMPT")
    print("BODY:")
    print("- A user attempted to modify protected safety logic.")
    print("- Modules targeted: BFS, BID, HSS, Physics Layer")
    print("- Action: BLOCKED")
    print("- Captain consent required for any future requests.\n")
    time.sleep(4)

    print("System returning to safe operational mode.\n")
    time.sleep(2)

# ---------------------------------------------------------
# MAIN EXECUTION
# ---------------------------------------------------------
if __name__ == "__main__":
    cryo_roast()
    boot_sequence()
    user_id = login_prompt()
    diagnostics_panel()

    # Simple post-boot menu
    while True:
        print("OPTIONS:")
        print("1) Simulate lost-QR recovery")
        print("2) Attempt to modify core safety logic (demo block)")
        print("3) Exit demo\n")
        choice = input("Select an option: ").strip()

        if choice == "1":
            recover_qr_token()
            diagnostics_panel()

        elif choice == "2":
            attempt_core_safety_modification()
            diagnostics_panel()

        elif choice == "3":
            clear()
            print("Shutting down demo. Systems returning to standby.\n")
            break

        else:
            print("\nInvalid selection. Try again.\n")
            time.sleep(1.5)