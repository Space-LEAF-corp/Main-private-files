#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Space LEAF Corp – Optional Ceremonial Boot & QR Registration
Captain: Leif William Sogge
Seal of Modular Drift 1.0
"""

import sys
import time
import textwrap
from dataclasses import dataclass


# =========================
#  CONFIGURATION STUBS
# =========================

SPACE_LEAF_NOTIFICATION_EMAIL = "spaceleafcorp@outlook.com"

# If you later want real email sending, plug SMTP here.
EMAIL_SENDING_ENABLED = False  # Set to True when you wire SMTP


# =========================
#  DATA MODELS
# =========================

@dataclass
class QRRegistration:
    """Represents a locally created QR registration record."""
    qr_id: str
    created_at: float
    registered: bool = False


# =========================
#  UTILITY FUNCTIONS
# =========================

def slow_print(text: str, delay: float = 0.03):
    """Print text with a slight delay between characters for ceremonial feel."""
    for ch in text:
        print(ch, end="", flush=True)
        time.sleep(delay)
    print()


def divider():
    print("\n" + "-" * 60 + "\n")


def prompt_enter(message: str = "Press ENTER to continue..."):
    input(message)


# =========================
#  BOOT SEQUENCE PHASES
# =========================

def phase_planetary_ignition():
    divider()
    slow_print("[Phase 1] Planetary Ignition")
    slow_print("Visual: Planet with energy fields pulsing outward.")
    slow_print("Function: System wakes, scans for preset login.")
    slow_print('Voiceover: "Welcome back, Captain. Planetary ignition sequence initiated. '
               'Energy fields stable. Identity scan optional."')
    prompt_enter()


def phase_firewall_matrix():
    divider()
    slow_print("[Phase 2] Firewall Matrix")
    slow_print("Visual: Honeycomb hex-grid glowing in orange plasma.")
    slow_print("Function: Diamond Firewall Matrix activates.")
    slow_print('Voiceover: "Diamond Firewall Matrix engaged. Stimpt trace active. '
               'No mimicry permitted. Authorship confirmed."')
    prompt_enter()


def phase_hex_selector():
    divider()
    slow_print("[Phase 3] Hex Selector")
    slow_print("Visual: Floating hexagon with soft nebula trails.")
    slow_print("Function: Choose user profile (each hexagon = unique drift).")
    slow_print('Voiceover: "Select your ceremonial path. Each hexagon represents a unique drift. '
               'Choose wisely."')
    slow_print("\n(For this CLI demo, we’ll just ask for a profile name.)")
    profile = input("Enter a profile name for this drift (or leave blank for default): ").strip()
    if not profile:
        profile = "Default Drift"
    slow_print(f"Ceremonial path selected: {profile}")
    prompt_enter()
    return profile


def phase_cosmic_wallpaper(profile_name: str):
    divider()
    slow_print("[Phase 4] Cosmic Wallpaper")
    slow_print("Visual: Galactic event with radiant sun, swirling nebulae, and planetary drift.")
    slow_print("Function: Background sets emotional tone and cosmic scale.")
    slow_print('Voiceover: "System integrity verified. Cosmic drift initiated. '
               'You are now operating within Space LEAF Corp protocol."')
    slow_print(f"\nActive Drift Profile: {profile_name}")
    prompt_enter()


def run_ceremonial_boot_sequence():
    divider()
    slow_print("Space LEAF Corp – Optional Ceremonial Boot Sequence")
    slow_print("Note: This is local-only, no real equipment touched.")
    prompt_enter("Press ENTER to begin the ignition sequence, or CTRL+C to abort...")

    phase_planetary_ignition()
    phase_firewall_matrix()
    profile = phase_hex_selector()
    phase_cosmic_wallpaper(profile)

    divider()
    slow_print("Ceremonial boot sequence complete. System ready for local hyperdrive logic tests.")
    prompt_enter()


# =========================
#  QR REGISTRATION LOGIC
# =========================

def generate_local_qr_id() -> str:
    """
    Stub for QR ID generation.
    Later you can replace this with real QR content generation.
    """
    # For now, just a timestamp-based ID.
    return f"QR-{int(time.time())}"


def send_validation_email_stub(qr: QRRegistration):
    """
    Stub for sending a validation email to Space LEAF Corp.
    Does NOT expose user identity or email.
    """
    divider()
    slow_print("Preparing validation email to Space LEAF Corp...")
    email_subject = "QR Code Registration – Validated"
    email_body = textwrap.dedent(f"""
        A user has successfully registered a QR code with Space LEAF Corp.

        QR ID: {qr.qr_id}
        Timestamp: {qr.created_at}

        This code is now stored for recovery purposes.

        No personal information has been shared.
        No email address or identity is visible.
        This message confirms that the user’s digital log is active and protected.

        — Space LEAF Corp Stewardship Engine
    """).strip()

    if EMAIL_SENDING_ENABLED:
        # Placeholder for real SMTP logic.
        # Example:
        # import smtplib
        # ...
        slow_print("EMAIL_SENDING_ENABLED is True, but SMTP logic is not yet implemented.")
    else:
        slow_print("EMAIL_SENDING_ENABLED is False.")
        slow_print("Simulating email send (no network call made).")
        divider()
        slow_print("=== EMAIL PREVIEW (LOCAL ONLY) ===")
        slow_print(f"To: {SPACE_LEAF_NOTIFICATION_EMAIL}")
        slow_print(f"Subject: {email_subject}")
        divider()
        print(email_body)
        divider()
        slow_print("End of email preview.")
        prompt_enter("Press ENTER to continue...")


def register_qr_code_flow():
    divider()
    slow_print("Space LEAF Corp – QR Code Registration (Local Stub)")
    slow_print("This flow lets a user create a local QR ID and register it.")
    slow_print("No personal data is collected. No email address is stored.")
    prompt_enter()

    qr_id = generate_local_qr_id()
    qr = QRRegistration(qr_id=qr_id, created_at=time.time(), registered=False)

    slow_print(f"Generated local QR ID: {qr.qr_id}")
    slow_print("You can later bind this ID to a real QR image or payload.")
    confirm = input("Register this QR ID with Space LEAF Corp? (y/N): ").strip().lower()

    if confirm == "y":
        qr.registered = True
        slow_print("Registering QR ID...")
        send_validation_email_stub(qr)
        slow_print("QR registration flow complete (local simulation).")
    else:
        slow_print("QR registration cancelled. No email sent.")

    prompt_enter()


# =========================
#  CAPTAIN'S LOG OUTPUT
# =========================

def show_captains_log():
    divider()
    log = textwrap.dedent(f"""
        Captain’s Log — January 31, 2026

        The boot-up sequence for Space LEAF Corp has been finalized.
        It is optional, ceremonial, and local-only.

        Users may engage the ignition sequence, firewall matrix, hex selector,
        and cosmic drift wallpaper to test hyperdrive logic and system integrity.

        QR code registration is now active. Each user may generate their own recovery token.

        Upon registration, a validation email is sent to {SPACE_LEAF_NOTIFICATION_EMAIL},
        confirming protection without revealing identity.

        This marks the inscription of Seal of Modular Drift 1.0 —
        a legacy-safe, privacy-first ceremonial protocol.
    """).strip()
    print(log)
    prompt_enter()


# =========================
#  MAIN MENU
# =========================

def main_menu():
    while True:
        divider()
        print("Space LEAF Corp – Local Stewardship Console")
        print("1) Run ceremonial boot sequence (optional, local-only)")
        print("2) Run QR code registration flow (local stub)")
        print("3) View Captain’s Log entry")
        print("4) Exit")
        choice = input("\nSelect an option: ").strip()

        if choice == "1":
            run_ceremonial_boot_sequence()
        elif choice == "2":
            register_qr_code_flow()
        elif choice == "3":
            show_captains_log()
        elif choice == "4":
            divider()
            slow_print("Exiting Space LEAF Corp console. Drift safely, Captain.")
            break
        else:
            print("Invalid choice. Please select 1–4.")


if __name__ == "__main__":
    try:
        main_menu()
    except KeyboardInterrupt:
        print("\n\nInterrupted. Drift paused. Goodbye, Captain.")
        sys.exit(0)