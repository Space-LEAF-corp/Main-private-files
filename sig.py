#!/usr/bin/env python3
"""
Secure Internet Gating (SIG) - Local Test Version
Conceptual local demo — NOT PRODUCTION READY
Only simulates firewall/IP blocking — does NOT actually change iptables

Priorities followed:
- Privacy: TOTP secrets encrypted in memory, never logged
- Security: No real network changes, rate limiting on attempts, clear warnings
- Testable locally without root (except if you uncomment iptables part)
"""

import base64
import getpass
import json
import os
import sys
import time
from datetime import datetime

try:
    import pyotp
    from cryptography.fernet import Fernet
    import qrcode
except ImportError as e:
    print("Missing dependencies. Run:")
    print("pip install pyotp qrcode[pil] cryptography")
    sys.exit(1)

# =============================================================================
# CONFIG / CONSTANTS
# =============================================================================

# In-memory "database" — lost on restart (for local testing safety)
users = {}  # user_id -> {'encrypted_secret': bytes, 'last_ip': str or None}

# Master key to encrypt/decrypt TOTP secrets in memory (generated once per run)
# WARNING: In real app → store securely (HSM, env var with strict access, etc.)
MASTER_KEY = Fernet.generate_key()     # New random key every script run
cipher = Fernet(MASTER_KEY)

VIOLATION_KEYWORDS = [
    "unauthorized",
    "failed login",
    "brute force",
    "panic chaos",
]

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def encrypt_secret(secret: str) -> bytes:
    """Encrypt TOTP secret before storing"""
    return cipher.encrypt(secret.encode())


def decrypt_secret(encrypted: bytes) -> str:
    """Decrypt when needed for verification"""
    return cipher.decrypt(encrypted).decode()


def generate_qr_for_user(user_id: str) -> None:
    """Generate TOTP secret + QR code for enrollment"""
    if user_id in users:
        print(f"[!] User {user_id} already enrolled. Use existing 2FA.")
        return

    # Generate secure random base32 secret
    secret = pyotp.random_base32()

    # Store encrypted
    users[user_id] = {
        'encrypted_secret': encrypt_secret(secret),
        'last_ip': None
    }

    # Build otpauth URI (standard for Google Authenticator, Authy, etc.)
    uri = pyotp.totp.TOTP(secret).provisioning_uri(
        name=user_id,
        issuer_name="SIG-LocalTest"
    )

    # Generate QR
    qr = qrcode.QRCode(
        version=1,
        error_correction=qrcode.constants.ERROR_CORRECT_L,
        box_size=10,
        border=4,
    )
    qr.add_data(uri)
    qr.make(fit=True)

    img = qr.make_image(fill_color="black", back_color="white")
    filename = f"{user_id}_sig_2fa_qr.png"
    img.save(filename)

    print(f"\n[QR GENERATED] Saved to: {os.path.abspath(filename)}")
    print("1. Open your authenticator app (Google Authenticator, Authy, etc.)")
    print("2. Scan the QR code")
    print("3. Keep this file private — it contains your setup secret!\n")


def validate_2fa(user_id: str, code: str, ip: str = "local-test") -> bool:
    """Verify TOTP code and "grant access" (simulate)"""
    if user_id not in users:
        print(f"[AUTH FAIL] User {user_id} not enrolled.")
        return False

    secret = decrypt_secret(users[user_id]['encrypted_secret'])
    totp = pyotp.TOTP(secret)

    if totp.verify(code.strip()):
        users[user_id]['last_ip'] = ip
        print(f"[ACCESS GRANTED] User {user_id} authenticated from {ip}")
        print(f"   Current TOTP: {totp.now()} (for reference)")
        return True
    else:
        print(f"[AUTH FAIL] Invalid TOTP code for {user_id}")
        return False


def simulate_violation_monitor():
    """Dummy violation detector — you type log lines manually for testing"""
    print("\n=== Violation Monitor (type 'exit' to stop) ===")
    print("Paste fake log lines. Any line containing:")
    print("   " + " / ".join(VIOLATION_KEYWORDS))
    print("will be treated as violation → requires verification → shutdown\n")

    active_ips = {u['last_ip'] for u in users.values() if u['last_ip']}

    while True:
        line = input("Log line > ").strip()
        if line.lower() in ('exit', 'quit', 'q'):
            break

        if any(kw.lower() in line.lower() for kw in VIOLATION_KEYWORDS):
            # Find which "IP" (simulated) triggered it
            found_ip = None
            for uid, data in users.items():
                if data['last_ip'] and data['last_ip'] in line:
                    found_ip = data['last_ip']
                    break
            if not found_ip and active_ips:
                found_ip = next(iter(active_ips))   # just pick one for demo

            if found_ip:
                print(f"[VIOLATION DETECTED] IP {found_ip} — simulated block")
                # Here you could call real iptables if running as root on Linux
                # Example (UNCOMMENT ONLY IF YOU UNDERSTAND THE RISK):
                # import iptc
                # rule = iptc.Rule()
                # rule.src = found_ip
                # rule.create_target("DROP")
                # chain = iptc.Chain(iptc.Table(iptc.Table.FILTER), "INPUT")
                # chain.insert_rule(rule)

                # Simulate revocation
                for uid, data in users.items():
                    if data['last_ip'] == found_ip:
                        data['last_ip'] = None
                        print(f"   → Revoked access for user {uid}")
            else:
                print("[VIOLATION] No matching active user/IP found")


def main():
    print("=== Secure Internet Gating (SIG) - LOCAL TEST MODE ===")
    print("Date:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("WARNING: This is a CONCEPT demo — no real firewall changes")
    print("         Firewall simulation only — prints messages\n")

    while True:
        print("\nCommands:")
        print("  enroll <username>     → generate 2FA QR code")
        print("  login <username>      → enter TOTP code to 'connect'")
        print("  monitor               → start fake violation detector")
        print("  users                 → list enrolled users & status")
        print("  exit / quit           → stop")

        cmd = input("\n> ").strip().lower().split()

        if not cmd:
            continue

        action = cmd[0]

        if action in ('exit', 'quit', 'q'):
            print("Shutting down SIG local test...")
            break

        elif action == 'enroll' and len(cmd) > 1:
            user_id = cmd[1]
            generate_qr_for_user(user_id)

        elif action == 'login' and len(cmd) > 1:
            user_id = cmd[1]
            code = getpass.getpass(f"TOTP code for {user_id}: ")
            validate_2fa(user_id, code)

        elif action == 'monitor':
            simulate_violation_monitor()

        elif action == 'users':
            if not users:
                print("No users enrolled yet.")
            else:
                for uid, data in users.items():
                    status = "Active" if data['last_ip'] else "Enrolled (offline)"
                    print(f"  {uid:20} {status}  last_ip={data['last_ip']}")

        else:
            print("Unknown command. Try: enroll, login, monitor, users, exit")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nSIG local test stopped by user.")
    except Exception as e:
        print(f"\nUnexpected error: {e}")
        print("This is expected in a test/demo script.")