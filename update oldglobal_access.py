"""
global_access.py - Enhanced version with Space LEAF Corp gatekeeper + emergency path

- Read:         default for authenticated users
- Modify:       editor/admin roles + confirmation
- Emergency:    ONLY space_leaf_admin role + separate emergency passphrase
- All critical actions heavily logged
"""

import hashlib
import json
import logging
import os
import secrets
import time
from getpass import getpass
from typing import Optional, Dict, Any

# ────────────────────────────────────────────────
# CONFIGURATION ───────────────────────────────────
# In real deployment: move secrets to env / vault / HSM
# ────────────────────────────────────────────────

LOG_FILE = "global_access_audit.log"
TOKEN_VALIDITY_SECONDS = 3600 * 24 * 7          # 7 days
SECRET_SALT = secrets.token_hex(32)             # Rotate in prod!

# Emergency passphrase hash (change immediately in real use!)
EMERGENCY_PASSPHRASE_HASH = hashlib.sha256(
    ("LEAF-OMEGA-EMERGENCY-2026-GATEKEEPER-!!!" + SECRET_SALT).encode()
).hexdigest()

# Simulated user database (passwords MUST be hashed!)
USERS: Dict[str, Dict[str, Any]] = {
    "readonly": {
        "password_hash": hashlib.sha256(("correct-horse-battery-staple" + SECRET_SALT).encode()).hexdigest(),
        "roles": ["readonly"],
    },
    "editor": {
        "password_hash": hashlib.sha256(("editor-very-secret-2026" + SECRET_SALT).encode()).hexdigest(),
        "roles": ["readonly", "editor"],
    },
    "admin": {
        "password_hash": hashlib.sha256(("change-me-immediately-!!!" + SECRET_SALT).encode()).hexdigest(),
        "roles": ["readonly", "editor", "admin"],
    },
    "space_leaf_corp": {
        "password_hash": hashlib.sha256(("space-leaf-alpha-secure-2026" + SECRET_SALT).encode()).hexdigest(),
        "roles": ["readonly", "editor", "admin", "space_leaf_admin"],
    }
}

# In-memory token store (use Redis / DB in production)
active_tokens: Dict[str, Dict] = {}

# ────────────────────────────────────────────────
# LOGGING SETUP
# ────────────────────────────────────────────────
logging.basicConfig(
    filename=LOG_FILE,
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | user=%(user)-16s | roles=%(roles)s | action=%(action)-20s | emerg=%(emergency)s | ip=%(ip)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

def log_event(level: str, message: str, extra: Dict = None):
    extra = extra or {}
    extra.setdefault("ip", "127.0.0.1")           # Replace with real IP in web context
    getattr(logging, level.lower())(message, extra=extra)


# ────────────────────────────────────────────────
# CORE AUTH FUNCTIONS
# ────────────────────────────────────────────────

def hash_password(password: str) -> str:
    return hashlib.sha256((password + SECRET_SALT).encode()).hexdigest()


def authenticate_user(username: str, password: str) -> Optional[str]:
    user = USERS.get(username.lower())
    if not user:
        return None
    if user["password_hash"] != hash_password(password):
        return None
    return username.lower()


def get_user_roles(username: str) -> list:
    return USERS.get(username.lower(), {}).get("roles", ["readonly"])


def generate_access_token(username: str, roles: list) -> Dict:
    token_value = secrets.token_urlsafe(48)
    expires = int(time.time()) + TOKEN_VALIDITY_SECONDS

    token_data = {
        "token": token_value,
        "user": username,
        "roles": roles,
        "expires": expires,
        "issued": int(time.time())
    }

    active_tokens[token_value] = token_data
    return token_data


def validate_token(token_str: str) -> Optional[Dict]:
    token_data = active_tokens.get(token_str)
    if not token_data:
        return None
    if time.time() > token_data["expires"]:
        del active_tokens[token_str]
        return None
    return token_data


def verify_emergency_passphrase(phrase: str) -> bool:
    return hash_password(phrase) == EMERGENCY_PASSPHRASE_HASH


# ────────────────────────────────────────────────
# PERMISSION CHECKS
# ────────────────────────────────────────────────

def can_read(token_data: Dict) -> bool:
    return bool(token_data and time.time() < token_data.get("expires", 0))


def can_modify(token_data: Dict, require_confirmation: bool = True) -> bool:
    if not can_read(token_data):
        return False

    roles = token_data.get("roles", [])
    if not ("editor" in roles or "admin" in roles):
        return False

    if require_confirmation:
        print("\n" + "="*60)
        print("*** STANDARD MODIFY ACCESS REQUEST ***")
        print("This action can affect production data/systems.")
        confirm = getpass("Enter 'CONFIRM-MODIFY' to proceed: ").strip()
        if confirm != "CONFIRM-MODIFY":
            log_event("warning", "Standard modify rejected — confirmation failed",
                      {"user": token_data["user"], "roles": ",".join(roles),
                       "action": "modify_deny_confirm", "emergency": False})
            return False

    return True


def can_emergency_modify(token_data: Dict) -> bool:
    if not can_read(token_data):
        return False

    roles = token_data.get("roles", [])
    if "space_leaf_admin" not in roles:
        print("Emergency modify restricted to Space LEAF Corp gatekeeper only.")
        return False

    print("\n" + "="*70)
    print("*** EMERGENCY MODIFY PROTOCOL — SPACE LEAF CORP GATEKEEPER ONLY ***")
    print("Bypasses normal controls. Use ONLY in genuine emergencies.")
    print("All attempts are permanently logged.")

    phrase = getpass("Emergency passphrase: ")
    if not verify_emergency_passphrase(phrase):
        log_event("critical", "EMERGENCY PASSPHRASE FAILURE — possible breach attempt",
                  {"user": token_data["user"], "roles": ",".join(roles),
                   "action": "emergency_auth_fail", "emergency": True})
        print("Emergency authentication FAILED. Action denied. Incident logged.")
        return False

    confirm = getpass("Final confirmation — type 'EMERGENCY-LEAF-AUTHORIZED': ").strip()
    if confirm != "EMERGENCY-LEAF-AUTHORIZED":
        log_event("warning", "Emergency modify aborted — final confirmation missing",
                  {"user": token_data["user"], "roles": ",".join(roles),
                   "action": "emergency_deny_confirm", "emergency": True})
        print("Emergency action aborted.")
        return False

    log_event("critical", "EMERGENCY MODIFY ACCESS GRANTED — SPACE LEAF CORP",
              {"user": token_data["user"], "roles": ",".join(roles),
               "action": "emergency_granted", "emergency": True})

    print("\n*** EMERGENCY MODIFY AUTHORIZED BY SPACE LEAF CORP GATEKEEPER ***")
    return True


# ────────────────────────────────────────────────
# MAIN INTERACTIVE LOOP
# ────────────────────────────────────────────────

def main():
    print("═" * 70)
    print("GLOBAL ACCESS CONTROL — 2026 Edition")
    print("Read-Only Default • Space LEAF Corp Gatekeeper Enforced")
    print("Emergency protocol restricted to Space LEAF Corp ONLY")
    print("═" * 70, "\n")

    # ── Login ─────────────────────────────────────
    username = input("Username: ").strip().lower()
    password = getpass("Password: ")

    auth_user = authenticate_user(username, password)
    if not auth_user:
        print("\nAuthentication failed.")
        log_event("error", "Login failed", {"user": username, "action": "login_fail", "emergency": False})
        return

    roles = get_user_roles(auth_user)
    print(f"\nAuthenticated as {auth_user}  |  Roles: {', '.join(roles)}")

    token_data = generate_access_token(auth_user, roles)

    print("\nAccess Token (keep secure!):")
    print(json.dumps(token_data, indent=2))
    print(f"Valid until: {time.ctime(token_data['expires'])}\n")

    # ── Token-based command loop ──────────────────
    while True:
        print("─" * 50)
        print("Available actions:  read  |  modify  |  emergency-modify  |  quit")
        action = input("Action > ").strip().lower()

        if action in ("q", "quit", "exit"):
            print("Session ended.")
            break

        if action == "read":
            if can_read(token_data):
                print("Read access granted (default for authenticated users).")
                log_event("info", "Read access granted",
                          {"user": auth_user, "roles": ",".join(roles), "action": "read_grant", "emergency": False})
            else:
                print("Token invalid or expired.")

        elif action == "modify":
            if can_modify(token_data):
                print("Standard modify access GRANTED.")
                log_event("info", "Standard modify granted",
                          {"user": auth_user, "roles": ",".join(roles), "action": "modify_grant", "emergency": False})
                # ── Example placeholder for actual modify logic ──
                print("  → Would now execute protected write operation...")
                print("  → Example: updating config / database / file / firewall rule")
            else:
                print("Standard modify DENIED.")

        elif action == "emergency-modify":
            if can_emergency_modify(token_data):
                print("EMERGENCY MODIFY ACCESS GRANTED — critical audit logged.")
                # ── Example placeholder for emergency logic ──
                print("  → Executing emergency override / backdoor / kill-switch removal...")
                print("  → All actions are being audited under emergency protocol.")
            else:
                print("Emergency modify DENIED.")

        else:
            print("Unknown action. Try: read / modify / emergency-modify / quit")

    # Cleanup (in real app you'd do this on logout / timeout)
    if token_data["token"] in active_tokens:
        del active_tokens[token_data["token"]]


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nSession terminated by user.")
    except Exception as e:
        log_event("error", f"Unexpected error: {str(e)}", {"action": "runtime_error", "emergency": False})
        print("An unexpected error occurred. Check audit log.")
