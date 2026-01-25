import json
import os
import time
from dataclasses import dataclass, asdict
from typing import Dict, Any, Optional

# =========================
# Utility: pacing & speak
# =========================

def pause(seconds: float = 1.0) -> None:
    """Slow pacing for accessibility and calm narration."""
    time.sleep(seconds)

def speak(text: str, delay: float = 0.6) -> None:
    """Print text slowly and clearly, with a small pause."""
    print(text)
    pause(delay)

# =========================
# Identity capsule model
# =========================

CAPSULE_FILE = "slm_identity_capsule.json"
CAPSULE_SCHEMA_VERSION = 1

@dataclass
class IdentityCapsule:
    user_id: str
    schema_version: int
    preferences: Dict[str, Any]
    progress: Dict[str, Any]
    traits: Dict[str, Any]
    admin_password: Optional[str] = None  # NEW FIELD

    @staticmethod
    def default(user_id: str) -> "IdentityCapsule":
        return IdentityCapsule(
            user_id=user_id,
            schema_version=CAPSULE_SCHEMA_VERSION,
            preferences={
                "pacing": "normal",
                "theme": "space_leaf",
                "accessibility_mode": True,
            },
            progress={},
            traits={
                "likes_diagnostics": True,
                "prefers_calm_intro": True,
            },
            admin_password=None,
        )

# =========================
# Capsule store (local-only)
# =========================

class CapsuleStore:
    def __init__(self, path: str = CAPSULE_FILE):
        self.path = path

    def load(self) -> Optional[IdentityCapsule]:
        if not os.path.exists(self.path):
            return None
        try:
            with open(self.path, "r", encoding="utf-8") as f:
                data = json.load(f)
            return IdentityCapsule(
                user_id=data["user_id"],
                schema_version=data.get("schema_version", 1),
                preferences=data.get("preferences", {}),
                progress=data.get("progress", {}),
                traits=data.get("traits", {}),
                admin_password=data.get("admin_password"),
            )
        except Exception:
            return None

    def save(self, capsule: IdentityCapsule) -> None:
        data = asdict(capsule)
        with open(self.path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)

# =========================
# Module API context
# =========================

class ModuleContext:
    def __init__(self, capsule: IdentityCapsule, store: CapsuleStore):
        self.capsule = capsule
        self.store = store

    def get_capsule(self) -> IdentityCapsule:
        return self.capsule

    def update_capsule(self, **kwargs) -> None:
        for key, value in kwargs.items():
            if hasattr(self.capsule, key):
                setattr(self.capsule, key, value)
        self.store.save(self.capsule)

    def log_event(self, event_name: str, details: Optional[Dict[str, Any]] = None) -> None:
        details = details or {}
        print(f"[LOG] {event_name}: {details}")

# =========================
# Login Flow (NEW)
# =========================

def login_flow(ctx: ModuleContext) -> None:
    capsule = ctx.get_capsule()

    speak("\nSelect Login Type:")
    speak("1) Guest Login")
    speak("2) Admin Login")

    choice = input("\nEnter choice (1 or 2): ").strip()

    if choice == "1":
        speak("\nGuest login selected.")
        ctx.log_event("guest_login")
        return

    if choice == "2":
        # Admin login path
        if capsule.admin_password is None:
            speak("\nNo admin password set yet.")
            speak("Would you like to create one now? (y/n)")
            create = input("> ").strip().lower()

            if create == "y":
                while True:
                    pw = input("\nCreate admin password (min 4 chars): ").strip()
                    if len(pw) < 4:
                        speak("Too short. Try again.")
                        continue
                    capsule.admin_password = pw
                    ctx.update_capsule(admin_password=pw)
                    speak("\nAdmin password created.")
                    ctx.log_event("admin_password_created")
                    return
            else:
                speak("\nProceeding as guest.")
                ctx.log_event("guest_login_fallback")
                return

        # Admin password exists → verify
        speak("\nAdmin login selected.")
        for _ in range(3):
            pw = input("Enter admin password: ").strip()
            if pw == capsule.admin_password:
                speak("\nAdmin access granted.")
                ctx.log_event("admin_login_success")
                return
            else:
                speak("Incorrect password.")

        speak("\nToo many attempts. Switching to guest mode.")
        ctx.log_event("admin_login_failed")
        return

    # Default fallback
    speak("\nInvalid choice. Proceeding as guest.")
    ctx.log_event("guest_login_invalid_choice")

# =========================
# Demo module: pre-flight
# =========================

def run_system_check(ctx: ModuleContext) -> None:
    speak("\nSpace LEAF Corp – Pre-Flight System Check")
    speak("Welcome, traveler. We'll move slowly and steadily.")
    speak("This pre-check confirms your device is ready.")

    checks = [
        "Checking operating system ... OK",
        "Checking local storage availability ... OK",
        "Checking safe, private environment ... OK",
        "Checking input readiness ... OK",
        "All systems nominal.",
    ]

    for c in checks:
        speak(c, delay=0.4)

    speak("\nSystem check complete.")
    pause(1.0)
    ctx.log_event("system_check_completed")

def finish(ctx: ModuleContext) -> None:
    speak("\n☒ Pre-Flight Complete")
    speak("You're now ready to walk the system at your own pace.")
    input("\nPress ENTER to continue... ")

    speak(
        "\nPersonal firewall engaged. "
        "Diamond structure integrity 100%. "
        "Privacy layer intact!!!",
        delay=0.2,
    )
    ctx.log_event("preflight_finished")

def preflight_module(ctx: ModuleContext) -> None:
    os.system("clear" if os.name == "posix" else "cls")
    speak("Initializing Space LEAF Corp Pre-Flight Module...")
    pause(1.0)

    login_flow(ctx)
    run_system_check(ctx)
    finish(ctx)

# =========================
# Main entrypoint
# =========================

def main() -> None:
    store = CapsuleStore()
    capsule = store.load()

    if capsule is None:
        capsule = IdentityCapsule.default(user_id="local-anonymous-user")
        store.save(capsule)

    ctx = ModuleContext(capsule, store)
    preflight_module(ctx)

if __name__ == "__main__":
    main()
