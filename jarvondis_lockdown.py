from typing import Set, Dict, Any

class LockdownPolicy:
    def __init__(
        self,
        owner_id: str,
        admin_secret: str,  # Keep it, but make it optional or required as needed
        allowlist_clients: Set[str],
    ):
        self.owner_id = owner_id
        self.admin_secret = admin_secret
        self.allowlist_clients = allowlist_clients

class AdministrativeLockdown:
    def __init__(self, policy: LockdownPolicy):
        self.policy = policy

    def respond(self, client_id: str, user_agent: str, command: str) -> Dict[str, Any]:
        if client_id not in self.policy.allowlist_clients:
            return {
                "status": "blocked",
                "reason": "privacy_lockdown",
                "message": "Client not on allowlist",
            }

        # Allowlisted → minimal response (expand as needed)
        return {
            "status": "ok",
            "mode": "minimal_lockdown",
            "owner": self.policy.owner_id,
            "command": command,
        }

# === Your fixed usage ===
# Define the secret properly (never hardcode in real code!)
admin_secret = "your-super-secret-admin-key-here"  # Or load from env/var

policy = LockdownPolicy(
    owner_id="leif.w.sogge",
    admin_secret=admin_secret,
    allowlist_clients={"dimitri", "captains-log"}
)

guard = AdministrativeLockdown(policy)

# Test allowlisted client
print(guard.respond("dimitri", "Dimitri/1.0", "status"))
# → {'status': 'ok', 'mode': 'minimal_lockdown', ...}

# Test non-allowlisted client
print(guard.respond("miko", "Miko/1.0", "status"))
# → {'status': 'blocked', 'reason': 'privacy_lockdown', ...}
