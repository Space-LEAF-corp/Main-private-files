# main.py
import json
import os
from jarvondis_lockdown import LockdownPolicy, AdministrativeLockdown

def load_policy(path: str) -> LockdownPolicy:
    try:
        with open(path, "r") as f:
            cfg = json.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"Policy configuration file not found: {path}")
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in policy file {path}: {e}")
    
    # Validate required fields
    if "owner_id" not in cfg:
        raise ValueError("Policy configuration missing required field: owner_id")
    
    admin_secret = os.environ.get("JARVONDIS_ADMIN_SECRET")
    if not admin_secret:
        raise ValueError("JARVONDIS_ADMIN_SECRET environment variable not set")
    
    return LockdownPolicy(
        owner_id=cfg["owner_id"],
        admin_secret=admin_secret,
        lockdown_active=cfg.get("lockdown_active", True),
        allowlist_clients=set(cfg.get("allowlist_clients", [])),
        audit_enabled=cfg.get("audit_enabled", True),
        max_skew_seconds=cfg.get("max_skew_seconds", 60)
    )

if __name__ == "__main__":
    policy = load_policy("jarvondis_policy.json")
    guard = AdministrativeLockdown(policy)

    # Example request
    client_id = "captains-log"
    user_agent = "Jarvondis/1.0 (ceremonial-client)"
    prompt = "Status check."

    result = guard.respond(client_id, user_agent, prompt)
    print(result)
