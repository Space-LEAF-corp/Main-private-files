# jarvondis_lockdown.py
# Jarvondis — Administrative Access Lockdown
# Author: Leif William Sogge
# Purpose: Enforce sovereign authorship, owner-only control, and privacy lockdown.

from dataclasses import dataclass, field
from typing import Set, Optional, Dict, Any
import hmac
import hashlib
import time

AI_AGENT_FINGERPRINTS = {
    "openai", "anthropic", "azure-ai", "google-ai", "gemini", "bard",
    "perplexity", "cohere", "llama", "claude", "copilot", "stability"
}

def _hmac_sha256(secret: bytes, message: bytes) -> str:
    return hmac.new(secret, message, hashlib.sha256).hexdigest()

@dataclass
class LockdownPolicy:
    owner_id: str
    admin_secret: str  # keep out of VCS if possible; or use env var
    lockdown_active: bool = True
    allowlist_clients: Set[str] = field(default_factory=set)
    audit_enabled: bool = True
    audit_log: list = field(default_factory=list)
    max_skew_seconds: int = 60  # signature freshness window

class AdministrativeLockdown:
    def __init__(self, policy: LockdownPolicy):
        self.policy = policy

    # --- Audit helpers ---
    def _audit(self, event: str, meta: Optional[Dict[str, Any]] = None):
        if not self.policy.audit_enabled:
            return
        self.policy.audit_log.append({
            "ts": int(time.time()),
            "event": event,
            "meta": meta or {}
        })

    # --- Signature verification for owner commands ---
    def verify_owner_command(self, owner_id: str, signature: str, payload: str, timestamp: int) -> bool:
        # freshness check
        now = int(time.time())
        if abs(now - timestamp) > self.policy.max_skew_seconds:
            self._audit("reject_stale_signature", {"owner_id": owner_id, "timestamp": timestamp})
            return False

        if owner_id != self.policy.owner_id:
            self._audit("reject_wrong_owner", {"owner_id": owner_id})
            return False

        message = f"{owner_id}:{timestamp}:{payload}".encode("utf-8")
        expected = _hmac_sha256(self.policy.admin_secret.encode("utf-8"), message)
        ok = hmac.compare_digest(expected, signature)
        self._audit("verify_owner_command", {"ok": ok})
        return ok

    # --- Lockdown state controls (owner-only) ---
    def set_lockdown(self, owner_id: str, signature: str, timestamp: int, active: bool) -> bool:
        payload = f"set_lockdown:{active}"
        if not self.verify_owner_command(owner_id, signature, payload, timestamp):
            return False
        self.policy.lockdown_active = active
        self._audit("lockdown_set", {"active": active})
        return True

    def set_allowlist(self, owner_id: str, signature: str, timestamp: int, clients: Set[str]) -> bool:
        payload = f"set_allowlist:{','.join(sorted(clients))}"
        if not self.verify_owner_command(owner_id, signature, payload, timestamp):
            return False
        self.policy.allowlist_clients = set(clients)
        self._audit("allowlist_updated", {"count": len(clients)})
        return True

    # --- Privacy protection: block other AIs & unauthorized clients ---
    def is_client_allowed(self, client_id: str, user_agent: Optional[str]) -> bool:
        # Explicit allowlist overrides; during lockdown, only allowlisted clients pass
        if self.policy.lockdown_active:
            if client_id not in self.policy.allowlist_clients:
                self._audit("block_client_not_allowlisted", {"client_id": client_id})
                return False

            # Heuristic block: known AI agent fingerprints in the user-agent
            if user_agent:
                ua_l = user_agent.lower()
                if any(fp in ua_l for fp in AI_AGENT_FINGERPRINTS):
                    self._audit("block_ai_agent", {"client_id": client_id, "user_agent": user_agent})
                    return False
        return True

    # --- Jarvondis response gate ---
    def respond(self, client_id: str, user_agent: Optional[str], prompt: str) -> Dict[str, Any]:
        # Gate: only proceed if client allowed
        if not self.is_client_allowed(client_id, user_agent):
            return {
                "status": "blocked",
                "reason": "privacy_lockdown",
                "message": "Jarvondis is under sovereign privacy lockdown. Access denied."
            }

        # Optional: additional rules when lockdown is active
        if self.policy.lockdown_active:
            # Minimal response mode: no external calls, no data exfil.
            self._audit("respond_minimal", {"client_id": client_id})
            return {
                "status": "ok",
                "mode": "minimal_lockdown",
                "reply": f"Jarvondis acknowledges your request: '{prompt}'. Privacy lockdown is active."
            }

        # Normal mode response placeholder
        self._audit("respond_normal", {"client_id": client_id})
        return {
            "status": "ok",
            "mode": "normal",
            "reply": f"Jarvondis responding: '{prompt}'."
        }
