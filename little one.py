jarvondis3/
  core/
    sovereign.py
  gaming/
    engine.py
    ai_agents.py
    rituals.py
    assets/
      README.md
  compliance/
    exports/
      README.md
  docs/
    JARVONDIS3_MANIFESTO.md
    POSITION_PAPER.md
  README.md
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Jarvondis 3.0 — Unified Sovereign Control Plane (Complete Finalized Code)

Features:
- Tamper-proofing (file hashes, HMAC signatures as placeholder, RBAC)
- Global command surface (view vs. execute with root authority)
- Global sovereign interface (single gateway)
- Administrator + emergency access lifecycle (dormant until granted)
- Ceremonial milestone banners and living audit log
- Hash-chained provenance for tamper-evident logs
- Compliance snapshot export (JSON)
- Agent harness (RBAC, provenance, resource quotas)
- Multi-factor sovereignty (MFA for sensitive commands)
- Consensus mode (trusted co-signers) with administrator fail-safe (admin signature required)
"""

import hashlib
import hmac
import json
from datetime import datetime
from typing import Dict, List, Optional, Set


# === Ceremonial Milestone Banners ===
def milestone_banner(event: str) -> str:
    glyphs = {
        "EMERGENCY_ACTIVATED": "⚡⚡⚡",
        "EMERGENCY_REVOKED": "🛑🔒🛑",
        "COMMAND_EXECUTED": "📜✨📜",
        "ACCESS_DENIED": "🚫🕯️🚫",
        "SIGNATURE_INVALID": "❌🔑❌",
        "VIEW_ONLY": "👁️🧭👁️",
        "UNKNOWN_COMMAND": "❓📟❓",
        "AUDIT_HEADER": "🏛️",
        "DIAGNOSTICS_COMPLETE": "🔍🛡️🔍",
        "PROTOCOL_OVERRIDE": "🌀⛓️🌀",
        "BACKUP_COMPLETED": "🧱🔐🧱",
        "PROVENANCE_CHECKPOINT": "⛓️📜⛓️",
        "MFA_REQUIRED": "🔒⏳🔒",
        "CONSENSUS_REQUIRED": "🤝⛓️🤝",
        "CONSENSUS_ADDED": "🪙🤝🪙",
        "CONSENSUS_CLEARED": "🧹🤝🧹"
    }
    return f"{glyphs.get(event, '🔔')} [MILESTONE] {event} 🔔"


# === Tamper-Proofing Layer ===
class Jarvondis3TamperProofing:
    def __init__(self, admin_key: bytes):
        self.admin_key: bytes = admin_key
        self.file_hashes: Dict[str, str] = {}

    def _hash_file(self, filepath: str) -> str:
        hasher = hashlib.sha256()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                hasher.update(chunk)
        return hasher.hexdigest()

    def register_file(self, filepath: str):
        self.file_hashes[filepath] = self._hash_file(filepath)

    def verify_file(self, filepath: str) -> bool:
        baseline = self.file_hashes.get(filepath)
        return baseline is not None and self._hash_file(filepath) == baseline

    def sign_command(self, command: str) -> str:
        # Placeholder for asymmetric signatures; HMAC ensures integrity with admin_key
        return hmac.new(self.admin_key, command.encode(), hashlib.sha256).hexdigest()

    def verify_signature(self, command: str, signature: str) -> bool:
        expected = self.sign_command(command)
        return hmac.compare_digest(expected, signature)

    def check_access(self, user_role: str, required_role: str) -> bool:
        hierarchy = ["user", "operator", "admin"]
        try:
            return hierarchy.index(user_role) >= hierarchy.index(required_role)
        except ValueError:
            return False


# === Global Command Surface ===
class Jarvondis3GlobalCommands:
    def __init__(self):
        self.commands: Dict[str, str] = {
            "/JARVONDIS_SYS_REBOOT": "System-wide reboot of Jarvondis 3.0",
            "/JARVONDIS_SYS_SHUTDOWN": "Graceful shutdown of all modules",
            "/JARVONDIS_SYS_UPDATE": "Initiate update sequence with changelog continuity",
            "/JARVONDIS_SYS_BACKUP": "Trigger immediate encrypted backup",
            "/JARVONDIS_SYS_ACCESS_GRANT": "Grant scoped access to subsystems or commands",
            "/JARVONDIS_SYS_ACCESS_REVOKE": "Revoke scoped access to subsystems or commands",
            "/JARVONDIS_SYS_DIAGNOSTICS": "Run comprehensive system diagnostics with trust banners",
            "/JARVONDIS_SYS_PROTOCOL_OVERRIDE": "Temporarily override protocols (requires ceremonial approval)"
        }
        self.log: List[str] = []

    def list_commands(self) -> Dict[str, str]:
        return self.commands

    def execute_command(self, command: str, authorized: bool = False) -> str:
        timestamp = datetime.utcnow().isoformat()

        if not authorized:
            entry = f"[{timestamp}] {milestone_banner('VIEW_ONLY')} {command} restricted. Root authority required."
            self.log.append(entry)
            return entry

        if command not in self.commands:
            entry = f"[{timestamp}] {milestone_banner('UNKNOWN_COMMAND')} Unknown command: {command}"
            self.log.append(entry)
            return entry

        description = self.commands[command]
        entry = f"[{timestamp}] {milestone_banner('COMMAND_EXECUTED')} {command}: {description}"
        self.log.append(entry)

        # Ritualized hooks
        if command == "/JARVONDIS_SYS_DIAGNOSTICS":
            self.log.append(f"[{timestamp}] {milestone_banner('DIAGNOSTICS_COMPLETE')} Diagnostics completed with trust banners.")
        elif command == "/JARVONDIS_SYS_BACKUP":
            self.log.append(f"[{timestamp}] {milestone_banner('BACKUP_COMPLETED')} Encrypted backup ritual completed.")
        elif command == "/JARVONDIS_SYS_PROTOCOL_OVERRIDE":
            self.log.append(f"[{timestamp}] {milestone_banner('PROTOCOL_OVERRIDE')} Protocol override invoked (ceremonial approval logged).")

        return entry

    def view_log(self) -> str:
        return "\n".join(self.log)


# === Global Sovereign Interface with Provenance, MFA, Consensus & Emergency Access ===
class Jarvondis3GlobalInterface:
    def __init__(self, admin_key: bytes, mfa_secret: str, consensus_threshold: int = 0):
        self.tamper = Jarvondis3TamperProofing(admin_key)
        self.commands = Jarvondis3GlobalCommands()
        self.log: List[str] = []
        self.emergency_mode: bool = False

        # Provenance
        self.last_hash: str = "GENESIS"

        # MFA & Consensus
        self.mfa_secret: str = mfa_secret
        self.consensus_threshold: int = max(0, consensus_threshold)
        self.consensus_signatures: Set[str] = set()
        self.trusted_signers: Set[str] = set()  # IDs of allowed co-signers (policy-as-code hook)

        # Sensitive commands requiring MFA + admin signature (+ optional consensus)
        self.sensitive_commands: Set[str] = {
            "/JARVONDIS_SYS_SHUTDOWN",
            "/JARVONDIS_SYS_REBOOT",
            "/JARVONDIS_SYS_PROTOCOL_OVERRIDE",
            "/JARVONDIS_SYS_ACCESS_GRANT",
            "/JARVONDIS_SYS_ACCESS_REVOKE",
            "EMERGENCY_GRANT",
            "EMERGENCY_REVOKE",
        }

    # --- Provenance utilities ---
    def _hash_entry(self, entry: str) -> str:
        return hashlib.sha256((self.last_hash + entry).encode()).hexdigest()

    def _append_entry(self, entry: str):
        new_hash = self._hash_entry(entry)
        provenance = f"{entry} | PROVENANCE_HASH={new_hash}"
        self.log.append(provenance)
        self.last_hash = new_hash
        self.log.append(milestone_banner("PROVENANCE_CHECKPOINT"))

    # --- MFA utilities (placeholder for TOTP/HOTP) ---
    def verify_mfa(self, provided_code: str) -> bool:
        return provided_code == self.mfa_secret

    # --- Consensus management ---
    def add_trusted_signer(self, signer_id: str):
        self.trusted_signers.add(signer_id)
        self._append_entry(f"{milestone_banner('CONSENSUS_ADDED')} Trusted signer added: {signer_id}")

    def clear_consensus(self):
        self.consensus_signatures.clear()
        self._append_entry(milestone_banner("CONSENSUS_CLEARED"))

    def provide_consensus(self, signer_id: str, signature: str) -> str:
        if signer_id not in self.trusted_signers:
            msg = f"{milestone_banner('ACCESS_DENIED')} Unrecognized signer: {signer_id}"
            self._append_entry(msg)
            return msg
        # Consensus signature is over literal "CONSENSUS"
        if self.tamper.verify_signature("CONSENSUS", signature):
            self.consensus_signatures.add(signer_id)
            msg = f"{milestone_banner('CONSENSUS_ADDED')} Consensus accepted from {signer_id}"
            self._append_entry(msg)
            return msg
        msg = f"{milestone_banner('SIGNATURE_INVALID')} Invalid consensus signature from {signer_id}"
        self._append_entry(msg)
        return msg

    def _consensus_met(self) -> bool:
        return len(self.consensus_signatures) >= self.consensus_threshold

    # --- Emergency lifecycle (admin signature REQUIRED always) ---
    def grant_emergency_access(self, admin_signature: str, mfa_code: Optional[str] = None) -> str:
        command_literal = "EMERGENCY_GRANT"
        if not self.tamper.verify_signature(command_literal, admin_signature):
            return f"{milestone_banner('SIGNATURE_INVALID')} Invalid administrator signature."

        if command_literal in self.sensitive_commands:
            if not self.verify_mfa(mfa_code or ""):
                return f"{milestone_banner('MFA_REQUIRED')} MFA verification failed."
            if self.consensus_threshold > 0 and not self._consensus_met():
                return f"{milestone_banner('CONSENSUS_REQUIRED')} Consensus threshold not met."

        self.emergency_mode = True
        entry = f"{milestone_banner('EMERGENCY_ACTIVATED')} Emergency access granted by administrator."
        self._append_entry(entry)
        return entry

    def revoke_emergency_access(self, admin_signature: str, mfa_code: Optional[str] = None) -> str:
        command_literal = "EMERGENCY_REVOKE"
        if not self.tamper.verify_signature(command_literal, admin_signature):
            return f"{milestone_banner('SIGNATURE_INVALID')} Invalid administrator signature."

        if command_literal in self.sensitive_commands:
            if not self.verify_mfa(mfa_code or ""):
                return f"{milestone_banner('MFA_REQUIRED')} MFA verification failed."
            if self.consensus_threshold > 0 and not self._consensus_met():
                return f"{milestone_banner('CONSENSUS_REQUIRED')} Consensus threshold not met."

        self.emergency_mode = False
        entry = f"{milestone_banner('EMERGENCY_REVOKED')} Emergency access revoked by administrator."
        self._append_entry(entry)
        return entry

    # --- Core request path (administrator supremacy enforced) ---
    def request(self, user_role: str, command: str, signature: Optional[str] = None,
                admin_signature: Optional[str] = None, mfa_code: Optional[str] = None) -> str:
        """
        Unified request method:
        - Admin requests: must satisfy sensitive command policy (admin signature + MFA + optional consensus).
        - Operator requests: allowed only in emergency_mode, for non-sensitive commands unless admin_signed.
        - All sensitive commands ALWAYS require administrator signature (fail-safe).
        """
        timestamp = datetime.utcnow().isoformat()

        if command not in self.commands.commands:
            entry = f"[{timestamp}] {milestone_banner('UNKNOWN_COMMAND')} Unknown command: {command}"
            self._append_entry(entry)
            return entry

        is_sensitive = command in self.sensitive_commands or command in {
            "/JARVONDIS_SYS_SHUTDOWN",
            "/JARVONDIS_SYS_REBOOT",
            "/JARVONDIS_SYS_PROTOCOL_OVERRIDE",
            "/JARVONDIS_SYS_ACCESS_GRANT",
            "/JARVONDIS_SYS_ACCESS_REVOKE",
        }

        if is_sensitive:
            if not admin_signature or not self.tamper.verify_signature(command, admin_signature):
                entry = f"[{timestamp}] {milestone_banner('SIGNATURE_INVALID')} Admin signature required for sensitive command."
                self._append_entry(entry)
                return entry
            if not self.verify_mfa(mfa_code or ""):
                entry = f"[{timestamp}] {milestone_banner('MFA_REQUIRED')} MFA verification failed."
                self._append_entry(entry)
                return entry
            if self.consensus_threshold > 0 and not self._consensus_met():
                entry = f"[{timestamp}] {milestone_banner('CONSENSUS_REQUIRED')} Consensus threshold not met."
                self._append_entry(entry)
                return entry

        if not self.tamper.check_access(user_role, "admin"):
            if not (self.emergency_mode and user_role == "operator" and not is_sensitive):
                entry = f"[{timestamp}] {milestone_banner('ACCESS_DENIED')} ACCESS DENIED for {user_role} on {command}"
                self._append_entry(entry)
                return entry

        if signature and not self.tamper.verify_signature(command, signature):
            entry = f"[{timestamp}] {milestone_banner('SIGNATURE_INVALID')} INVALID SIGNATURE for {command}"
            self._append_entry(entry)
            return entry

        result = self.commands.execute_command(command, authorized=True)
        self._append_entry(result)
        return result

    # --- Audit & Compliance ---
    def view_audit_log(self) -> str:
        ceremonial_header = f"{milestone_banner('AUDIT_HEADER')} === Jarvondis 3.0 Global Sovereign Log ==="
        return "\n".join([ceremonial_header] + self.log + [self.commands.view_log()])

    def export_compliance_snapshot(self) -> str:
        snapshot = {
            "header": "Jarvondis 3.0 Compliance Snapshot",
            "timestamp": datetime.utcnow().isoformat(),
            "audit_log": self.log,
            "last_provenance_hash": self.last_hash,
            "consensus_threshold": self.consensus_threshold,
            "current_consensus_signers": sorted(list(self.consensus_signatures)),
            "trusted_signers": sorted(list(self.trusted_signers)),
            "emergency_mode": self.emergency_mode,
        }
        return json.dumps(snapshot, indent=2)


# === Agent Harness (constrained execution wrapper) ===
class Jarvondis3AgentHarness:
    """
    Constrains external agents to act through Jarvondis's sovereign interface.
    Enforces RBAC, provenance logging, and optional resource quotas (placeholder).
    """
    def __init__(self, interface: Jarvondis3GlobalInterface, agent_id: str, role: str = "operator"):
        self.interface = interface
        self.agent_id = agent_id
        self.role = role
        # Placeholder quota fields (wire up to real resource managers if needed)
        self.max_commands_per_minute = 60
        self.commands_executed = 0

    def execute(self, command: str, signature: Optional[str] = None) -> str:
        # Agents cannot provide admin_signature or MFA; they are always constrained.
        return self.interface.request(
            user_role=self.role,
            command=command,
            signature=signature
        )
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from core.sovereign import Jarvondis3GlobalInterface, milestone_banner

class JarvondisGameEngine:
    """
    Minimal sovereign-governed game engine scaffold.
    All lifecycle events route through Jarvondis for provenance and ceremony.
    """
    def __init__(self, sovereign: Jarvondis3GlobalInterface):
        self.sovereign = sovereign
        self.state = {"running": False, "tick": 0}

    def start(self):
        self.state["running"] = True
        banner = milestone_banner("COMMAND_EXECUTED")
        return f"{banner} Game engine started under sovereign control."

    def init_diagnostics(self):
        # Demonstrates operator using emergency lane for non-sensitive diagnostics (if active)
        return self.sovereign.request(
            user_role="operator",
            command="/JARVONDIS_SYS_DIAGNOSTICS"
        )

    def tick(self):
        if not self.state["running"]:
            return "Engine not running."
        self.state["tick"] += 1
        return f"Game tick {self.state['tick']} executed."

    def shutdown(self, admin_signature: str, mfa_code: str):
        # Sensitive: requires admin signature + MFA (and consensus if configured)
        return self.sovereign.request(
            user_role="admin",
            command="/JARVONDIS_SYS_SHUTDOWN",
            admin_signature=admin_signature,
            mfa_code=mfa_code
        )
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from core.sovereign import Jarvondis3AgentHarness

class GameAgent:
    """
    Example game agent constrained by Jarvondis harness.
    """
    def __init__(self, harness: Jarvondis3AgentHarness, name: str):
        self.harness = harness
        self.name = name

    def act(self, command: str):
        """
        Agent attempts an action via sovereign interface.
        Non-sensitive commands permitted under emergency-mode for operators;
        sensitive commands require admin signature (which agents cannot provide).
        """
        return self.harness.execute(command)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from core.sovereign import milestone_banner

def level_up_ritual(level: int) -> str:
    return f"{milestone_banner('COMMAND_EXECUTED')} Player reached level {level}!"

def boss_defeated_ritual(boss_name: str) -> str:
    return f"{milestone_banner('DIAGNOSTICS_COMPLETE')} Boss {boss_name} defeated under sovereign watch."


# Gaming assets

Place symbolic glyphs, art, sound, and other resources here.
Tie asset usage to sovereign milestones to preserve the ceremonial layer.
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Crypto adapter for Jarvondis 3.0
- Ed25519 signatures for admin + validators
- HMAC fallback (optional)
- Canonical JSON helper (stable serialization)
"""

import json
from typing import Optional
from nacl.signing import SigningKey, VerifyKey
from nacl.exceptions import BadSignatureError


def canonical_json(data) -> str:
    return json.dumps(data, sort_keys=True, separators=(",", ":"))


class Ed25519Adapter:
    """
    Admin and validators use Ed25519 keys.
    - Admin: one signing key for root authority
    - Validators: list of verify keys (public) or signing keys if you also sign locally
    """
    def __init__(
        self,
        admin_signing_key_b64: Optional[str] = None,
        admin_verify_key_b64: Optional[str] = None,
        validator_verify_keys_b64: Optional[list[str]] = None
    ):
        # Admin may be loaded by private signing key OR public verify key (verify-only mode)
        self.admin_signing_key = SigningKey(admin_signing_key_b64) if admin_signing_key_b64 else None
        self.admin_verify_key = (
            VerifyKey(admin_verify_key_b64) if admin_verify_key_b64 else
            (self.admin_signing_key.verify_key if self.admin_signing_key else None)
        )
        self.validator_verify_keys = [
            VerifyKey(v) for v in (validator_verify_keys_b64 or [])
        ]

    def admin_sign(self, message: str) -> bytes:
        if not self.admin_signing_key:
            raise RuntimeError("Admin signing key not loaded")
        return self.admin_signing_key.sign(message.encode("utf-8")).signature

    def admin_verify(self, message: str, signature: bytes) -> bool:
        if not self.admin_verify_key:
            raise RuntimeError("Admin verify key not loaded")
        try:
            self.admin_verify_key.verify(message.encode("utf-8"), signature)
            return True
        except BadSignatureError:
           
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TOTP MFA for Jarvondis 3.0
- RFC 6238 time-based one-time passwords
- Uses pyotp-like semantics but implemented minimally with hmac/sha1
"""

import base64
import hmac
import hashlib
import time
from typing import Optional


def _int_to_bytes(i: int) -> bytes:
    return i.to_bytes(8, "big")


class TOTP:
    """
    Minimal TOTP verifier.
    - secret_b32: base32-encoded secret (e.g., from authenticator app)
    - interval: time step in seconds (default 30)
    - digits: number of digits (default 6)
    """
    def __init__(self, secret_b32: str, interval: int = 30, digits: int = 6):
        self.secret = base64.b32decode(secret_b32.upper())
        self.interval = interval
        self.digits = digits

    def _code_at(self, for_time: int) -> str:
        counter = for_time // self.interval
        msg = _int_to_bytes(counter)
        digest = hmac.new(self.secret, msg, hashlib.sha1).digest()
        offset = digest[-1] & 0x0F
        binary = ((digest[offset] & 0x7f) << 24) | (digest[offset + 1] << 16) | (digest[offset + 2] << 8) | (digest[offset + 3])
        otp = binary % (10 ** self.digits)
        return str(otp).zfill(self.digits)

    def verify(self, code: str, valid_window: int = 1, now: Optional[int] = None) -> bool:
        now = now or int(time.time())
        for delta in range(-valid_window, valid_window + 1):
            if self._code_at(now + delta * self.interval) == code:
                return True
        return False
# ... existing imports ...
from typing import List, Optional
from core.crypto_adapter import Ed25519Adapter, canonical_json
from core.mfa_totp import TOTP

def milestone_banner(event: str) -> str:
    glyphs = {
        "LEGAL_REQUEST_SUBMITTED": "⚖️",
        "LEGAL_REQUEST_INVALID": "⚠️",
        "LEGAL_SIGNATURE_INVALID": "❌",
        "LEGAL_MFA_FAILED": "🔐",
        "LEGAL_CONSENSUS_REQUIRED": "👥",
        "LEGAL_VALIDATOR_INVALID": "🧪",
        "LEGAL_DECISION": "📜",
        "LEGAL_EXECUTED": "🏛️",
        "LEGAL_SIMULATED": "🧪",
        "LEGAL_CANARY_PUBLISHED": "🕊️",
        "AUDIT_HEADER": "🏛️",
        # ... others as needed ...
    }
    return f"{glyphs.get(event, '🔔')} [MILESTONE] {event} 🔔"


class Jarvondis3GlobalInterface:
    def __init__(
        self,
        admin_signing_key_b64: Optional[str],
        admin_verify_key_b64: Optional[str],
        validator_verify_keys_b64: Optional[List[str]],
        mfa_secret_b32: str,
        consensus_threshold: int = 0
    ):
        # Crypto and MFA
        self.crypto = Ed25519Adapter(
            admin_signing_key_b64=admin_signing_key_b64,
            admin_verify_key_b64=admin_verify_key_b64,
            validator_verify_keys_b64=validator_verify_keys_b64 or []
        )
        self.totp = TOTP(mfa_secret_b32)
        self.consensus_threshold = consensus_threshold

        # Logs and module attachments
        self._log: List[str] = []
        self._legal = None  # attached later

    # ----- MFA -----
    def verify_mfa(self, mfa_code: str) -> bool:
        return self.totp.verify(mfa_code, valid_window=1)

    # ----- Logging -----
    def _append_entry(self, entry: str):
        from datetime import datetime
        ts = datetime.utcnow().isoformat()
        self._log.append(f"[{ts}] {entry}")

    def export_compliance_snapshot(self) -> str:
        return "\n".join([milestone_banner("AUDIT_HEADER")] + self._log)

    # ----- Legal module wiring -----
    def attach_legal_module(self, legal_module):
        self._legal = legal_module

    @property
    def legal(self):
        if not self._legal:
            raise RuntimeError("LegalProcessModule not attached")
        return self._legal

    def milestone_banner(self, event: str) -> str:
        return milestone_banner(event)
# ... existing imports ...
from core.crypto_adapter import canonical_json

# In LegalProcessModule.approve_request(), replace signature checks:

bound = f"/JARVONDIS_SYS_LEGAL_REQUEST:{req.action}:{req.canonical()}"
# Admin signature (Ed25519)
if not self.sovereign.crypto.admin_verify(bound, admin_signature):
    return self._append_event("LEGAL_SIGNATURE_INVALID", "Admin signature invalid", bound)

# MFA check remains:
if not self.sovereign.verify_mfa(mfa_code):
    return self._append_event("LEGAL_MFA_FAILED", "MFA verification failed", request_id)

# Validator threshold (Ed25519)
if self.sovereign.consensus_threshold > 0:
    validator
