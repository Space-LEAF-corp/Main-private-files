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
        "CONSENSUS_CLEARED": "🧹🤝🧹",
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
    def __init__(
        self,
        admin_signing_key_b64: Optional[str] = None,
        admin_verify_key_b64: Optional[str] = None,
        validator_verify_keys_b64: Optional[List[str]] = None,
        mfa_secret_b32: Optional[str] = None,
        consensus_threshold: int = 0,
        admin_key: Optional[bytes] = None
    ):
        # For backward compatibility with HMAC-based tamper proofing
        if admin_key is None:
            # Generate a random admin key if not provided
            import secrets
            admin_key = secrets.token_bytes(32)
        self.tamper = Jarvondis3TamperProofing(admin_key)
        self.commands = Jarvondis3GlobalCommands()
        self.log: List[str] = []
        self.emergency_mode: bool = False

        # Provenance
        self.last_hash: str = "GENESIS"

        # MFA & Consensus
        self.mfa_secret_b32: Optional[str] = mfa_secret_b32
        self.consensus_threshold: int = max(0, consensus_threshold)
        self.consensus_signatures: Set[str] = set()
        self.trusted_signers: Set[str] = set()

        # Crypto adapter (Ed25519) - optional
        self.crypto = None
        self.totp = None
        
        if admin_signing_key_b64 or admin_verify_key_b64:
            from .crypto_adapter import Ed25519Adapter
            self.crypto = Ed25519Adapter(
                admin_signing_key_b64=admin_signing_key_b64,
                admin_verify_key_b64=admin_verify_key_b64,
                validator_verify_keys_b64=validator_verify_keys_b64 or []
            )
        
        if mfa_secret_b32:
            from .mfa_totp import TOTP
            self.totp = TOTP(mfa_secret_b32)

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

        # Legal module attachment
        self._legal = None

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
        if self.totp:
            return self.totp.verify(provided_code)
        # Fallback for testing
        return provided_code == self.mfa_secret_b32

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
            if not mfa_code:
                return f"{milestone_banner('MFA_REQUIRED')} MFA code required for sensitive command."
            if not self.verify_mfa(mfa_code):
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
            if not mfa_code:
                return f"{milestone_banner('MFA_REQUIRED')} MFA code required for sensitive command."
            if not self.verify_mfa(mfa_code):
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

        is_sensitive = command in self.sensitive_commands

        if is_sensitive:
            if not admin_signature or not self.tamper.verify_signature(command, admin_signature):
                entry = f"[{timestamp}] {milestone_banner('SIGNATURE_INVALID')} Admin signature required for sensitive command."
                self._append_entry(entry)
                return entry
            if not mfa_code:
                entry = f"[{timestamp}] {milestone_banner('MFA_REQUIRED')} MFA code required for sensitive command."
                self._append_entry(entry)
                return entry
            if not self.verify_mfa(mfa_code):
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

    # --- Legal module wiring ---
    def attach_legal_module(self, legal_module):
        self._legal = legal_module

    @property
    def legal(self):
        if not self._legal:
            raise RuntimeError("LegalProcessModule not attached")
        return self._legal

    def milestone_banner(self, event: str) -> str:
        return milestone_banner(event)


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
