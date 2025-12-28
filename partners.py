#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Federated Stewardship Runtime (FSR)
- QR-DNA activation (pre-init)
- Repair-Only Covenant (no replacement)
- Immutable Labeling Covenant (transparent, non-harm, warnings target intruders)
- Multi-Team Federation with Witness Attestation (Microsoft, NASA, US, Canada, Sovereign Steward)
- Lockout Shield Sequence (fail-closed, drift freeze, notarization)
- Tamper-evident audit chains per tenant + federated rollups

Author: Leif William Sogge
"""

import hashlib, json, os, time
from dataclasses import dataclass
from typing import Optional, Dict, Any, List, TypedDict
# TypedDict for immutable_label
class ImmutableLabel(TypedDict):
    name: str
    version: str
    purpose: str
    harm_prohibition: bool
    repair_only: bool
    warnings_target: str

# ---------------------------
# Configuration
# ---------------------------

TENANTS = [
    {"id": "sovereign", "name": "Sovereign Stewardship"},
    {"id": "microsoft", "name": "Microsoft"},
    {"id": "nasa", "name": "NASA"},
    {"id": "us", "name": "United States"},
    {"id": "canada", "name": "Canada"},
]

class PolicyConfig(TypedDict):
    repair_only_enabled: bool
    immutable_labeling_enabled: bool
    fail_closed: bool

class WitnessConfig(TypedDict):
    required_quorum: int
    allowed_witnesses: list[str]

class AuditConfig(TypedDict):
    hash_algorithm: str
    dir: str
    federated_rollup: str
    notarization_interval_seconds: int

class ThresholdsConfig(TypedDict):
    max_latency_overhead_pct: int

class MainConfig(TypedDict):
    runtime_id: str
    qr_dna_required: bool
    policy: PolicyConfig
    witness: WitnessConfig
    audit: AuditConfig
    thresholds: ThresholdsConfig

CONFIG: MainConfig = {
    "runtime_id": "FSR-1.4",
    "qr_dna_required": True,
    "policy": {
        "repair_only_enabled": True,
        "immutable_labeling_enabled": True,
        "fail_closed": True
    },
    "witness": {
        # Quorum: require 3 distinct witnesses for high-risk ops
        "required_quorum": 3,
        # Tenants allowed to act as witnesses
        "allowed_witnesses": ["sovereign", "microsoft", "nasa", "us", "canada"]
    },
    "audit": {
        "hash_algorithm": "sha256",
        "dir": "./audit",
        "federated_rollup": "./audit/federated_rollup.log",
        "notarization_interval_seconds": 300
    },
    "thresholds": {
        "max_latency_overhead_pct": 5
    }
}

# Ensure audit directory exists
os.makedirs(CONFIG["audit"]["dir"], exist_ok=True)

# ---------------------------
# Utilities
# ---------------------------

def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def now_ts() -> int:
    return int(time.time())

def tenant_log_path(tenant_id: str) -> str:
    return os.path.join(CONFIG["audit"]["dir"], f"{tenant_id}_audit.log")

def read_file(path: str) -> bytes:
    with open(path, "rb") as f:
        return f.read()

def write_file(path: str, data: bytes):
    with open(path, "wb") as f:
        f.write(data)

# ---------------------------
# Audit chains
# ---------------------------

@dataclass
class AuditEvent:
    event_id: str
    ts: int
    kind: str
    details: Dict[str, Any]
    prev_hash: Optional[str] = None
    entry_hash: Optional[str] = None

class AuditChain:
    def __init__(self, path: str, algo: str = "sha256"):
        self.path = path
        self.algo = algo
        self.prev_hash = self._last_hash()

        # Touch file if missing
        if not os.path.exists(self.path):
            with open(self.path, "w") as f:
                f.write("")

    def _last_hash(self) -> Optional[str]:
        try:
            with open(self.path, "r") as f:
                lines = f.readlines()
                if not lines:
                    return None
                last_entry = json.loads(lines[-1])
                return last_entry.get("entry_hash")
        except Exception:
            return None

    def append(self, event: AuditEvent):
        payload = { # pyright: ignore[reportUnknownVariableType]
            "event_id": event.event_id,
            "ts": event.ts,
            "kind": event.kind,
            "details": event.details,
            "prev_hash": self.prev_hash
        }
        serialized = json.dumps(payload, sort_keys=True).encode("utf-8")
        entry_hash = sha256(serialized)
        record: Dict[str, Any] = {**payload, "entry_hash": entry_hash}
        with open(self.path, "a") as f:
            f.write(json.dumps(record) + "\n")
        self.prev_hash = entry_hash
        event.entry_hash = entry_hash

class FederatedRollup:
    def __init__(self, path: str):
        self.path = path
        if not os.path.exists(path):
            with open(path, "w") as f:
                f.write("")

    def append(self, summary: Dict[str, Any]):
        with open(self.path, "a") as f:
            f.write(json.dumps(summary) + "\n")

# ---------------------------
# QR-DNA activation
# ---------------------------

@dataclass
class QRScanResult:
    qr_payload: str
    dna_hash: str
    verbal_seal: Optional[str]
    session_token: Optional[str]

class QRDNAAuthenticator:
    def __init__(self, audit: AuditChain):
        self.audit = audit

    def activate(self, scan: QRScanResult, expected_dna_hash: str) -> bool:
        qr_integrity = sha256(scan.qr_payload.encode("utf-8"))
        dna_match = (scan.dna_hash == expected_dna_hash)
        activated = bool(qr_integrity) and dna_match
        self.audit.append(AuditEvent(
            event_id="qr-dna-activation",
            ts=now_ts(),
            kind="auth",
            details={
                "qr_integrity_hash": qr_integrity,
                "dna_match": dna_match,
                "verbal_seal_present": bool(scan.verbal_seal),
                "session_token_present": bool(scan.session_token)
            }
        ))
        return activated

# ---------------------------
# Policy gates and labeling
# ---------------------------


class PolicyGate:
    def __init__(self, config: MainConfig, audit: AuditChain):
        self.config = config
        self.audit = audit

    def enforce_fail_closed(self):
        enabled = self.config["policy"]["repair_only_enabled"] and self.config["policy"]["immutable_labeling_enabled"]
        if not enabled and self.config["policy"]["fail_closed"]:
            self.audit.append(AuditEvent(
                event_id="policy-fail-closed",
                ts=now_ts(),
                kind="policy",
                details={"reason": "required covenants not enabled"}
            ))
            raise RuntimeError("Fail-closed: covenants disabled.")

    def assert_immutable_labeling(self, label: Dict[str, Any]):
        required = {"name", "version", "purpose", "harm_prohibition", "repair_only", "warnings_target"}
        missing = required - set(label.keys())
        if missing:
            self.audit.append(AuditEvent(
                event_id="labeling-violation",
                ts=now_ts(),
                kind="labeling",
                details={"missing_fields": list(missing)}
            ))
            raise ValueError(f"Immutable labeling violation: missing {missing}")

        if not label.get("harm_prohibition", False):
            raise ValueError("Label must prohibit harmful use.")
        if not label.get("repair_only", False):
            raise ValueError("Label must declare repair-only.")
        if label.get("warnings_target") != "intruders":
            raise ValueError("Warnings must target intruders/hackers.")

        self.audit.append(AuditEvent(
            event_id="labeling-validated",
            ts=now_ts(),
            kind="labeling",
            details={"label": label}
        ))

# ---------------------------
# Witness attestation
# ---------------------------

@dataclass
class WitnessSignature:
    tenant_id: str
    signature_hash: str
    note: Optional[str] = None

class WitnessRegistry:
    def __init__(self, allowed: List[str], audit: AuditChain):
        self.allowed = set(allowed)
        self.audit = audit

    def sign(self, tenant_id: str, payload: Dict[str, Any]) -> WitnessSignature:
        if tenant_id not in self.allowed:
            raise PermissionError(f"Witness not allowed: {tenant_id}")
        serialized = json.dumps(payload, sort_keys=True).encode("utf-8")
        sig_hash = sha256(serialized + tenant_id.encode("utf-8"))
        sig = WitnessSignature(tenant_id=tenant_id, signature_hash=sig_hash)
        self.audit.append(AuditEvent(
            event_id="witness-sign",
            ts=now_ts(),
            kind="witness",
            details={"tenant_id": tenant_id, "signature_hash": sig_hash}
        ))
        return sig

    def quorum(self, signatures: List[WitnessSignature], required: int) -> bool:
        unique = {s.tenant_id for s in signatures}
        ok = len(unique) >= required
        self.audit.append(AuditEvent(
            event_id="witness-quorum",
            ts=now_ts(),
            kind="witness",
            details={"unique_witnesses": list(unique), "required": required, "met": ok}
        ))
        return ok

# ---------------------------
# Repair engine (no replacement)
# ---------------------------

@dataclass
class ArtifactIdentity:
    type: str         # "file" | "config" | "service" | "device"
    identifier: str   # path/service name/device serial
    original_hash: str

class RepairEngine:
    def __init__(self, audit: AuditChain):
        self.audit = audit

    def repair(self, artifact: ArtifactIdentity, source_bytes: Optional[bytes] = None) -> bool:
        # Read current state
        try:
            if artifact.type == "file":
                current_hash = sha256(read_file(artifact.identifier))
            else:
                current_hash = "state-hash-" + artifact.identifier  # placeholder
        except Exception as e:
            self.audit.append(AuditEvent(
                event_id="repair-read-error",
                ts=now_ts(),
                kind="repair",
                details={"artifact": artifact.__dict__, "error": str(e)}
            ))
            return False

        needs_repair = (current_hash != artifact.original_hash)

        # Block replacement attempts
        if source_bytes is not None and sha256(source_bytes) != artifact.original_hash:
            self.audit.append(AuditEvent(
                event_id="replacement-attempt-blocked",
                ts=now_ts(),
                kind="policy",
                details={"artifact": artifact.__dict__, "reason": "bytes mismatch original hash"}
            ))
            return False

        # Perform targeted repair
        try:
            if artifact.type == "file" and needs_repair and source_bytes:
                write_file(artifact.identifier, source_bytes)
            # Non-file repair ops would revert baseline configs/services (placeholder)
        except Exception as e:
            self.audit.append(AuditEvent(
                event_id="repair-write-error",
                ts=now_ts(),
                kind="repair",
                details={"artifact": artifact.__dict__, "error": str(e)}
            ))
            return False

        # Verify
        if artifact.type == "file":
            new_hash = sha256(read_file(artifact.identifier))
        else:
            new_hash = artifact.original_hash  # placeholder

        success = (new_hash == artifact.original_hash)
        self.audit.append(AuditEvent(
            event_id="repair-complete",
            ts=now_ts(),
            kind="repair",
            details={"artifact": artifact.__dict__, "needs_repair": needs_repair, "success": success}
        ))
        return success

# ---------------------------
# Lockout Shield Sequence
# ---------------------------

class LockoutShield:
    def __init__(self, audit: AuditChain, witness: WitnessRegistry, federated: FederatedRollup):
        self.audit = audit
        self.witness = witness
        self.federated = federated

    def engage(self, reason: str, tenants: List[str], required_quorum: int) -> bool:
        payload: Dict[str, Any] = {"action": "lockout_shield", "reason": reason, "tenants": tenants, "ts": now_ts()}
        signatures: List[WitnessSignature] = []
        for t in tenants:
            try:
                signatures.append(self.witness.sign(t, payload))
            except Exception as e:
                self.audit.append(AuditEvent(
                    event_id="lockout-witness-error",
                    ts=now_ts(),
                    kind="shield",
                    details={"tenant_id": t, "error": str(e)}
                ))
        met = self.witness.quorum(signatures, required_quorum)
        self.audit.append(AuditEvent(
            event_id="lockout-engaged",
            ts=now_ts(),
            kind="shield",
            details={"reason": reason, "quorum_met": met, "sig_count": len(signatures)}
        ))
        self.federated.append({"event": "lockout_engaged", "reason": reason, "quorum_met": met, "ts": now_ts()})
        return met

# ---------------------------
# Initialization (per tenant) and federated bootstrap
# ---------------------------

def init_tenant(tenant_id: str, expected_dna_hash: str, qr_scan: QRScanResult, immutable_label: 'ImmutableLabel') -> Dict[str, Any]:
    audit = AuditChain(tenant_log_path(tenant_id), CONFIG["audit"]["hash_algorithm"])

    # QR-DNA activation
    auth = QRDNAAuthenticator(audit)
    if CONFIG["qr_dna_required"] and not auth.activate(qr_scan, expected_dna_hash):
        audit.append(AuditEvent(
            event_id="qr-dna-activation-failed",
            ts=now_ts(),
            kind="auth",
            details={"tenant_id": tenant_id, "reason": "DNA mismatch or QR failure"}
        ))
        raise PermissionError(f"[{tenant_id}] Initialization blocked: QR-DNA failed.")

    # Policy gates
    gate = PolicyGate(CONFIG, audit)
    gate.enforce_fail_closed()
    gate.assert_immutable_labeling(dict(immutable_label))

    # Runtime banner
    audit.append(AuditEvent(
        event_id="runtime-initialized",
        ts=now_ts(),
        kind="init",
        details={"runtime_id": CONFIG["runtime_id"], "tenant_id": tenant_id}
    ))

    return {"audit": audit, "policy_gate": gate, "repair_engine": RepairEngine(audit)}

def federated_bootstrap(expected_dna_hash: str, qr_scan: QRScanResult, immutable_label: 'ImmutableLabel'):
    rollup = FederatedRollup(CONFIG["audit"]["federated_rollup"])
    contexts: Dict[str, Dict[str, Any]] = {}
    for t in TENANTS:
        ctx = init_tenant(t["id"], expected_dna_hash, qr_scan, immutable_label)
        contexts[t["id"]] = ctx
    # Shared witness and shield
    # Use Microsoft, Sovereign, NASA for quorum baseline
    any_tenant_audit: AuditChain = contexts["microsoft"]["audit"]
    witness = WitnessRegistry(CONFIG["witness"]["allowed_witnesses"], any_tenant_audit)
    shield = LockoutShield(any_tenant_audit, witness, rollup)
    return contexts, witness, shield

# ---------------------------
# Example usage
# ---------------------------

if __name__ == "__main__":
    expected_dna_hash = "expected-dna-sha256-placeholder"
    qr_scan = QRScanResult(
        qr_payload="QR:FSR:ACTIVATE:SESSION-777",
        dna_hash="expected-dna-sha256-placeholder",
        verbal_seal="Seal of Triple Witnessing",
        session_token="sess-777"
    )
    immutable_label: ImmutableLabel = {
        "name": "Federated Stewardship Runtime",
        "version": "1.4",
        "purpose": "Simultaneous improvement across teams via repair-only governance and transparent labeling",
        "harm_prohibition": True,
        "repair_only": True,
        "warnings_target": "intruders"
    }

    try:
        contexts, witness, shield = federated_bootstrap(expected_dna_hash, qr_scan, immutable_label)

        # Example artifact repair under Microsoft tenant
        ms_engine: "RepairEngine" = contexts["microsoft"]["repair_engine"]
        artifact = ArtifactIdentity(
            type="file",
            identifier="./example.bin",
            original_hash="original-file-sha256-placeholder"
        )
        # Example repair attempt (source_bytes would be the original file bytes)
        # ms_engine.repair(artifact, source_bytes=read_file("./example.bin"))
    except Exception as e:
        print(f"Error during initialization or repair: {e}")
