"""MJ Protocol: Mirror Junction — Lineage-Safe Dual-Layer Pivot Node.

This module provides:
1. Resonance types and ceremonial seals
2. Seal integrity verification
3. Satellite broadcast coordination
4. Alexandria archive logging
"""
from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass, asdict
from enum import Enum
from typing import Dict, Optional


class ResonanceType(Enum):
    """Ceremonial seal types for MJ protocol."""
    FIREWALL = "Seal of Firewall Resonance 1.0"
    RESILIENCE = "Seal of Ridiculous Resilience 1.0"
    REMINDER = "Seal of Gentle Reminder 1.0"
    MIRROR_JUNCTION = "Seal of Mirror Junction 1.0"
    UNKNOWN = "Unclassified"


@dataclass
class SealEntry:
    """A ceremonial seal entry in the Captain's Log."""
    seal_type: ResonanceType
    timestamp: float
    author: str
    description: str
    commit_hash: Optional[str] = None
    integrity_hash: str = ""

    def compute_integrity(self) -> str:
        """Compute SHA256 integrity hash for this seal."""
        data = json.dumps({
            "seal": self.seal_type.value,
            "timestamp": self.timestamp,
            "author": self.author,
            "description": self.description,
            "commit": self.commit_hash,
        }, sort_keys=True)
        return hashlib.sha256(data.encode()).hexdigest()

    def verify(self) -> bool:
        """Verify seal integrity."""
        if not self.integrity_hash:
            return False
        computed = self.compute_integrity()
        return computed == self.integrity_hash

    def to_dict(self) -> Dict:
        return {
            "seal_type": self.seal_type.name,  # Store as string name
            "timestamp": self.timestamp,
            "author": self.author,
            "description": self.description,
            "commit_hash": self.commit_hash,
            "integrity_hash": self.integrity_hash,
        }


class CeremonialsManager:
    """Manages ceremonial seals and resonance tagging."""

    def __init__(self, archive_file: str = "alexandria_of_joy.json"):
        self.archive_file = archive_file
        self._seals: list[SealEntry] = []
        self._load_archive()

    def _load_archive(self) -> None:
        """Load seals from Alexandria archive."""
        try:
            with open(self.archive_file, "r", encoding="utf-8") as f:
                data = json.load(f)
            self._seals = [
                SealEntry(
                    seal_type=ResonanceType[entry.get("seal_type", "UNKNOWN")],
                    timestamp=entry.get("timestamp", time.time()),
                    author=entry.get("author", "unknown"),
                    description=entry.get("description", ""),
                    commit_hash=entry.get("commit_hash"),
                    integrity_hash=entry.get("integrity_hash", ""),
                )
                for entry in data.get("seals", [])
            ]
        except Exception:
            self._seals = []

    def _save_archive(self) -> None:
        """Save seals to Alexandria archive."""
        data = {
            "version": "1.0",
            "timestamp": time.time(),
            "seals": [s.to_dict() for s in self._seals],
        }
        with open(self.archive_file, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)

    def inscribe_seal(
        self,
        seal_type: ResonanceType,
        author: str,
        description: str,
        commit_hash: Optional[str] = None,
    ) -> SealEntry:
        """Inscribe a new ceremonial seal."""
        seal = SealEntry(
            seal_type=seal_type,
            timestamp=time.time(),
            author=author,
            description=description,
            commit_hash=commit_hash,
        )
        seal.integrity_hash = seal.compute_integrity()
        self._seals.append(seal)
        self._save_archive()
        return seal

    def verify_seal(self, seal_index: int) -> bool:
        """Verify a seal by index."""
        if 0 <= seal_index < len(self._seals):
            return self._seals[seal_index].verify()
        return False

    def list_seals(self) -> list[SealEntry]:
        """List all inscribed seals."""
        return list(self._seals)

    def get_seals_by_type(self, seal_type: ResonanceType) -> list[SealEntry]:
        """Get all seals of a specific type."""
        return [s for s in self._seals if s.seal_type == seal_type]

    def heritage_chain(self) -> Dict:
        """Generate inheritance chain for future captains."""
        verified = sum(1 for s in self._seals if s.verify())
        return {
            "total_seals": len(self._seals),
            "verified_seals": verified,
            "integrity_rate": verified / len(self._seals) if self._seals else 0,
            "seals_by_type": {
                rt.name: len(self.get_seals_by_type(rt))
                for rt in ResonanceType
            },
            "timestamp": time.time(),
        }


class MJLocalLayer:
    """Local use layer of MJ protocol (ground station)."""

    def __init__(self, auth_manager, firewall):
        self.auth = auth_manager
        self.firewall = firewall
        self.ceremonials = CeremonialsManager()

    def register_and_seal(self, user_id: str, password: str, dna: str) -> Dict:
        """Register user and inscribe seal."""
        try:
            self.auth.register(user_id, password, dna)
            seal = self.ceremonials.inscribe_seal(
                ResonanceType.FIREWALL,
                author=user_id,
                description=f"User {user_id} registered and bound to lineage",
            )
            return {
                "status": "ok",
                "seal": seal.to_dict(),
                "message": f"User {user_id} registered with seal integrity",
            }
        except ValueError as e:
            return {"status": "error", "message": str(e)}

    def login_and_check(self, user_id: str, password: str, dna: str) -> Dict:
        """Login and verify against seals."""
        res1 = self.firewall.login_step1(user_id, password, dna)
        if res1.get("status") == "challenge":
            # Log the attempt
            seal = self.ceremonials.inscribe_seal(
                ResonanceType.REMINDER,
                author=user_id,
                description=f"Login attempt by {user_id} — OTP challenge issued",
            )
            return {
                "status": "challenge",
                "otp_token": res1.get("otp_token"),
                "seal_log": seal.to_dict(),
            }
        return res1

    def complete_login_and_verify(self, user_id: str, otp: str, commit_hash: Optional[str] = None) -> Dict:
        """Complete login and verify seal chain."""
        res2 = self.firewall.login_step2(user_id, otp)
        if res2.get("status") == "ok":
            seal = self.ceremonials.inscribe_seal(
                ResonanceType.MIRROR_JUNCTION,
                author=user_id,
                description=f"User {user_id} authenticated — mirror junction active",
                commit_hash=commit_hash,
            )
            heritage = self.ceremonials.heritage_chain()
            return {
                "status": "ok",
                "session_token": res2.get("session_token"),
                "seal_verification": seal.to_dict(),
                "heritage_integrity": heritage,
            }
        return res2


class MJSatelliteLayer:
    """Satellite layer of MJ protocol (orbital broadcast) — stub for activation."""

    def __init__(self):
        self.broadcast_queue: list[Dict] = []
        self.connected = False

    def activate(self) -> Dict:
        """Activate satellite connection."""
        self.connected = True
        return {"status": "ok", "message": "Satellite layer activated"}

    def broadcast_seal(self, seal: SealEntry) -> Dict:
        """Broadcast a seal to distributed nodes."""
        if not self.connected:
            return {"status": "error", "message": "Satellite not connected"}
        self.broadcast_queue.append(seal.to_dict())
        return {"status": "ok", "message": f"Seal broadcast queued"}

    def status(self) -> Dict:
        return {
            "connected": self.connected,
            "queue_length": len(self.broadcast_queue),
            "timestamp": time.time(),
        }
