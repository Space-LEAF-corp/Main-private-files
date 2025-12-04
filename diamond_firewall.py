"""Diamond Firewall — cleaned and importable version.

This file is a cleaned copy of the original `# diamond_firewall.py` (renamed
to remove the leading #). It implements a minimal lineage-safe firewall
simulation with thread-safety and simple unit tests at the bottom.
"""
from __future__ import annotations
import hashlib
import logging
import threading
import time
import os
from typing import Dict, List, Optional

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DiamondFirewall")


class AccessDeniedError(Exception):
    pass


class DiamondFirewall:
    """
    Minimal lineage-safe firewall simulation:
      - mirror_layer stores attacker traps (fake hashes)
      - diamond_layer holds renewal fingerprints for self-renewal cycles
      - dna_qr_check enforces collective consent
    Thread-safe for concurrent intrusion checks.
    """

    def __init__(
        self,
        guardians: List[str],
        captains: List[str],
        un_consent: bool,
        required_captains: int = 4,
        hash_algo: str = "sha256",
    ):
        self.guardians = list(guardians)
        self.captains = list(captains)
        self.un_consent = bool(un_consent)
        self.required_captains = int(required_captains)
        self.hash_algo = hash_algo
        self.mirror_layer: Dict[str, str] = {}
        self.diamond_layer: Dict[str, str] = {}
        self.active = True
        self._lock = threading.RLock()

    def _make_hash(self, data: str) -> str:
        h = hashlib.new(self.hash_algo)
        h.update(data.encode("utf-8"))
        return h.hexdigest()

    def dna_qr_check(self, dna_code: str, required_captains: Optional[int] = None) -> bool:
        """
        Verify DNA QR code and collective consent.
        dna_code must begin with "LINEAGE_SAFE" to be considered valid.
        """
        required = self.required_captains if required_captains is None else required_captains
        with self._lock:
            ok = (
                isinstance(dna_code, str)
                and dna_code.startswith("LINEAGE_SAFE")
                and len(self.captains) >= required
                and self.un_consent is True
                and len(self.guardians) > 0
            )
            logger.debug("dna_qr_check result=%s", ok)
            return ok

    def intrusion_attempt(self, attacker_id: str, payload: str) -> str:
        """
        Redirect attacker into mirror layer (phantom system).
        Returns a summary message. Trapped values are fake hashes.
        """
        if not self.active:
            logger.warning("Firewall inactive; intrusion attempt ignored.")
            return "Firewall inactive."

        salt = f"{time.time_ns()}:{os.getpid()}"
        fake_hash = self._make_hash(payload + salt)
        with self._lock:
            self.mirror_layer[attacker_id] = fake_hash
            logger.info("Attacker %s trapped in mirror layer.", attacker_id)
        return f"Attacker {attacker_id} trapped in mirror layer."

    def collapse_cycle(self, purge_mirror: bool = True) -> str:
        """
        Firewall collapse + inversion cycle (self-renewal).
        Optionally purges mirror_layer (simulated cleanup).
        """
        if not self.active:
            logger.info("Collapse cycle requested while inactive — re-activating.")
            self.active = True

        timestamp = time.time()
        renewal_hash = self._make_hash(str(timestamp))
        with self._lock:
            self.diamond_layer["renewal"] = renewal_hash
            if purge_mirror:
                self.mirror_layer.clear()
                logger.info("Mirror layer purged during collapse cycle.")
        logger.info("Diamond firewall collapsed inward — shield renewed.")
        return "Diamond firewall collapsed inward, anomalies purged, shield renewed."

    def access_request(self, dna_code: str) -> str:
        """
        Attempt to access true firewall. Raises AccessDeniedError on refusal.
        Returns success message on grant.
        """
        if self.dna_qr_check(dna_code):
            logger.info("Access granted: collective consent confirmed.")
            return "Access granted: collective consent confirmed."
        logger.warning("Access denied: lineage-safe consent required.")
        raise AccessDeniedError("Access denied: lineage-safe consent required.")

    def is_trapped(self, attacker_id: str) -> Optional[str]:
        with self._lock:
            return self.mirror_layer.get(attacker_id)

    def __repr__(self) -> str:
        return (
            f"<DiamondFirewall guardians={len(self.guardians)} "
            f"captains={len(self.captains)} active={self.active}>"
        )


# Basic unit tests
if __name__ == "__main__":
    import unittest

    class TestDiamondFirewall(unittest.TestCase):
        def setUp(self):
            self.guardians = ["Miko", "JD", "Guardian_A", "Guardian_B"]
            self.captains = ["Captain_1", "Captain_2", "Captain_3", "Captain_4"]
            self.firewall = DiamondFirewall(self.guardians, self.captains, True)

        def test_dna_qr_check_pass(self):
            self.assertTrue(self.firewall.dna_qr_check("LINEAGE_SAFE_QR123"))

        def test_access_request_granted(self):
            self.assertEqual(
                self.firewall.access_request("LINEAGE_SAFE_QR123"),
                "Access granted: collective consent confirmed.",
            )

        def test_access_request_denied_raises(self):
            with self.assertRaises(AccessDeniedError):
                self.firewall.access_request("INVALID_QR")

        def test_intrusion_and_trap(self):
            res = self.firewall.intrusion_attempt("Hacker_007", "malware_payload")
            self.assertIn("trapped in mirror layer", res)
            self.assertIsNotNone(self.firewall.is_trapped("Hacker_007"))

        def test_collapse_cycle_purges_mirror(self):
            self.firewall.intrusion_attempt("Hacker_007", "payload")
            self.assertIn("Hacker_007", self.firewall.mirror_layer)
            self.firewall.collapse_cycle()
            self.assertNotIn("Hacker_007", self.firewall.mirror_layer)

    unittest.main()
