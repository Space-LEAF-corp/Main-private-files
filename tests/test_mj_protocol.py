"""Unit tests for MJ Protocol."""

import os
import tempfile
import unittest
import time
from typing import Any, Dict

from mj_protocol import (
    CeremonialsManager,
    MJLocalLayer,
    MJSatelliteLayer,
    ResonanceType,
    SealEntry,
)
from secured_firewall import SecuredFirewall


class TestResonanceType(unittest.TestCase):
    def test_resonance_types_exist(self):
        self.assertEqual(ResonanceType.FIREWALL.value, "Seal of Firewall Resonance 1.0")
        self.assertEqual(ResonanceType.RESILIENCE.value, "Seal of Ridiculous Resilience 1.0")
        self.assertEqual(ResonanceType.REMINDER.value, "Seal of Gentle Reminder 1.0")
        self.assertEqual(ResonanceType.MIRROR_JUNCTION.value, "Seal of Mirror Junction 1.0")


class TestSealEntry(unittest.TestCase):
    def test_seal_integrity(self):
        seal = SealEntry(
            seal_type=ResonanceType.FIREWALL,
            timestamp=time.time(),
            author="Leif",
            description="Test seal",
        )
        integrity = seal.compute_integrity()
        seal.integrity_hash = integrity
        self.assertTrue(seal.verify())

    def test_seal_tampering_detected(self):
        seal = SealEntry(
            seal_type=ResonanceType.FIREWALL,
            timestamp=time.time(),
            author="Leif",
            description="Test seal",
        )
        integrity = seal.compute_integrity()
        seal.integrity_hash = integrity
        seal.author = "tampered"  # Tamper with the seal
        self.assertFalse(seal.verify())


class TestCeremonialsManager(unittest.TestCase):
    def setUp(self):
        fd, self.archive_file = tempfile.mkstemp(suffix=".json")
        os.close(fd)
        try:
            os.remove(self.archive_file)
        except Exception:
            pass

    def tearDown(self):
        try:
            os.remove(self.archive_file)
        except Exception:
            pass

    def test_inscribe_seal(self):
        cm = CeremonialsManager(self.archive_file)
        seal = cm.inscribe_seal(
            ResonanceType.FIREWALL,
            "Guardian",
            "Test inscription",
        )
        self.assertIsNotNone(seal.integrity_hash)
        self.assertTrue(cm.verify_seal(0))

    def test_heritage_chain(self):
        cm = CeremonialsManager(self.archive_file)
        cm.inscribe_seal(ResonanceType.FIREWALL, "Guardian", "Seal 1")
        cm.inscribe_seal(ResonanceType.RESILIENCE, "Guardian", "Seal 2")
        heritage: Dict[str, Any] = cm.heritage_chain()  # type: ignore
        self.assertEqual(heritage["total_seals"], 2)
        self.assertEqual(heritage["verified_seals"], 2)


class TestMJLocalLayer(unittest.TestCase):
    def setUp(self):
        fd, self.auth_file = tempfile.mkstemp(suffix=".json")
        os.close(fd)
        try:
            os.remove(self.auth_file)
        except Exception:
            pass
        fd, self.archive_file = tempfile.mkstemp(suffix=".json")
        os.close(fd)
        try:
            os.remove(self.archive_file)
        except Exception:
            pass

        self.firewall = SecuredFirewall(users_file=self.auth_file)
        self.mj_local = MJLocalLayer(self.firewall.auth_manager, self.firewall)
        self.mj_local.ceremonials.archive_file = self.archive_file

    def tearDown(self):
        try:
            os.remove(self.auth_file)
        except Exception:
            pass
        try:
            os.remove(self.archive_file)
        except Exception:
            pass

    def test_register_and_seal(self):
        res: Dict[str, Any] = self.mj_local.register_and_seal(  # type: ignore[reportUnknownMemberType]
            "alice", "password123", "LINEAGE_SAFE_ALICE_001"
        )
        assert isinstance(res, dict)
        self.assertEqual(res.get("status"), "ok")
        self.assertIn("seal", res)

    def test_login_and_check(self):
        self.mj_local.register_and_seal("bob", "password123", "LINEAGE_SAFE_BOB_001")  # type: ignore[reportUnknownMemberType]
        res: Dict[str, Any] = self.mj_local.login_and_check("bob", "password123", "LINEAGE_SAFE_BOB_001")  # type: ignore[no-untyped-call, assignment, arg-type]
        self.assertEqual(res.get("status"), "challenge")
        self.assertIn("otp_token", res)


class TestMJSatelliteLayer(unittest.TestCase):
    def test_satellite_activation(self):
        sat = MJSatelliteLayer()
        res: Dict[str, Any] = sat.activate()  # type: ignore[assignment, call-arg, unknown-member, reportUnknownMemberType, reportUnknownMemberType, type-arg, reportUnknownMemberType, type-unknown]
        self.assertEqual(res.get("status"), "ok")
        self.assertTrue(sat.connected)

    def test_satellite_status(self):
        sat = MJSatelliteLayer()
        status: Dict[str, Any] = sat.status()  # type: ignore
        self.assertFalse(status.get("connected"))
        sat.activate() # pyright: ignore[reportUnknownMemberType]
        status = sat.status()  # type: ignore
        self.assertTrue(status.get("connected"))


if __name__ == "__main__":
    unittest.main()
