"""Tests for lockdown keys and progressive intruder notifications."""
import unittest
from diamond_firewall import DiamondFirewall


class TestLockdownAndNotifications(unittest.TestCase):
    def setUp(self):
        self.guardians = ["G1", "G2"]
        self.captains = ["C1", "C2", "C3", "C4"]
        self.firewall = DiamondFirewall(self.guardians, self.captains, True)
        # set keys
        self.firewall.set_lockdown_keys(primary="SPACELEAF_MAIN_KEY", secondary_un="UN_SECONDARY_KEY")

    def test_activate_lockdown_primary(self):
        ok = self.firewall.activate_lockdown_with_primary("SPACELEAF_MAIN_KEY")
        self.assertTrue(ok)
        with self.assertRaises(Exception):
            self.firewall.access_request(dna_code="LINEAGE_SAFE_1")

    def test_activate_lockdown_un(self):
        # deactivate first
        self.firewall.deactivate_lockdown("SPACELEAF_MAIN_KEY")
        ok = self.firewall.activate_lockdown_with_un("UN_SECONDARY_KEY")
        self.assertTrue(ok)
        with self.assertRaises(Exception):
            self.firewall.access_request(dna_code="LINEAGE_SAFE_1")

    def test_deactivate_lockdown_with_un_or_primary(self):
        self.firewall.activate_lockdown_with_primary("SPACELEAF_MAIN_KEY")
        self.assertTrue(self.firewall.deactivate_lockdown("UN_SECONDARY_KEY"))
        # reactivate then deactivate with primary
        self.firewall.activate_lockdown_with_un("UN_SECONDARY_KEY")
        self.assertTrue(self.firewall.deactivate_lockdown("SPACELEAF_MAIN_KEY"))

    def test_progressive_warnings_and_escalation(self):
        attacker = "Hacker_X"
        # First attempt -> first warning
        msg1 = self.firewall.intrusion_attempt(attacker, "payload1")
        self.assertIn("WARNING 1", msg1)
        # Second attempt -> second warning
        msg2 = self.firewall.intrusion_attempt(attacker, "payload2")
        self.assertIn("WARNING 2", msg2)
        # Third attempt -> final warning + escalation
        msg3 = self.firewall.intrusion_attempt(attacker, "payload3")
        self.assertIn("WARNING 3", msg3)
        # After escalation, attacker should be marked escalated
        self.assertEqual(self.firewall.is_trapped(attacker), "ESCALATED_TO_AUTHORITIES")

    def test_notifier_callback_invoked(self):
        from typing import List, Tuple, Any
        called: List[Tuple[str, Any, int]] = []

        def notifier(attacker_id: str, message: Any, step: int):
            called.append((attacker_id, message, step))

        setattr(self.firewall, "notifier_callback", notifier)
        attacker = "Hacker_Y"
        self.firewall.intrusion_attempt(attacker, "p1")
        self.firewall.intrusion_attempt(attacker, "p2")
        self.firewall.intrusion_attempt(attacker, "p3")
        # notifier should have been called 3 times
        self.assertEqual(len(called), 3)
        self.assertEqual(called[-1][2], 3)


if __name__ == "__main__":
    unittest.main()
