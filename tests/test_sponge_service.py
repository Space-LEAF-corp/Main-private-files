"""Tests for HoneySponge service (defensive honeypot/sponger)."""
import unittest
from diamond_firewall import DiamondFirewall


class TestSpongeService(unittest.TestCase):
    def setUp(self):
        guardians = ["G1", "G2"]
        captains = ["C1", "C2", "C3", "C4"]
        self.firewall = DiamondFirewall(guardians, captains, True)
        self.notifications = []

        def cb(attacker, message, step):
            self.notifications.append((attacker, message, step))

        self.firewall.notifier_callback = cb
        # Set recipients to None (no real HTTP calls)
        self.firewall.set_notification_recipients(space_leaf=None, un=None, president_of_the_united_states=None)

    def test_run_sponge_cycle_soft_harden(self):
        # Simulate an attacker in mirror_layer but with low warnings
        self.firewall.mirror_layer["att1"] = "fakehash"
        self.firewall.intruder_warnings["att1"] = 1
        summary = self.firewall.run_sponge_cycle()
        self.assertEqual(summary["processed"], 1)
        # Defensive policy hardened should trigger a notifier callback
        self.assertTrue(any("Defensive policy hardened" in n[1] for n in self.notifications))

    def test_run_sponge_cycle_escalate_and_block(self):
        # Simulate attacker with many warnings
        self.firewall.mirror_layer["badguy"] = "fakehash"
        self.firewall.intruder_warnings["badguy"] = 5
        summary = self.firewall.run_sponge_cycle()
        self.assertIn("badguy", summary["blocked"]) or self.assertIn("badguy", summary["escalated"]) 
        self.assertEqual(self.firewall.mirror_layer.get("badguy"), "ESCALATED_TO_AUTHORITIES")
        self.assertIn("badguy", self.firewall.blocklist)


if __name__ == "__main__":
    unittest.main()
