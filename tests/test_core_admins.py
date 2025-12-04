"""Tests for core admin immutable behavior (Leif William Sogge)."""
import unittest
from diamond_firewall import DiamondFirewall


class TestCoreAdmins(unittest.TestCase):
    def setUp(self):
        guardians = ["Leif William Sogge", "Guardian_A", "Guardian_B"]
        captains = ["C1", "C2", "C3", "C4"]
        self.firewall = DiamondFirewall(guardians, captains, True)
        # Ensure Leif is a core admin
        self.assertTrue(self.firewall.is_core_admin("Leif William Sogge"))

    def test_leif_cannot_update_core_without_consent(self):
        # Attempt to set lockdown keys without guardian consent (other than Leif)
        ok = self.firewall.update_core_config(
            user="Leif William Sogge",
            key="lockdown_keys",
            value={"primary": "NEW_PRIMARY", "secondary_un": "NEW_UN"},
            consenting_guardians=[],
        )
        self.assertFalse(ok)

    def test_leif_can_update_core_with_guardian_consent(self):
        ok = self.firewall.update_core_config(
            user="Leif William Sogge",
            key="lockdown_keys",
            value={"primary": "NEW_PRIMARY", "secondary_un": "NEW_UN"},
            consenting_guardians=["Guardian_A"],
        )
        self.assertTrue(ok)
        # Verify keys were applied
        self.assertEqual(self.firewall.lockdown_keys["primary"], "NEW_PRIMARY")
        self.assertEqual(self.firewall.lockdown_keys["secondary_un"], "NEW_UN")

    def test_non_core_cannot_update_core(self):
        ok = self.firewall.update_core_config(
            user="RandomUser",
            key="lockdown_keys",
            value={"primary": "X", "secondary_un": "Y"},
            consenting_guardians=["Guardian_A"],
        )
        self.assertFalse(ok)

    def test_update_notification_recipients_with_consent(self):
        ok = self.firewall.update_core_config(
            user="Leif William Sogge",
            key="notification_recipients",
            value={"space_leaf": "http://sl.local/hook", "un": None, "president_of_the_united_states": None},
            consenting_guardians=["Guardian_B"],
        )
        self.assertTrue(ok)
        self.assertEqual(self.firewall.notification_recipients["space_leaf"], "http://sl.local/hook")


if __name__ == "__main__":
    unittest.main()
