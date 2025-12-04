"""Tests for webhook configuration and UN-limited lockdown access."""
import unittest
from diamond_firewall import DiamondFirewall


class TestWebhookAndUNLimited(unittest.TestCase):
    def setUp(self):
        self.guardians = ["G1", "G2"]
        self.captains = ["C1", "C2", "C3", "C4"]
        self.firewall = DiamondFirewall(self.guardians, self.captains, True)
        self.firewall.set_lockdown_keys(primary="SPACELEAF_MAIN_KEY", secondary_un="UN_KEY")

    def test_set_notifier_webhook(self):
        # Set webhook URL (won't actually be called in test)
        self.firewall.set_notifier_webhook("http://example.local/hook", {"X-API-KEY": "abc"})
        self.assertEqual(self.firewall.notifier_webhook_url, "http://example.local/hook")

    def test_un_limited_access_during_lockdown(self):
        # Activate lockdown via UN key
        self.firewall.activate_lockdown_with_un("UN_KEY")
        # Attempt access with UN key should result in limited access message
        res = self.firewall.access_request(network="internet", metadata={"admin_key": "UN_KEY"})
        self.assertIn("limited access", res)

    def test_primary_key_full_access_during_lockdown(self):
        self.firewall.activate_lockdown_with_primary("SPACELEAF_MAIN_KEY")
        res = self.firewall.access_request(network="internet", metadata={"admin_key": "SPACELEAF_MAIN_KEY"})
        self.assertIn("Access granted by primary lockdown key", res)


if __name__ == "__main__":
    unittest.main()
