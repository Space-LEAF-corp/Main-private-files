"""Unit tests for network-aware DiamondFirewall policies."""
import unittest
from diamond_firewall import DiamondFirewall, AccessDeniedError


class TestNetworkFirewall(unittest.TestCase):
    def setUp(self):
        self.guardians = ["G1", "G2"]
        self.captains = ["C1", "C2", "C3", "C4"]
        self.firewall = DiamondFirewall(self.guardians, self.captains, True)

    def test_default_internet_access_requires_qr(self):
        with self.assertRaises(AccessDeniedError):
            self.firewall.access_request()

    def test_wifi_trusted_ssid_allows_without_qr(self):
        # Configure trusted SSID
        self.firewall.set_network_policy("wifi", {"trusted_ssids": ["HomeNet"] , "require_qr": False})
        res = self.firewall.access_request(network="wifi", metadata={"ssid": "HomeNet"})
        self.assertIn("Access granted", res)

    def test_wifi_untrusted_requires_qr(self):
        # Ensure wifi requires QR by default
        with self.assertRaises(AccessDeniedError):
            self.firewall.access_request(network="wifi")

    def test_cellular_allowed_carrier(self):
        # Allow a carrier
        self.firewall.set_network_policy("cellular", {"allowed_carriers": ["TelcoX"], "uplink_strict": True})
        # Provide carrier metadata and dna_code
        res = self.firewall.access_request(dna_code="LINEAGE_SAFE_ABC", network="cellular", metadata={"carrier": "TelcoX"})
        self.assertIn("Access granted", res)

    def test_cellular_denied_unknown_carrier(self):
        self.firewall.set_network_policy("cellular", {"allowed_carriers": ["TelcoX"], "uplink_strict": True})
        with self.assertRaises(AccessDeniedError):
            self.firewall.access_request(dna_code="LINEAGE_SAFE_ABC", network="cellular", metadata={"carrier": "BadCarrier"})

    def test_satellite_allowed(self):
        self.firewall.set_network_policy("satellite", {"allowed_satellites": ["Sat-9"], "uplink_strict": True})
        res = self.firewall.access_request(dna_code="LINEAGE_SAFE_1", network="satellite", metadata={"satellite_id": "Sat-9"})
        self.assertIn("Access granted", res)

    def test_satellite_denied_without_allowed_list(self):
        # By default, strict uplink will block unless allowed list empty and not provided
        self.firewall.set_network_policy("satellite", {"uplink_strict": True, "allowed_satellites": ["Sat-9"]})
        with self.assertRaises(AccessDeniedError):
            self.firewall.access_request(dna_code="LINEAGE_SAFE_1", network="satellite", metadata={"satellite_id": "Unknown"})

if __name__ == "__main__":
    unittest.main()
