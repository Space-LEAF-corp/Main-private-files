"""Unit tests for Science Lab module."""
import unittest

from science_lab import (
    CryptoLab,
    NetworkSecurityLab,
    AuthenticationProtocolLab,
    SecurityMetricsLab,
    ScienceLab
)


class TestCryptoLab(unittest.TestCase):
    """Tests for CryptoLab class."""
    
    def test_demo_hashing(self):
        """Test hashing demonstration."""
        result = CryptoLab.demo_hashing("test data")
        
        self.assertEqual(result["experiment"], "Hashing Algorithms")
        self.assertEqual(result["input"], "test data")
        self.assertIn("algorithms", result)
        self.assertIn("SHA-256", result["algorithms"])
        self.assertIn("SHA-512", result["algorithms"])
        self.assertIn("MD5", result["algorithms"])
        
    def test_demo_pbkdf2(self):
        """Test PBKDF2 demonstration."""
        result = CryptoLab.demo_pbkdf2("TestPassword")
        
        self.assertEqual(result["experiment"], "PBKDF2 Password Hashing")
        self.assertEqual(result["password"], "TestPassword")
        self.assertIn("salt_hex", result)
        self.assertIn("derived_key", result)
        self.assertIn("computation_time_ms", result)
        self.assertEqual(result["iterations"], 100000)
        
    def test_demo_hmac(self):
        """Test HMAC demonstration."""
        result = CryptoLab.demo_hmac("test message", "secret")
        
        self.assertEqual(result["experiment"], "HMAC - Message Authentication")
        self.assertEqual(result["message"], "test message")
        self.assertEqual(result["key"], "secret")
        self.assertIn("hmac", result)
        self.assertEqual(result["verification"], "Valid")
        
    def test_demo_token_generation(self):
        """Test token generation demonstration."""
        result = CryptoLab.demo_token_generation()
        
        self.assertEqual(result["experiment"], "Secure Token Generation")
        self.assertIn("tokens", result)
        self.assertIn("16_bytes", result["tokens"])
        self.assertIn("32_bytes", result["tokens"])
        self.assertIn("64_bytes", result["tokens"])
        
        # Verify token lengths
        self.assertEqual(len(result["tokens"]["16_bytes"]), 32)  # 16 bytes = 32 hex chars
        self.assertEqual(len(result["tokens"]["32_bytes"]), 64)  # 32 bytes = 64 hex chars
        
    def test_analyze_password_strength_weak(self):
        """Test password strength analysis for weak password."""
        result = CryptoLab.analyze_password_strength("abc")
        
        self.assertEqual(result["experiment"], "Password Strength Analysis")
        self.assertEqual(result["password"], "abc")
        self.assertEqual(result["length"], 3)
        self.assertEqual(result["strength"], "Weak")
        
    def test_analyze_password_strength_strong(self):
        """Test password strength analysis for strong password."""
        result = CryptoLab.analyze_password_strength("MyStr0ng!P@ssw0rd")
        
        self.assertEqual(result["length"], 17)
        self.assertTrue(result["characteristics"]["uppercase"])
        self.assertTrue(result["characteristics"]["lowercase"])
        self.assertTrue(result["characteristics"]["digits"])
        self.assertTrue(result["characteristics"]["special_chars"])
        self.assertEqual(result["strength"], "Strong")


class TestNetworkSecurityLab(unittest.TestCase):
    """Tests for NetworkSecurityLab class."""
    
    def test_analyze_firewall_config(self):
        """Test firewall configuration analysis."""
        guardians = ["G1", "G2", "G3"]
        captains = ["C1", "C2", "C3", "C4"]
        
        result = NetworkSecurityLab.analyze_firewall_config(guardians, captains)
        
        self.assertEqual(result["experiment"], "Firewall Configuration Analysis")
        self.assertEqual(result["guardians"]["count"], 3)
        self.assertEqual(result["captains"]["count"], 4)
        self.assertEqual(result["guardians"]["redundancy"], "High")
        self.assertEqual(result["captains"]["redundancy"], "High")
        
    def test_analyze_firewall_config_low_redundancy(self):
        """Test firewall config with low redundancy."""
        guardians = ["G1"]
        captains = ["C1", "C2"]
        
        result = NetworkSecurityLab.analyze_firewall_config(guardians, captains)
        
        self.assertEqual(result["guardians"]["redundancy"], "Low")
        self.assertEqual(result["captains"]["redundancy"], "Low")
        
    def test_simulate_intrusion_scenarios(self):
        """Test intrusion scenario simulation."""
        result = NetworkSecurityLab.simulate_intrusion_scenarios()
        
        self.assertEqual(result["experiment"], "Intrusion Scenario Simulation")
        self.assertIn("scenarios", result)
        self.assertGreater(len(result["scenarios"]), 0)
        
        # Check scenarios have required fields
        for scenario in result["scenarios"]:
            self.assertIn("scenario", scenario)
            self.assertIn("description", scenario)
            self.assertIn("detection", scenario)
            self.assertIn("mitigation", scenario)
            
    def test_measure_authentication_latency(self):
        """Test authentication latency measurement."""
        result = NetworkSecurityLab.measure_authentication_latency()
        
        self.assertEqual(result["experiment"], "Authentication Latency Measurement")
        self.assertIn("operations", result)
        self.assertIn("registration", result["operations"])
        self.assertIn("login_step1", result["operations"])
        self.assertIn("login_step2", result["operations"])
        self.assertIn("total", result["operations"])


class TestAuthenticationProtocolLab(unittest.TestCase):
    """Tests for AuthenticationProtocolLab class."""
    
    def test_explain_two_factor_auth(self):
        """Test two-factor authentication explanation."""
        result = AuthenticationProtocolLab.explain_two_factor_auth()
        
        self.assertEqual(result["experiment"], "Two-Factor Authentication (2FA) Explanation")
        self.assertIn("factors", result)
        self.assertIn("something_you_know", result["factors"])
        self.assertIn("something_you_have", result["factors"])
        self.assertIn("system_implementation", result)
        self.assertIn("flow", result)
        
    def test_explain_session_management(self):
        """Test session management explanation."""
        result = AuthenticationProtocolLab.explain_session_management()
        
        self.assertEqual(result["experiment"], "Session Management Explained")
        self.assertIn("token_lifecycle", result)
        self.assertIn("security_measures", result)
        self.assertIn("best_practices", result)
        
    def test_analyze_qr_dna_binding(self):
        """Test QR-DNA binding analysis."""
        result = AuthenticationProtocolLab.analyze_qr_dna_binding()
        
        self.assertEqual(result["experiment"], "QR-DNA Binding Analysis")
        self.assertEqual(result["format"], "LINEAGE_SAFE_<IDENTITY>_<SEQUENCE>")
        self.assertEqual(result["example"], "LINEAGE_SAFE_ALICE_001")
        self.assertIn("security_properties", result)
        self.assertIn("attack_resistance", result)


class TestSecurityMetricsLab(unittest.TestCase):
    """Tests for SecurityMetricsLab class."""
    
    def test_calculate_entropy_simple(self):
        """Test entropy calculation for simple password."""
        result = SecurityMetricsLab.calculate_entropy("abc")
        
        self.assertEqual(result["experiment"], "Password Entropy Analysis")
        self.assertEqual(result["password"], "abc")
        self.assertEqual(result["length"], 3)
        self.assertIn("entropy_bits", result)
        self.assertIn("security_rating", result)
        
    def test_calculate_entropy_complex(self):
        """Test entropy calculation for complex password."""
        result = SecurityMetricsLab.calculate_entropy("MyC0mpl3x!P@ss")
        
        self.assertEqual(result["length"], 14)
        self.assertGreater(result["entropy_bits"], 50)
        self.assertIn(result["security_rating"], ["Fair", "Good", "Excellent"])
        
    def test_compare_hashing_performance(self):
        """Test hashing performance comparison."""
        result = SecurityMetricsLab.compare_hashing_performance()
        
        self.assertEqual(result["experiment"], "Hashing Algorithm Performance Comparison")
        self.assertIn("results", result)
        self.assertIn("sha256", result["results"])
        self.assertIn("md5", result["results"])
        self.assertIn("fastest", result)
        
        # Verify all algorithms have required metrics
        for algo, metrics in result["results"].items():
            self.assertIn("total_time_ms", metrics)
            self.assertIn("avg_time_us", metrics)


class TestScienceLab(unittest.TestCase):
    """Tests for main ScienceLab interface."""
    
    def setUp(self):
        """Set up science lab."""
        self.lab = ScienceLab()
        
    def test_initialization(self):
        """Test science lab initializes correctly."""
        self.assertIsNotNone(self.lab.crypto)
        self.assertIsNotNone(self.lab.network)
        self.assertIsNotNone(self.lab.auth)
        self.assertIsNotNone(self.lab.metrics)
        
    def test_run_all_experiments(self):
        """Test running all experiments."""
        results = self.lab.run_all_experiments()
        
        self.assertEqual(results["science_lab"], "Complete Experiment Suite")
        self.assertIn("experiments", results)
        self.assertIn("cryptography", results["experiments"])
        self.assertIn("network_security", results["experiments"])
        self.assertIn("authentication", results["experiments"])
        self.assertIn("metrics", results["experiments"])
        
    def test_get_experiment_catalog(self):
        """Test getting experiment catalog."""
        catalog = self.lab.get_experiment_catalog()
        
        self.assertEqual(catalog["catalog"], "Science Lab Experiments")
        self.assertIn("categories", catalog)
        self.assertIn("Cryptography", catalog["categories"])
        self.assertIn("Network Security", catalog["categories"])
        self.assertIn("Authentication Protocols", catalog["categories"])
        self.assertIn("Security Metrics", catalog["categories"])
        
    def test_crypto_experiments(self):
        """Test cryptography experiments are accessible."""
        hashing = self.lab.crypto.demo_hashing()
        self.assertIn("algorithms", hashing)
        
        pbkdf2 = self.lab.crypto.demo_pbkdf2()
        self.assertIn("derived_key", pbkdf2)
        
    def test_network_experiments(self):
        """Test network security experiments are accessible."""
        firewall = self.lab.network.analyze_firewall_config(["G1"], ["C1"])
        self.assertIn("guardians", firewall)
        
        scenarios = self.lab.network.simulate_intrusion_scenarios()
        self.assertIn("scenarios", scenarios)
        
    def test_auth_experiments(self):
        """Test authentication experiments are accessible."""
        two_fa = self.lab.auth.explain_two_factor_auth()
        self.assertIn("factors", two_fa)
        
        sessions = self.lab.auth.explain_session_management()
        self.assertIn("token_lifecycle", sessions)
        
    def test_metrics_experiments(self):
        """Test metrics experiments are accessible."""
        entropy = self.lab.metrics.calculate_entropy("test123")
        self.assertIn("entropy_bits", entropy)
        
        perf = self.lab.metrics.compare_hashing_performance()
        self.assertIn("results", perf)


if __name__ == "__main__":
    unittest.main()
