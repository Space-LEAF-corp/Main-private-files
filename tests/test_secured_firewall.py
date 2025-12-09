"""Unit tests for SecuredFirewall integration."""
import os
import tempfile
import unittest

from secured_firewall import SecuredFirewall


class TestSecuredFirewall(unittest.TestCase):
    def setUp(self):
        fd, self.auth_file = tempfile.mkstemp(suffix=".json")
        os.close(fd)
        try:
            os.remove(self.auth_file)
        except Exception:
            pass
        self.sf = SecuredFirewall(users_file=self.auth_file)

    def tearDown(self):
        try:
            os.remove(self.auth_file)
        except Exception:
            pass

    def test_full_auth_and_access_flow(self):
        # Register
        reg = self.sf.register_user("bob", "securePass123", "LINEAGE_SAFE_BOB_001")
        self.assertEqual(reg.get("status"), "ok")

        # Login step 1
        login1 = self.sf.login_step1("bob", "securePass123", "LINEAGE_SAFE_BOB_001")
        self.assertEqual(login1.get("status"), "challenge")
        otp = login1.get("otp_token")
        self.assertIsNotNone(otp)

        # Login step 2
        login2 = self.sf.login_step2("bob", otp)
        self.assertEqual(login2.get("status"), "ok")
        session = login2.get("session_token")
        self.assertIsNotNone(session)

        # Access firewall
        access = self.sf.access_firewall(session)
        self.assertEqual(access.get("status"), "ok")

    def test_intrusion_attempt(self):
        # First 3 attempts give warnings
        res1 = self.sf.intrude_attempt("hacker", "malicious_payload")
        self.assertEqual(res1.get("status"), "trapped")
        self.assertIn("WARNING 1", res1.get("message"))
        
        res2 = self.sf.intrude_attempt("hacker", "malicious_payload")
        self.assertEqual(res2.get("status"), "trapped")
        self.assertIn("WARNING 2", res2.get("message"))
        
        res3 = self.sf.intrude_attempt("hacker", "malicious_payload")
        self.assertEqual(res3.get("status"), "trapped")
        self.assertIn("WARNING 3", res3.get("message"))
        
        # Fourth attempt traps in mirror layer
        res4 = self.sf.intrude_attempt("hacker", "malicious_payload")
        self.assertEqual(res4.get("status"), "trapped")
        self.assertIn("trapped in mirror layer", res4.get("message"))


if __name__ == "__main__":
    unittest.main()
