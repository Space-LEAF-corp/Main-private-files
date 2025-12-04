import os
import tempfile
import unittest

from auth import AuthManager


class TestAuthManager(unittest.TestCase):
    def setUp(self):
        fd, self.path = tempfile.mkstemp(suffix=".json")
        os.close(fd)
        try:
            os.remove(self.path)
        except Exception:
            pass
        self.am = AuthManager(users_file=self.path)

    def tearDown(self):
        try:
            os.remove(self.path)
        except Exception:
            pass

    def test_register_and_login_flow(self):
        user = "alice"
        pwd = "s3cureP@ssw0rd"
        dna = "LINEAGE_SAFE_ALICE_001"
        self.assertTrue(self.am.register(user, pwd, dna))

        # attempt login with correct password and dna -> challenge with OTP
        res = self.am.start_login(user, pwd, dna_code=dna)
        self.assertEqual(res.get("status"), "challenge")
        otp = res.get("otp")
        self.assertIsNotNone(otp)

        # complete login with OTP
        res2 = self.am.complete_login(user, otp)
        self.assertEqual(res2.get("status"), "ok")
        self.assertIn("session_token", res2)

    def test_wrong_password(self):
        user = "bob"
        pwd = "anotherPass1"
        dna = "LINEAGE_SAFE_BOB_001"
        self.am.register(user, pwd, dna)
        res = self.am.start_login(user, "badpass", dna_code=dna)
        self.assertEqual(res.get("status"), "denied")

    def test_dna_bind_on_first_login(self):
        user = "carol"
        pwd = "strongpassword"
        dna = "LINEAGE_SAFE_CAROL_001"
        # register by creating user record directly by mimicking register
        self.am.register(user, pwd, dna)
        # start login without dna should be denied because dna required (already bound)
        res = self.am.start_login(user, pwd, dna_code=None)
        self.assertEqual(res.get("status"), "denied")


if __name__ == "__main__":
    unittest.main()
