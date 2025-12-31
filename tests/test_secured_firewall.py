"""Unit tests for SecuredFirewall integration."""
import os
import tempfile




import unittest
from typing import Any, Dict


from typing import Any, Dict
from secured_firewall import SecuredFirewall

# Patch for type checking: add type hints if missing in SecuredFirewall
def _patch_sf_annotations():
    # Get unbound functions from the class, not bound methods from the instance
    funcs: list[tuple[str, dict[str, Any]]] = [
        ('register_user', {'user_id': str, 'password': str, 'dna_code': str, 'return': Dict[str, Any]}),
        ('login_step1', {'user_id': str, 'password': str, 'dna_code': str, 'return': Dict[str, Any]}),
        ('login_step2', {'user_id': str, 'otp_token': str, 'return': Dict[str, Any]}),
        ('access_firewall', {
            'session_token': str,
            'network': str,
            'metadata': 'dict[str, Any] | None',
            'return': Dict[str, Any]
        }),
        ('intrude_attempt', {'user_id': str, 'payload': str, 'return': Dict[str, Any]}),
    ]
    for fname, ann in funcs:
        func = getattr(SecuredFirewall, fname, None)
        if func is not None and not hasattr(func, '__annotations__'):
            func.__annotations__ = ann
        elif func is not None and not func.__annotations__:
            func.__annotations__ = ann
_patch_sf_annotations()


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

        from typing import Any, Dict, Callable
        register_user: Callable[[str, str, str], Dict[str, Any]] = self.sf.register_user  # type: ignore
        reg: Dict[str, Any] = register_user("bob", "securePass123", "LINEAGE_SAFE_BOB_001")  # type: ignore[typeddict-item]
        self.assertEqual(reg.get("status"), "ok")

        # Login step 1
        from typing import Callable
        login_step1: Callable[[str, str, str], Dict[str, Any]] = self.sf.login_step1  # type: ignore
        login1: Dict[str, Any] = login_step1("bob", "securePass123", "LINEAGE_SAFE_BOB_001")  # type: ignore
        self.assertEqual(login1.get("status"), "challenge")
        from typing import Optional
        otp: Optional[str] = login1.get("otp_token")
        self.assertIsNotNone(otp)
        assert isinstance(otp, str)  # type: ignore

        # Login step 2
        from typing import Dict, Any
        from typing import Dict, Any
        login2: Dict[str, Any] = self.sf.login_step2("bob", otp)  # type: ignore
        self.assertEqual(login2.get("status"), "ok")
        session = login2.get("session_token")
        self.assertIsNotNone(session)

        # Access firewall
        assert isinstance(session, str)  # Ensure session is str for type checker
        from typing import Dict, Any, Callable
        access_firewall: Callable[[str, str, dict[str, Any] | None], Dict[str, Any]] = self.sf.access_firewall  # type: ignore
        access: Dict[str, Any] = access_firewall(session_token=session, network="internet", metadata=None)  # type: ignore
        self.assertEqual(access.get("status"), "ok")

    def test_intrusion_attempt(self):
        # First 3 attempts give warnings
        res1: Dict[str, Any] = self.sf.intrude_attempt("hacker", "malicious_payload")  # type: ignore[call-arg]
        status1: str = res1.get("status", "")
        message1: str = res1.get("message", "")
        self.assertEqual(status1, "trapped")
        self.assertIn("WARNING 1", message1)

        res2: Dict[str, Any] = self.sf.intrude_attempt("hacker", "malicious_payload")  # type: ignore
        status2: str = res2.get("status", "")
        message2: str = res2.get("message", "")
        self.assertEqual(status2, "trapped")
        self.assertIn("WARNING 2", message2)

        res3: Dict[str, Any] = self.sf.intrude_attempt("hacker", "malicious_payload")  # type: ignore
        status3: str = res3.get("status", "")
        message3: str = res3.get("message", "")
        self.assertEqual(status3, "trapped")
        self.assertIn("WARNING 3", message3)

        # Fourth attempt traps in mirror layer
        res4: Dict[str, Any] = self.sf.intrude_attempt("hacker", "malicious_payload")  # type: ignore
        status4: str = res4.get("status", "")
        message4: str = res4.get("message", "")
        self.assertEqual(status4, "trapped")
        self.assertIn("trapped in mirror layer", message4)


if __name__ == "__main__":
    unittest.main()
