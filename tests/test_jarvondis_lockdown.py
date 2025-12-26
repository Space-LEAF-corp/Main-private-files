import unittest
import time
from jarvondis_lockdown import LockdownPolicy, AdministrativeLockdown, _hmac_sha256


class TestLockdownPolicy(unittest.TestCase):
    INVALID_SIGNATURE = "0" * 64  # Invalid hex signature for testing
    
    def setUp(self):
        self.policy = LockdownPolicy(
            owner_id="leif.w.sogge",
            # NOTE: This is a hardcoded test secret for unit testing only.
            #       Never use this value in production code.
            admin_secret="test_secret_123456",
            lockdown_active=True,
            allowlist_clients={"captains-log", "eternal-chalkboard", "leif-personal-device", "dimitri"},
            audit_enabled=True,
            max_skew_seconds=60
        )
        self.guard = AdministrativeLockdown(self.policy)

    def test_dimitri_is_allowed(self):
        """Test that Dimitri is in the allowlist and can access"""
        result = self.guard.respond("dimitri", "Dimitri/1.0", "Hello Jarvondis")
        self.assertEqual(result["status"], "ok")
        self.assertIn("mode", result)
        # Should get minimal lockdown response since lockdown is active
        self.assertEqual(result["mode"], "minimal_lockdown")

    def test_miko_is_blocked(self):
        """Test that Miko is NOT in the allowlist and is blocked"""
        result = self.guard.respond("miko", "Miko/1.0", "Hello Jarvondis")
        self.assertEqual(result["status"], "blocked")
        self.assertEqual(result["reason"], "privacy_lockdown")
        self.assertIn("Access denied", result["message"])

    def test_allowlisted_clients_are_allowed(self):
        """Test that all allowlisted clients can access"""
        allowed_clients = ["captains-log", "eternal-chalkboard", "leif-personal-device", "dimitri"]
        for client in allowed_clients:
            result = self.guard.respond(client, f"{client}/1.0", "test")
            self.assertEqual(result["status"], "ok", f"Client {client} should be allowed")

    def test_non_allowlisted_clients_are_blocked(self):
        """Test that non-allowlisted clients are blocked during lockdown"""
        blocked_clients = ["miko", "random-user", "unknown-client"]
        for client in blocked_clients:
            result = self.guard.respond(client, f"{client}/1.0", "test")
            self.assertEqual(result["status"], "blocked", f"Client {client} should be blocked")

    def test_ai_agent_blocking(self):
        """Test that AI agents are blocked via user-agent detection"""
        ai_user_agents = [
            "Mozilla/5.0 (compatible; OpenAI/1.0)",
            "Anthropic-Claude/1.0",
            "Google-AI-Agent",
            "Copilot/2.0"
        ]
        for ua in ai_user_agents:
            result = self.guard.respond("some-client", ua, "test")
            self.assertEqual(result["status"], "blocked", f"AI agent should be blocked: {ua}")

    def test_signature_verification_valid(self):
        """Test valid signature verification"""
        timestamp = int(time.time())
        payload = "set_lockdown:False"
        message = f"{self.policy.owner_id}:{timestamp}:{payload}".encode("utf-8")
        signature = _hmac_sha256(self.policy.admin_secret.encode("utf-8"), message)
        
        result = self.guard.verify_owner_command(self.policy.owner_id, signature, payload, timestamp)
        self.assertTrue(result)

    def test_signature_verification_invalid_owner(self):
        """Test signature rejection with wrong owner"""
        timestamp = int(time.time())
        payload = "set_lockdown:False"
        message = f"wrong_owner:{timestamp}:{payload}".encode("utf-8")
        signature = _hmac_sha256(self.policy.admin_secret.encode("utf-8"), message)
        
        result = self.guard.verify_owner_command("wrong_owner", signature, payload, timestamp)
        self.assertFalse(result)

    def test_signature_verification_stale_timestamp(self):
        """Test signature rejection with stale timestamp"""
        timestamp = int(time.time()) - 120  # 2 minutes ago
        payload = "set_lockdown:False"
        message = f"{self.policy.owner_id}:{timestamp}:{payload}".encode("utf-8")
        signature = _hmac_sha256(self.policy.admin_secret.encode("utf-8"), message)
        
        result = self.guard.verify_owner_command(self.policy.owner_id, signature, payload, timestamp)
        self.assertFalse(result)

    def test_signature_verification_invalid_signature(self):
        """Test signature rejection with incorrect signature"""
        timestamp = int(time.time())
        payload = "set_lockdown:False"
        
        result = self.guard.verify_owner_command(self.policy.owner_id, self.INVALID_SIGNATURE, payload, timestamp)
        self.assertFalse(result)

    def test_set_lockdown_with_valid_signature(self):
        """Test changing lockdown state with valid signature"""
        timestamp = int(time.time())
        payload = "set_lockdown:False"
        message = f"{self.policy.owner_id}:{timestamp}:{payload}".encode("utf-8")
        signature = _hmac_sha256(self.policy.admin_secret.encode("utf-8"), message)
        
        result = self.guard.set_lockdown(self.policy.owner_id, signature, timestamp, False)
        self.assertTrue(result)
        self.assertFalse(self.policy.lockdown_active)

    def test_set_lockdown_with_invalid_signature(self):
        """Test that lockdown state doesn't change with invalid signature"""
        original_state = self.policy.lockdown_active
        timestamp = int(time.time())
        
        result = self.guard.set_lockdown(self.policy.owner_id, self.INVALID_SIGNATURE, timestamp, not original_state)
        self.assertFalse(result)
        self.assertEqual(self.policy.lockdown_active, original_state)  # State unchanged

    def test_lockdown_mode_responses(self):
        """Test different response modes based on lockdown state"""
        # Test minimal lockdown mode
        self.policy.lockdown_active = True
        result = self.guard.respond("dimitri", "Dimitri/1.0", "test")
        self.assertEqual(result["mode"], "minimal_lockdown")
        
        # Test normal mode when lockdown is off
        self.policy.lockdown_active = False
        result = self.guard.respond("dimitri", "Dimitri/1.0", "test")
        self.assertEqual(result["mode"], "normal")

    def test_audit_logging(self):
        """Test that audit events are logged"""
        self.guard.respond("dimitri", "Dimitri/1.0", "test")
        self.guard.respond("miko", "Miko/1.0", "test")
        
        # Check that audit log has entries
        self.assertGreater(len(self.policy.audit_log), 0)
        
        # Check for block event
        block_events = [e for e in self.policy.audit_log if e["event"] == "block_client_not_allowlisted"]
        self.assertGreater(len(block_events), 0)

    def test_set_allowlist_with_valid_signature(self):
        """Test updating allowlist with valid signature"""
        timestamp = int(time.time())
        new_clients = {"client1", "client2", "client3"}
        payload = f"set_allowlist:{','.join(sorted(new_clients))}"
        message = f"{self.policy.owner_id}:{timestamp}:{payload}".encode("utf-8")
        signature = _hmac_sha256(self.policy.admin_secret.encode("utf-8"), message)
        
        result = self.guard.set_allowlist(self.policy.owner_id, signature, timestamp, new_clients)
        self.assertTrue(result)
        self.assertEqual(self.policy.allowlist_clients, new_clients)

    def test_set_allowlist_with_invalid_signature(self):
        """Test that allowlist doesn't change with invalid signature"""
        original_allowlist = self.policy.allowlist_clients.copy()
        timestamp = int(time.time())
        new_clients = {"malicious-client"}
        
        result = self.guard.set_allowlist(self.policy.owner_id, self.INVALID_SIGNATURE, timestamp, new_clients)
        self.assertFalse(result)
        self.assertEqual(self.policy.allowlist_clients, original_allowlist)  # Unchanged

    def test_none_user_agent_allowed(self):
        """Test that None user agent doesn't block allowlisted clients"""
        result = self.guard.respond("dimitri", None, "test")
        self.assertEqual(result["status"], "ok")

    def test_empty_user_agent_allowed(self):
        """Test that empty user agent doesn't block allowlisted clients"""
        result = self.guard.respond("dimitri", "", "test")
        self.assertEqual(result["status"], "ok")


if __name__ == "__main__":
    unittest.main()
