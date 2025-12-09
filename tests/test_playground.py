"""Unit tests for Playground module."""
import os
import tempfile
import unittest

from playground import Playground, PlaygroundSession


class TestPlaygroundSession(unittest.TestCase):
    """Tests for PlaygroundSession class."""
    
    def test_session_creation(self):
        """Test creating a playground session."""
        session = PlaygroundSession("test_session_1")
        self.assertEqual(session.session_id, "test_session_1")
        self.assertEqual(len(session.commands_history), 0)
        
    def test_log_command(self):
        """Test logging commands in session."""
        session = PlaygroundSession("test_session_2")
        session.log_command("test_command", {"result": "success"})
        
        history = session.get_history()
        self.assertEqual(len(history), 1)
        self.assertEqual(history[0]["command"], "test_command")
        self.assertEqual(history[0]["result"]["result"], "success")


class TestPlayground(unittest.TestCase):
    """Tests for Playground class."""
    
    def setUp(self):
        """Set up test playground."""
        self.playground = Playground(isolated=True)
        
    def tearDown(self):
        """Clean up playground."""
        self.playground.cleanup()
        
    def test_playground_initialization(self):
        """Test playground initializes correctly."""
        self.assertIsNotNone(self.playground.secured_firewall)
        self.assertTrue(self.playground.isolated)
        
    def test_session_management(self):
        """Test creating and retrieving sessions."""
        session = self.playground.create_session("session1")
        self.assertIsNotNone(session)
        self.assertEqual(session.session_id, "session1")
        
        retrieved = self.playground.get_session("session1")
        self.assertEqual(retrieved.session_id, "session1")
        
        non_existent = self.playground.get_session("nonexistent")
        self.assertIsNone(non_existent)
        
    def test_demo_registration(self):
        """Test registration demo."""
        result = self.playground.demo_registration("test_user_1")
        
        self.assertEqual(result["demo"], "registration")
        self.assertIn("user_id", result)
        self.assertIn("password", result)
        self.assertIn("dna_code", result)
        self.assertIn("result", result)
        self.assertEqual(result["result"]["status"], "ok")
        
    def test_demo_login(self):
        """Test login demo."""
        # First register
        self.playground.demo_registration("test_user_2")
        
        # Then login
        result = self.playground.demo_login("test_user_2")
        
        self.assertEqual(result["demo"], "login")
        self.assertIn("step1_status", result)
        self.assertIn("step2_status", result)
        self.assertEqual(result["step2_status"], "ok")
        self.assertIsNotNone(result.get("session_token"))
        
    def test_demo_firewall_access(self):
        """Test firewall access demo."""
        # Register first
        self.playground.demo_registration("test_user_3")
        
        # Access firewall
        result = self.playground.demo_firewall_access("test_user_3")
        
        self.assertEqual(result["demo"], "firewall_access")
        self.assertIn("access_result", result)
        self.assertEqual(result["firewall_status"], "ok")
        
    def test_demo_intrusion_detection(self):
        """Test intrusion detection demo."""
        result = self.playground.demo_intrusion_detection()
        
        self.assertEqual(result["demo"], "intrusion_detection")
        self.assertIn("result", result)
        self.assertEqual(result["result"]["status"], "trapped")
        
    def test_tutorial_flow(self):
        """Test complete tutorial flow."""
        # Start tutorial
        start = self.playground.start_tutorial()
        self.assertEqual(start["step"], 1)
        self.assertEqual(self.playground.tutorial_step, 1)
        
        # Step 1: Registration
        step1 = self.playground.tutorial_next()
        self.assertEqual(step1["demo"], "registration")
        self.assertEqual(self.playground.tutorial_step, 2)
        
        # Step 2: Login
        step2 = self.playground.tutorial_next()
        self.assertEqual(step2["demo"], "login")
        self.assertEqual(self.playground.tutorial_step, 3)
        
        # Step 3: Firewall Access
        step3 = self.playground.tutorial_next()
        self.assertEqual(step3["demo"], "firewall_access")
        self.assertEqual(self.playground.tutorial_step, 4)
        
        # Step 4: Intrusion Detection
        step4 = self.playground.tutorial_next()
        self.assertEqual(step4["demo"], "intrusion_detection")
        self.assertEqual(self.playground.tutorial_step, 0)
        self.assertIn("congratulations", step4)
        
    def test_tutorial_not_started(self):
        """Test tutorial_next without starting tutorial."""
        result = self.playground.tutorial_next()
        self.assertIn("error", result)
        
    def test_custom_experiment(self):
        """Test custom experiment."""
        result = self.playground.custom_experiment(
            user_id="custom_user",
            password="CustomPass123!",
            dna_code="LINEAGE_SAFE_CUSTOM_001",
            test_intrusion=True
        )
        
        self.assertEqual(result["experiment"], "custom")
        self.assertIn("steps", result)
        self.assertGreater(len(result["steps"]), 0)
        
        # Check all steps completed
        actions = [step["action"] for step in result["steps"]]
        self.assertIn("register", actions)
        self.assertIn("login_step1", actions)
        self.assertIn("login_step2", actions)
        self.assertIn("firewall_access", actions)
        self.assertIn("intrusion_test", actions)
        
    def test_custom_experiment_without_intrusion(self):
        """Test custom experiment without intrusion test."""
        result = self.playground.custom_experiment(
            user_id="custom_user_2",
            password="CustomPass456!",
            dna_code="LINEAGE_SAFE_CUSTOM2_001",
            test_intrusion=False
        )
        
        actions = [step["action"] for step in result["steps"]]
        self.assertNotIn("intrusion_test", actions)
        
    def test_get_help(self):
        """Test help information."""
        help_info = self.playground.get_help()
        
        self.assertIn("playground_help", help_info)
        self.assertIn("demos", help_info["playground_help"])
        self.assertIn("tutorial", help_info["playground_help"])
        self.assertIn("custom", help_info["playground_help"])
        self.assertIn("utility", help_info["playground_help"])
        
    def test_isolated_mode(self):
        """Test that isolated mode creates temporary files."""
        playground = Playground(isolated=True)
        # File is created when first used, not at initialization
        # But we can verify it's a temp file path
        self.assertTrue(playground.auth_file.startswith("/tmp") or "playground_" in playground.auth_file)
        # Trigger file creation by registering a user
        playground.demo_registration("temp_user")
        self.assertTrue(os.path.exists(playground.auth_file))
        playground.cleanup()
        self.assertFalse(os.path.exists(playground.auth_file))
        
    def test_non_isolated_mode(self):
        """Test non-isolated mode uses named file."""
        playground = Playground(isolated=False)
        self.assertEqual(playground.auth_file, "playground_users.json")
        # Clean up any created file
        try:
            os.remove("playground_users.json")
        except Exception:
            pass


if __name__ == "__main__":
    unittest.main()
