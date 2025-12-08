#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for Jarvondis 3.0 components"""

import unittest
import tempfile
import os
from jarvondis3.core import (
    Jarvondis3GlobalInterface,
    Jarvondis3AgentHarness,
    milestone_banner
)
from jarvondis3.gaming.engine import JarvondisGameEngine
from jarvondis3.gaming.ai_agents import GameAgent
from jarvondis3.gaming.rituals import level_up_ritual, boss_defeated_ritual


class TestJarvondis3Core(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures"""
        self.interface = Jarvondis3GlobalInterface(
            mfa_secret_b32="JBSWY3DPEHPK3PXP",
            consensus_threshold=0,
            admin_key=b"test_admin_key_12345"
        )

    def test_milestone_banner(self):
        """Test milestone banner generation"""
        banner = milestone_banner("COMMAND_EXECUTED")
        self.assertIn("MILESTONE", banner)
        self.assertIn("COMMAND_EXECUTED", banner)

    def test_list_commands(self):
        """Test command listing"""
        commands = self.interface.commands.list_commands()
        self.assertIn("/JARVONDIS_SYS_DIAGNOSTICS", commands)
        self.assertIn("/JARVONDIS_SYS_SHUTDOWN", commands)
        self.assertGreater(len(commands), 0)

    def test_admin_command_with_signature(self):
        """Test admin command execution with signature"""
        command = "/JARVONDIS_SYS_DIAGNOSTICS"
        signature = self.interface.tamper.sign_command(command)
        
        result = self.interface.request(
            user_role="admin",
            command=command,
            admin_signature=signature,
            mfa_code="JBSWY3DPEHPK3PXP"  # Use secret as code for testing
        )
        
        self.assertIn("DIAGNOSTICS", result)

    def test_operator_access_denied_without_emergency(self):
        """Test that operators are denied access without emergency mode"""
        command = "/JARVONDIS_SYS_DIAGNOSTICS"
        
        result = self.interface.request(
            user_role="operator",
            command=command
        )
        
        self.assertIn("ACCESS_DENIED", result)

    def test_emergency_mode_lifecycle(self):
        """Test emergency mode activation and revocation"""
        # Grant emergency access
        grant_sig = self.interface.tamper.sign_command("EMERGENCY_GRANT")
        result = self.interface.grant_emergency_access(grant_sig)
        self.assertIn("EMERGENCY_ACTIVATED", result)
        self.assertTrue(self.interface.emergency_mode)
        
        # Revoke emergency access
        revoke_sig = self.interface.tamper.sign_command("EMERGENCY_REVOKE")
        result = self.interface.revoke_emergency_access(revoke_sig)
        self.assertIn("EMERGENCY_REVOKED", result)
        self.assertFalse(self.interface.emergency_mode)

    def test_provenance_logging(self):
        """Test provenance hash chaining"""
        initial_hash = self.interface.last_hash
        self.assertEqual(initial_hash, "GENESIS")
        
        # Execute a command to trigger provenance
        command = "/JARVONDIS_SYS_UPDATE"
        signature = self.interface.tamper.sign_command(command)
        self.interface.request(
            user_role="admin",
            command=command,
            admin_signature=signature,
            mfa_code="JBSWY3DPEHPK3PXP"
        )
        
        # Check that hash changed
        self.assertNotEqual(self.interface.last_hash, "GENESIS")

    def test_audit_log_export(self):
        """Test audit log viewing"""
        command = "/JARVONDIS_SYS_BACKUP"
        signature = self.interface.tamper.sign_command(command)
        self.interface.request(
            user_role="admin",
            command=command,
            admin_signature=signature,
            mfa_code="JBSWY3DPEHPK3PXP"
        )
        
        audit_log = self.interface.view_audit_log()
        self.assertIn("AUDIT_HEADER", audit_log)
        self.assertIn("BACKUP", audit_log)

    def test_compliance_snapshot(self):
        """Test compliance snapshot export"""
        snapshot = self.interface.export_compliance_snapshot()
        self.assertIn("Jarvondis 3.0 Compliance Snapshot", snapshot)
        self.assertIn("audit_log", snapshot)
        self.assertIn("emergency_mode", snapshot)

    def test_consensus_management(self):
        """Test consensus signer management"""
        self.interface.add_trusted_signer("validator_001")
        self.assertIn("validator_001", self.interface.trusted_signers)
        
        # Provide consensus
        sig = self.interface.tamper.sign_command("CONSENSUS")
        result = self.interface.provide_consensus("validator_001", sig)
        self.assertIn("CONSENSUS_ADDED", result)
        
        # Clear consensus
        self.interface.clear_consensus()
        self.assertEqual(len(self.interface.consensus_signatures), 0)


class TestJarvondis3Gaming(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures"""
        self.interface = Jarvondis3GlobalInterface(
            mfa_secret_b32="JBSWY3DPEHPK3PXP",
            admin_key=b"test_admin_key_12345"
        )
        self.engine = JarvondisGameEngine(self.interface)

    def test_game_engine_start(self):
        """Test game engine startup"""
        result = self.engine.start()
        self.assertIn("Game engine started", result)
        self.assertTrue(self.engine.state["running"])

    def test_game_engine_tick(self):
        """Test game tick execution"""
        self.engine.start()
        result = self.engine.tick()
        self.assertIn("tick", result.lower())
        self.assertEqual(self.engine.state["tick"], 1)

    def test_game_engine_not_running(self):
        """Test tick when engine not running"""
        result = self.engine.tick()
        self.assertIn("not running", result)

    def test_level_up_ritual(self):
        """Test level up ritual"""
        result = level_up_ritual(10)
        self.assertIn("level 10", result)
        self.assertIn("MILESTONE", result)

    def test_boss_defeated_ritual(self):
        """Test boss defeated ritual"""
        result = boss_defeated_ritual("Dragon Lord")
        self.assertIn("Dragon Lord", result)
        self.assertIn("defeated", result)


class TestJarvondis3AgentHarness(unittest.TestCase):
    def setUp(self):
        """Set up test fixtures"""
        self.interface = Jarvondis3GlobalInterface(
            mfa_secret_b32="JBSWY3DPEHPK3PXP",
            admin_key=b"test_admin_key_12345"
        )
        self.harness = Jarvondis3AgentHarness(
            self.interface,
            agent_id="test_agent",
            role="operator"
        )
        self.agent = GameAgent(self.harness, "TestBot")

    def test_agent_creation(self):
        """Test agent harness creation"""
        self.assertEqual(self.harness.agent_id, "test_agent")
        self.assertEqual(self.harness.role, "operator")
        self.assertEqual(self.agent.name, "TestBot")

    def test_agent_denied_without_emergency(self):
        """Test agent access denied without emergency mode"""
        result = self.agent.act("/JARVONDIS_SYS_DIAGNOSTICS")
        self.assertIn("ACCESS_DENIED", result)

    def test_agent_access_with_emergency(self):
        """Test agent can access non-sensitive commands in emergency mode"""
        # Grant emergency access
        grant_sig = self.interface.tamper.sign_command("EMERGENCY_GRANT")
        self.interface.grant_emergency_access(grant_sig)
        
        # Agent can now access non-sensitive commands
        result = self.agent.act("/JARVONDIS_SYS_DIAGNOSTICS")
        self.assertIn("DIAGNOSTICS", result)


class TestJarvondis3Crypto(unittest.TestCase):
    def test_with_ed25519(self):
        """Test Ed25519 signature integration"""
        try:
            from nacl.signing import SigningKey
            import nacl.encoding
            
            # Generate admin key pair
            admin_key = SigningKey.generate()
            admin_signing_key_b64 = admin_key.encode(encoder=nacl.encoding.Base64Encoder).decode()
            
            # Initialize with Ed25519
            interface = Jarvondis3GlobalInterface(
                admin_signing_key_b64=admin_signing_key_b64,
                mfa_secret_b32="JBSWY3DPEHPK3PXP"
            )
            
            # Verify crypto adapter is initialized
            self.assertIsNotNone(interface.crypto)
            self.assertIsNotNone(interface.crypto.admin_signing_key)
            self.assertIsNotNone(interface.crypto.admin_verify_key)
            
        except ImportError:
            self.skipTest("PyNaCl not installed")

    def test_totp_initialization(self):
        """Test TOTP initialization"""
        interface = Jarvondis3GlobalInterface(
            mfa_secret_b32="JBSWY3DPEHPK3PXP"
        )
        
        self.assertIsNotNone(interface.totp)
        
        # Test TOTP code generation
        code = interface.totp.now()
        self.assertEqual(len(code), 6)
        self.assertTrue(code.isdigit())


class TestJarvondis3TamperProofing(unittest.TestCase):
    def test_file_integrity(self):
        """Test file integrity monitoring"""
        interface = Jarvondis3GlobalInterface(
            admin_key=b"test_key"
        )
        
        # Create a temporary file
        fd, path = tempfile.mkstemp(suffix=".txt")
        try:
            with os.fdopen(fd, 'w') as f:
                f.write("test content")
            
            # Register and verify
            interface.tamper.register_file(path)
            self.assertTrue(interface.tamper.verify_file(path))
            
            # Modify file
            with open(path, 'w') as f:
                f.write("modified content")
            
            # Verification should fail
            self.assertFalse(interface.tamper.verify_file(path))
            
        finally:
            try:
                os.remove(path)
            except Exception:
                pass

    def test_command_signing(self):
        """Test command signing and verification"""
        interface = Jarvondis3GlobalInterface(
            admin_key=b"test_key"
        )
        
        command = "/JARVONDIS_SYS_DIAGNOSTICS"
        signature = interface.tamper.sign_command(command)
        
        self.assertTrue(interface.tamper.verify_signature(command, signature))
        self.assertFalse(interface.tamper.verify_signature(command, "invalid_sig"))

    def test_rbac_hierarchy(self):
        """Test RBAC role hierarchy"""
        interface = Jarvondis3GlobalInterface(
            admin_key=b"test_key"
        )
        
        # Admin can access admin-level commands
        self.assertTrue(interface.tamper.check_access("admin", "admin"))
        self.assertTrue(interface.tamper.check_access("admin", "operator"))
        self.assertTrue(interface.tamper.check_access("admin", "user"))
        
        # Operator can access operator-level but not admin
        self.assertFalse(interface.tamper.check_access("operator", "admin"))
        self.assertTrue(interface.tamper.check_access("operator", "operator"))
        self.assertTrue(interface.tamper.check_access("operator", "user"))
        
        # User can only access user-level
        self.assertFalse(interface.tamper.check_access("user", "admin"))
        self.assertFalse(interface.tamper.check_access("user", "operator"))
        self.assertTrue(interface.tamper.check_access("user", "user"))


if __name__ == "__main__":
    unittest.main()
