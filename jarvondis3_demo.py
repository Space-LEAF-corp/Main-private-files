#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Jarvondis 3.0 Demonstration Script

This script demonstrates the key features of Jarvondis 3.0:
- Tamper-proofing and command signing
- Global command surface
- Emergency mode
- Provenance logging
- Gaming integration
- Agent harness
"""

from jarvondis3.core import Jarvondis3GlobalInterface, Jarvondis3AgentHarness
from jarvondis3.gaming.engine import JarvondisGameEngine
from jarvondis3.gaming.ai_agents import GameAgent
from jarvondis3.gaming.rituals import level_up_ritual, boss_defeated_ritual


def print_section(title):
    """Print a section header"""
    print(f"\n{'=' * 60}")
    print(f"  {title}")
    print(f"{'=' * 60}\n")


def main():
    print_section("Jarvondis 3.0 Demonstration")
    
    # Initialize the sovereign interface
    print("Initializing Jarvondis 3.0 Sovereign Interface...")
    interface = Jarvondis3GlobalInterface(
        mfa_secret_b32="JBSWY3DPEHPK3PXP",
        consensus_threshold=0,
        admin_key=b"demo_admin_key_12345"
    )
    print("✓ Sovereign Interface initialized")
    
    # List available commands
    print_section("Available Global Commands")
    commands = interface.commands.list_commands()
    for cmd, desc in commands.items():
        print(f"  {cmd}")
        print(f"    → {desc}")
    
    # Execute a command with proper authorization
    print_section("Executing Authorized Command (Non-Sensitive)")
    command = "/JARVONDIS_SYS_UPDATE"
    signature = interface.tamper.sign_command(command)
    print(f"Command: {command}")
    print(f"Signature: {signature[:32]}...")
    
    result = interface.request(
        user_role="admin",
        command=command,
        admin_signature=signature
    )
    print(f"Result: {result}")
    
    # Demonstrate access denial
    print_section("Access Control - Operator Without Emergency Mode")
    result = interface.request(
        user_role="operator",
        command="/JARVONDIS_SYS_UPDATE"
    )
    print(f"Result: {result[:100]}...")
    
    # Emergency mode demonstration
    print_section("Emergency Mode Lifecycle")
    print("Granting emergency access...")
    grant_sig = interface.tamper.sign_command("EMERGENCY_GRANT")
    mfa_code = interface.totp.now()
    print(f"Using TOTP MFA code: {mfa_code}")
    result = interface.grant_emergency_access(grant_sig, mfa_code)
    print(f"✓ {result[:80]}...")
    
    print("\nOperator can now access non-sensitive commands:")
    result = interface.request(
        user_role="operator",
        command="/JARVONDIS_SYS_UPDATE"
    )
    print(f"Result: {result[:100]}...")
    
    print("\nRevoking emergency access...")
    revoke_sig = interface.tamper.sign_command("EMERGENCY_REVOKE")
    mfa_code = interface.totp.now()
    result = interface.revoke_emergency_access(revoke_sig, mfa_code)
    print(f"✓ {result[:80]}...")
    
    # Provenance demonstration
    print_section("Provenance Hash Chain")
    print(f"Initial hash: {interface.last_hash}")
    
    command = "/JARVONDIS_SYS_UPDATE"
    signature = interface.tamper.sign_command(command)
    interface.request(
        user_role="admin",
        command=command,
        admin_signature=signature
    )
    
    print(f"Hash after command: {interface.last_hash}")
    print("✓ Provenance chain maintained")
    
    # Demonstrate sensitive command with MFA
    print_section("Sensitive Command with MFA")
    command = "/JARVONDIS_SYS_BACKUP"
    signature = interface.tamper.sign_command(command)
    mfa_code = interface.totp.now()
    print(f"Command: {command}")
    print(f"MFA Code: {mfa_code}")
    result = interface.request(
        user_role="admin",
        command=command,
        admin_signature=signature,
        mfa_code=mfa_code
    )
    print(f"Result: {result[:100]}...")
    
    # Gaming integration
    print_section("Gaming Engine Integration")
    engine = JarvondisGameEngine(interface)
    
    print("Starting game engine...")
    result = engine.start()
    print(f"✓ {result}")
    
    print("\nExecuting game ticks...")
    for i in range(3):
        result = engine.tick()
        print(f"  {result}")
    
    print("\nCeremonial rituals:")
    print(f"  {level_up_ritual(5)}")
    print(f"  {boss_defeated_ritual('Shadow Dragon')}")
    
    # Agent harness demonstration
    print_section("Agent Harness (Constrained Execution)")
    
    # Grant emergency access for agent demo
    grant_sig = interface.tamper.sign_command("EMERGENCY_GRANT")
    mfa_code = interface.totp.now()
    interface.grant_emergency_access(grant_sig, mfa_code)
    
    harness = Jarvondis3AgentHarness(interface, agent_id="demo_agent", role="operator")
    agent = GameAgent(harness, "BotAlpha")
    
    print(f"Agent: {agent.name}")
    print(f"Agent ID: {harness.agent_id}")
    print(f"Agent Role: {harness.role}")
    
    print("\nAgent attempting action:")
    result = agent.act("/JARVONDIS_SYS_UPDATE")
    print(f"Result: {result[:100]}...")
    
    # Compliance snapshot
    print_section("Compliance Snapshot Export")
    snapshot = interface.export_compliance_snapshot()
    print("Compliance snapshot generated:")
    print(snapshot[:500] + "...")
    
    # Audit log
    print_section("Audit Log")
    audit_log = interface.view_audit_log()
    lines = audit_log.split('\n')
    print(f"Total log entries: {len(lines)}")
    print(f"\nFirst 5 entries:")
    for line in lines[:5]:
        print(f"  {line[:100]}...")
    
    print_section("Demo Complete")
    print("✓ All Jarvondis 3.0 features demonstrated successfully!")
    print("\nKey Features Shown:")
    print("  • Tamper-proofing and command signing")
    print("  • Global command surface with RBAC")
    print("  • Emergency mode lifecycle")
    print("  • Provenance hash chaining")
    print("  • Gaming engine integration")
    print("  • Agent harness with constraints")
    print("  • Compliance snapshot export")
    print("  • Comprehensive audit logging")


if __name__ == "__main__":
    main()
