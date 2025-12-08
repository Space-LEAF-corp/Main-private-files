# Jarvondis 3.0 - Unified Sovereign Control Plane

Jarvondis 3.0 is a comprehensive sovereign control plane that combines cryptographic security, multi-factor authentication, consensus mechanisms, and ceremonial governance.

## Features

- **Tamper-Proofing**: File hashes, HMAC signatures, and RBAC
- **Global Command Surface**: Unified interface for all system operations
- **Sovereign Interface**: Single gateway with comprehensive audit logging
- **Emergency Access**: Administrator-controlled lifecycle with provenance logging
- **Ceremonial Milestones**: Living audit log with ceremonial banners
- **Cryptographic Security**: Ed25519 signatures for administrators and validators
- **Multi-Factor Authentication**: TOTP-based MFA for sensitive operations
- **Consensus Mechanism**: Multi-party approval with administrator fail-safe
- **Agent Harness**: Constrained execution wrapper for external agents
- **Provenance System**: Hash-chained audit trail for tamper evidence
- **Compliance Export**: JSON snapshots for regulatory requirements

## Directory Structure

```
jarvondis3/
├── core/
│   ├── sovereign.py        # Core sovereign interface and components
│   ├── crypto_adapter.py   # Ed25519 cryptographic signatures
│   └── mfa_totp.py         # TOTP multi-factor authentication
├── gaming/
│   ├── engine.py           # Sovereign-governed game engine
│   ├── ai_agents.py        # Game agent harness
│   ├── rituals.py          # Ceremonial game events
│   └── assets/             # Game assets and resources
├── compliance/
│   └── exports/            # Compliance snapshot exports
└── docs/
    ├── JARVONDIS3_MANIFESTO.md  # Vision and principles
    └── POSITION_PAPER.md         # Technical design document
```

## Installation

```bash
pip install pynacl
```

## Quick Start

```python
from jarvondis3.core import Jarvondis3GlobalInterface

# Initialize with MFA secret (base32-encoded)
interface = Jarvondis3GlobalInterface(
    mfa_secret_b32="JBSWY3DPEHPK3PXP",
    consensus_threshold=0
)

# Execute a command (requires admin role)
result = interface.request(
    user_role="admin",
    command="/JARVONDIS_SYS_DIAGNOSTICS",
    admin_signature=interface.tamper.sign_command("/JARVONDIS_SYS_DIAGNOSTICS"),
    mfa_code="123456"  # Get from TOTP authenticator
)

print(result)

# View audit log
print(interface.view_audit_log())

# Export compliance snapshot
print(interface.export_compliance_snapshot())
```

## With Ed25519 Signatures

```python
from nacl.signing import SigningKey
import nacl.encoding

# Generate admin key pair
admin_key = SigningKey.generate()
admin_signing_key_b64 = admin_key.encode(encoder=nacl.encoding.Base64Encoder).decode()
admin_verify_key_b64 = admin_key.verify_key.encode(encoder=nacl.encoding.Base64Encoder).decode()

# Initialize with Ed25519
interface = Jarvondis3GlobalInterface(
    admin_signing_key_b64=admin_signing_key_b64,
    mfa_secret_b32="JBSWY3DPEHPK3PXP"
)

# Sign and execute command
command = "/JARVONDIS_SYS_DIAGNOSTICS"
signature = interface.crypto.admin_sign(command)

result = interface.request(
    user_role="admin",
    command=command,
    admin_signature=interface.tamper.sign_command(command),
    mfa_code="123456"
)
```

## Gaming Integration

```python
from jarvondis3.gaming.engine import JarvondisGameEngine
from jarvondis3.gaming.rituals import level_up_ritual

# Create game engine with sovereign control
engine = JarvondisGameEngine(interface)
engine.start()

# Execute game tick
engine.tick()

# Ceremonial events
print(level_up_ritual(5))
```

## Agent Harness

```python
from jarvondis3.core import Jarvondis3AgentHarness
from jarvondis3.gaming.ai_agents import GameAgent

# Create constrained agent
harness = Jarvondis3AgentHarness(interface, agent_id="agent_001", role="operator")
agent = GameAgent(harness, name="Bot Alpha")

# Agent attempts action (constrained by RBAC)
result = agent.act("/JARVONDIS_SYS_DIAGNOSTICS")
```

## Emergency Mode

```python
# Grant emergency access (requires admin signature + MFA)
admin_sig = interface.tamper.sign_command("EMERGENCY_GRANT")
result = interface.grant_emergency_access(admin_sig, mfa_code="123456")

# Operators can now access non-sensitive commands
# ...

# Revoke emergency access
admin_sig = interface.tamper.sign_command("EMERGENCY_REVOKE")
result = interface.revoke_emergency_access(admin_sig, mfa_code="123456")
```

## Consensus Mode

```python
# Initialize with consensus threshold
interface = Jarvondis3GlobalInterface(
    mfa_secret_b32="JBSWY3DPEHPK3PXP",
    consensus_threshold=2  # Requires 2 validator signatures
)

# Add trusted signers
interface.add_trusted_signer("validator_001")
interface.add_trusted_signer("validator_002")

# Provide consensus signatures
sig1 = interface.tamper.sign_command("CONSENSUS")
interface.provide_consensus("validator_001", sig1)

sig2 = interface.tamper.sign_command("CONSENSUS")
interface.provide_consensus("validator_002", sig2)

# Now sensitive commands can proceed (with admin signature + MFA)
```

## Documentation

- [Jarvondis 3.0 Manifesto](docs/JARVONDIS3_MANIFESTO.md) - Vision and principles
- [Position Paper](docs/POSITION_PAPER.md) - Technical design and security analysis

## Security

Jarvondis 3.0 implements multiple security layers:

1. **Cryptographic**: Ed25519 signatures for non-repudiation
2. **Multi-Factor**: TOTP-based MFA for sensitive operations
3. **Access Control**: RBAC with administrator supremacy
4. **Consensus**: Multi-party approval for critical decisions
5. **Provenance**: Hash-chained audit trail for tamper evidence
6. **Agent Constraints**: Prevents privilege escalation

See [POSITION_PAPER.md](docs/POSITION_PAPER.md) for detailed security analysis.

## License

See repository LICENSE file.
