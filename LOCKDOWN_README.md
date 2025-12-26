# Jarvondis Administrative Access Lockdown

## Overview

The Jarvondis Lockdown Protocol enforces sovereign authorship, owner-only control, and privacy lockdown for the Jarvondis system.

## Components

### Core Files

1. **jarvondis_lockdown.py** - Main lockdown implementation
   - `LockdownPolicy`: Configuration dataclass
   - `AdministrativeLockdown`: Main enforcement class with HMAC signature verification, allowlist checking, and AI blocking

2. **jarvondis_policy.json** - Configuration file
   - **IMPORTANT**: This file is gitignored to prevent secret leakage. Copy `jarvondis_policy.json.example` to `jarvondis_policy.json` and set `JARVONDIS_ADMIN_SECRET` environment variable.
   - Configure `owner_id`, `allowlist_clients`, and other settings

3. **main.py** - Example usage demonstrating policy loading and guard usage

4. **sign_owner_command.py** - Utility for generating signed owner commands

## Current Allowlist

The following clients have trusted access:
- `captains-log`
- `eternal-chalkboard`
- `leif-personal-device`
- `dimitri` ✅ (Restored to allowlist)

**Blocked clients:**
- `miko` ❌ (Not in allowlist)
- All AI agents detected via user-agent fingerprints

## Usage

### Basic Usage

```python
import json
import os
from jarvondis_lockdown import LockdownPolicy, AdministrativeLockdown

# Load secret from environment (recommended)
admin_secret = os.environ.get("JARVONDIS_ADMIN_SECRET")
if not admin_secret:
    raise ValueError("JARVONDIS_ADMIN_SECRET environment variable must be set")

# Load policy
with open("jarvondis_policy.json", "r") as f:
    cfg = json.load(f)
policy = LockdownPolicy(
    owner_id=cfg["owner_id"],
    admin_secret=admin_secret,
    lockdown_active=cfg.get("lockdown_active", True),
    allowlist_clients=set(cfg.get("allowlist_clients", []))
)

# Create guard
guard = AdministrativeLockdown(policy)

# Check access
result = guard.respond("dimitri", "Dimitri/1.0", "Hello")
print(result)
```

### Owner Commands

To change lockdown settings, use signed commands:

```python
import time
from sign_owner_command import make_signature

# Generate signature
sig, ts = make_signature("leif.w.sogge", "YOUR_SECRET", "set_lockdown:False")

# Apply command
guard.set_lockdown("leif.w.sogge", sig, ts, False)
```

## Security Features

1. **HMAC Signature Verification**: Owner commands require valid HMAC-SHA256 signatures
2. **Timestamp Validation**: Signatures expire after 60 seconds (configurable)
3. **Allowlist Enforcement**: During lockdown, only allowlisted clients can access
4. **AI Agent Blocking**: Heuristic detection of AI agents via user-agent strings
5. **Audit Logging**: All access attempts and policy changes are logged

## Testing

Run the test suite:

```bash
python -m unittest tests.test_jarvondis_lockdown -v
```

All 12 tests should pass, validating:
- Dimitri is allowed access
- Miko is blocked
- Signature verification works correctly
- AI agents are blocked
- Audit logging functions properly

## Protocol Update

As specified in the issue:
- ✅ **Dimitri is restored to the allowlist** (trusted access permitted)
- ❌ **Miko remains blocked** (no access granted)
- ✅ All other privacy clauses remain intact
