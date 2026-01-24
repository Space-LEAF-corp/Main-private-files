 updated license 1.0 as of 2026
Jarvondis (refactor)
====================

This workspace contains a refactored `Jarvondis` module and a small CLI runner.

Playground Quick Start
---------------------

- Install dependencies (optional, recommended):

```bash
python -m pip install -r requirements.txt
```

- Run interactively:

**Guided tutorials for learning authentication flows**
**Custom experiment support for advanced users**

Files added

- `jarvondis/jarvondis.py` — refactored class and persistence
- `auth/` — authentication module (QR-DNA, password, OTP multi-layer)
  - `auth/__init__.py` — exports AuthManager
  - `auth/manager.py` — core auth logic (PBKDF2 hashing, OTP, registration)
- `diamond_firewall.py` — lineage-safe firewall (thread-safe, mirror + diamond layers)
- `secured_firewall.py` — example integrating firewall + auth
- `cli.py` — small interactive CLI
- `tests/test_jarvondis.py` — unit tests (Jarvondis memory)
- `tests/test_auth.py` — unit tests (auth flows)
- `requirements.txt` — pandas, numpy

Authentication & Security
------------------------

**Auth module** (`auth/manager.py`)

- Supports registration with QR-DNA (prefix: `LINEAGE_SAFE`), password (min 8 chars), and user_id
- Multi-layer login: password → DNA bind → OTP challenge → session token
- Uses PBKDF2-HMAC-SHA256 for password hashing (no external crypto deps)
- OTP tokens expire after 5 minutes

**SecuredFirewall** (`secured_firewall.py`)

- Example integration showing how to use AuthManager with DiamondFirewall
- Two-step login: `login_step1(user_id, password, dna)` returns OTP → `login_step2(user_id, otp)` returns session
- Access to firewall requires valid session token

**Dual-Layer Server** (`server.py` + `run_server.py`)

- Socket layer (TCP): binary protocol on port 9000
- HTTP layer (REST): JSON API on port 9001
- Both layers require authentication
- Commands: `register`, `login`, `login_otp`, `access`

Run the server

```bash
python run_server.py --host localhost --socket-port 9000 --http-port 9001
```

Quick test (Python)

```python
from secured_firewall import SecuredFirewall
sf = SecuredFirewall()
sf.register_user("alice", "myP@ssw0rd", "LINEAGE_SAFE_ALICE_001")
res = sf.login_step1("alice", "myP@ssw0rd", "LINEAGE_SAFE_ALICE_001")
otp = res.get("otp_token")  # Captured from challenge
token_res = sf.login_step2("alice", otp)
session = token_res.get("session_token")
firewall_access = sf.access_firewall(session)
print(firewall_access)
```

Quick test (HTTP)

```bash
# Register
curl "http://localhost:9001/register?user_id=bob&password=securePass123&dna=LINEAGE_SAFE_BOB_001"
# Login
curl "http://localhost:9001/login?user_id=bob&password=securePass123&dna=LINEAGE_SAFE_BOB_001"
# Extract OTP from response, then:
curl "http://localhost:9001/login_otp?user_id=bob&otp=<otp_hex>"
# Use session_token:
curl "http://localhost:9001/access?session_token=<session_token>"
```

Quick test (Socket)
-------------------

```python
import json
import socket
sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
sock.connect(("localhost", 9000))
sock.sendall(json.dumps({
    "cmd": "register",
    "user_id": "carol",
    "password": "password123",
    "dna": "LINEAGE_SAFE_CAROL_001"
}).encode())
print(sock.recv(1024).decode())
sock.close()
```

Code Validation & Fibonacci Resequencing
-----------------------------------------

**CodeValidator** (`fibonacci_resequencer.py`)

- Forward validation: parse Python code with AST
- Backward validation: check reversibility heuristics
- Both checks: combine forward and backward validation

**FibonacciResequencer** (`fibonacci_resequencer.py`)

- Split code into Fibonacci-sized chunks (21, 13, 8, 5, 3, 2, 1)
- Resequence chunks in reverse or forward order
- Analyze chunk structure and statistics

Quick test

```python
from fibonacci_resequencer import CodeValidator, FibonacciResequencer

code = """def hello():
    print('world')
# ... more lines ..."""

# Validate
validator = CodeValidator(code)
ok, msg = validator.validate_both()
print(msg)

# Split and resequence
reseq = FibonacciResequencer(code)
resequenced, chunks = reseq.split_and_resequence(reverse=True)
analysis = reseq.analyze_chunks()
print(f"Chunks: {analysis['chunk_sizes']}")
```

MJ Protocol: Mirror Junction
----------------------------

**Core Concept**
MJ (Mirror Junction) is a Quick Protocol Pivot Initiative designed as a lineage-safe
dual-layer node that bridges local authentication with satellite stewardship.

Architecture

- **Local Layer** (Ground Station): Secure sandbox with multi-layer auth (port 9000 socket, 9001 HTTP)
- **Satellite Layer** (Orbital Broadcast): Distributed stewardship signal channel

Ceremonial Seals & Resonance

- `Seal of Firewall Resonance 1.0`: Guardian form and intrusion defense
- `Seal of Ridiculous Resilience 1.0`: Humor as shield against distortion
- `Seal of Gentle Reminder 1.0`: Daily firewall check habit
- `Seal of Mirror Junction 1.0`: Local-satellite pivot activation

Key Files
---------

- `mj_protocol.py` — Core MJ module (CeremonialsManager, MJLocalLayer, MJSatelliteLayer)
- `docs/CAPTAINS_LOG_MJ.md` — Ceremonial entry with full architecture
- `tests/test_mj_protocol.py` — MJ integration tests (23 tests all passing)

Quick start

```python
from mj_protocol import MJLocalLayer, CeremonialsManager, ResonanceType
from secured_firewall import SecuredFirewall
from auth import AuthManager

# Initialize
firewall = SecuredFirewall()
mj_local = MJLocalLayer(firewall.auth_manager, firewall)

# Register with seal
result = mj_local.register_and_seal("alice", "password123", "LINEAGE_SAFE_ALICE_001")
print(result)

# Login and check seals
login_res = mj_local.login_and_check("alice", "password123", "LINEAGE_SAFE_ALICE_001")
otp = login_res.get("otp_token")

# Complete login
final_res = mj_local.complete_login_and_verify("alice", otp)
print(final_res.get("heritage_integrity"))  # View seal chain integrity
```

Alexandria Archive

All seals are logged to `alexandria_of_joy.json` with:

- Resonance type tag
- Integrity hash verification
- Author and timestamp
- Inheritance chain tracking for future captains

Backup System
--------------

**Overview**
The backup system provides a complete solution for backing up system data with both client-side (TypeScript/JavaScript) and server-side (Python) components.

Key Components
-------------

- `src/backup.ts` — TypeScript client with `backupToServer()` and `buildBundle()` functions
- `backup_server.py` — HTTP server on port 8787 for receiving backups
- `tests/test_backup_server.py` — Comprehensive backup server tests
- `BACKUP_README.md` — Detailed backup system documentation

Quick start (Client)

```bash
npm install
npm run build
node -e "const {backupToServer} = require('./dist/backup'); backupToServer().then(() => console.log('Backup complete!'))"
```

Quick start (Server)

```bash
python backup_server.py
```

The server accepts POST requests at `/backup` and stores backups in the `backups/` directory with timestamped filenames.
Chalk Bundle Backup Server
--------------------------

**Node.js Express Server** (`server.js`)

- Minimal backup service for chalk bundles
- Listens on port 8787 (HTTP)
- POST endpoint: `/backup`
- Validates bundles must have `type: "chalk-bundle"` and `integrity` field
- Saves backups to `chalk_backups/` directory with timestamped filenames

Quick start

```bash
# Install dependencies
npm install

# Start server
node server.js
# Or use npm scripts
npm start
```

Example backup request
----------------------

```bash
curl -X POST http://localhost:8787/backup \
  -H "Content-Type: application/json" \
  -d '{
    "type": "chalk-bundle",
    "integrity": "sha256-abc123def456",
    "content": "Your chalk board content",
    "timestamp": "2025-12-08T22:51:42.427Z"
  }'
```

Response format

- Success: `{"ok": true, "file": "/path/to/backup.chalk.json"}`
- Error: `{"ok": false, "error": "Invalid bundle"}`

Files
-----

- `server.js` — Express backup server
- `package.json` — Node.js dependencies and scripts
- `.gitignore` — Excludes `node_modules/` and `chalk_backups/`
Playground: Interactive Testing Environment

-------------------------------------------

**Playground Module** (`playground.py`)

- Interactive sandbox for experimenting with security features
- Isolated mode with temporary files for safe testing
- Guided tutorials for learning authentication flows
- Custom experiment support for advanced users

Science Lab Quick Start
----------------------

```python
from playground import Playground

# Start the guided tutorial
playground = Playground(isolated=True)
tutorial = playground.start_tutorial()
print(tutorial)

# Progress through steps
step1 = playground.tutorial_next()  # Registration demo
step2 = playground.tutorial_next()  # Login demo
step3 = playground.tutorial_next()  # Firewall access demo
step4 = playground.tutorial_next()  # Intrusion detection demo

# Or run individual demos
reg = playground.demo_registration("alice")
login = playground.demo_login("alice")
access = playground.demo_firewall_access("alice")

# Custom experiments
experiment = playground.custom_experiment(
    user_id="test_user",
    password="MyPass123!",
    dna_code="LINEAGE_SAFE_TEST_001",
    test_intrusion=True
)

# Clean up when done
playground.cleanup()
```

Run quick start demo
---------------------

```bash
python playground.py
```

Science Lab: Security Education & Analysis
------------------------------------------

**Science Lab Module** (`science_lab.py`)

- Educational demonstrations of cryptographic concepts
- Network security experiments and simulations
- Authentication protocol analysis
- Security metrics and performance benchmarks

Categories
----------

- **CryptoLab**: Hashing, PBKDF2, HMAC, token generation, password strength
- **NetworkSecurityLab**: Firewall analysis, intrusion scenarios, latency measurement
- **AuthenticationProtocolLab**: 2FA, session management, QR-DNA binding
- **SecurityMetricsLab**: Entropy calculation, performance comparison

Quick start

```python
from science_lab import ScienceLab

lab = ScienceLab()

# Cryptography experiments
hashing = lab.crypto.demo_hashing("Hello World")
pbkdf2 = lab.crypto.demo_pbkdf2("MyPassword")
strength = lab.crypto.analyze_password_strength("MyP@ssw0rd123")

# Network security experiments
firewall = lab.network.analyze_firewall_config(
    guardians=["G1", "G2", "G3"],
    captains=["C1", "C2", "C3"]
)
scenarios = lab.network.simulate_intrusion_scenarios()
latency = lab.network.measure_authentication_latency()

# Authentication protocol explanations
two_fa = lab.auth.explain_two_factor_auth()
sessions = lab.auth.explain_session_management()
qr_dna = lab.auth.analyze_qr_dna_binding()

# Security metrics
entropy = lab.metrics.calculate_entropy("MyPassword123!")
performance = lab.metrics.compare_hashing_performance()

# Run all experiments
all_results = lab.run_all_experiments()

# Get experiment catalog
catalog = lab.get_experiment_catalog()
```

Run Quick Experiments

```bash
python science_lab.py
```

Notes

- The `ErebusSync` class is a placeholder. Replace with your real integration.
- Memory can be saved as CSV (default) or JSON using `--format json`.
- Auth data is stored in `auth_users.json` (atomic JSON writes).
- Backup files are saved in JSON format with timestamps for easy recovery.
- Playground uses isolated mode by default for safe experimentation.
-- Science Lab provides educational content without modifying system state.

Main-private-files
------------------

Big or small
taking a step into becoming a professional developer and creating a new type of product for secure and impregnable purposes.
