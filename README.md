Jarvondis (refactor)
=====================

This workspace contains a refactored `Jarvondis` module and a small CLI runner.

Quick start
-----------

- Install dependencies (optional, recommended):

```bash
python -m pip install -r requirements.txt
```

- Run interactively:

```bash
python cli.py run --memory-file jarvondis_memory.csv --format csv --tone witty
```

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

**Run the server**
```bash
python run_server.py --host localhost --socket-port 9000 --http-port 9001
```

**Quick test (Python)**
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

**Quick test (HTTP)**
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

**Quick test (Socket)**
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

**Quick test**
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

Notes
- The `ErebusSync` class is a placeholder. Replace with your real integration.
- Memory can be saved as CSV (default) or JSON using `--format json`.
- Auth data is stored in `auth_users.json` (atomic JSON writes).
# Main-private-files
Big or small
taking a step into becoming a professional developer and creating a new type of product for secure and impregnable purposes.