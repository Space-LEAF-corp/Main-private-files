"""
sign_owner_command.py - Generate HMAC-SHA256 signatures for Jarvondis owner commands

Usage:
    python sign_owner_command.py

Security:
    - Never commit real secrets to version control
    - Use environment variables or interactive input for production secrets
    - Generated signatures expire after 60 seconds (default freshness window)
"""
import time, hashlib, hmac

def hmac_sha256(secret: str, message: str) -> str:
    return hmac.new(secret.encode("utf-8"), message.encode("utf-8"), hashlib.sha256).hexdigest()

def make_signature(owner_id: str, secret: str, payload: str) -> tuple:
    ts = int(time.time())
    msg = f"{owner_id}:{ts}:{payload}"
    return hmac_sha256(secret, msg), ts

if __name__ == "__main__":
    owner_id = "leif.w.sogge"
    secret = "CHANGE_ME_TO_A_STRONG_SECRET"
    payload = "set_lockdown:False"
    sig, ts = make_signature(owner_id, secret, payload)
    print({"signature": sig, "timestamp": ts})
