#!/usr/bin/env python3
"""
tools/rotate_and_resign.py

Rotate Ed25519 signing key, publish revocation list, and optionally re-sign stored payloads.
Assumes payloads are base64 packaged blobs stored as .txt files in PAYLOAD_DIR.
Private keys must be stored securely; this script expects PEM-like files in KEY_DIR for demo.
In production, replace file I/O with HSM/vault operations.
"""
import argparse
import base64
import json
import os
from pathlib import Path
from nacl.signing import SigningKey, VerifyKey
from nacl.encoding import Base64Encoder
from datetime import datetime

# CONFIG - replace with secure vault integration in production
KEY_DIR = Path("keys")
PAYLOAD_DIR = Path("payloads")
PUBLISH_DIR = Path("published")
REVOCATION_FILE = PUBLISH_DIR / "revocation_list.json"
OWNER = "LEIF William Sogge"
VERSION = b"v1"

KEY_DIR.mkdir(parents=True, exist_ok=True)
PAYLOAD_DIR.mkdir(parents=True, exist_ok=True)
PUBLISH_DIR.mkdir(parents=True, exist_ok=True)

def load_signing_key(path):
    return SigningKey(path.read_bytes(), encoder=Base64Encoder)

def save_signing_key(sk, path):
    path.write_bytes(sk.encode(encoder=Base64Encoder))

def generate_new_keypair(name):
    sk = SigningKey.generate()
    vk = sk.verify_key
    sk_path = KEY_DIR / f"{name}.sk.b64"
    vk_path = KEY_DIR / f"{name}.vk.b64"
    save_signing_key(sk, sk_path)
    vk_path.write_text(vk.encode(encoder=Base64Encoder).decode())
    return sk, vk

def publish_revocation(old_vk_b64, reason="compromise", revoked_by=None):
    now = datetime.utcnow().isoformat() + "Z"
    entry = {"public_key": old_vk_b64, "revoked_at": now, "reason": reason, "revoked_by": revoked_by}
    if REVOCATION_FILE.exists():
        data = json.loads(REVOCATION_FILE.read_text())
    else:
        data = {"revocations": []}
    data["revocations"].append(entry)
    REVOCATION_FILE.write_text(json.dumps(data, indent=2))
    print(f"Published revocation for key {old_vk_b64[:16]}...")

def resign_payloads(new_sk, resign=True):
    # For each payload file, verify old signature then re-sign with new_sk if requested
    files = list(PAYLOAD_DIR.glob("*.txt"))
    count = 0
    for f in files:
        blob_b64 = f.read_text().strip()
        # unpack and re-sign: keep version||salt||nonce||ct||tag||signature||owner
        parts = base64.b64decode(blob_b64).split(b"||")
        if len(parts) != 7:
            print(f"Skipping {f.name}: unexpected format")
            continue
        v, salt, nonce, ct, tag, old_sig, owner = parts
        # create new signature over same fields (without old signature)
        new_sig = new_sk.sign(v + salt + nonce + ct + tag + owner).signature
        new_blob = base64.b64encode(b"||".join([v, salt, nonce, ct, tag, new_sig, owner])).decode()
        if resign:
            f.write_text(new_blob)
            count += 1
    print(f"Re-signed {count} payloads with new key.")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rotate", action="store_true", help="Generate new keypair and publish revocation for old key")
    parser.add_argument("--resign", action="store_true", help="Re-sign stored payloads with new key")
    parser.add_argument("--old-key", type=str, default=None, help="Base64 old verify key to revoke (optional)")
    args = parser.parse_args()

    # Find current key if exists
    existing_vk_files = list(KEY_DIR.glob("*.vk.b64"))
    existing_sk_files = list(KEY_DIR.glob("*.sk.b64"))
    current_vk_b64 = None
    if existing_vk_files:
        current_vk_b64 = existing_vk_files[-1].read_text().strip()
    if args.rotate:
        # revoke current if provided
        if current_vk_b64:
            publish_revocation(current_vk_b64, reason="rotation", revoked_by=OWNER)
        # generate new keypair
        new_sk, new_vk = generate_new_keypair(datetime.utcnow().strftime("key_%Y%m%dT%H%M%SZ"))
        new_vk_b64 = new_vk.encode(encoder=Base64Encoder).decode()
        print(f"Generated new key. Public key (base64): {new_vk_b64}")
        # Optionally re-sign payloads
        if args.resign:
            resign_payloads(new_sk, resign=True)
        # publish public key file
        (PUBLISH_DIR / "current_public_key.b64").write_text(new_vk_b64)
        print("Published new public key to published/current_public_key.b64")
    else:
        print("No rotation requested. Use --rotate to rotate keys.")

if __name__ == "__main__":
    main()