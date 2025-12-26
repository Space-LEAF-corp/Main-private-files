#!/usr/bin/env python3
import os
import secrets
import string
from datetime import datetime

ENV_FILE = ".env"
BACKUP_DIR = "secret_backups"

ALPHABET = string.ascii_letters + string.digits + "!@#$%^&*-_=+"

def generate_secret(length=64):
    return "".join(secrets.choice(ALPHABET) for _ in range(length))

def load_env(path):
    env = {}
    if not os.path.exists(path):
        return env
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            env[k] = v
    return env

def save_env(path, env):
    with open(path, "w") as f:
        for k, v in env.items():
            f.write(f"{k}={v}\n")

def backup_env(path):
    os.makedirs(BACKUP_DIR, exist_ok=True)
    timestamp = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    backup_path = os.path.join(BACKUP_DIR, f".env.{timestamp}.bak")
    with open(path, "r") as src, open(backup_path, "w") as dst:
        dst.write(src.read())
    print(f"Backed up {path} to {backup_path}")

def main():
    if not os.path.exists(ENV_FILE):
        print(f"{ENV_FILE} not found. Create it from .env.example first.")
        return

    print("Loading current env...")
    env = load_env(ENV_FILE)

    # Decide which keys to rotate (your rotation matrix)
    keys_to_rotate = [
        "ADMIN_INITIAL_PASSWORD",
        "DB_PASSWORD",
        "JWT_SECRET"
    ]

    # Backup current env
    backup_env(ENV_FILE)

    # Rotate
    for key in keys_to_rotate:
        old = env.get(key, "")
        new = generate_secret(32 if key.endswith("_PASSWORD") else 64)
        env[key] = new
        print(f"Rotated {key}: {len(old)} chars -> {len(new)} chars")

    # Save updated env
    save_env(ENV_FILE, env)
    print("Rotation complete. Remember to restart the app so changes take effect.")

if __name__ == "__main__":
    main()
