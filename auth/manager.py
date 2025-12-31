"""Authentication manager (clean module).
"""
from __future__ import annotations

import json
import os
import secrets
import tempfile
import time
from dataclasses import dataclass, asdict
from hashlib import pbkdf2_hmac
from typing import Dict, Optional


DEFAULT_USERS_FILE = "auth_users.json"
PBKDF2_ITER = 150_000


def _atomic_write(path: str, data: str) -> None:
    dir_ = os.path.dirname(os.path.abspath(path)) or None
    tmp = tempfile.NamedTemporaryFile("w", delete=False, dir=dir_, encoding="utf-8")
    try:
        with open(tmp.name, "w", encoding="utf-8") as f:
            f.write(data)
        os.replace(tmp.name, path)
    finally:
        if os.path.exists(tmp.name):
            try:
                os.remove(tmp.name)
            except Exception:
                pass


def _hash_password(password: str, salt: bytes) -> str:
    dk = pbkdf2_hmac("sha256", password.encode("utf-8"), salt, PBKDF2_ITER)
    return dk.hex()


def _generate_salt() -> bytes:
    return secrets.token_bytes(16)


@dataclass
class UserRecord:
    user_id: str
    password_hash: str
    salt_hex: str
    dna_bound: Optional[str] = None
    otp_secret: Optional[str] = None
    otp_expiry: Optional[float] = None

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


class AuthManager:
    def __init__(self, users_file: str = DEFAULT_USERS_FILE):
        self.users_file = users_file
        self._users: Dict[str, UserRecord] = {}
        self._load()

    def _load(self) -> None:
        if not os.path.exists(self.users_file):
            self._users = {}
            return
        try:
            with open(self.users_file, "r", encoding="utf-8") as f:
                data = json.load(f)
            self._users = {uid: UserRecord(**rec) for uid, rec in data.items()}
        except Exception:
            self._users = {}

    def save(self) -> None:
        data: Dict[str, Dict[str, object]] = {uid: rec.to_dict() for uid, rec in self._users.items()}
        _atomic_write(self.users_file, json.dumps(data, ensure_ascii=False, indent=2))

    def register(self, user_id: str, password: str, dna_code: str) -> bool:
        if not user_id:
            raise ValueError("user_id required")
        if len(password) < 8:
            raise ValueError("password must be at least 8 characters")
        if not dna_code.startswith("LINEAGE_SAFE"):
            raise ValueError("invalid dna_code")
        if user_id in self._users:
            raise ValueError("user already exists")

        salt = _generate_salt()
        password_hash = _hash_password(password, salt)
        rec = UserRecord(
            user_id=user_id,
            password_hash=password_hash,
            salt_hex=salt.hex(),
            dna_bound=dna_code,
        )
        self._users[user_id] = rec
        self.save()
        return True

    def _get_user(self, user_id: str) -> Optional[UserRecord]:
        return self._users.get(user_id)

    def start_login(self, user_id: str, password: str, dna_code: Optional[str] = None) -> Dict[str, str]:
        rec = self._get_user(user_id)
        if rec is None:
            return {"status": "denied", "reason": "unknown user"}

        salt = bytes.fromhex(rec.salt_hex)
        attempted = _hash_password(password, salt)
        if not secrets.compare_digest(attempted, rec.password_hash):
            return {"status": "denied", "reason": "invalid credentials"}

        if rec.dna_bound:
            if not dna_code or dna_code != rec.dna_bound:
                return {"status": "denied", "reason": "dna mismatch"}
        else:
            if not dna_code or not dna_code.startswith("LINEAGE_SAFE"):
                return {"status": "denied", "reason": "dna required for first-time bind"}
            rec.dna_bound = dna_code

        otp = secrets.token_hex(3)
        rec.otp_secret = otp
        rec.otp_expiry = time.time() + 300
        self.save()
        return {"status": "challenge", "method": "otp", "otp": otp}

    def complete_login(self, user_id: str, otp: str) -> Dict[str, str]:
        rec = self._get_user(user_id)
        if rec is None:
            return {"status": "denied", "reason": "unknown user"}
        if not rec.otp_secret or not rec.otp_expiry:
            return {"status": "denied", "reason": "no otp requested"}
        if time.time() > rec.otp_expiry:
            return {"status": "denied", "reason": "otp expired"}
        if not secrets.compare_digest(rec.otp_secret, otp):
            return {"status": "denied", "reason": "invalid otp"}
        rec.otp_secret = None
        rec.otp_expiry = None
        self.save()
        session_token = secrets.token_hex(16)
        return {"status": "ok", "session_token": session_token}

    def reset_password(self, user_id: str, new_password: str) -> bool:
        rec = self._get_user(user_id)
        if rec is None:
            raise ValueError("unknown user")
        if not new_password or len(new_password) < 8:
            raise ValueError("password must be at least 8 characters")
        salt = _generate_salt()
        rec.salt_hex = salt.hex()
        rec.password_hash = _hash_password(new_password, salt)
        self.save()
        return True
