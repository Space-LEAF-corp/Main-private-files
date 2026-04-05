import hashlib
import json
from dataclasses import dataclass
from typing import Dict

@dataclass
class UserAccount:
    username: str
    password_hash: str  # salted hash in a real system

@dataclass
class ProjectMetadata:
    owner: str
    integrity_hash: str

class SecurityManager:
    def __init__(self):
        self._users: Dict[str, UserAccount] = {}
        self._projects: Dict[str, ProjectMetadata] = {}

    def create_user(self, username: str, password: str) -> None:
        if username in self._users:
            raise ValueError("User already exists")
        pw_hash = self._hash_password(password)
        self._users[username] = UserAccount(username=username, password_hash=pw_hash)

    def authenticate(self, username: str, password: str) -> bool:
        user = self._users.get(username)
        if not user:
            return False
        return user.password_hash == self._hash_password(password)

    def register_project(self, project_id: str, owner: str, payload: dict[str, object]) -> None:
        integrity_hash = self._hash_payload(payload)
        self._projects[project_id] = ProjectMetadata(owner=owner, integrity_hash=integrity_hash)

    def verify_project(self, project_id: str, payload: dict[str, object]) -> bool:
        meta = self._projects.get(project_id)
        if not meta:
            return False
        return meta.integrity_hash == self._hash_payload(payload)

    @staticmethod
    def _hash_password(password: str) -> str:
        # Use PBKDF2-HMAC with SHA-256 for password hashing.
        # Note: In a real system, use a per-user random salt stored alongside the hash.
        salt = b"ceremonial-terrain-static-salt"
        iterations = 100_000
        dk = hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), salt, iterations)
        return dk.hex()

    @staticmethod
    def _hash_payload(payload: dict[str, object]) -> str:
        data = json.dumps(payload, sort_keys=True).encode("utf-8")
        return hashlib.sha256(data).hexdigest()
