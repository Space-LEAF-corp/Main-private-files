import hashlib
import hmac
import os
from datetime import datetime

class Jarvondis3TamperProofing:
    def __init__(self, admin_key: bytes):
        self.admin_key = admin_key
        self.file_hashes = {}

    # 1. Hash-Based Integrity Verification
    def register_file(self, filepath: str):
        self.file_hashes[filepath] = self._hash_file(filepath)

    def verify_file(self, filepath: str) -> bool:
        current_hash = self._hash_file(filepath)
        return current_hash == self.file_hashes.get(filepath)

    def _hash_file(self, filepath: str) -> str:
        hasher = hashlib.sha256()
        with open(filepath, "rb") as f:
            hasher.update(f.read())
        return hasher.hexdigest()

    # 2. Digital Signatures (mocked with HMAC)
    def sign_command(self, command: str) -> str:
        return hmac.new(self.admin_key, command.encode(), hashlib.sha256).hexdigest()

    def verify_signature(self, command: str, signature: str) -> bool:
        expected = self.sign_command(command)
        return hmac.compare_digest(expected, signature)

    # 3. Role-Based Access Control (simplified)
    def check_access(self, user_role: str, required_role: str) -> bool:
        hierarchy = ["user", "operator", "admin"]
        return hierarchy.index(user_role) >= hierarchy.index(required_role)

    # 4. Intrusion Detection & Prevention (mock)
    def detect_intrusion(self, event: str) -> bool:
        suspicious_patterns = ["unauthorized", "bruteforce", "tamper"]
        return any(p in event.lower() for p in suspicious_patterns)

    def prevent_intrusion(self, event: str):
        if self.detect_intrusion(event):
            print(f"[ALERT] Intrusion detected: {event}")
            # rollback, lockout, or alert logic here


class Jarvondis3GlobalCommands:
    def __init__(self):
        # Sovereign command namespace
        self.commands = {
            "/JARVONDIS_SYS_REBOOT": "System-wide reboot of Jarvondis 3.0",
            "/JARVONDIS_SYS_SHUTDOWN": "Graceful shutdown of all modules",
            "/JARVONDIS_SYS_UPDATE": "Initiate update sequence with changelog continuity",
            "/JARVONDIS_SYS_BACKUP": "Trigger immediate encrypted backup",
            "/JARVONDIS_SYS_ACCESS_GRANT": "Grant scoped access to subsystems or commands",
            "/JARVONDIS_SYS_ACCESS_REVOKE": "Revoke scoped access to subsystems or commands",
            "/JARVONDIS_SYS_DIAGNOSTICS": "Run comprehensive system diagnostics with trust banners",
            "/JARVONDIS_SYS_PROTOCOL_OVERRIDE": "Temporarily override protocols (requires ceremonial approval)"
        }
        self.log = []

    def list_commands(self):
        return self.commands

    def execute_command(self, command: str, authorized: bool = False):
        timestamp = datetime.utcnow().isoformat()
        if not authorized:
            entry = f"[{timestamp}] [VIEW ONLY] {command} is restricted. Root authority required."
            self.log.append(entry)
            return entry
        if command in self.commands:
            entry = f"[{timestamp}] [EXECUTED] {command}: {self.commands[command]}"
            self.log.append(entry)
            return entry
        entry = f"[{timestamp}] [ERROR] Unknown command: {command}"
        self.log.append(entry)
        return entry

    def view_log(self):
        return "\n".join(self.log)


# === Example ceremonial usage ===
if __name__ == "__main__":
    # Initialize with a mock admin key
    tp = Jarvondis3TamperProofing(admin_key=b'sovereign_secret_key')
    gc = Jarvondis3GlobalCommands()

    # Viewing mode
    print("Available Commands:", gc.list_commands())

    # Attempt execution without authority
    print(gc.execute_command("/JARVONDIS_SYS_UPDATE"))

    # Authorized execution
    print(gc.execute_command("/JARVONDIS_SYS_UPDATE", authorized=True))

    # View ceremonial log
    print("\n--- Command Log ---")
    print(gc.view_log())
