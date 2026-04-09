
PPS-OS-System
Pass-Phrase Signature Operating System
Captain Leif Edition — Local-First Demo
"""

from dataclasses import dataclass
from datetime import datetime
import hashlib
import hmac
import secrets


# ---------------------------------------------------------
# 1. Core PPS Engine
# ---------------------------------------------------------

@dataclass
class PassPhraseSignature:
    phrase: str
    created_at: datetime = datetime.utcnow()

    def derive_key(self) -> str:
        """
        Converts the pass-phrase signature into a stable,
        high-entropy cryptographic key using SHA-256 + salt.
        """
        salt = b"PPS-OS-Salt-v1"  # static for demo; can be rotated later
        digest = hmac.new(salt, self.phrase.encode("utf-8"), hashlib.sha256).hexdigest()
        return digest


# ---------------------------------------------------------
# 2. PPS-OS Identity Token
# ---------------------------------------------------------

@dataclass
class PPSIdentityToken:
    user_id: str
    signature_key: str
    issued_at: datetime = datetime.utcnow()

    def verify(self, phrase: str) -> bool:
        """
        Re-derive the key from the provided phrase and compare.
        """
        expected = PassPhraseSignature(phrase).derive_key()
        return secrets.compare_digest(expected, self.signature_key)


# ---------------------------------------------------------
# 3. PPS-OS System Controller
# ---------------------------------------------------------

class PPSOS:
    def __init__(self):
        self.registered_tokens = {}

    def register_user(self, user_id: str, phrase: str) -> PPSIdentityToken:
        """
        Creates a new PPS identity token for a user.
        """
        sig = PassPhraseSignature(phrase)
        key = sig.derive_key()
        token = PPSIdentityToken(user_id=user_id, signature_key=key)
        self.registered_tokens[user_id] = token
        return token

    def authenticate(self, user_id: str, phrase: str) -> bool:
        """
        Validates a pass-phrase signature for a user.
        """
        token = self.registered_tokens.get(user_id)
        if not token:
            return False
        return token.verify(phrase)


# ---------------------------------------------------------
# 4. Demo Run (safe to execute locally)
# ---------------------------------------------------------

if __name__ == "__main__":
    print("⭐ PPS-OS-System — Local Test ⭐\n")

    pps = PPSOS()

    # Captain registers a pass-phrase signature
    captain_phrase = "SSK, GO 💿."
    captain_id = "captain_leif"

    print("Registering PPS signature...")
    token = pps.register_user(captain_id, captain_phrase)
    print("Signature Key:", token.signature_key)
    print("Issued At:", token.issued_at)
    print("\n---\n")

    # Test authentication
    print("Testing authentication with correct phrase...")
    result = pps.authenticate(captain_id, "SSK, GO 💿.")
    print("Auth Result:", result)

    print("\nTesting authentication with incorrect phrase...")
    result = pps.authenticate(captain_id, "SSK, STOP 💿.")
    print("Auth Result:", result)

    print("\n⭐ PPS-OS-System Test Complete ⭐")
