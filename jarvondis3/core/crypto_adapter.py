#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Crypto adapter for Jarvondis 3.0
- Ed25519 signatures for admin + validators
- HMAC fallback (optional)
- Canonical JSON helper (stable serialization)
"""

import json
from typing import Optional, List
from nacl.signing import SigningKey, VerifyKey
from nacl.exceptions import BadSignatureError
import nacl.encoding


def canonical_json(data) -> str:
    return json.dumps(data, sort_keys=True, separators=(",", ":"))


class Ed25519Adapter:
    """
    Admin and validators use Ed25519 keys.
    - Admin: one signing key for root authority
    - Validators: list of verify keys (public) or signing keys if you also sign locally
    """
    def __init__(
        self,
        admin_signing_key_b64: Optional[str] = None,
        admin_verify_key_b64: Optional[str] = None,
        validator_verify_keys_b64: Optional[List[str]] = None
    ):
        # Admin may be loaded by private signing key OR public verify key (verify-only mode)
        self.admin_signing_key = None
        if admin_signing_key_b64:
            self.admin_signing_key = SigningKey(admin_signing_key_b64.encode(), encoder=nacl.encoding.Base64Encoder)
        
        self.admin_verify_key = None
        if admin_verify_key_b64:
            self.admin_verify_key = VerifyKey(admin_verify_key_b64.encode(), encoder=nacl.encoding.Base64Encoder)
        elif self.admin_signing_key:
            self.admin_verify_key = self.admin_signing_key.verify_key
            
        self.validator_verify_keys = []
        for v in (validator_verify_keys_b64 or []):
            self.validator_verify_keys.append(
                VerifyKey(v.encode(), encoder=nacl.encoding.Base64Encoder)
            )

    def admin_sign(self, message: str) -> bytes:
        if not self.admin_signing_key:
            raise RuntimeError("Admin signing key not loaded")
        return self.admin_signing_key.sign(message.encode("utf-8")).signature

    def admin_verify(self, message: str, signature: bytes) -> bool:
        if not self.admin_verify_key:
            raise RuntimeError("Admin verify key not loaded")
        try:
            self.admin_verify_key.verify(message.encode("utf-8"), signature)
            return True
        except BadSignatureError:
            return False

    def validator_verify(self, message: str, signature: bytes, validator_index: int) -> bool:
        if validator_index < 0 or validator_index >= len(self.validator_verify_keys):
            return False
        try:
            self.validator_verify_keys[validator_index].verify(message.encode("utf-8"), signature)
            return True
        except BadSignatureError:
            return False
