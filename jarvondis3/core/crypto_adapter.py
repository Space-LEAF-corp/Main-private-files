#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Crypto adapter for Jarvondis 3.0
- Ed25519 signatures for admin + validators
- HMAC fallback (optional)
- Canonical JSON helper (stable serialization)

This adapter exposes a small, stable verification API that accepts signatures
as bytes, base64 strings, or hex strings. Main verification entry points:
- verify_with_pubkey(pubkey_b64: str, message: str, signature: Union[bytes,str]) -> bool
- verify_signature(pubkey_b64: str, message: str, signature: Union[bytes,str]) -> bool  # alias
- verify_admin_signature(message: str, signature: Union[bytes,str]) -> bool
- admin_sign(message: str) -> bytes
- admin_sign_base64(message: str) -> str
"""
import base64
import binascii
import json
from typing import Optional, List, Union
from nacl.signing import SigningKey, VerifyKey
from nacl.exceptions import BadSignatureError
import nacl.encoding


def canonical_json(data) -> str:
    return json.dumps(data, sort_keys=True, separators=(",", ":"))


class Ed25519Adapter:
    """
    Admin and validators use Ed25519 keys.

    Exposed stable verification API:
    - verify_with_pubkey(pubkey_b64, message, signature) -> bool
      pubkey_b64: base64-encoded ed25519 public key
      signature: bytes | base64-string | hex-string

    - verify_admin_signature(message, signature) -> bool
      Uses configured admin_verify_key if available.

    Signature encoding:
    - Accepts bytes
    - Accepts base64 string
    - Accepts hex string
    """

    def __init__(
        self,
        admin_signing_key_b64: Optional[str] = None,
        admin_verify_key_b64: Optional[str] = None,
        validator_verify_keys_b64: Optional[List[str]] = None
    ):
        # Keep the raw b64 values for callers that expect them
        self.admin_signing_key_b64 = admin_signing_key_b64
        self.admin_verify_key_b64 = admin_verify_key_b64

        # Admin may be loaded by private signing key OR public verify key (verify-only mode)
        self.admin_signing_key: Optional[SigningKey] = None
        if admin_signing_key_b64:
            # SigningKey expects raw bytes or a Base64Encoder-wrapped string
            # We use the Base64Encoder to construct from base64 text
            self.admin_signing_key = SigningKey(admin_signing_key_b64.encode(), encoder=nacl.encoding.Base64Encoder)

        self.admin_verify_key: Optional[VerifyKey] = None
        if admin_verify_key_b64:
            self.admin_verify_key = VerifyKey(admin_verify_key_b64.encode(), encoder=nacl.encoding.Base64Encoder)
        elif self.admin_signing_key:
            self.admin_verify_key = self.admin_signing_key.verify_key
            # also set admin_verify_key_b64 for callers
            try:
                self.admin_verify_key_b64 = self.admin_verify_key.encode(encoder=nacl.encoding.Base64Encoder).decode()
            except Exception:
                # keep previous value if encoding fails
                pass

        # Validator verify keys list (VerifyKey objects)
        self.validator_verify_keys: List[VerifyKey] = []
        for v in (validator_verify_keys_b64 or []):
            self.validator_verify_keys.append(
                VerifyKey(v.encode(), encoder=nacl.encoding.Base64Encoder)
            )

    # --- Helpers for signature encoding normalization ---
    def _decode_signature(self, signature: Union[bytes, str]) -> bytes:
        """
        Decode signature provided as bytes, base64 string, or hex string into raw bytes.
        """
        if isinstance(signature, bytes):
            return signature
        if isinstance(signature, str):
            # Try base64 first
            try:
                return base64.b64decode(signature, validate=True)
            except (binascii.Error, ValueError):
                pass
            # Try hex
            try:
                return bytes.fromhex(signature)
            except (ValueError, TypeError):
                pass
            # Last attempt: try raw utf-8 bytes (not recommended)
            return signature.encode("utf-8")
        raise TypeError("Unsupported signature type. Expected bytes or str.")

    def _decode_pubkey_b64(self, pubkey_b64: str) -> VerifyKey:
        """
        Construct a VerifyKey from a base64-encoded public key string.
        """
        if not isinstance(pubkey_b64, str):
            raise TypeError("pubkey_b64 must be a base64-encoded string")
        return VerifyKey(pubkey_b64.encode(), encoder=nacl.encoding.Base64Encoder)

    # --- Admin sign/verify helpers ---
    def admin_sign(self, message: str) -> bytes:
        """
        Return raw signature bytes. Raises RuntimeError if admin signing key not configured.
        """
        if not self.admin_signing_key:
            raise RuntimeError("Admin signing key not loaded")
        signed = self.admin_signing_key.sign(message.encode("utf-8"))
        return signed.signature

    def admin_sign_base64(self, message: str) -> str:
        """
        Return signature encoded as base64 string (convenience wrapper).
        """
        sig = self.admin_sign(message)
        return base64.b64encode(sig).decode()

    def admin_verify(self, message: str, signature: Union[bytes, str]) -> bool:
        """
        Verify signature with the configured admin_verify_key.
        Accepts signature as bytes, base64 string, or hex string.
        """
        if not self.admin_verify_key:
            raise RuntimeError("Admin verify key not loaded")
        sig_bytes = self._decode_signature(signature)
        try:
            self.admin_verify_key.verify(message.encode("utf-8"), sig_bytes)
            return True
        except BadSignatureError:
            return False

    # expose a stable alias expected by sovereign.py
    def verify_admin_signature(self, message: str, signature: Union[bytes, str]) -> bool:
        return self.admin_verify(message, signature)

    # --- Generic verification API ---
    def verify_with_pubkey(self, pubkey_b64: str, message: str, signature: Union[bytes, str]) -> bool:
        """
        Verify an Ed25519 signature using an explicit base64-encoded public key.
        Returns True on success, False otherwise.
        Accepts signature as bytes, base64 string, or hex string.
        """
        try:
            vk = self._decode_pubkey_b64(pubkey_b64)
        except Exception:
            return False
        sig_bytes = self._decode_signature(signature)
        try:
            vk.verify(message.encode("utf-8"), sig_bytes)
            return True
        except BadSignatureError:
            return False

    # alias name used by sovereign.py and other callers
    def verify_signature(self, pubkey_b64: str, message: str, signature: Union[bytes, str]) -> bool:
        return self.verify_with_pubkey(pubkey_b64, message, signature)

    # --- Validator helpers (by index) ---
    def validator_verify(self, message: str, signature: Union[bytes, str], validator_index: int) -> bool:
        """
        Verify a validator signature by its index in the configured validator list.
        This keeps backward compatibility with earlier API.
        """
        if validator_index < 0 or validator_index >= len(self.validator_verify_keys):
            return False
        sig_bytes = self._decode_signature(signature)
        try:
            self.validator_verify_keys[validator_index].verify(message.encode("utf-8"), sig_bytes)
            return True
        except BadSignatureError:
            return False
