#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TOTP MFA for Jarvondis 3.0
- RFC 6238 time-based one-time passwords
- Uses pyotp-like semantics but implemented minimally with hmac/sha1
"""

import base64
import hmac
import hashlib
import time
from typing import Optional


def _int_to_bytes(i: int) -> bytes:
    return i.to_bytes(8, "big")


class TOTP:
    """
    Minimal TOTP verifier.
    - secret_b32: base32-encoded secret (e.g., from authenticator app)
    - interval: time step in seconds (default 30)
    - digits: number of digits (default 6)
    """
    def __init__(self, secret_b32: str, interval: int = 30, digits: int = 6):
        self.secret = base64.b32decode(secret_b32.upper())
        self.interval = interval
        self.digits = digits

    def _code_at(self, for_time: int) -> str:
        counter = for_time // self.interval
        msg = _int_to_bytes(counter)
        digest = hmac.new(self.secret, msg, hashlib.sha1).digest()
        offset = digest[-1] & 0x0F
        binary = ((digest[offset] & 0x7f) << 24) | (digest[offset + 1] << 16) | (digest[offset + 2] << 8) | (digest[offset + 3])
        otp = binary % (10 ** self.digits)
        return str(otp).zfill(self.digits)

    def verify(self, code: str, valid_window: int = 1, now: Optional[int] = None) -> bool:
        now = now or int(time.time())
        for delta in range(-valid_window, valid_window + 1):
            if self._code_at(now + delta * self.interval) == code:
                return True
        return False

    def now(self) -> str:
        """Generate current TOTP code"""
        return self._code_at(int(time.time()))
