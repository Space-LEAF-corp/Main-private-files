from __future__ import annotations
import time
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Callable
import hashlib

# Your organic symbol definitions (circle counts + special roles)
class OrganicSymbol:
    # Circle counts and properties — exactly as you defined
    ZERO  = {"value": 0, "circles": 0, "unique": True, "binary_fold": 0}
    ONE   = {"value": 1, "circles": 0, "unique": True, "binary_fold": 1}  # or single dot
    TWO   = {"value": 2, "circles": 2}
    THREE = {"value": 3, "circles": 3}
    FOUR_TO_EIGHT = {
        "values": [4,5,6,7,8],
        "circles": 3,                  # shared [3°] clover glyph
        "compressed": True,
        "requires_context": True
    }
    NINE  = {"value": 9, "circles": 9, "role": "master_key", "privileged": True}
    TEN   = {"value": 10, "circles": 10, "role": "master_key", "privileged": True}

@dataclass
class SymbolSequence:
    """A validated sequence of organic symbols — your "key"."""
    symbols: List[dict]  # list of symbol definitions (e.g., OrganicSymbol.NINE)
    hash: str = field(init=False)
    
    def __post_init__(self):
        self.hash = hashlib.sha3_256(str([s["value"] for s in self.symbols]).encode()).hexdigest()

class TimeEngineSecurity:
    """Security layer using your organic multinary symbols."""
    
    def __init__(self, master_sequence: Optional[SymbolSequence] = None):
        self.authorized = False
        self.master_sequence = master_sequence  # the true key only you hold
        self.failed_attempts = 0
        self.lockout_threshold = 3
    
    def authenticate(self, candidate: List[dict]) -> bool:
        """Authenticate using a sequence of organic symbols."""
        if self.lockout_threshold > 0 and self.failed_attempts >= self.lockout_threshold:
            return False  # permanent lockout after too many fails
        
        # Must contain at least one master key (9 or 10)
        has_master = any(s in (OrganicSymbol.NINE, OrganicSymbol.TEN) for s in candidate)
        if not has_master:
            self.failed_attempts += 1
            return False
        
        # If we have a stored master sequence, it must match exactly
        if self.master_sequence:
            candidate_seq = SymbolSequence(candidate)
            match = candidate_seq.hash == self.master_sequence.hash
            if not match:
                self.failed_attempts += 1
            return match
        
        # First valid master-key sequence becomes the root of trust
        self.master_sequence = SymbolSequence(candidate)
        self.authorized = True
        return True
    
    def require_auth(self, func: Callable) -> Callable:
        """Decorator: only allow operation if authenticated with symbols."""
        def wrapper(engine: 'SecureTimeEngine', *args, **kwargs):
            if not engine.security.authorized:
                raise PermissionError("TimeEngine locked: valid organic symbol sequence with master key required.")
            return func(engine, *args, **kwargs)
        return wrapper

class SecureTimeEngine(TimeEngine):
    """TimeEngine with organic symbol security layer."""
    
    def __init__(self, master_symbol_key: Optional[List[dict]] = None):
        super().__init__()
        self.security = TimeEngineSecurity(
            master_sequence=SymbolSequence(master_symbol_key) if master_symbol_key else None
        )
    
    def authenticate_with_symbols(self, symbol_sequence: List[dict]) -> Dict:
        """Public authentication method."""
        success = self.security.authenticate(symbol_sequence)
        return {
            "authenticated": success,
            "status": "authorized" if success else "denied",
            "failed_attempts": self.security.failed_attempts
        }
    
    # Override critical methods with auth decorator
    @TimeEngineSecurity.require_auth
    def add_item(self, *args, **kwargs):
        return super().add_item(*args, **kwargs)
    
    @TimeEngineSecurity.require_auth
    def review_item(self, *args, **kwargs):
        return super().review_item(*args, **kwargs)
    
    @TimeEngineSecurity.require_auth
    def get_next_reviews(self, *args, **kwargs):
        return super().get_next_reviews(*args, **kwargs)
    
    @TimeEngineSecurity.require_auth
    def learning_analytics(self, *args, **kwargs):
        return super().learning_analytics(*args, **kwargs)
    
    @TimeEngineSecurity.require_auth
    def export_state(self, *args, **kwargs):
        return super().export_state(*args, **kwargs)

# Example of locking it with your true key
# You define this once, privately — e.g., [NINE, THREE, TEN, ONE, NINE]
TRUE_MASTER_KEY = [
    OrganicSymbol.NINE,
    OrganicSymbol.THREE,
    OrganicSymbol.TEN,
    OrganicSymbol.ONE,
    OrganicSymbol.NINE
]

engine = SecureTimeEngine(master_symbol_key=TRUE_MASTER_KEY)

# Anyone else trying to use it without the exact symbol sequence (with 9/10) gets blocked