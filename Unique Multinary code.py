from enum import Enum
from typing import List

class OrganicDigit(Enum):
    ZERO = {
        "value": 0,
        "circles": 0,                  # completely unique: no circles
        "description": "empty / void / absence",
        "binary_fold": 0,
        "special": "binary heartbeat"
    }
    
    ONE = {
        "value": 1,
        "circles": 0,                  # completely unique: single point or minimal mark
        "description": "single dot / primal mark",
        "binary_fold": 1,
        "special": "binary heartbeat"
    }
    
    TWO = {"value": 2, "circles": 2, "description": "two circles, paired"}
    THREE = {"value": 3, "circles": 3, "description": "three circles, triangle"}
    FOUR_TO_EIGHT = {
        "values": [4,5,6,7,8],
        "circles": 3,                  # shared [3°] glyph — the clover/knot
        "description": "three-circle clover [3°], disambiguated by context/master key",
        "compressed": True
    }
    
    NINE = {
        "value": 9,
        "circles": 9,                  # or specific pattern — you define exact arrangement
        "description": "master key — nine circles",
        "role": "privileged mode / authentication / escape"
    }
    
    TEN = {
        "value": 10,
        "circles": 10,                 # or ring + center, etc.
        "description": "master key — ten circles",
        "role": "fold control / dimension shift / trust gate"
    }

class OrganicCode:
    def __init__(self):
        self.stream = []  # list of OrganicDigit entries
    
    def encode(self, message: str) -> List:
        # We'll build this together based on your full rules
        pass
    
    def decode(self, stream) -> str:
        pass
    
    def fold_to_binary(self):
        # Collapse to pure 0/1 stream where possible
        pass
    
    def requires_master_key(self, stream) -> bool:
        # Detect if 9 or 10 present
        pass