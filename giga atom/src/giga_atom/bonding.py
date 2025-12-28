"""
bonding.py — Bonding heuristics and partner scanning
Space Leaf Corp — Internal Use Only
"""

import numpy as np
import pandas as pd
from .simulations import PHYSICAL_SHELLS

def valence_from_Z(Z: int) -> int:
    remaining = Z
    for cap in PHYSICAL_SHELLS:
        if remaining <= cap:
            return remaining
        remaining -= cap
    return remaining

def bonding_score(val_giga: int, val_elem: int) -> float:
    if val_giga <= 0 or val_elem <= 0:
        return 0.0
    base = min(val_giga, val_elem)
    mismatch = abs(val_elem - val_giga) / max(1, val_elem + val_giga)
    raw = base * (1 - mismatch)
    return raw / base

def bonding_scan(val_giga: int = 1, maxZ: int = 118):
    rows = []
    for Z in range(1, maxZ + 1):
        val = valence_from_Z(Z)
        score = bonding_score(val_giga, val)
        rows.append({"Z": Z, "valence": val, "bonding_score": score})
    return pd.DataFrame(rows)

def combined_stability_for_partner(val_giga: int, val_elem: int, days: int = 3650, d: float = 2.0):
    mismatch = abs(val_giga - val_elem) / max(1, val_elem + val_giga)
    t = np.arange(0, days + 1)
    stability = np.exp(-d * mismatch * t / 365.0)
    return t, stability
