"""
simulations.py — Stability, stress tests, hardening, Monte Carlo
Space Leaf Corp — Internal Use Only
"""

import numpy as np
from typing import List, Dict, Tuple
from .model import ShellState, compute_excess_fraction, compute_lambda

PHYSICAL_SHELLS = [2, 8, 18, 32, 50, 72, 98, 128]
TOTAL_ELECTRONS = sum(PHYSICAL_SHELLS)

def baseline_counts() -> List[int]:
    return PHYSICAL_SHELLS.copy()

def apply_overload(base: List[int], overloads: Dict[int, float]) -> List[int]:
    counts = base.copy()
    for idx, frac in overloads.items():
        i = idx - 1
        counts[i] = int(round(counts[i] * (1 + frac)))
    return counts

def compute_states(counts: List[int], k: float = 3.0) -> List[ShellState]:
    states = []
    for i, (phys, cnt) in enumerate(zip(PHYSICAL_SHELLS, counts), start=1):
        excess = compute_excess_fraction(cnt, phys)
        lam = compute_lambda(excess, k)
        states.append(ShellState(i, phys, cnt, excess, lam))
    return states

def stability_time_series(states: List[ShellState], days: int = 3650):
    t = np.arange(0, days + 1)
    per_shell = np.zeros((len(states), len(t)))

    for i, s in enumerate(states):
        if s.excess_fraction == 0:
            per_shell[i, :] = 1.0
        else:
            per_shell[i, :] = np.exp(-s.lambda_daily * t)

    overall = per_shell.mean(axis=0)
    return t, overall, per_shell

def run_scenario(name: str, counts: List[int], k: float = 3.0, days: int = 3650):
    states = compute_states(counts, k)
    t, overall, per_shell = stability_time_series(states, days)

    def first_below(series, thresh):
        idx = np.where(series < thresh)[0]
        return int(idx[0]) if idx.size > 0 else None

    return {
        "name": name,
        "states": states,
        "counts": counts,
        "t": t,
        "overall": overall,
        "per_shell": per_shell,
        "day0_overall": float(overall[0]),
        "dayN_overall": float(overall[-1]),
        "first_below_0.9": first_below(overall, 0.9),
        "first_below_0.75": first_below(overall, 0.75),
        "first_below_0.5": first_below(overall, 0.5),
    }
