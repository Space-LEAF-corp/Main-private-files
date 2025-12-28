"""
model.py — Core data structures for the Giga Atom
Space Leaf Corp — Internal Use Only
"""

from dataclasses import dataclass

@dataclass
class ShellState:
    index: int
    physical_max: int
    count: int
    excess_fraction: float
    lambda_daily: float

def compute_excess_fraction(count: int, physical_max: int) -> float:
    """Return excess fraction above physical shell capacity."""
    if count <= physical_max:
        return 0.0
    return (count - physical_max) / physical_max

def compute_lambda(excess_fraction: float, k: float = 3.0) -> float:
    """Daily decay constant."""
    return (k * excess_fraction) / 365.0
