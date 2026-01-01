"""
Ceremonial Terrain v1 – Map Generator Prototype

This package implements:
- Terrain generation (surface + ocean floor)
- Tectonic + volcanic evolution over time
- Basic security scaffolding for projects
"""

from .config import SimulationConfig
from .models import WorldState
from .simulation import run_365_day_simulation, validate_forward_backward

__all__ = [
    "SimulationConfig",
    "WorldState",
    "run_365_day_simulation",
    "validate_forward_backward",
]
