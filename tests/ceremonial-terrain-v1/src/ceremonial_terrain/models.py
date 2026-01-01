from dataclasses import dataclass, field
from typing import List, Tuple, Dict
import numpy as np

@dataclass
class Plate:
    id: int
    velocity: Tuple[float, float]  # (vx, vy)
    is_oceanic: bool

@dataclass
class VolcanoEvent:
    day: int
    position: Tuple[int, int]
    strength: float

@dataclass
class WorldState:
    day: int
    # Elevations: surface and seafloor
    surface_height: np.ndarray  # shape: (H, W)
    seafloor_height: np.ndarray  # shape: (H, W)

    # Tectonic plates
    plates: List[Plate]

    # Volcanic history
    volcano_events: List[VolcanoEvent] = field(default_factory=list)

Snapshot = WorldState  # alias for clearer semantics
