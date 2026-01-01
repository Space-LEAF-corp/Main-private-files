from typing import List, Tuple
import numpy as np
from .models import Plate

def initialize_plates(width: int, height: int, plate_count: int, rng: np.random.Generator) -> List[Plate]:
    plates = []
    for pid in range(plate_count):
        angle = rng.uniform(0, 2 * np.pi)
        speed = rng.uniform(0.001, 0.01)
        vx, vy = speed * np.cos(angle), speed * np.sin(angle)
        is_oceanic = bool(rng.integers(0, 2))
        plates.append(Plate(id=pid, velocity=(vx, vy), is_oceanic=is_oceanic))
    return plates

def update_surface_for_tectonics(surface: np.ndarray, plates: List[Plate], day: int) -> np.ndarray:
    """
    Placeholder tectonic deformation:
    - Slight uplift or subsidence based on plate ids and pseudo boundaries.
    This is where a real plate-boundary computation would go.
    """
    h, w = surface.shape
    deformed = surface.copy()
    # Simple pattern: uplift along a diagonal "boundary"
    boundary_idx = (np.arange(h)[:, None] + np.arange(w)[None, :]) % 20
    uplift_mask = boundary_idx < 2
    subsidence_mask = boundary_idx > 17

    uplift_amount = 0.0005
    subsidence_amount = 0.0003

    deformed[uplift_mask] += uplift_amount
    deformed[subsidence_mask] -= subsidence_amount

    # Clamp back to [0, 1]
    deformed -= deformed.min()
    maxv = max(deformed.max(), 1e-9)
    deformed /= maxv
    return deformed
