from typing import List
import numpy as np
from .models import VolcanoEvent, WorldState

def update_volcanic_activity(state: WorldState, rng: np.random.Generator) -> WorldState:
    """
    Simple volcanic model:
    - Low probability random eruptions in low or mid elevation zones.
    - Eruptions raise local elevation slightly.
    """
    surface = state.surface_height.copy()
    h, w = surface.shape

    # Probability of eruption per day
    eruption_prob = 0.01

    if rng.uniform() < eruption_prob:
        i = int(rng.integers(0, h))
        j = int(rng.integers(0, w))
        strength = float(rng.uniform(0.01, 0.05))
        _apply_eruption(surface, i, j, strength)
        state.volcano_events.append(VolcanoEvent(day=state.day, position=(i, j), strength=strength))

    # Renormalize
    surface -= surface.min()
    maxv = max(surface.max(), 1e-9)
    surface /= maxv
    state.surface_height = surface
    return state

def _apply_eruption(surface: np.ndarray, i: int, j: int, strength: float, radius: int = 3) -> None:
    h, w = surface.shape
    for di in range(-radius, radius+1):
        for dj in range(-radius, radius+1):
            ii, jj = i + di, j + dj
            if 0 <= ii < h and 0 <= jj < w:
                dist = (di**2 + dj**2) ** 0.5
                if dist <= radius:
                    surface[ii, jj] += strength * (1 - dist / (radius + 1e-9))
