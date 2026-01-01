

import numpy as np
from numpy.random import Generator

from .models import VolcanoEvent, WorldState




def update_volcanic_activity(state: WorldState, rng: np.random.Generator) -> WorldState:
    """
    Simple volcanic model:
    - Low probability random eruptions in low or mid elevation zones.
    - Eruptions raise local elevation slightly.
    """
    # Ensure surface is a numpy ndarray for type checkers
    surface: np.ndarray = state.surface_height.copy()  # type: ignore[attr-defined]
    h_raw, w_raw = surface.shape  # type: ignore[attr-defined]
    h = int(h_raw) # pyright: ignore[reportUnknownArgumentType]
    w = int(w_raw) # pyright: ignore[reportUnknownArgumentType]

    # Probability of eruption per day
    eruption_prob = 0.01



    if rng.uniform() < eruption_prob: # pyright: ignore[reportUnknownMemberType]
        i = int(rng.integers(low=0, high=h, dtype=int)) # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
        j = int(rng.integers(low=0, high=w, dtype=int)) # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
        low, high = 0.01, 0.05
        strength = float(rng.uniform(low, high)) # pyright: ignore[reportUnknownArgumentType]
        _apply_eruption(surface, i, j, strength) # pyright: ignore[reportUnknownArgumentType]
        state.volcano_events.append(VolcanoEvent(day=state.day, position=(i, j), strength=strength))

    # Renormalize
    surface = np.asarray(surface)  # type: ignore[attr-defined]
    minv = float(np.min(surface)) # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
    surface = surface - minv # pyright: ignore[reportUnknownVariableType]
    maxval = float(np.max(surface)) # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
    maxv = max(1e-9, maxval)
    surface = surface / maxv # pyright: ignore[reportUnknownVariableType]
    state.surface_height = surface
    return state

def _apply_eruption(surface: np.ndarray, i: int, j: int, strength: float, radius: int = 3) -> None: # pyright: ignore[reportUnknownParameterType]
    h, w = surface.shape # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    for di in range(-radius, radius+1):
        for dj in range(-radius, radius+1):
            ii, jj = i + di, j + dj
            if 0 <= ii < h and 0 <= jj < w:
                dist = (di**2 + dj**2) ** 0.5
                if dist <= radius:
                    surface[ii, jj] += strength * (1 - dist / (radius + 1e-9))
