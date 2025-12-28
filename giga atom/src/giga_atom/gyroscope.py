"""
gyroscope.py — Gyroscopic stability model
Space Leaf Corp — Internal Use Only
"""

from typing import Optional
import numpy as np # pyright: ignore[reportMissingImports]

def gyroscopic_S_time_series( # pyright: ignore[reportUnknownParameterType]
    counts: list[float],
    days: int = 3650,
    omega0: float = 10.0,
    c: float = 0.5,
    impulse_every: 'Optional[int]' = None,
    impulse_drop: float = 0.05
):
    radii = [(i + 1) ** 2 for i in range(len(counts))]
    masses = [cnt * 1e-3 for cnt in counts]
    I_n = [(2/3) * m * (r**2) for m, r in zip(masses, radii)]
    I_total = sum(I_n)

    t = np.arange(0, days + 1) # pyright: ignore[reportUnknownVariableType, reportUnknownMemberType]
    omega: np.ndarray = omega0 * np.exp(-c * t / days) # pyright: ignore[reportUnknownVariableType, reportUnknownMemberType]

    if impulse_every is not None:
        try:
            impulse_every_int = int(impulse_every)
        except (TypeError, ValueError):
            raise ValueError("impulse_every must be an integer or convertible to int")
        for day in range(impulse_every_int, days + 1, impulse_every_int):
            omega[day:] *= (1 - impulse_drop)

    S = omega / omega0 # pyright: ignore[reportUnknownVariableType]
    return t, S, I_total # pyright: ignore[reportUnknownVariableType]
