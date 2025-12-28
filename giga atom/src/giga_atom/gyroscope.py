"""
gyroscope.py — Gyroscopic stability model
Space Leaf Corp — Internal Use Only
"""

import numpy as np

def gyroscopic_S_time_series(
    counts,
    days: int = 3650,
    omega0: float = 10.0,
    c: float = 0.5,
    impulse_every=None,
    impulse_drop=0.05
):
    radii = [(i + 1) ** 2 for i in range(len(counts))]
    masses = [cnt * 1e-3 for cnt in counts]
    I_n = [(2/3) * m * (r**2) for m, r in zip(masses, radii)]
    I_total = sum(I_n)

    t = np.arange(0, days + 1)
    omega = omega0 * np.exp(-c * t / days)

    if impulse_every:
        for day in range(impulse_every, days + 1, impulse_every):
            omega[day:] *= (1 - impulse_drop)

    S = omega / omega0
    return t, S, I_total
