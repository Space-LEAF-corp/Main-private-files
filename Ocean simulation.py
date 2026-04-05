#!/usr/bin/env python3
"""
Simple ocean salinity dilution simulator for hydrogen engine water discharge.

Compares a small-scale engine vs a large industrial-scale engine over 365 days.
"""

from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class EngineScenario:
    name: str
    daily_discharge_m3: float  # m^3/day of water from hydrogen combustion


@dataclass
class OceanBox:
    volume_m3: float          # Local ocean box volume
    salinity_initial: float   # PSU
    salinity_background: float  # PSU (outside water)
    flushing_fraction: float  # fraction of box exchanged per day (0-1)


def simulate_salinity(
    ocean: OceanBox,
    engine: EngineScenario,
    discharge_salinity_psu: float = 0.0,
    days: int = 365
) -> List[float]:
    """
    Run a time-step simulation of salinity in a local ocean box.

    Returns a list of salinity values (PSU) for each day, including day 0.
    """
    S = ocean.salinity_initial
    V = ocean.volume_m3
    S_bg = ocean.salinity_background
    f = ocean.flushing_fraction
    Q = engine.daily_discharge_m3
    S_in = discharge_salinity_psu

    history = [S]

    for _ in range(days):
        # Step 1: add discharge water
        S_prime = (S * V + S_in * Q) / (V + Q)

        # Step 2: flushing/mixing with background
        S_next = (1.0 - f) * S_prime + f * S_bg

        history.append(S_next)
        S = S_next

    return history


def compare_scenarios(
    ocean: OceanBox,
    small_engine: EngineScenario,
    large_engine: EngineScenario,
    days: int = 365
) -> Tuple[List[float], List[float]]:
    """
    Simulate both small and large engine scenarios and return their salinity histories.
    """
    small_hist = simulate_salinity(ocean, small_engine, days=days)
    large_hist = simulate_salinity(ocean, large_engine, days=days)
    return small_hist, large_hist


def summarize_results(
    salinity_history: List[float],
    S0: float
) -> dict:
    """
    Compute summary metrics for a salinity time series.
    """
    S_final = salinity_history[-1]
    delta_S = S_final - S0
    dilution_ratio = (S0 - S_final) / S0 if S0 != 0 else 0.0
    min_S = min(salinity_history)
    max_S = max(salinity_history)

    return {
        "S_final": S_final,
        "delta_S": delta_S,
        "dilution_ratio": dilution_ratio,
        "min_S": min_S,
        "max_S": max_S,
    }


def main():
    # --- Define the ocean box ---
    # Example: a small coastal cell or bay segment
    ocean = OceanBox(
        volume_m3=1e9,          # 1 billion m^3
        salinity_initial=35.0,  # PSU
        salinity_background=35.0,
        flushing_fraction=0.05  # 5% of the box exchanged per day
    )

    # --- Define engine scenarios ---
    # These are example values; tweak them as you like.
    small_engine = EngineScenario(
        name="Small hydrogen engine",
        daily_discharge_m3=0.1  # 0.1 m^3/day (~100 kg/day of water)
    )

    large_engine = EngineScenario(
        name="Large industrial hydrogen plant",
        daily_discharge_m3=100.0  # 100 m^3/day (~100,000 kg/day of water)
    )

    days = 365

    # --- Run simulations ---
    small_hist, large_hist = compare_scenarios(ocean, small_engine, large_engine, days=days)

    # --- Summaries ---
    S0 = ocean.salinity_initial
    small_summary = summarize_results(small_hist, S0)
    large_summary = summarize_results(large_hist, S0)

    print("=== Ocean Salinity Dilution Simulation ===")
    print(f"Days simulated: {days}")
    print(f"Initial salinity: {S0:.4f} PSU")
    print()

    print(f"Scenario: {small_engine.name}")
    for k, v in small_summary.items():
        print(f"  {k}: {v:.6f}")
    print()

    print(f"Scenario: {large_engine.name}")
    for k, v in large_summary.items():
        print(f"  {k}: {v:.6f}")
    print()

    # Optional: show a few sample days
    sample_days = [0, 1, 7, 30, 180, 365]
    print("\nSample salinity values (PSU):")
    print("Day | Small_engine | Large_engine")
    for d in sample_days:
        if d <= days:
            print(f"{d:3d} | {small_hist[d]:12.6f} | {large_hist[d]:12.6f}")


if __name__ == "__main__":
    main()