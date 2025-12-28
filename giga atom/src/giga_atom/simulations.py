# src/giga_atom/simulations.py
"""
Giga Atom simulation suite
- Baseline stability (2n^2 shells)
- Stress scenarios and hardening
- Gyroscopic tests
- Bonding scan across Z=1..118
- Monte Carlo robustness trials

Author: Leif William Sogge (adapted for Space Leaf Corp private repo)
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import random

# -------------------------
# Model definitions
# -------------------------
PHYSICAL_SHELLS = [2, 8, 18, 32, 50, 72, 98, 128]  # 2n^2 rule for shells 1..8
TOTAL_ELECTRONS = sum(PHYSICAL_SHELLS)

@dataclass
class ShellState:
    index: int
    physical_max: int
    count: int
    excess_fraction: float = 0.0
    lambda_daily: float = 0.0

def compute_excess_and_lambda(counts: List[int], k: float = 3.0) -> List[ShellState]:
    """Compute excess fraction and daily lambda for each shell."""
    states = []
    for i, (phys, cnt) in enumerate(zip(PHYSICAL_SHELLS, counts), start=1):
        if cnt <= phys:
            excess = 0.0
        else:
            excess = (cnt - phys) / phys
        # lambda per day (using 365-day normalization)
        lam = k * excess / 365.0
        states.append(ShellState(index=i, physical_max=phys, count=cnt,
                                 excess_fraction=excess, lambda_daily=lam))
    return states

def stability_time_series(states: List[ShellState], days: int = 3650) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute per-shell and overall stability time series.
    Per-shell stability(t) = exp(-lambda * t)
    Overall stability(t) = mean(per-shell stabilities)
    """
    t = np.arange(0, days + 1)  # inclusive
    per_shell = np.zeros((len(states), len(t)))
    for i, s in enumerate(states):
        if s.excess_fraction == 0:
            per_shell[i, :] = 1.0
        else:
            per_shell[i, :] = np.exp(-s.lambda_daily * t)
    overall = per_shell.mean(axis=0)
    return t, overall, per_shell

# -------------------------
# Scenario runners
# -------------------------
def baseline_counts() -> List[int]:
    return PHYSICAL_SHELLS.copy()

def apply_overload(base_counts: List[int], overloads: Dict[int, float]) -> List[int]:
    """
    overloads: mapping shell_index (1-based) -> fractional overload (e.g., 0.5 for +50%)
    Returns new counts (rounded to int)
    """
    counts = base_counts.copy()
    for idx, frac in overloads.items():
        i = idx - 1
        counts[i] = int(round(counts[i] * (1.0 + frac)))
    return counts

def run_scenario(name: str, counts: List[int], k: float = 3.0, days: int = 3650):
    states = compute_excess_and_lambda(counts, k=k)
    t, overall, per_shell = stability_time_series(states, days=days)
    # threshold crossing days
    def first_below(series, thresh):
        idx = np.where(series < thresh)[0]
        return int(idx[0]) if idx.size > 0 else None
    cross_09 = first_below(overall, 0.9)
    cross_075 = first_below(overall, 0.75)
    cross_05 = first_below(overall, 0.5)
    summary = {
        "name": name,
        "counts": counts,
        "states": states,
        "day0_overall": float(overall[0]),
        "dayN_overall": float(overall[-1]),
        "first_below_0.9": cross_09,
        "first_below_0.75": cross_075,
        "first_below_0.5": cross_05,
        "t": t,
        "overall": overall,
        "per_shell": per_shell
    }
    return summary

# -------------------------
# Gyroscopic model
# -------------------------
def gyroscopic_S_time_series(counts: List[int], days: int = 3650,
                             omega0: float = 10.0, c: float = 0.5,
                             impulse_every: Optional[int] = None, impulse_drop: float = 0.05):
    """
    Compute S(t) = omega(t)/omega0 with damping and optional periodic impulses.
    c is the dimensionless damping applied over the full period (higher -> faster decay).
    If impulse_every is set (days), apply instantaneous multiplicative drop at those days.
    """
    # compute moment of inertia (not used in normalized S but kept for completeness)
    radii = [(n+1)**2 for n in range(len(counts))]  # n from 0..7 -> shell index n+1
    masses = [cnt * 1e-3 for cnt in counts]
    I_n = [(2.0/3.0) * m * (r**2) for m, r in zip(masses, radii)]
    I_total = sum(I_n)
    t = np.arange(0, days + 1)
    # continuous damping factor per day: exp(-c * t / days)
    omega = omega0 * np.exp(-c * t / days)
    if impulse_every:
        for day in range(impulse_every, days + 1, impulse_every):
            omega[day:] *= (1.0 - impulse_drop)
    S = omega / omega0
    return t, S, I_total

# -------------------------
# Bonding heuristics
# -------------------------
def valence_from_Z(Z: int) -> int:
    """
    Compute valence as electrons in outermost occupied shell using 2n^2 filling.
    This is a simplified heuristic: fill shells sequentially with PHYSICAL_SHELLS capacities.
    """
    remaining = Z
    last_filled = 0
    for cap in PHYSICAL_SHELLS:
        if remaining <= cap:
            last_filled = remaining
            break
        remaining -= cap
    else:
        # if Z > sum(PHYSICAL_SHELLS), approximate valence as remainder
        last_filled = remaining
    return int(last_filled)

def bonding_score(val_giga: int, val_elem: int) -> float:
    """Heuristic bonding score normalized to 0..1."""
    if val_giga <= 0 or val_elem <= 0:
        return 0.0
    base = min(val_giga, val_elem)
    mismatch = abs(val_elem - val_giga) / max(1, val_elem + val_giga)
    raw = base * (1.0 - mismatch)
    # normalize by max possible (which is base when mismatch=0)
    norm = raw / base if base > 0 else 0.0
    return float(norm)

def bonding_scan(val_giga: int = 1, maxZ: int = 118):
    rows = []
    for Z in range(1, maxZ + 1):
        val = valence_from_Z(Z)
        score = bonding_score(val_giga, val)
        rows.append({"Z": Z, "valence": val, "bonding_score": score})
    df = pd.DataFrame(rows)
    return df

def combined_stability_for_partner(val_giga: int, val_elem: int, days: int = 3650, d: float = 2.0):
    mismatch = abs(val_giga - val_elem) / max(1, val_elem + val_giga)
    t = np.arange(0, days + 1)
    stability = np.exp(-d * mismatch * t / 365.0)
    return t, stability

# -------------------------
# Monte Carlo harness
# -------------------------
def monte_carlo_trials(n_trials: int = 200, days: int = 3650, seed: int = 42):
    random.seed(seed)
    np.random.seed(seed)
    results = []
    stress_scenarios = {
        "A": {8: 0.5},
        "B": {7: 0.3, 8: 0.8},
        "C": {5: 0.2, 6: 0.2, 7: 0.2, 8: 0.2}
    }
    for trial in range(n_trials):
        scenario_key = random.choice(list(stress_scenarios.keys()))
        base = baseline_counts()
        overloads = stress_scenarios[scenario_key].copy()
        # apply random noise to each overload fraction
        noisy_overloads = {}
        for idx, frac in overloads.items():
            noise = np.random.normal(0.0, 0.1)
            noise = np.clip(noise, -0.2, 0.5)
            noisy_overloads[idx] = frac * (1.0 + noise)
        counts = apply_overload(base, noisy_overloads)
        # vary k by ±20%
        k = 3.0 * (1.0 + np.random.uniform(-0.2, 0.2))
        summary = run_scenario(f"MC_{trial}_{scenario_key}", counts, k=k, days=days)
        # record crossing days
        results.append({
            "trial": trial,
            "scenario": scenario_key,
            "k": k,
            "counts": counts,
            "day0": summary["day0_overall"],
            "dayN": summary["dayN_overall"],
            "first_below_0.9": summary["first_below_0.9"],
            "first_below_0.75": summary["first_below_0.75"],
            "first_below_0.5": summary["first_below_0.5"]
        })
    df = pd.DataFrame(results)
    return df

# -------------------------
# Plot helpers
# -------------------------
def plot_overall(t: np.ndarray, series: Dict[str, np.ndarray], title: str, filename: Optional[str] = None):
    plt.figure(figsize=(10, 5))
    for label, arr in series.items():
        plt.plot(t, arr, label=label)
    plt.xlabel("Days")
    plt.ylabel("Overall Stability")
    plt.title(title)
    plt.legend()
    plt.grid(True)
    if filename:
        plt.savefig(filename, dpi=150)
    plt.show()

def plot_gyroscopic(t: np.ndarray, S_dict: Dict[str, np.ndarray], filename: Optional[str] = None):
    plt.figure(figsize=(10, 5))
    for label, S in S_dict.items():
        plt.plot(t, S, label=label)
    plt.xlabel("Days")
    plt.ylabel("Normalized Angular Momentum S(t)")
    plt.title("Gyroscopic Stability Over Time")
    plt.legend()
    plt.grid(True)
    if filename:
        plt.savefig(filename, dpi=150)
    plt.show()

# -------------------------
# Example run (main)
# -------------------------
def main_run():
    days = 3650
    # Baseline
    base = baseline_counts()
    baseline = run_scenario("Baseline", base, k=3.0, days=days)

    # Stress scenarios
    A_counts = apply_overload(base, {8: 0.5})
    B_counts = apply_overload(base, {7: 0.3, 8: 0.8})
    C_counts = apply_overload(base, {5: 0.2, 6: 0.2, 7: 0.2, 8: 0.2})

    A = run_scenario("Stress_A_shell8_+50%", A_counts, k=3.0, days=days)
    B = run_scenario("Stress_B_7+30_8+80", B_counts, k=3.0, days=days)
    C = run_scenario("Stress_C_5-8_+20%", C_counts, k=3.0, days=days)

    # Hardening
    A_H1 = run_scenario("A_H1_k=1.5", A_counts, k=1.5, days=days)
    A_H2_counts = apply_overload(base, {8: 0.5 * 0.6})  # reduce overload by 40%
    A_H2 = run_scenario("A_H2_reduce_overload40%", A_H2_counts, k=3.0, days=days)

    # Gyroscopic cases
    t_g0, S_g0, I0 = gyroscopic_S_time_series(base, days=days, c=0.5, impulse_every=None)
    t_g1, S_g1, I1 = gyroscopic_S_time_series(base, days=days, c=1.5, impulse_every=None)
    t_g2, S_g2, I2 = gyroscopic_S_time_series(base, days=days, c=0.5, impulse_every=365, impulse_drop=0.05)

    # Bonding scan
    bonding_df = bonding_scan(val_giga=1, maxZ=118)
    # filter non-radioactive (Z <= 83)
    bonding_nonradio = bonding_df[bonding_df["Z"] <= 83].copy()
    top10 = bonding_nonradio.sort_values("bonding_score", ascending=False).head(10)
    top5 = top10.head(5)

    # Combined stability for top5 partners
    combined = {}
    for _, row in top5.iterrows():
        Z = int(row["Z"])
        val = int(row["valence"])
        t_c, stab = combined_stability_for_partner(1, val, days=days, d=2.0)
        combined[f"Z{Z}_val{val}"] = stab

    # Monte Carlo
    mc_df = monte_carlo_trials(n_trials=200, days=days, seed=42)

    # Print concise numeric summaries
    print("=== SUMMARY ===")
    print(f"Total electrons (baseline): {TOTAL_ELECTRONS}")
    print(f"Baseline Day0 overall: {baseline['day0_overall']:.3f}, Day{days} overall: {baseline['dayN_overall']:.3f}")
    for s in [A, B, C]:
        print(f"{s['name']}: Day0={s['day0_overall']:.3f}, Day{days}={s['dayN_overall']:.3f}, first_below_0.5={s['first_below_0.5']}")

    # Export CSV-style tables
    # Per-scenario per-shell Day3650 stabilities
    def export_scenario_table(summary):
        states = summary["states"]
        t = summary["t"]
        per_shell = summary["per_shell"]
        last = per_shell[:, -1]
        rows = []
        for s, val in zip(states, last):
            rows.append({
                "scenario": summary["name"],
                "shell": s.index,
                "physical_max": s.physical_max,
                "count": s.count,
                "excess_fraction": s.excess_fraction,
                "lambda_daily": s.lambda_daily,
                "day3650_stability": float(val)
            })
        return pd.DataFrame(rows)

    df_list = []
    for summary in [baseline, A, B, C, A_H1, A_H2]:
        df_list.append(export_scenario_table(summary))
    per_shell_df = pd.concat(df_list, ignore_index=True)

    # Top partners table
    top_partners_df = top10.copy()

    # Monte Carlo summary
    mc_summary = {
        "median_first_below_0.5": int(mc_df["first_below_0.5"].dropna().median()) if mc_df["first_below_0.5"].notna().any() else None,
        "pct_survive_above_0.5": float((mc_df["first_below_0.5"].isna().sum()) / len(mc_df))
    }

    # Save CSVs locally (caller can upload to repo)
    per_shell_df.to_csv("per_shell_summaries.csv", index=False)
    top_partners_df.to_csv("top_partners.csv", index=False)
    mc_df.to_csv("monte_carlo_trials.csv", index=False)

    # Plots
    plot_overall(baseline["t"], {"Baseline": baseline["overall"]}, "Baseline Overall Stability", filename="baseline_overall.png")
    plot_overall(A["t"], {"Stress A": A["overall"], "A_H1_k1.5": A_H1["overall"], "A_H2_reduced": A_H2["overall"]}, "Stress A and Hardened Variants", filename="stressA_compare.png")
    plot_gyroscopic(t_g0, {"G0_c0.5": S_g0, "G1_c1.5": S_g1, "G2_impulses": S_g2}, filename="gyroscopic.png")

    # Bonding score plot
    plt.figure(figsize=(10,4))
    plt.plot(bonding_df["Z"], bonding_df["bonding_score"], label="bonding_score")
    plt.xlabel("Atomic Number Z")
    plt.ylabel("Bonding Score (normalized)")
    plt.title("Bonding Score vs Atomic Number (Giga Atom valence=1)")
    plt.grid(True)
    plt.savefig("bonding_score_vs_Z.png", dpi=150)
    plt.show()

    # Combined stability for top5 partners
    plt.figure(figsize=(10,5))
    for label, arr in combined.items():
        plt.plot(t_c, arr, label=label)
    plt.xlabel("Days")
    plt.ylabel("Combined Stability")
    plt.title("Combined Stability for Top 5 Partners")
    plt.legend()
    plt.grid(True)
    plt.savefig("combined_top5.png", dpi=150)
    plt.show()

    # Monte Carlo survival curve (fraction above 0.5 over time)
    # approximate by checking each trial's day of first below 0.5
    survival = []
    for day in range(0, days+1, 10):
        alive = (mc_df["first_below_0.5"].isna()) | (mc_df["first_below_0.5"] > day)
        survival.append(alive.mean())
    plt.figure(figsize=(10,5))
    plt.plot(np.arange(0, days+1, 10), survival)
    plt.xlabel("Days")
    plt.ylabel("Fraction surviving above 0.5")
    plt.title("Monte Carlo Survival Curve (above 0.5)")
    plt.grid(True)
    plt.savefig("mc_survival_curve.png", dpi=150)
    plt.show()

    # Return key objects for programmatic use
    return {
        "baseline": baseline,
        "A": A, "B": B, "C": C,
        "A_H1": A_H1, "A_H2": A_H2,
        "gyroscopic": {"g0": (t_g0, S_g0), "g1": (t_g1, S_g1), "g2": (t_g2, S_g2)},
        "bonding_df": bonding_df,
        "top_partners": top_partners_df,
        "monte_carlo": mc_df,
        "per_shell_df": per_shell_df,
        "mc_summary": mc_summary
    }

if __name__ == "__main__":
    results = main_run()
