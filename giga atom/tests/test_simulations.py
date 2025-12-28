"""
tests/test_simulations.py
Space Leaf Corp — Internal Use Only
All rights reserved under the Space Leaf Corp Proprietary License 1.0.

These tests verify:
- Baseline stability behavior
- Stress scenario behavior
- Hardening scenario improvements
- Gyroscopic model consistency
- Bonding scan validity
- Monte Carlo harness integrity
"""

import numpy as np
from giga_atom.simulations import (
    baseline_counts,
    apply_overload,
    run_scenario,
    gyroscopic_S_time_series,
    bonding_scan,
    combined_stability_for_partner,
    monte_carlo_trials,
    TOTAL_ELECTRONS
)

# ---------------------------------------------------------
# Baseline Tests
# ---------------------------------------------------------

def test_baseline_counts_correct():
    """Baseline should match the physical 2n^2 shell model."""
    counts = baseline_counts()
    assert counts == [2, 8, 18, 32, 50, 72, 98, 128]

def test_baseline_stability_constant():
    """Baseline stability should remain 1.0 for all 3650 days."""
    summary = run_scenario("Baseline", baseline_counts(), k=3.0, days=3650)
    assert summary["day0_overall"] == 1.0
    assert summary["dayN_overall"] == 1.0
    assert summary["first_below_0.5"] is None

# ---------------------------------------------------------
# Stress Scenario Tests
# ---------------------------------------------------------

def test_stress_A_reduces_stability():
    """Stress A (shell 8 +50%) should reduce long-term stability."""
    base = baseline_counts()
    counts = apply_overload(base, {8: 0.5})
    summary = run_scenario("Stress A", counts, k=3.0, days=3650)
    assert summary["dayN_overall"] < 1.0
    assert summary["first_below_0.5"] is not None

def test_stress_B_more_severe_than_A():
    """Stress B should degrade stability faster than Stress A."""
    base = baseline_counts()
    A_counts = apply_overload(base, {8: 0.5})
    B_counts = apply_overload(base, {7: 0.3, 8: 0.8})

    A = run_scenario("A", A_counts, k=3.0, days=3650)
    B = run_scenario("B", B_counts, k=3.0, days=3650)

    assert B["dayN_overall"] < A["dayN_overall"]

# ---------------------------------------------------------
# Hardening Tests
# ---------------------------------------------------------

def test_hardening_improves_stress_A():
    """Hardening should improve Stress A stability."""
    base = baseline_counts()
    A_counts = apply_overload(base, {8: 0.5})

    A = run_scenario("A", A_counts, k=3.0, days=3650)
    A_H1 = run_scenario("A_H1", A_counts, k=1.5, days=3650)

    assert A_H1["dayN_overall"] > A["dayN_overall"]

# ---------------------------------------------------------
# Gyroscopic Tests
# ---------------------------------------------------------

def test_gyroscopic_S_monotonic():
    """Gyroscopic S(t) should monotonically decrease under damping."""
    base = baseline_counts()
    t, S, _ = gyroscopic_S_time_series(base, days=3650, c=0.5)
    assert S[0] == 1.0
    assert S[-1] < 1.0
    assert np.all(np.diff(S) <= 1e-9)  # allow tiny float noise

# ---------------------------------------------------------
# Bonding Tests
# ---------------------------------------------------------

def test_bonding_scan_valid():
    """Bonding scan should return 118 rows and valid scores."""
    df = bonding_scan(val_giga=1)
    assert len(df) == 118
    assert df["bonding_score"].between(0, 1).all()

def test_combined_stability_monotonic():
    """Combined stability with a partner should decay monotonically."""
    t, stab = combined_stability_for_partner(1, 2, days=3650)
    assert stab[0] == 1.0
    assert stab[-1] <= 1.0
    assert np.all(np.diff(stab) <= 1e-9)

# ---------------------------------------------------------
# Monte Carlo Tests
# ---------------------------------------------------------

def test_monte_carlo_runs():
    """Monte Carlo should produce the correct number of trials."""
    df = monte_carlo_trials(n_trials=20, days=3650, seed=123)
    assert len(df) == 20
    assert "first_below_0.5" in df.columns
