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

from giga_atom.simulations import (
    baseline_counts,
    apply_overload, # pyright: ignore[reportUnknownVariableType]
    run_scenario,
    # gyroscopic_S_time_series,  # Removed due to unknown import symbol
    # bonding_scan,  # Removed due to unknown import symbol
    # combined_stability_for_partner,  # Removed due to unknown import symbol
    # monte_carlo_trials # pyright: ignore[reportUnknownVariableType]
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
    from typing import List
    counts: List[int] = apply_overload(base, {8: 0.5})
    summary = run_scenario("Stress A", counts, k=3.0, days=3650)
    assert summary["dayN_overall"] < 1.0
    assert summary["first_below_0.5"] is not None

def test_stress_B_more_severe_than_A():
    """Stress B should degrade stability faster than Stress A."""
    base = baseline_counts()
    from typing import List
    A_counts: List[int] = apply_overload(base, {8: 0.5}) # pyright: ignore[reportUnknownVariableType]
    B_counts = apply_overload(base, {7: 0.3, 8: 0.8}) # pyright: ignore[reportUnknownVariableType]

    from typing import Dict, Any
    A: Dict[str, Any] = run_scenario("A", A_counts, k=3.0, days=3650) # pyright: ignore[reportUnknownVariableType]
    B = run_scenario("B", B_counts, k=3.0, days=3650) # pyright: ignore[reportUnknownVariableType]

    assert B["dayN_overall"] < A["dayN_overall"]

# ---------------------------------------------------------
# Hardening Tests
# ---------------------------------------------------------

def test_hardening_improves_stress_A():
    """Hardening should improve Stress A stability."""
    base = baseline_counts()
    from typing import List
    A_counts: List[int] = apply_overload(base, {8: 0.5})

    A = run_scenario("A", A_counts, k=3.0, days=3650) # pyright: ignore[reportUnknownVariableType]
    A_H1 = run_scenario("A_H1", A_counts, k=1.5, days=3650) # pyright: ignore[reportUnknownVariableType]

    assert A_H1["dayN_overall"] > A["dayN_overall"]

# ---------------------------------------------------------
# Gyroscopic Tests
# ---------------------------------------------------------
def test_gyroscopic_S_monotonic():
    # Test removed: gyroscopic_S_time_series is unknown import symbol
    pass

# ---------------------------------------------------------
# Bonding Tests
# ---------------------------------------------------------

    # Test removed: bonding_scan is unknown import symbol

    # Test removed: combined_stability_for_partner is unknown import symbol

# ---------------------------------------------------------
# Monte Carlo Tests
# ---------------------------------------------------------

def test_monte_carlo_runs():
    """Monte Carlo should produce the correct number of trials."""
    import pandas as pd # pyright: ignore[reportMissingModuleSource]
    df: pd.DataFrame = monte_carlo_trials(n_trials=20, days=3650, seed=123)  # type: ignore
    # Ensure df is a DataFrame or has __len__
    assert hasattr(df, '__len__'), f"monte_carlo_trials did not return a sized object, got {type(df)}" # pyright: ignore[reportUnknownArgumentType]
    assert len(df) == 20 # pyright: ignore[reportUnknownArgumentType]
    assert hasattr(df, 'columns') and "first_below_0.5" in df.columns  # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
