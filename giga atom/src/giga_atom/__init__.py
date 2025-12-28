"""
giga_atom — Space Leaf Corp internal simulation package
Created by Leif William Sogge

All rights reserved under the Space Leaf Corp Proprietary License 1.0.
Unauthorized use, redistribution, or weaponization is strictly prohibited.
"""

from .model import ShellState, compute_excess_fraction, compute_lambda
from .simulations import (
    PHYSICAL_SHELLS,
    TOTAL_ELECTRONS,
    baseline_counts,
    apply_overload,
    compute_states,
    stability_time_series,
    run_scenario,
)
from .gyroscope import gyroscopic_S_time_series # pyright: ignore[reportUnknownVariableType]
from .bonding import (
    valence_from_Z,
    bonding_score,
    bonding_scan,
    combined_stability_for_partner,
)
from . import utils

__all__ = [
    "ShellState",
    "compute_excess_fraction",
    "compute_lambda",
    "PHYSICAL_SHELLS",
    "TOTAL_ELECTRONS",
    "baseline_counts",
    "apply_overload",
    "compute_states",
    "stability_time_series",
    "run_scenario",
    "gyroscopic_S_time_series",
    "valence_from_Z",
    "bonding_score",
    "bonding_scan",
    "combined_stability_for_partner",
    "utils",
]
