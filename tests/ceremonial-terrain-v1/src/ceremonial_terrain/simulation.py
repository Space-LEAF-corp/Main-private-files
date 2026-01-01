from typing import Dict, List
import numpy as np

from .config import SimulationConfig
from .models import WorldState, Snapshot
from .terrain import generate_initial_heightmap, generate_seafloor_from_surface
from .tectonics import initialize_plates, update_surface_for_tectonics
from .volcanoes import update_volcanic_activity
from .ocean import update_seafloor_for_tectonics_and_volcanoes

def run_365_day_simulation(config: SimulationConfig) -> Dict[int, Snapshot]:
    rng = np.random.default_rng(config.seed)

    surface = generate_initial_heightmap(config.width, config.height, rng)
    seafloor = generate_seafloor_from_surface(surface)
    plates = initialize_plates(config.width, config.height, config.plate_count, rng)

    snapshots: Dict[int, Snapshot] = {}
    state = WorldState(day=0, surface_height=surface, seafloor_height=seafloor, plates=plates)

    snapshots[0] = _copy_state(state)

    for day in range(1, config.days + 1):
        state.day = day
        state.surface_height = update_surface_for_tectonics(state.surface_height, state.plates, day)
        state = update_volcanic_activity(state, rng)
        state.seafloor_height = update_seafloor_for_tectonics_and_volcanoes(
            state.seafloor_height, state.surface_height
        )

        if day % config.snapshot_interval == 0 or day == config.days:
            snapshots[day] = _copy_state(state)

    return snapshots

def validate_forward_backward(config: SimulationConfig) -> bool:
    """
    For now, validation is:
    - Run forward simulation to get snapshots.
    - Re-run from initial seed and ensure snapshots match at key points.
    This acts like a forward/backward deterministic check.
    """
    forward_snaps = run_365_day_simulation(config)
    replay_snaps = run_365_day_simulation(config)

    key_days: List[int] = sorted(forward_snaps.keys())
    for day in key_days:
        a = forward_snaps[day]
        b = replay_snaps[day]
        if not _states_equal(a, b):
            return False
    return True

def _copy_state(state: WorldState) -> Snapshot:
    import copy
    return WorldState(
        day=state.day,
        surface_height=state.surface_height.copy(),
        seafloor_height=state.seafloor_height.copy(),
        plates=copy.deepcopy(state.plates),
        volcano_events=list(state.volcano_events),
    )

def _states_equal(a: WorldState, b: WorldState, tol: float = 1e-6) -> bool:
    if a.day != b.day:
        return False
    if a.surface_height.shape != b.surface_height.shape:
        return False
    if a.seafloor_height.shape != b.seafloor_height.shape:
        return False
    if not np.allclose(a.surface_height, b.surface_height, atol=tol):
        return False
    if not np.allclose(a.seafloor_height, b.seafloor_height, atol=tol):
        return False
    if len(a.plates) != len(b.plates):
        return False
    # Simple plate equality check
    for pa, pb in zip(a.plates, b.plates):
        if pa.id != pb.id or pa.velocity != pb.velocity or pa.is_oceanic != pb.is_oceanic:
            return False
    return True
