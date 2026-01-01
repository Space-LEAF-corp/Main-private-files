
import numpy as np
from typing import Dict, List

from .config import SimulationConfig
from .models import WorldState, Snapshot

from .terrain import generate_initial_heightmap, generate_seafloor_from_surface  # type: ignore
from .tectonics import initialize_plates, update_surface_for_tectonics
from .volcanoes import update_volcanic_activity # pyright: ignore[reportUnknownVariableType]
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from numpy import ndarray
    from typing import Callable
    update_seafloor_for_tectonics_and_volcanoes: Callable[["ndarray", "ndarray"], "ndarray"]
else:
    from .ocean import update_seafloor_for_tectonics_and_volcanoes

# Type hint for update_seafloor_for_tectonics_and_volcanoes
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from numpy import ndarray
    def update_seafloor_for_tectonics_and_volcanoes(seafloor: 'ndarray', surface: 'ndarray') -> 'ndarray': ... # pyright: ignore[reportUnknownParameterType]


def run_365_day_simulation(config: SimulationConfig) -> Dict[int, Snapshot]:
    rng = np.random.default_rng(config.seed)

    surface: np.ndarray = generate_initial_heightmap(config.width, config.height, rng) # pyright: ignore[reportUnknownArgumentType]
    seafloor: np.ndarray = generate_seafloor_from_surface(surface)  # type: ignore[reportUnknownVariableType]
    plates = initialize_plates(config.width, config.height, config.plate_count, rng)  # type: ignore[call-arg]

    snapshots: Dict[int, Snapshot] = {}
    state = WorldState(day=0, surface_height=surface, seafloor_height=seafloor, plates=plates)  # type: ignore[reportUnknownMemberType]

    snapshots[0] = _copy_state(state)

    for day in range(1, config.days + 1):
        state.day = day
        state.surface_height = update_surface_for_tectonics(
            state.surface_height, state.plates, day  # type: ignore[arg-type]
        )
        state = update_volcanic_activity(state, rng) # pyright: ignore[reportUnknownArgumentType]
        state.seafloor_height = update_seafloor_for_tectonics_and_volcanoes(
            state.seafloor_height, state.surface_height # pyright: ignore[reportUnknownArgumentType]
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

import copy
from typing import cast
from numpy import ndarray

def _copy_state(state: WorldState) -> Snapshot:
    # Ensure surface_height and seafloor_height are np.ndarray for type checkers
    surface_height: np.ndarray = np.array(state.surface_height, copy=True)
    seafloor_height: np.ndarray = np.array(state.seafloor_height, copy=True)
    # Handle volcano_events if not present
    volcano_events = list(getattr(state, 'volcano_events', []))
    return WorldState(
        day=state.day,
        surface_height=surface_height,
        seafloor_height=seafloor_height,
        plates=copy.deepcopy(state.plates),
        volcano_events=volcano_events,
    )

def _states_equal(a: WorldState, b: WorldState, tol: float = 1e-6) -> bool:
    if a.day != b.day:
        return False
    # Ensure surface_height and seafloor_height are np.ndarray for type checkers
    a_surface = cast(np.ndarray, a.surface_height)
    b_surface = cast(np.ndarray, b.surface_height)
    a_seafloor = cast(np.ndarray, a.seafloor_height)
    b_seafloor = cast(np.ndarray, b.seafloor_height)
    if a_surface.shape != b_surface.shape:  # type: ignore[attr-defined]
        return False
    if a_seafloor.shape != b_seafloor.shape:  # type: ignore[attr-defined]
        return False
    if not np.allclose(a_surface, b_surface, atol=tol):  # type: ignore[attr-defined]
        return False
    if not np.allclose(a_seafloor, b_seafloor, atol=tol):  # type: ignore[attr-defined]
        return False
    if len(a.plates) != len(b.plates):
        return False
    # Simple plate equality check
    for pa, pb in zip(a.plates, b.plates):
        if pa.id != pb.id or pa.velocity != pb.velocity or pa.is_oceanic != pb.is_oceanic:
            return False
    return True
