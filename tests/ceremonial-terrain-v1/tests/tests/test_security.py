# tests/test_simulation.py


from typing import Any
from ceremonial_terrain.config import SimulationConfig  # type: ignore
from ceremonial_terrain.simulation import validate_forward_backward  # type: ignore
# from ceremonial_terrain.simulation import run_365_day_simulation

# Explicit type annotation for run_365_day_simulation to resolve unknown type error
from typing import Callable, Any, List
# from ceremonial_terrain.config import SimulationConfig
run_365_day_simulation: Callable[[SimulationConfig], List[Any]]

# Type hint for run_365_day_simulation to resolve unknown type error
from typing import Callable, Any, List
run_365_day_simulation: Callable[[SimulationConfig], List[Any]]

def test_run_simulation():
    cfg: SimulationConfig = SimulationConfig(seed=123, days=30) # pyright: ignore[reportUnknownVariableType]
    snaps: list[Any] = run_365_day_simulation(cfg) # pyright: ignore[reportUnknownArgumentType]
    assert len(snaps) > 0

def test_validation():
    cfg: SimulationConfig = SimulationConfig(seed=123, days=30) # pyright: ignore[reportUnknownVariableType]
    assert validate_forward_backward(cfg)
