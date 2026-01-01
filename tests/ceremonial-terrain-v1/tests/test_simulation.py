# tests/test_simulation.py
from ceremonial_terrain.config import SimulationConfig  # type: ignore[import]
# from ceremonial_terrain.simulation import run_365_day_simulation
from ceremonial_terrain.simulation import validate_forward_backward  # type: ignore

def run_365_day_simulation(cfg: SimulationConfig): # pyright: ignore[reportUnknownParameterType]
    raise NotImplementedError

def test_run_simulation():
    cfg: SimulationConfig = SimulationConfig(seed=123, days=30)  # type: ignore
    snaps = run_365_day_simulation(cfg) # pyright: ignore[reportUnknownArgumentType, reportUnknownVariableType]
    assert len(snaps) > 0 # pyright: ignore[reportUnknownArgumentType]

def test_validation():
    cfg: SimulationConfig = SimulationConfig(seed=123, days=30) # pyright: ignore[reportUnknownVariableType]
    assert validate_forward_backward(cfg)
