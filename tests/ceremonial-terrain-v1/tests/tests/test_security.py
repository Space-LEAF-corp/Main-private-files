# tests/test_simulation.py
from ceremonial_terrain.config import SimulationConfig
from ceremonial_terrain.simulation import run_365_day_simulation, validate_forward_backward

def test_run_simulation():
    cfg = SimulationConfig(seed=123, days=30)
    snaps = run_365_day_simulation(cfg)
    assert len(snaps) > 0

def test_validation():
    cfg = SimulationConfig(seed=123, days=30)
    assert validate_forward_backward(cfg)
