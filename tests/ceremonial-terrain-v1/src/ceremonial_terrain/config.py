from dataclasses import dataclass

@dataclass
class SimulationConfig:
    seed: int = 42
    days: int = 365
    width: int = 128          # grid width
    height: int = 128         # grid height
    plate_count: int = 6
    snapshot_interval: int = 5  # store state every N days
