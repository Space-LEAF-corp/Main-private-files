import numpy as np # pyright: ignore[reportMissingImports]
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np # pyright: ignore[reportMissingImports]
from dataclasses import dataclass, field
from typing import List, Tuple

if TYPE_CHECKING:
    import numpy as np # pyright: ignore[reportMissingImports]

@dataclass
class Plate:
    id: int
    velocity: Tuple[float, float]  # (vx, vy)
    is_oceanic: bool

@dataclass
class VolcanoEvent:
    day: int
    position: Tuple[int, int]
    strength: float


@dataclass
class WorldState:
    day: int
    # Elevations: surface and seafloor
    surface_height: 'np.ndarray'  # pyright: ignore[reportUnknownMemberType] # shape: (H, W)
    seafloor_height: 'np.ndarray'  # pyright: ignore[reportUnknownMemberType] # shape: (H, W)

    # Tectonic plates
    plates: List[Plate]

    # Volcanic history
    volcano_events: 'List[VolcanoEvent]' = field(default_factory=list) # pyright: ignore[reportUnknownVariableType]

    def __post_init__(self):
        global np
        try:
            import numpy as np # pyright: ignore[reportMissingImports]
        except ImportError:
            raise ImportError("numpy is required for WorldState but could not be imported.")

Snapshot = WorldState  # alias for clearer semantics
