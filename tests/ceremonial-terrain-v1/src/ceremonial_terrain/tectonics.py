



from typing import List
try:
    import numpy as np # pyright: ignore[reportMissingImports]
except ImportError:
    np = None
import math

from typing import Any

# Always use Any for NDArray to avoid import issues
NDArray = Any

    # Removed empty try block that caused syntax error
from .models import Plate
from typing import Any


def initialize_plates(width: int, height: int, plate_count: int, rng: Any) -> List[Plate]:
    plates: List[Plate] = []
    for pid in range(plate_count):
        angle = float(rng.uniform(0, 2 * np.pi)) # pyright: ignore[reportUnknownMemberType, reportOptionalMemberAccess, reportUnknownArgumentType]
        speed: float = float(rng.uniform(0.001, 0.01)) # pyright: ignore[reportUnknownArgumentType]
        if np is not None:
            vx: float = float(speed * float(np.cos(float(angle)))) # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
            vy: float = float(speed * float(np.sin(float(angle)))) # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
        else:
            vx: float = float(speed * math.cos(angle))
            vy: float = float(speed * math.sin(angle))
        is_oceanic = bool(rng.integers(0, 2)) # pyright: ignore[reportUnknownArgumentType]
        plates.append(Plate(id=pid, velocity=(vx, vy), is_oceanic=is_oceanic))
    return plates

def update_surface_for_tectonics(surface: "NDArray | Any", plates: List[Plate], day: int) -> "NDArray | Any": # pyright: ignore[reportInvalidTypeForm, reportUnknownParameterType]
    """
    Placeholder tectonic deformation:
    - Slight uplift or subsidence based on plate ids and pseudo boundaries.
    This is where a real plate-boundary computation would go.
    """
    h, w = surface.shape # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    deformed: Any = surface.copy() # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    # Simple pattern: uplift along a diagonal "boundary"
    if np is not None:
        boundary_idx = (np.arange(h)[:, None] + np.arange(w)[None, :]) % 20 # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    else:
        # Fallback using Python lists
        boundary_idx = [[(i + j) % 20 for j in range(w)] for i in range(h)] # pyright: ignore[reportUnknownArgumentType]
    if np is not None:
        uplift_mask = boundary_idx < 2 # pyright: ignore[reportOperatorIssue, reportUnknownVariableType]
        subsidence_mask = boundary_idx > 17 # pyright: ignore[reportOperatorIssue, reportUnknownVariableType]
    else:
        # For Python lists, create boolean masks manually
        from typing import cast
        boundary_idx_list: List[List[int]] = cast(List[List[int]], boundary_idx)
        uplift_mask: List[List[bool]] = [[cell < 2 for cell in row] for row in boundary_idx_list]
        subsidence_mask: List[List[bool]] = [[cell > 17 for cell in row] for row in boundary_idx_list]

    uplift_amount = 0.0005
    subsidence_amount = 0.0003

    if np is not None:
        deformed[uplift_mask] += uplift_amount
        deformed[subsidence_mask] -= subsidence_amount
        # Clamp back to [0, 1]
        deformed -= np.min(deformed) # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
        maxv = max(deformed.max(), 1e-9) # pyright: ignore[reportUnknownMemberType, reportUnknownArgumentType]
    else:
        # fallback for non-numpy arrays
        # Explicitly type deformed as List[List[float]] for type safety
        from typing import cast
        deformed_list: List[List[float]] = cast(List[List[float]], deformed)
        # Apply uplift and subsidence using the boolean masks
        for row, uplift_row, subs_row in zip(deformed_list, uplift_mask, subsidence_mask): # pyright: ignore[reportUnknownArgumentType]
            for j, (uplift, subs) in enumerate(zip(uplift_row, subs_row)): # pyright: ignore[reportUnknownVariableType, reportUnknownArgumentType]
                if uplift:
                    row[j] += uplift_amount
                if subs:
                    row[j] -= subsidence_amount
        minv = min([min(row) for row in deformed_list])
        deformed_list = [[v - minv for v in row] for row in deformed_list]
        maxv = max([max(row) for row in deformed_list] + [1e-9])
        deformed_list = [[v / maxv for v in row] for row in deformed_list]
        deformed = deformed_list
    return deformed # pyright: ignore[reportUnknownVariableType]
