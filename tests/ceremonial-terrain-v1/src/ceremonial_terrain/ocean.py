try:
    import numpy as np  # type: ignore
except ImportError:
    raise ImportError("numpy is required for ceremonial_terrain.ocean module. Please install numpy.")

def update_seafloor_for_tectonics_and_volcanoes(seafloor: np.ndarray, surface: np.ndarray) -> np.ndarray: # pyright: ignore[reportUnknownMemberType, reportUnknownParameterType]
    """
    For now, we approximate the seafloor as:
    - Some fraction of the surface tectonic deformation
    - Deeper where surface is low
    """
    # type: (np.ndarray, np.ndarray) -> np.ndarray
    assert isinstance(seafloor, np.ndarray), "seafloor must be a numpy ndarray" # pyright: ignore[reportUnknownMemberType]
    sea_level = 0.4
    ocean_mask: np.ndarray = surface < sea_level # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    new_seafloor: np.ndarray = seafloor.copy()  # type: ignore[attr-defined]
    new_seafloor[ocean_mask] = -(sea_level - surface[ocean_mask]) * 1.5
    new_seafloor[~ocean_mask] = 0.0
    return new_seafloor  # pyright: ignore[reportUnknownVariableType] # type: np.ndarray
