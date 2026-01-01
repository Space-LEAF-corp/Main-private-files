import numpy as np

def update_seafloor_for_tectonics_and_volcanoes(seafloor: np.ndarray, surface: np.ndarray) -> np.ndarray:
    """
    For now, we approximate the seafloor as:
    - Some fraction of the surface tectonic deformation
    - Deeper where surface is low
    """
    sea_level = 0.4
    ocean_mask = surface < sea_level
    new_seafloor = seafloor.copy()
    new_seafloor[ocean_mask] = -(sea_level - surface[ocean_mask]) * 1.5
    new_seafloor[~ocean_mask] = 0.0
    return new_seafloor
