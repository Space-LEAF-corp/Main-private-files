import numpy as np  # type: ignore

def generate_initial_heightmap(width: int, height: int, rng: np.random.Generator) -> np.ndarray:  # pyright: ignore[reportUnknownParameterType, reportUnknownMemberType]
    """
    Generate a simple procedural heightmap.
    This is a placeholder for:
    - noise-based terrain
    - or image-based tea/coffee residue maps
    """
    base: np.ndarray = rng.normal(loc=0.0, scale=1.0, size=(height, width)) # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    # Smooth with a simple blur kernel
    kernel: np.ndarray = np.array([[1, 2, 1], # pyright: ignore[reportUnknownVariableType] # pyright: ignore[reportUnknownMemberType] # pyright: ignore[reportUnknownMemberType] # type: ignore
                                   [2, 4, 2],
                                   [1, 2, 1]], dtype=float)
    kernel /= kernel.sum() # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    for _ in range(3):
        base = _convolve2d(base, kernel) # pyright: ignore[reportUnknownVariableType, reportUnknownArgumentType]
    # Normalize to [0, 1]
    base -= base.min() # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    base /= np.maximum(base.max(), 1e-9) # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    return base # pyright: ignore[reportUnknownVariableType]

def generate_seafloor_from_surface(surface: np.ndarray, sea_level: float = 0.4) -> np.ndarray: # pyright: ignore[reportUnknownParameterType, reportUnknownMemberType]
    """
    Simple bathymetry: below sea_level, define depth; above sea_level, keep shallow.
    """

    seafloor = surface.copy() # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    # Make ocean deeper where surface is low
    ocean_mask = surface < sea_level # pyright: ignore[reportUnknownVariableType]
    seafloor[ocean_mask] = -(sea_level - surface[ocean_mask]) # pyright: ignore[reportUndefinedVariable]
    seafloor[~ocean_mask] = 0.0 # pyright: ignore[reportUndefinedVariable]
    return seafloor # pyright: ignore[reportUnknownVariableType]

def _convolve2d(arr: np.ndarray, kernel: np.ndarray) -> np.ndarray: # pyright: ignore[reportUnknownMemberType, reportUnknownParameterType]
    h, w = arr.shape # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    kh, kw = kernel.shape # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    pad_h, pad_w = kh // 2, kw // 2 # pyright: ignore[reportUnknownVariableType]
    padded = np.pad(arr, ((pad_h, pad_h), (pad_w, pad_w)), mode="edge") # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    out = np.zeros_like(arr) # pyright: ignore[reportUnknownVariableType, reportUnknownMemberType]
    for i in range(h): # pyright: ignore[reportUnknownArgumentType]
        for j in range(w): # pyright: ignore[reportUnknownArgumentType]
            region = padded[i:i+kh, j:j+kw] # pyright: ignore[reportUnknownVariableType]
            out[i, j] = np.sum(region * kernel) # pyright: ignore[reportUnknownMemberType]
    return out # pyright: ignore[reportUnknownVariableType]
