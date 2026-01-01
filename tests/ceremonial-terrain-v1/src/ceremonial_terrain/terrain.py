import numpy as np

def generate_initial_heightmap(width: int, height: int, rng: np.random.Generator) -> np.ndarray:
    """
    Generate a simple procedural heightmap.
    This is a placeholder for:
    - noise-based terrain
    - or image-based tea/coffee residue maps
    """
    base = rng.normal(loc=0.0, scale=1.0, size=(height, width))
    # Smooth with a simple blur kernel
    kernel = np.array([[1, 2, 1],
                       [2, 4, 2],
                       [1, 2, 1]], dtype=float)
    kernel /= kernel.sum()
    for _ in range(3):
        base = _convolve2d(base, kernel) # pyright: ignore[reportUnknownArgumentType]
    # Normalize to [0, 1]
    base -= base.min()
    base /= np.maximum(base.max(), 1e-9)
    return base

def generate_seafloor_from_surface(surface: np.ndarray, sea_level: float = 0.4) -> np.ndarray:
    """
    Simple bathymetry: below sea_level, define depth; above sea_level, keep shallow.
    """
    seafloor = surface.copy()
    ocean_mask = surface < sea_level
    # Make ocean deeper where surface is low
    seafloor[ocean_mask] = -(sea_level - surface[ocean_mask])
    seafloor[~ocean_mask] = 0.0
    return seafloor

def _convolve2d(arr: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    h, w = arr.shape # pyright: ignore[reportUnknownVariableType]
    kh, kw = kernel.shape
    pad_h, pad_w = kh // 2, kw // 2
    padded = np.pad(arr, ((pad_h, pad_h), (pad_w, pad_w)), mode="edge")
    out = np.zeros_like(arr)
    for i in range(h):
        for j in range(w):
            region = padded[i:i+kh, j:j+kw]
            out[i, j] = np.sum(region * kernel)
    return out
