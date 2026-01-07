"""
hawkeye_timeline_prototype.py

Single-file Hawkeye Timeline Prototype:
- SensorSample model
- Fibonacci decoding
- Mandala mapping
- Star anchors + Vulture Stone
- DNA-style reverse/complement
- Orchestrator + simple demo
"""

from dataclasses import dataclass
from typing import List, Dict, Tuple
import numpy as np


# ---------- Core data models ----------

@dataclass
class SensorSample:
    time: float                  # seconds
    angle: float                 # radians
    wavelength_band: str         # e.g. "VIS", "IR"
    intensity: float             # signal strength
    meta: Dict[str, float]       # extra info if needed


@dataclass
class FibonacciBlock:
    block_id: int
    fib_size: int
    samples: List[SensorSample]
    mean_angle: float
    spectrum: Dict[str, float]
    stability_score: float


@dataclass
class MandalaCell:
    ring_index: int
    sector_index: int
    samples: List[FibonacciBlock]
    aggregated_intensity: float
    coherence: float


@dataclass
class StarAnchor:
    name: str
    central_angle: float   # radians
    angular_width: float   # radians
    epoch_tag: str


@dataclass
class VultureStoneAnchor:
    name: str
    associated_constellations: List[str]
    epoch_start_bce: int
    epoch_end_bce: int
    warning_weight: float


@dataclass
class HawkeyeConfig:
    num_rings: int
    num_sectors: int
    fibonacci_cycle_len: int
    mandala_radius_scale: float
    star_anchors: List[StarAnchor]
    vulture_anchor: VultureStoneAnchor


# ---------- Fibonacci decoding ----------

def fibonacci_sequence(n: int) -> List[int]:
    fib = [1, 1]
    while len(fib) < n:
        fib.append(fib[-1] + fib[-2])
    return fib[:n]


def fibonacci_remap_indices(num_samples: int, cycle_len: int) -> List[int]:
    fib = fibonacci_sequence(cycle_len + 1)
    remap = []
    for n in range(num_samples):
        f_n = fib[n % cycle_len]
        j_n = n + f_n
        remap.append(j_n)
    return remap


def build_fibonacci_blocks(
    samples: List[SensorSample],
    cycle_len: int
) -> List[FibonacciBlock]:
    if not samples:
        return []

    num_samples = len(samples)
    remapped = fibonacci_remap_indices(num_samples, cycle_len)
    indexed = list(zip(remapped, samples))
    indexed.sort(key=lambda x: x[0])

    fib = fibonacci_sequence(20)
    blocks: List[FibonacciBlock] = []
    cursor = 0
    block_id = 0

    for f_size in fib:
        if cursor + f_size > num_samples:
            break
        block_samples = [s for (_, s) in indexed[cursor:cursor + f_size]]
        mean_angle = float(np.mean([s.angle for s in block_samples]))

        spectrum: Dict[str, List[float]] = {}
        for s in block_samples:
            spectrum.setdefault(s.wavelength_band, []).append(s.intensity)
        spectrum_mean = {band: float(np.mean(vals)) for band, vals in spectrum.items()}

        intensities = [s.intensity for s in block_samples]
        var = np.var(intensities) if len(intensities) > 1 else 0.0
        stability = float(1.0 / (1.0 + var))

        blocks.append(
            FibonacciBlock(
                block_id=block_id,
                fib_size=f_size,
                samples=block_samples,
                mean_angle=mean_angle,
                spectrum=spectrum_mean,
                stability_score=stability,
            )
        )
        cursor += f_size
        block_id += 1

    return blocks


# ---------- Mandala mapping ----------

def angle_to_sector(angle: float, num_sectors: int) -> int:
    angle = (angle + 2 * np.pi) % (2 * np.pi)
    sector_width = 2 * np.pi / num_sectors
    return int(angle // sector_width)


def build_mandala(
    blocks: List[FibonacciBlock],
    config: HawkeyeConfig
) -> List[MandalaCell]:
    cells: Dict[Tuple[int, int], MandalaCell] = {}

    for block in blocks:
        ring_index = min(block.fib_size, config.num_rings - 1)
        sector_index = angle_to_sector(block.mean_angle, config.num_sectors)
        key = (ring_index, sector_index)

        if key not in cells:
            cells[key] = MandalaCell(
                ring_index=ring_index,
                sector_index=sector_index,
                samples=[],
                aggregated_intensity=0.0,
                coherence=0.0,
            )
        cell = cells[key]
        cell.samples.append(block)
        cell.aggregated_intensity += sum(block.spectrum.values())

    if cells:
        max_intensity = max(c.aggregated_intensity for c in cells.values()) or 1.0
    else:
        max_intensity = 1.0

    for cell in cells.values():
        cell.coherence = cell.aggregated_intensity / max_intensity

    return list(cells.values())


# ---------- Star anchors + Vulture Stone ----------

def tag_cells_with_star_anchors(
    cells: List[MandalaCell],
    config: HawkeyeConfig
) -> Dict[Tuple[int, int], List[str]]:
    sector_angle_width = 2 * np.pi / config.num_sectors
    cell_tags: Dict[Tuple[int, int], List[str]] = {}

    for cell in cells:
        sector_center = (cell.sector_index + 0.5) * sector_angle_width - np.pi
        tags: List[str] = []
        for anchor in config.star_anchors:
            delta = abs(((sector_center - anchor.central_angle + np.pi) % (2 * np.pi)) - np.pi)
            if delta <= anchor.angular_width / 2:
                tags.append(anchor.name)
        cell_tags[(cell.ring_index, cell.sector_index)] = tags

    return cell_tags


def apply_vulture_warning_weight(
    cells: List[MandalaCell],
    config: HawkeyeConfig,
    cell_tags: Dict[Tuple[int, int], List[str]]
) -> None:
    vulture = config.vulture_anchor
    warning_constellations = set(vulture.associated_constellations)

    for cell in cells:
        tags = cell_tags.get((cell.ring_index, cell.sector_index), [])
        if any(t in warning_constellations for t in tags):
            cell.coherence *= (1.0 + vulture.warning_weight)


# ---------- DNA-style reverse/complement ----------

def build_cell_lookup(cells: List[MandalaCell]) -> Dict[Tuple[int, int], MandalaCell]:
    return {(c.ring_index, c.sector_index): c for c in cells}


def dna_reverse_complement(
    cells: List[MandalaCell],
    config: HawkeyeConfig
) -> List[MandalaCell]:
    complement_cells: List[MandalaCell] = []

    for cell in cells:
        comp_ring = (config.num_rings - 1) - cell.ring_index
        comp_sector = (cell.sector_index + config.num_sectors // 2) % config.num_sectors
        comp_coherence = 1.0 - cell.coherence

        complement_cells.append(
            MandalaCell(
                ring_index=comp_ring,
                sector_index=comp_sector,
                samples=cell.samples,
                aggregated_intensity=cell.aggregated_intensity,
                coherence=comp_coherence,
            )
        )

    return complement_cells


def overlay_original_and_complement(
    original: List[MandalaCell],
    complement: List[MandalaCell]
) -> Dict[Tuple[int, int], Dict[str, object]]:
    orig_lookup = build_cell_lookup(original)
    comp_lookup = build_cell_lookup(complement)
    result: Dict[Tuple[int, int], Dict[str, object]] = {}

    for key, orig_cell in orig_lookup.items():
        comp_cell = comp_lookup.get(key)
        if comp_cell is None:
            continue
        agreement = 1.0 - abs(orig_cell.coherence - comp_cell.coherence)
        result[key] = {
            "orig": orig_cell,
            "comp": comp_cell,
            "agreement": agreement,
        }

    return result


# ---------- Hawkeye Timeline orchestrator ----------

class HawkeyeTimeline:
    def __init__(self, config: HawkeyeConfig):
        self.config = config

    def process_samples(self, samples: List[SensorSample]) -> Dict[str, object]:
        blocks = build_fibonacci_blocks(samples, self.config.fibonacci_cycle_len)
        mandala_cells = build_mandala(blocks, self.config)
        cell_tags = tag_cells_with_star_anchors(mandala_cells, self.config)
        apply_vulture_warning_weight(mandala_cells, self.config, cell_tags)
        complement_cells = dna_reverse_complement(mandala_cells, self.config)
        overlay = overlay_original_and_complement(mandala_cells, complement_cells)

        return {
            "blocks": blocks,
            "mandala_cells": mandala_cells,
            "cell_tags": cell_tags,
            "complement_cells": complement_cells,
            "overlay": overlay,
        }


# ---------- Simple demo ----------

def example_run():
    # Define some rough star anchors (angles are just placeholders)
    star_anchors = [
        StarAnchor(
            name="Scorpius",
            central_angle=np.deg2rad(240),
            angular_width=np.deg2rad(30),
            epoch_tag="archaic",
        ),
        StarAnchor(
            name="Cygnus",
            central_angle=np.deg2rad(300),
            angular_width=np.deg2rad(30),
            epoch_tag="archaic",
        ),
    ]

    vulture_anchor = VultureStoneAnchor(
        name="Vulture Stone",
        associated_constellations=["Scorpius", "Cygnus"],
        epoch_start_bce=9600,
        epoch_end_bce=8000,
        warning_weight=0.5,
    )

    config = HawkeyeConfig(
        num_rings=16,
        num_sectors=32,
        fibonacci_cycle_len=8,
        mandala_radius_scale=1.0,
        star_anchors=star_anchors,
        vulture_anchor=vulture_anchor,
    )

    # Fake sensor data
    samples: List[SensorSample] = []
    for n in range(500):
        angle = float(np.random.uniform(-np.pi, np.pi))
        intensity = float(np.random.random())
        band = np.random.choice(["VIS", "IR"])
        samples.append(
            SensorSample(
                time=n * 0.1,
                angle=angle,
                wavelength_band=band,
                intensity=intensity,
                meta={"ring_id": 0},
            )
        )

    engine = HawkeyeTimeline(config)
    results = engine.process_samples(samples)

    print("Blocks:", len(results["blocks"]))
    print("Mandala cells:", len(results["mandala_cells"]))
    print("Overlay cells:", len(results["overlay"]))

    # Example: show top 5 highest-agreement cells
    overlay = results["overlay"]
    sorted_cells = sorted(
        overlay.items(),
        key=lambda kv: kv[1]["agreement"],
        reverse=True
    )[:5]

    print("\nTop 5 high-agreement cells (ring, sector, agreement):")
    for (ring, sector), data in sorted_cells:
        print(ring, sector, f"{data['agreement']:.3f}")


if __name__ == "__main__":
    example_run()
