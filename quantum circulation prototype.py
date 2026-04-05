from dataclasses import dataclass, asdict
from typing import List, Dict, Any, Tuple
import hashlib
import random

# -----------------------------
# 1. Core data structures
# -----------------------------

@dataclass
class AtomicNode:
    element: str          # e.g. "C", "Si"
    role: str             # e.g. "core-bridge", "stabilizer-node"
    state_level: int      # e.g. 1–9
    position_label: Tuple[float, float, float]  # (x, y, z) in nm

    def to_dict(self) -> Dict[str, Any]:
        return {
            "element": self.element,
            "role": self.role,
            "state_level": str(self.state_level),
            "position_label": f"({self.position_label[0]}, "
                              f"{self.position_label[1]}, "
                              f"{self.position_label[2]}) nm"
        }


@dataclass
class AtomicFingerprint:
    satellite_id: str
    nodes: List[AtomicNode]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "satellite_id": self.satellite_id,
            "nodes": [n.to_dict() for n in self.nodes]
        }


# -----------------------------
# 2. Deterministic RNG helpers
# -----------------------------

def make_seed(*parts: str) -> int:
    """Create a deterministic integer seed from multiple string parts."""
    joined = "|".join(parts)
    h = hashlib.sha256(joined.encode("utf-8")).hexdigest()
    # Use a slice of the hash as an integer seed
    return int(h[:16], 16)


def make_rng(*parts: str) -> random.Random:
    return random.Random(make_seed(*parts))


# -----------------------------
# 3. Configuration spaces
# -----------------------------

ELEMENTS = ["C", "Si", "Fe", "O", "N"]
ROLES = [
    "core-bridge",
    "stabilizer-node",
    "sensor-node",
    "shield-node",
    "tunnel-gate"
]


def pick_element(rng: random.Random) -> str:
    return rng.choice(ELEMENTS)


def pick_role(rng: random.Random) -> str:
    return rng.choice(ROLES)


def generate_position_pattern(rng: random.Random) -> Tuple[float, float, float]:
    """
    Very simple placeholder:
    small offsets in nm around origin.
    """
    x = round(rng.uniform(-0.5, 0.5), 3)
    y = round(rng.uniform(-0.5, 0.5), 3)
    z = round(rng.uniform(-0.5, 0.5), 3)
    return (x, y, z)


# -----------------------------
# 4. Fingerprint generator
# -----------------------------

def generate_atomic_fingerprint(
    satellite_id: str,
    launch_ts: str,
    orbit_class: str,
    mission_type: str,
    num_nodes: int = 2
) -> AtomicFingerprint:
    """
    Deterministically generate an atomic fingerprint for a satellite.
    """
    rng = make_rng(satellite_id, launch_ts, orbit_class, mission_type)

    nodes: List[AtomicNode] = []
    for i in range(num_nodes):
        element = pick_element(rng)
        role = pick_role(rng)
        state_level = rng.randint(1, 9)
        position = generate_position_pattern(rng)

        node = AtomicNode(
            element=element,
            role=role,
            state_level=state_level,
            position_label=position
        )
        nodes.append(node)

    return AtomicFingerprint(satellite_id=satellite_id, nodes=nodes)


# -----------------------------
# 5. Example: tiny atomic circuit
# -----------------------------

def example_circuit() -> AtomicFingerprint:
    """
    Hard-coded example similar to your console output:
    C as core-bridge, Si as stabilizer-node.
    """
    nodes = [
        AtomicNode(
            element="C",
            role="core-bridge",
            state_level=2,
            position_label=(0.0, 0.0, 0.0)
        ),
        AtomicNode(
            element="Si",
            role="stabilizer-node",
            state_level=5,
            position_label=(0.2, 0.1, 0.0)
        )
    ]
    return AtomicFingerprint(satellite_id="SAT-EXAMPLE-001", nodes=nodes)


# -----------------------------
# 6. Demo usage
# -----------------------------

if __name__ == "__main__":
    # Deterministic fingerprint from metadata
    fp = generate_atomic_fingerprint(
        satellite_id="SAT-001",
        launch_ts="2026-02-02T05:00:00Z",
        orbit_class="LEO",
        mission_type="imaging",
        num_nodes=4
    )
    print("Generated fingerprint:")
    for node in fp.nodes:
        print(node.to_dict())

    # Hard-coded quantum-ish circuit example
    print("\nExample circuit:")
    circuit_fp = example_circuit()
    for node in circuit_fp.nodes:
        print(node.to_dict())