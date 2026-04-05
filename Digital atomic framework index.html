from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
from enum import Enum, IntEnum
import uuid


class AtomicState(IntEnum):
    """
    Multi-level state (1–5) to encode
    refinement / role / interaction level.
    You can rename these labels for kids.
    """
    RAW_OBSERVED      = 1  # just seen
    IDENTIFIED        = 2  # matched to element
    CONTEXTUALIZED    = 3  # environment known (bonding/neighbors)
    FUNCTION_ASSIGNED = 4  # role in a structure (e.g. 'support', 'conductor')
    INTEGRATED_MODEL  = 5  # part of a validated Krystal Core pattern


@dataclass
class AtomicIdentity:
    """
    Digital identity for a single atomic site.
    """
    id: str                      # unique digital ID
    element_symbol: str          # e.g. "C", "Si", "O"
    atomic_number: int           # periodic table link
    position: Tuple[float, float, float]  # (x, y, z) in nanometers
    state: AtomicState           # 1–5 refinement / role state
    role_label: str = ""         # human/kid-friendly role (e.g. "bridge", "shield", "memory-node")
    metadata: Dict[str, str] = field(default_factory=dict)  # extra tags (sample, chamber, etc.)


@dataclass
class DigitalAtomicFrameworkIndex:
    """
    Core container: a map from digital IDs to atomic identities.
    This is your 'digital atomic framework index'.
    """
    atoms: Dict[str, AtomicIdentity] = field(default_factory=dict)

    def add_atom(
        self,
        element_symbol: str,
        atomic_number: int,
        position: Tuple[float, float, float],
        state: AtomicState = AtomicState.RAW_OBSERVED,
        role_label: str = "",
        metadata: Optional[Dict[str, str]] = None,
    ) -> str:
        """
        Add a new atomic identity to the index and return its ID.
        """
        atom_id = str(uuid.uuid4())
        atom = AtomicIdentity(
            id=atom_id,
            element_symbol=element_symbol,
            atomic_number=atomic_number,
            position=position,
            state=state,
            role_label=role_label,
            metadata=metadata or {},
        )
        self.atoms[atom_id] = atom
        return atom_id

    def update_state(self, atom_id: str, new_state: AtomicState, new_role_label: Optional[str] = None):
        """
        Refine an atom's state (e.g. from 1 → 3 → 5) as the microscope
        and analysis pipeline learn more.
        """
        if atom_id not in self.atoms:
            raise KeyError(f"No atom with id {atom_id}")
        self.atoms[atom_id].state = new_state
        if new_role_label is not None:
            self.atoms[atom_id].role_label = new_role_label

    def find_by_element(self, element_symbol: str) -> List[AtomicIdentity]:
        """
        Get all atoms of a given element (e.g. all 'Si' atoms).
        """
        return [
            atom for atom in self.atoms.values()
            if atom.element_symbol == element_symbol
        ]

    def find_by_state(self, state: AtomicState) -> List[AtomicIdentity]:
        """
        Get all atoms at a given refinement/role state.
        """
        return [
            atom for atom in self.atoms.values()
            if atom.state == state
        ]

    def export_for_kids_view(self) -> List[Dict[str, str]]:
        """
        Flatten the index into a kid-friendly view:
        only show symbol, role, and a simple position label.
        """
        simplified = []
        for atom in self.atoms.values():
            x, y, z = atom.position
            simplified.append({
                "element": atom.element_symbol,
                "role": atom.role_label or "unassigned",
                "state_level": str(int(atom.state)),
                "position_label": f"({x:.1f}, {y:.1f}, {z:.1f}) nm",
            })
        return simplified


# --- Example usage ---

if __name__ == "__main__":
    index = DigitalAtomicFrameworkIndex()

    # Add a few atoms as if they came from the tunneling microscope map
    c_id = index.add_atom(
        element_symbol="C",
        atomic_number=6,
        position=(0.0, 0.0, 0.0),
        state=AtomicState.IDENTIFIED,
        role_label="core-bridge",
        metadata={"chamber": "memory", "sample": "KrystalCorePrototype01"},
    )

    si_id = index.add_atom(
        element_symbol="Si",
        atomic_number=14,
        position=(0.2, 0.1, 0.0),
        state=AtomicState.CONTEXTUALIZED,
        role_label="signal-path",
        metadata={"chamber": "interface", "sample": "KrystalCorePrototype01"},
    )

    # Refine one atom as the model improves
    index.update_state(si_id, AtomicState.INTEGRATED_MODEL, new_role_label="stabilizer-node")

    # Get all silicon atoms
    silicon_atoms = index.find_by_element("Si")

    # Export a kid-friendly view
    kids_view = index.export_for_kids_view()
    for entry in kids_view:
        print(entry)