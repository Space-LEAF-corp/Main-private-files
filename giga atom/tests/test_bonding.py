# tests/test_bonding.py
# type: ignore
from giga_atom.simulations import valence_from_Z, bonding_score

def test_valence_small():
    assert valence_from_Z(1) == 1
    assert valence_from_Z(2) == 2

def test_bonding_score_symmetry():
    assert bonding_score(1, 1) == 1.0
    assert 0.0 <= bonding_score(1, 2) <= 1.0
