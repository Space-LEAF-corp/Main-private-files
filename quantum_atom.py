# quantum_atom.py
# Builds a toy "quantum atom" using hydrogen-like and oscillator models,
# with infinite light amplification embedded coherently.

from dataclasses import dataclass, field
from typing import List, Dict, Optional
from math import pi

from .infinite_light import (
    InfiniteLightAmplifier, SpectralScaler, H_BAR, M_E,
    hydrogen_energy_ev, energy_harmonic_level, omega_from_period,
    bohr_radius_Z, C
)

@dataclass
class QuantumAtom:
    """
    Toy quantum atom constructed from:
    - Hydrogen-like energy levels (Bohr model)
    - Harmonic oscillator shell (didactic vibrational mode)
    - Infinite light amplifier to scale spectral content safely

    Parameters:
    - Z: nuclear charge (default: 1 for hydrogen-like)
    - base_period_T: optional period to set an oscillator ω; if None, use ω_e (electron Compton analogue)
    - u: operational speed for Γ(c), must be < c
    - levels: number of bound levels to compute
    - scale_hydrogen: whether to scale hydrogen energies by Γ(c) (symbolic choice)
    """
    Z: int = 1
    base_period_T: Optional[float] = None
    u: float = 0.0
    levels: int = 5
    scale_hydrogen: bool = False

    amplifier: InfiniteLightAmplifier = field(init=False)
    scaler: SpectralScaler = field(init=False)

    def __post_init__(self):
        self.amplifier = InfiniteLightAmplifier(u=self.u, c=C)
        self.scaler = SpectralScaler(amp=self.amplifier)

    def oscillator_frequency(self) -> float:
        """
        Choose a didactic ω:
        - If base_period_T provided: ω = 2π / T
        - Else: use ω ≈ m_e c^2 / ħ (Compton-like angular frequency as a teaching anchor)
        """
        if self.base_period_T is not None:
            omega = omega_from_period(self.base_period_T)
        else:
            omega = (M_E * (C**2)) / H_BAR  # angular Compton frequency analogue
        return omega

    def build(self) -> Dict:
        """
        Construct the atom summary with:
        - Γ(c), oscillator ω and scaled ω_∞
        - Harmonic oscillator levels (J)
        - Hydrogen-like levels (eV), optionally Γ-scaled
        - Characteristic radius (Bohr), and scaled radius (symbolic)
        """
        gamma = self.amplifier.gamma()
        omega = self.oscillator_frequency()
        omega_inf = self.scaler.scale_frequency(omega)

        # Harmonic levels (Joules)
        ho_levels = [energy_harmonic_level(n, omega) for n in range(self.levels)]
        ho_levels_inf = [energy_harmonic_level(n, omega_inf) for n in range(self.levels)]

        # Hydrogen-like levels (eV)
        h_levels = [hydrogen_energy_ev(n) for n in range(1, self.levels + 1)]
        if self.scale_hydrogen:
            h_levels_inf = [self.scaler.scale_energy(E) for E in h_levels]
        else:
            h_levels_inf = h_levels[:]  # no scaling (physical hydrogen energies are not Γ-scaled; symbolic only)

        # Characteristic radius
        a0 = bohr_radius_Z(self.Z)
        # Optional symbolic radius scaling (nonphysical, ceremonial): a0 / Γ to represent "tightening under light"
        a0_inf = a0 / gamma

        return {
            "Z": self.Z,
            "gamma": gamma,
            "oscillator": {
                "omega": omega,
                "omega_infinite_light": omega_inf,
                "levels_J": ho_levels,
                "levels_J_infinite_light": ho_levels_inf,
            },
            "hydrogen_like": {
                "levels_eV": h_levels,
                "levels_eV_infinite_light": h_levels_inf,
                "scaled": self.scale_hydrogen,
            },
            "radius": {
                "bohr_m": a0,
                "bohr_m_infinite_light_symbolic": a0_inf,
            },
        }
