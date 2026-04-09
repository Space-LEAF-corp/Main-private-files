# quantum_chamber_sim.py
#
# Pure software simulation of:
# - Excitation chamber
# - Energy buffer (collector)
# - Containment shell
# - Guardian controller
#
# No hardware control. No real fusion. Just logic.

import random
from enum import Enum, auto


class ChamberState(Enum):
    IDLE = auto()
    EXCITING = auto()
    DAMPING = auto()
    LOCKDOWN = auto()


class ContainmentState(Enum):
    NOMINAL = auto()
    STRESSED = auto()
    CRITICAL = auto()
    FAILED = auto()


class ExcitationChamber:
    def __init__(self):
        self.state = ChamberState.IDLE
        self.excitation_level = 0.0  # arbitrary units
        self.excitation_rate = 0.0   # change per tick

    def excite(self, rate: float):
        if self.state not in {ChamberState.LOCKDOWN}:
            self.state = ChamberState.EXCITING
            self.excitation_rate = rate

    def damp(self, rate: float):
        if self.state not in {ChamberState.LOCKDOWN}:
            self.state = ChamberState.DAMPING
            self.excitation_rate = -abs(rate)

    def tick(self):
        # update excitation level
        self.excitation_level += self.excitation_rate

        # clamp to non-negative
        if self.excitation_level < 0:
            self.excitation_level = 0

        # if nothing happening, go idle
        if abs(self.excitation_rate) < 1e-6 and self.excitation_level == 0:
            self.state = ChamberState.IDLE

    def lockdown(self):
        self.state = ChamberState.LOCKDOWN
        self.excitation_rate = 0.0

    def snapshot(self):
        return {
            "state": self.state.name,
            "excitation_level": self.excitation_level,
            "excitation_rate": self.excitation_rate,
        }


class EnergyBuffer:
    def __init__(self, capacity: float):
        self.capacity = capacity
        self.stored = 0.0

    def absorb(self, amount: float):
        # absorb up to capacity
        absorbed = min(amount, self.capacity - self.stored)
        self.stored += absorbed
        return absorbed

    def bleed_off(self, amount: float):
        # simulate redistribution / safe discharge
        released = min(amount, self.stored)
        self.stored -= released
        return released

    def snapshot(self):
        return {
            "stored": self.stored,
            "capacity": self.capacity,
            "fill_ratio": self.stored / self.capacity if self.capacity > 0 else 0.0,
        }


class ContainmentShell:
    def __init__(self, stress_threshold: float, critical_threshold: float):
        self.state = ContainmentState.NOMINAL
        self.stress = 0.0
        self.stress_threshold = stress_threshold
        self.critical_threshold = critical_threshold

    def apply_stress(self, amount: float):
        self.stress += amount
        if self.stress >= self.critical_threshold:
            self.state = ContainmentState.CRITICAL
        elif self.stress >= self.stress_threshold:
            self.state = ContainmentState.STRESSED
        else:
            self.state = ContainmentState.NOMINAL

    def relieve_stress(self, factor: float = 0.1):
        # passive relaxation
        self.stress *= (1 - factor)
        if self.stress < 0:
            self.stress = 0
        if self.stress == 0:
            self.state = ContainmentState.NOMINAL

    def snapshot(self):
        return {
            "state": self.state.name,
            "stress": self.stress,
            "stress_threshold": self.stress_threshold,
            "critical_threshold": self.critical_threshold,
        }


class GuardianController:
    """
    Guardian logic:
    - Monitors chamber, buffer, containment
    - Decides when to:
      * increase excitation
      * shunt energy to buffer
      * bleed off buffer
      * trigger lockdown
    """

    def __init__(
        self,
        chamber: ExcitationChamber,
        buffer: EnergyBuffer,
        shell: ContainmentShell,
    ):
        self.chamber = chamber
        self.buffer = buffer
        self.shell = shell

        # arbitrary thresholds for this sim
        self.max_excitation = 100.0
        self.buffer_high = 0.8  # 80% full
        self.buffer_low = 0.2   # 20% full

    def tick(self):
        # 1. Chamber evolution
        self.chamber.tick()

        # 2. Convert excitation into "excess energy" + stress
        excess_energy = max(0.0, self.chamber.excitation_level - self.max_excitation)
        if excess_energy > 0:
            # apply stress proportional to excess
            self.shell.apply_stress(excess_energy * 0.1)

        # 3. Try to absorb some excitation into buffer
        absorbed = self.buffer.absorb(self.chamber.excitation_level * 0.05)
        # reduce excitation slightly when absorbed
        self.chamber.excitation_level -= absorbed * 0.5
        if self.chamber.excitation_level < 0:
            self.chamber.excitation_level = 0

        # 4. Passive containment relaxation
        self.shell.relieve_stress(factor=0.05)

        # 5. Guardian decisions
        self._guardian_logic()

    def _guardian_logic(self):
        # If containment is critical -> lockdown
        if self.shell.state == ContainmentState.CRITICAL:
            self.chamber.lockdown()
            # bleed buffer aggressively
            self.buffer.bleed_off(self.buffer.stored * 0.5)
            return

        # If buffer too full -> bleed off gently
        if self.buffer.snapshot()["fill_ratio"] > self.buffer_high:
            self.buffer.bleed_off(self.buffer.stored * 0.2)

        # If buffer low and containment nominal -> allow some excitation
        if (
            self.buffer.snapshot()["fill_ratio"] < self.buffer_low
            and self.shell.state == ContainmentState.NOMINAL
            and self.chamber.state != ChamberState.LOCKDOWN
        ):
            # random small excitation
            self.chamber.excite(rate=random.uniform(0.1, 0.5))
        else:
            # otherwise, damp
            if self.chamber.state == ChamberState.EXCITING:
                self.chamber.damp(rate=0.2)

    def snapshot(self):
        return {
            "chamber": self.chamber.snapshot(),
            "buffer": self.buffer.snapshot(),
            "shell": self.shell.snapshot(),
        }


def run_sim(steps: int = 100):
    chamber = ExcitationChamber()
    buffer = EnergyBuffer(capacity=500.0)
    shell = ContainmentShell(stress_threshold=50.0, critical_threshold=150.0)
    guardian = GuardianController(chamber, buffer, shell)

    history = []
    for t in range(steps):
        guardian.tick()
        snap = guardian.snapshot()
        snap["t"] = t
        history.append(snap)

    return history


if __name__ == "__main__":
    hist = run_sim(steps=50)
    for h in hist:
        print(
            f"t={h['t']:02d} | "
            f"Chamber: {h['chamber']} | "
            f"Buffer: {h['buffer']} | "
            f"Shell: {h['shell']}"
        )
