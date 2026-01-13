"""
Supersonic Echo Lock & Dormant Hyperlane Module
-----------------------------------------------

Concepts implemented (purely as math / logic):

1. Supersonic Echo Lock (convoy following)
   - Distance lock via echo-time measurement
   - Relative velocity lock
   - Works with any generic propulsion / power source

2. Autonomous Convoy Lane
   - Lane path defined in time
   - Leader and follower roles
   - Follower uses echo lock to stay in formation

3. Dormant Hyperlane Timing Layer (Hyperdrive Framework)
   - t = ((1 - f) * d) / beta
   - Used ONLY as a mission-planning / timing abstraction
   - Marked dormant and inactive; does NOT control hardware

This file is intentionally:
- Energy-source agnostic (no reactor design, no fusion details)
- Mechanism-agnostic (no spacetime / warp implementation)
- Navigation- and ethics-focused (readiness gates, safety-first)

Safe to store, share, and use as conceptual scaffolding.
"""

from dataclasses import dataclass
from typing import Callable, Tuple
import math


# ============================================================
# 1. Core data structures
# ============================================================

@dataclass
class Vec3:
    """Simple 3D vector for positions and velocities."""
    x: float
    y: float
    z: float

    def __add__(self, other: "Vec3") -> "Vec3":
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: "Vec3") -> "Vec3":
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, s: float) -> "Vec3":
        return Vec3(self.x * s, self.y * s, self.z * s)

    def norm(self) -> float:
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)


@dataclass
class ConvoyConfig:
    """Configuration for Supersonic Echo Lock convoy behavior."""
    d_target: float          # target separation distance (m)
    eps_d: float             # allowable distance error (m)
    eps_v: float             # allowable relative speed error (m/s)
    c_signal: float          # signal speed (m/s) (sound in air, EM in space)
    v_min: float             # minimum allowed convoy speed magnitude (m/s)
    v_max: float             # maximum allowed convoy speed magnitude (m/s)


@dataclass
class ConvoyState:
    """State of a leader-follower convoy pair at a given instant."""
    t: float                 # current time (s)
    x_leader: Vec3           # leader position
    v_leader: Vec3           # leader velocity
    x_follower: Vec3         # follower position
    v_follower: Vec3         # follower velocity


@dataclass
class HyperlaneParams:
    """
    Parameters for the dormant hyperlane timing layer.

    t = ((1 - f) * d) / beta

    NOTE: This is a mission-planning abstraction only.
    It does NOT control propulsion or spacetime. Dormant until explicitly activated.
    """
    f: float                 # collapse fraction [0, 1)
    beta: float              # effective speed factor (>0)
    dormant: bool = True     # MUST remain True until future tech & safety proven


@dataclass
class ReadinessContext:
    """
    Context used to decide if any advanced mode is allowed.

    All fields are abstract indicators (0.0–1.0 or booleans).
    They do NOT connect to any real hardware in this module.
    """
    tech_readiness: float        # e.g., TRL-like metric [0, 1]
    safety_certainty: float      # confidence in safety [0, 1]
    ethical_clearance: float     # ethical / governance readiness [0, 1]
    crew_onboard: bool           # whether humans are on the craft
    emergency_only: bool         # whether we are in an actual emergency mode


# ============================================================
# 2. Supersonic Echo Lock (convoy navigation)
# ============================================================

def echo_round_trip_time(cfg: ConvoyConfig, state: ConvoyState) -> float:
    """
    Compute the round-trip echo time between follower and leader.

    follower -> leader -> follower
    Δt_echo = 2 * distance / c_signal
    """
    r = state.x_follower - state.x_leader
    d = r.norm()
    return 2.0 * d / max(cfg.c_signal, 1e-9)


def distance_lock_error(cfg: ConvoyConfig, state: ConvoyState) -> float:
    """
    Distance lock error: |d - d_target|.
    """
    r = state.x_follower - state.x_leader
    d = r.norm()
    return abs(d - cfg.d_target)


def relative_speed(state: ConvoyState) -> float:
    """
    Magnitude of relative velocity between follower and leader.
    """
    v_rel = state.v_follower - state.v_leader
    return v_rel.norm()


def is_echo_lock_satisfied(cfg: ConvoyConfig, state: ConvoyState) -> bool:
    """
    Check if the convoy meets the Echo Lock constraints:

    1) Distance near target
    2) Relative speed near zero
    3) Both vehicles within safe speed band
    """
    d_err = distance_lock_error(cfg, state)
    v_rel_mag = relative_speed(state)

    vL = state.v_leader.norm()
    vF = state.v_follower.norm()

    distance_ok = d_err <= cfg.eps_d
    rel_speed_ok = v_rel_mag <= cfg.eps_v

    # Safe band: v_min <= speed <= v_max
    leader_speed_ok = (cfg.v_min <= vL <= cfg.v_max)
    follower_speed_ok = (cfg.v_min <= vF <= cfg.v_max)

    return distance_ok and rel_speed_ok and leader_speed_ok and follower_speed_ok


def follower_velocity_update(
    cfg: ConvoyConfig,
    state: ConvoyState,
    k_p: float = 0.5,
    k_v: float = 0.5
) -> Vec3:
    """
    Simple proportional controller for follower velocity
    to approach distance + relative velocity lock.

    This is conceptual only and does not touch any propulsion model.
    """
    # Distance error direction
    r = state.x_follower - state.x_leader
    d = r.norm()
    if d < 1e-9:
        # Avoid division by zero; return current velocity
        return state.v_follower

    # Unit vector from leader to follower
    u = Vec3(r.x / d, r.y / d, r.z / d)

    # Distance error
    d_err = d - cfg.d_target

    # Relative velocity
    v_rel = state.v_follower - state.v_leader

    # Control law (conceptual)
    # Move opposite direction of distance error and relative velocity
    v_correction = u * (-k_p * d_err) + v_rel * (-k_v)

    # Target follower velocity
    v_new = state.v_leader + v_correction

    # Optionally clip to safe band magnitude
    v_mag = v_new.norm()
    if v_mag > cfg.v_max:
        scale = cfg.v_max / max(v_mag, 1e-9)
        v_new = v_new * scale
    elif v_mag < cfg.v_min:
        # Scale up to minimum if not zero
        if v_mag > 1e-9:
            scale = cfg.v_min / v_mag
            v_new = v_new * scale

    return v_new


# ============================================================
# 3. Dormant Hyperlane Timing (hyperdrive framework)
# ============================================================

def hyperlane_travel_time(params: HyperlaneParams, distance: float) -> float:
    """
    Hyperlane timing function (mission-planning abstraction):

        t = ((1 - f) * d) / beta

    - distance: scalar path length (e.g., meters or normalized units)
    - f: collapse fraction (0 <= f < 1)
    - beta: effective speed factor (>0)

    NOTE:
    - This DOES NOT perform any real spacetime manipulation.
    - It is a dormant planning model only.
    """
    f = max(0.0, min(params.f, 0.999999))   # clamp for safety of abstraction
    beta = max(params.beta, 1e-9)
    d_eff = (1.0 - f) * distance
    return d_eff / beta


def hyperlane_is_active(params: HyperlaneParams, readiness: ReadinessContext) -> bool:
    """
    Determine if the hyperlane timing layer is allowed to be ACTIVE.

    By design, this should return False in current era. The logic is here
    to make the "dormant" nature explicit and auditable.

    Conditions (example, conservative):
    - params.dormant must be False (explicitly flipped by future, verified tech)
    - tech_readiness, safety_certainty, ethical_clearance all near 1.0
    - emergency_only is respected: if True, activation only in real emergencies
    - crew_onboard may tighten requirements further

    In our current use, we keep params.dormant = True,
    so this function always returns False.
    """
    if params.dormant:
        return False

    # Thresholds are placeholders for future governance decisions.
    tech_ok = readiness.tech_readiness >= 0.99
    safety_ok = readiness.safety_certainty >= 0.99
    ethics_ok = readiness.ethical_clearance >= 0.99

    if not (tech_ok and safety_ok and ethics_ok):
        return False

    # Example: only allow if emergency-only mode is satisfied
    if readiness.emergency_only is False and readiness.crew_onboard:
        # Example policy: no activation for routine crewed travel
        return False

    return True


# ============================================================
# 4. Unified "Autonomous Convoy Hyperlane" interface
# ============================================================

@dataclass
class LaneWaypoint:
    """A simple waypoint in spacetime for the convoy lane."""
    t: float       # target time (s)
    pos: Vec3      # target position (arbitrary frame)


@dataclass
class HyperlanePlan:
    """
    A high-level plan for a route:

    - total_distance: scalar distance along the path
    - baseline_time: time without any collapse factor (f=0, beta=1)
    - hyper_time: hypothetical hyperlane time, if ever allowed
    """
    total_distance: float
    baseline_time: float
    hyper_time: float


def plan_hyperlane_mission(
    distance: float,
    params: HyperlaneParams
) -> HyperlanePlan:
    """
    Plan a notional mission using:
    - baseline travel: t = d (normalized units; e.g., d=years at c=1)
    - hyperlane travel: t = ((1 - f) * d) / beta (abstract)

    This does not control any real vehicle.
    """
    baseline_time = distance  # assuming v=1 in normalized units
    hyper_time = hyperlane_travel_time(params, distance)
    return HyperlanePlan(
        total_distance=distance,
        baseline_time=baseline_time,
        hyper_time=hyper_time
    )


# ============================================================
# 5. Example usage (safe, conceptual)
# ============================================================

if __name__ == "__main__":
    # Example convoy configuration for a "space lane"
    cfg = ConvoyConfig(
        d_target=1000.0,          # 1 km separation
        eps_d=5.0,               # +/- 5 m tolerance
        eps_v=0.1,               # 0.1 m/s relative speed tolerance
        c_signal=3.0e8,          # speed of light (EM echo in space)
        v_min=100.0,             # min convoy speed (m/s) - arbitrary
        v_max=3.0e4              # max convoy speed (m/s) - arbitrary
    )

    # Simple initial state
    state = ConvoyState(
        t=0.0,
        x_leader=Vec3(0.0, 0.0, 0.0),
        v_leader=Vec3(1000.0, 0.0, 0.0),     # leader moving along +x
        x_follower=Vec3(-1200.0, 0.0, 0.0),  # follower starts 200 m too far
        v_follower=Vec3(900.0, 0.0, 0.0)     # follower slower than leader
    )

    # Update follower velocity toward Echo Lock
    v_new = follower_velocity_update(cfg, state)
    print("New follower velocity (conceptual):", v_new)

    # Check lock satisfaction
    lock_ok = is_echo_lock_satisfied(cfg, state)
    print("Initial Echo Lock satisfied?", lock_ok)

    # Dormant hyperlane params (stays inactive)
    hyper_params = HyperlaneParams(
        f=0.35,
        beta=1.0,
        dormant=True  # CRITICAL: remains True in current era
    )

    # Abstract mission planning example
    distance_normalized = 4.0  # e.g., "4 light-years at c=1"
    plan = plan_hyperlane_mission(distance_normalized, hyper_params)

    print("Hyperlane mission plan (conceptual):")
    print("  Distance (normalized):", plan.total_distance)
    print("  Baseline time:", plan.baseline_time)
    print("  Hyperlane time (dormant model):", plan.hyper_time)

    # Readiness (intentionally low, to keep hyperlane inactive)
    readiness = ReadinessContext(
        tech_readiness=0.1,
        safety_certainty=0.1,
        ethical_clearance=0.1,
        crew_onboard=False,
        emergency_only=False
    )

    print("Hyperlane active?", hyperlane_is_active(hyper_params, readiness))
