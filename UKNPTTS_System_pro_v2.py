UKNPTTS_System_pro_v2.py
# ================================================================
#  UKNPTTS — UK/UN Planetary Transportation Train System
#  CONFIG HEADER — PROFESSIONAL EDITION
#
#  Author: Leif William Sogge
#  Seal: UKNPTTS • Interplanetary Rail Stewardship Division
#
#  Mission Statement:
#  "To build safe, resilient, robot‑first transportation networks
#   that prepare humanity for life beyond Earth."
# ================================================================

from enum import Enum
import random
import math

# -----------------------------
# Core enums and base classes
# -----------------------------

class Role(Enum):
    ROBOT = "robot"
    HUMAN_REMOTE = "human_remote"
    HUMAN_ONSITE = "human_onsite"

class LocationType(Enum):
    EARTH = "Earth"
    MOON = "Moon"
    MARS = "Mars"
    ORBIT = "Orbit"
    TUNNEL = "Tunnel"
    HABITAT = "Habitat"

class ModuleState(Enum):
    IDLE = 0
    ALIGNING = 1
    LOCKED = 2
    HAZARD = 3
    OFFLINE = 4

class RailMode(Enum):
    CONVENTIONAL = 0
    MAGNETIC = 1

# -----------------------------
# Track and modular rail system
# -----------------------------

class TrackModule:
    """
    One modular, self-aligning rail segment.
    Can represent conventional or magnetic rail.
    """
    def __init__(self, module_id, position, target_position, rail_mode=RailMode.MAGNETIC):
        self.module_id = module_id
        self.position = position              # (x, y, z)
        self.target_position = target_position
        self.state = ModuleState.IDLE
        self.locked = False
        self.hazard_flag = False
        self.rail_mode = rail_mode

    def read_alignment_error(self):
        dx = self.target_position[0] - self.position[0]
        dy = self.target_position[1] - self.position[1]
        dz = self.target_position[2] - self.position[2]
        error = math.sqrt(dx*dx + dy*dy + dz*dz)
        return error, (dx, dy, dz)

    def read_hazard_sensors(self):
        # Placeholder: simulate rare hazard
        self.hazard_flag = (random.random() < 0.001)
        return self.hazard_flag

    def apply_actuators(self, delta):
        # Move a fraction toward target
        step_scale = 0.1
        self.position = (
            self.position[0] + delta[0] * step_scale,
            self.position[1] + delta[1] * step_scale,
            self.position[2] + delta[2] * step_scale,
        )

    def lock_in_place(self):
        self.locked = True
        self.state = ModuleState.LOCKED

    def unlock(self):
        self.locked = False
        self.state = ModuleState.IDLE

    def update(self, max_allowed_error):
        if self.read_hazard_sensors():
            self.state = ModuleState.HAZARD
            self.locked = True
            return

        if self.state == ModuleState.HAZARD:
            return

        error, delta = self.read_alignment_error()
        if error <= max_allowed_error:
            if not self.locked:
                self.lock_in_place()
            self.state = ModuleState.LOCKED
        else:
            self.state = ModuleState.ALIGNING
            self.unlock()
            self.apply_actuators(delta)


class TrackController:
    """
    Oversees a set of track modules, keeps them aligned and safe.
    """
    def __init__(self, modules, max_allowed_error=0.002):
        self.modules = modules
        self.max_allowed_error = max_allowed_error
        self.global_halt = False

    def update_all_modules(self):
        for m in self.modules:
            m.update(self.max_allowed_error)

    def check_global_safety(self):
        self.global_halt = any(m.state == ModuleState.HAZARD for m in self.modules)

    def step(self):
        self.update_all_modules()
        self.check_global_safety()
        return {
            "global_halt": self.global_halt,
            "modules": [
                {
                    "id": m.module_id,
                    "state": m.state.name,
                    "locked": m.locked,
                    "hazard": m.hazard_flag,
                    "rail_mode": m.rail_mode.name,
                }
                for m in self.modules
            ],
        }

# -----------------------------
# Safety spine: plasma + thermal power + comms
# -----------------------------

class ThermalHeatGenerator:
    """
    THG — Thermal Heat Generator
    Uses temperature gradient (plasma heating vs inner structure) to
    produce emergency power and heat flux data.
    """
    def __init__(self):
        self.last_power_output_w = 0.0
        self.last_heat_flux = 0.0

    def sense_heat(self, surface_temp, inner_temp):
        gradient = max(0.0, surface_temp - inner_temp)
        self.last_heat_flux = gradient
        return gradient

    def generate_power(self, gradient):
        # Conceptual: more gradient = more emergency power
        efficiency = 0.01  # placeholder
        self.last_power_output_w = gradient * efficiency
        return self.last_power_output_w


class PlasmaEmergencyPower:
    """
    PEP — Plasma Emergency Power
    Conceptual MHD-style power extraction from ionized flow.
    """
    def __init__(self):
        self.last_power_output_w = 0.0

    def generate_power(self, plasma_density, flow_speed):
        # Conceptual: power ~ density * speed^2 * small factor
        k = 0.0001
        self.last_power_output_w = k * plasma_density * (flow_speed ** 2)
        return self.last_power_output_w


class PlasmaAntennaNode:
    """
    Plasma antenna / 'lightning rod' node.
    Uses available emergency power to push minimal safety beacons
    through any available RF path.
    """
    def __init__(self):
        self.last_tx_success = False

    def push_minimal_packet(self, beacon, link_conditions, available_power_w):
        # Conceptual: if link_conditions and power are above thresholds, succeed
        min_power = 0.1
        if not link_conditions.get("path_open", False):
            self.last_tx_success = False
            return False
        if available_power_w < min_power:
            self.last_tx_success = False
            return False
        # If we reach here, assume packet gets out
        self.last_tx_success = True
        return True


class CrewSafetyPlasmaMonitor:
    """
    CSPM — Crew Safety Plasma Monitor
    Watches plasma, thermal data, and crew/vehicle health.
    Builds minimal safety beacons and uses PlasmaAntennaNode to transmit.
    """
    def __init__(self):
        self.crew_state = "UNKNOWN"   # OK / DEGRADED / CRITICAL
        self.plasma_state = "UNKNOWN" # NOMINAL / HOT / UNSTABLE
        self.thermal_state = "UNKNOWN"
        self.last_beacon = None

    def sense_plasma(self, rf_data, probe_data):
        # Placeholder: analyze RF reflections, sheath density, etc.
        # For now, simple thresholds on "density"
        density = probe_data.get("density", 0.0)
        if density < 1.0:
            self.plasma_state = "NOMINAL"
        elif density < 5.0:
            self.plasma_state = "HOT"
        else:
            self.plasma_state = "UNSTABLE"

    def sense_thermal(self, heat_flux):
        if heat_flux < 10.0:
            self.thermal_state = "COOL"
        elif heat_flux < 50.0:
            self.thermal_state = "WARM"
        else:
            self.thermal_state = "STRESSED"

    def sense_crew(self, sensors):
        # Placeholder: cabin pressure, temp, life support, G-load, etc.
        if sensors.get("life_support_ok", True) and sensors.get("cabin_pressure_ok", True):
            self.crew_state = "OK"
        else:
            self.crew_state = "DEGRADED"

    def build_beacon(self):
        return {
            "crew_state": self.crew_state,
            "plasma_state": self.plasma_state,
            "thermal_state": self.thermal_state,
        }

    def try_transmit_beacon(self, plasma_antenna, link_conditions, available_power_w):
        beacon = self.build_beacon()
        success = plasma_antenna.push_minimal_packet(beacon, link_conditions, available_power_w)
        if success:
            self.last_beacon = beacon
        return success

# -----------------------------
# Tunneling and construction robots
# -----------------------------

class TBMRobot:
    """
    Tunneling machine with robotic arms for track placement.
    Operated remotely during the 5-year robot-only phase.
    """
    def __init__(self, robot_id, location=LocationType.MOON):
        self.robot_id = robot_id
        self.location = location
        self.progress_m = 0.0
        self.operational = True

    def bore_step(self, meters=1.0):
        if not self.operational:
            return
        self.progress_m += meters

    def place_track_modules(self, track_controller, modules_per_step=1):
        for _ in range(modules_per_step):
            new_id = len(track_controller.modules)
            target_pos = (self.progress_m, 0.0, 0.0)
            start_pos = (
                target_pos[0] + random.uniform(-0.05, 0.05),
                target_pos[1] + random.uniform(-0.05, 0.05),
                target_pos[2] + random.uniform(-0.05, 0.05),
            )
            module = TrackModule(
                module_id=new_id,
                position=start_pos,
                target_position=target_pos,
                rail_mode=RailMode.MAGNETIC,
            )
            track_controller.modules.append(module)

    def step(self, track_controller):
        self.bore_step(meters=1.0)
        self.place_track_modules(track_controller, modules_per_step=1)

# -----------------------------
# Worker rail cars and habitats
# -----------------------------

class RailCarType(Enum):
    WORKER_TRANSPORT = 0
    BREAK_ROOM = 1
    BUNK_MEN = 2
    BUNK_WOMEN = 3

class RailCar:
    """
    Represents a single rail car: worker transport, break room, or bunk.
    """
    def __init__(self, car_id, car_type, capacity):
        self.car_id = car_id
        self.car_type = car_type
        self.capacity = capacity
        self.occupants = 0
        self.location = LocationType.TUNNEL

    def board(self, count):
        self.occupants = min(self.capacity, self.occupants + count)

    def disembark(self, count):
        self.occupants = max(0, self.occupants - count)


class HabitatTrain:
    """
    Worker habitat consist:
    - Women’s bunk car
    - Men’s bunk car
    - Shared break-room car
    - Worker transport car
    """
    def __init__(self):
        self.cars = []
        self._build_default_consist()

    def _build_default_consist(self):
        self.cars.append(RailCar(0, RailCarType.BUNK_WOMEN, capacity=20))
        self.cars.append(RailCar(1, RailCarType.BUNK_MEN, capacity=20))
        self.cars.append(RailCar(2, RailCarType.BREAK_ROOM, capacity=30))
        self.cars.append(RailCar(3, RailCarType.WORKER_TRANSPORT, capacity=40))

    def status(self):
        return [
            {
                "car_id": c.car_id,
                "type": c.car_type.name,
                "capacity": c.capacity,
                "occupants": c.occupants,
            }
            for c in self.cars
        ]

# -----------------------------
# Mission control and 5-year remote phase
# -----------------------------

class MissionPhase(Enum):
    YEAR_1_SURVEY = 1
    YEAR_2_TUNNEL_START = 2
    YEAR_3_NETWORK_EXPAND = 3
    YEAR_4_HABITAT_BUILD = 4
    YEAR_5_FULL_TEST = 5

class MissionControl:
    """
    High-level controller for the 5-year remote-only UKNPTTS build-out.
    Robots only on-site; humans remote.
    Includes safety spine: THG + PEP + CSPM + PlasmaAntennaNode.
    """
    def __init__(self):
        self.phase = MissionPhase.YEAR_1_SURVEY
        self.track_controller = TrackController(modules=[])
        self.tbm_robot = TBMRobot(robot_id="TBM-01", location=LocationType.MOON)
        self.habitat_train = HabitatTrain()

        # Safety spine
        self.thg = ThermalHeatGenerator()
        self.pep = PlasmaEmergencyPower()
        self.cspm = CrewSafetyPlasmaMonitor()
        self.plasma_antenna = PlasmaAntennaNode()

    def advance_phase_if_ready(self):
        if self.phase == MissionPhase.YEAR_1_SURVEY and len(self.track_controller.modules) > 0:
            self.phase = MissionPhase.YEAR_2_TUNNEL_START
        elif self.phase == MissionPhase.YEAR_2_TUNNEL_START and self.tbm_robot.progress_m > 1000:
            self.phase = MissionPhase.YEAR_3_NETWORK_EXPAND
        elif self.phase == MissionPhase.YEAR_3_NETWORK_EXPAND and self.tbm_robot.progress_m > 5000:
            self.phase = MissionPhase.YEAR_4_HABITAT_BUILD
        elif self.phase == MissionPhase.YEAR_4_HABITAT_BUILD:
            self.phase = MissionPhase.YEAR_5_FULL_TEST

    def simulate_safety_spine(self):
        # Conceptual environment values
        surface_temp = random.uniform(200.0, 2000.0)
        inner_temp = random.uniform(200.0, 400.0)
        plasma_density = random.uniform(0.0, 10.0)
        flow_speed = random.uniform(0.0, 8000.0)

        # THG: thermal gradient → power + heat flux
        gradient = self.thg.sense_heat(surface_temp, inner_temp)
        thg_power = self.thg.generate_power(gradient)

        # PEP: plasma flow → power
        pep_power = self.pep.generate_power(plasma_density, flow_speed)

        # Total emergency power available to safety systems
        emergency_power = thg_power + pep_power

        # CSPM: sense plasma + thermal + crew
        rf_data = {}  # placeholder
        probe_data = {"density": plasma_density}
        self.cspm.sense_plasma(rf_data, probe_data)
        self.cspm.sense_thermal(self.thg.last_heat_flux)
        crew_sensors = {"life_support_ok": True, "cabin_pressure_ok": True}
        self.cspm.sense_crew(crew_sensors)

        # Try to transmit minimal beacon
        link_conditions = {"path_open": True}  # conceptual
        tx_success = self.cspm.try_transmit_beacon(
            self.plasma_antenna,
            link_conditions,
            emergency_power_w=emergency_power
        )

        return {
            "emergency_power_w": emergency_power,
            "tx_success": tx_success,
            "crew_state": self.cspm.crew_state,
            "plasma_state": self.cspm.plasma_state,
            "thermal_state": self.cspm.thermal_state,
        }

    def simulate_step(self):
        # Construction / tunneling
        self.tbm_robot.step(self.track_controller)
        track_status = self.track_controller.step()
        self.advance_phase_if_ready()
        habitat_status = self.habitat_train.status()

        # Safety spine
        safety_status = self.simulate_safety_spine()

        snapshot = {
            "phase": self.phase.name,
            "tbm_progress_m": self.tbm_robot.progress_m,
            "track_global_halt": track_status["global_halt"],
            "num_track_modules": len(self.track_controller.modules),
            "habitat_train": habitat_status,
            "safety": safety_status,
        }
        return snapshot


if __name__ == "__main__":
    mc = MissionControl()
    for step in range(2000):
        snapshot = mc.simulate_step()
        if step % 200 == 0:
            print(f"Step {step}")
            print(f"  Phase: {snapshot['phase']}")
            print(f"  TBM progress (m): {snapshot['tbm_progress_m']:.1f}")
            print(f"  Track modules: {snapshot['num_track_modules']}")
            print(f"  Track global halt: {snapshot['track_global_halt']}")
            print(f"  Habitat cars: {snapshot['habitat_train']}")
            print(f"  Safety: {snapshot['safety']}")
            print("-" * 40)
