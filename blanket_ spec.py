"""
blanket_spec.py
Space LEAF Corp – Nano-Fiber Thermal Blanket System Spec

This module encodes the full specification of the nano-fiber thermal blanket:
- Nano-thread heating grid
- Fleece thermal conduction layer
- Dual solar + thermal harvesting
- Magnetic inlay strip + gold blade zipper bus
- Fingerprint-gated disengagement
- Child supervised usage mode
- Repair capsules with "nano-ants" and authenticated repair protocols
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional


# =========================
#  ENUMS & CONSTANTS
# =========================

class UsageMode(Enum):
    CHILD_SUPERVISED = "CHILD_SUPERVISED_USAGE"
    ADULT_PERSONAL = "ADULT_PERSONAL_USE"
    SERVICE = "SERVICE_MODE"


class AgeBand(Enum):
    TODDLER = "TODDLER"
    CHILD = "CHILD"
    TEEN = "TEEN"
    ADULT = "ADULT"


class SensitivityLevel(Enum):
    NORMAL = "NORMAL"
    HIGH = "HIGH_SENSITIVITY"


class ThermalEventType(Enum):
    NORMAL = "NORMAL"
    OVERHEAT_RISK = "OVERHEAT_RISK"
    HARD_SHUTDOWN = "HARD_SHUTDOWN"


class RepairEventType(Enum):
    GRID_BREAK = "GRID_BREAK"
    CONTACT_CORROSION = "CONTACT_CORROSION"
    MODULE_REPLACEMENT = "MODULE_REPLACEMENT"
    DIAGNOSTIC_ONLY = "DIAGNOSTIC_ONLY"


# =========================
#  CORE MATERIALS & LAYERS
# =========================

@dataclass
class NanoThreadGrid:
    """Nano-fiber heating and sensing grid."""
    material: str = "carbon-nanotube-infused polymer fiber"
    max_power_watts: float = 40.0
    operating_voltage: float = 5.0  # low-voltage USB-C compatible
    sensor_density_per_sq_ft: int = 16  # thermal sensors per square foot
    redundancy_factor: float = 1.5  # extra paths to avoid single-point failure


@dataclass
class FleeceThermalLayer:
    """Fleece layer for thermal conduction and comfort."""
    material: str = "recycled microfleece"
    thickness_mm: float = 2.5
    thermal_conductivity_w_mk: float = 0.035
    purpose: str = (
        "Capture and redistribute heat from nano-thread grid and body, "
        "providing comfort and mild thermal harvesting."
    )


@dataclass
class SafetyMembraneLayer:
    """Moisture and electrical isolation layer."""
    material: str = "breathable waterproof membrane"
    moisture_barrier_rating: str = "IPX4-equivalent"
    dielectric_strength_v_per_mm: float = 20.0


@dataclass
class OuterShellLayer:
    """Outer protective and aesthetic layer."""
    material: str = "durable woven polyester-cotton blend"
    abrasion_resistance_cycles: int = 20000
    wash_cycle_rating: int = 200


@dataclass
class InnerComfortLayer:
    """Skin-contact layer."""
    material: str = "hypoallergenic brushed microfiber"
    softness_rating: int = 9  # 1–10 subjective scale


@dataclass
class BlanketLayerStack:
    """Full blanket layer architecture."""
    outer_shell: OuterShellLayer = field(default_factory=OuterShellLayer)
    thermal_buffer: FleeceThermalLayer = field(default_factory=FleeceThermalLayer)
    nano_grid: NanoThreadGrid = field(default_factory=NanoThreadGrid)
    safety_membrane: SafetyMembraneLayer = field(default_factory=SafetyMembraneLayer)
    inner_comfort: InnerComfortLayer = field(default_factory=InnerComfortLayer)


# =========================
#  POWER & HARVESTING
# =========================

@dataclass
class ThermalHarvestingModule:
    """Harvests small amounts of energy from temperature gradients."""
    enabled: bool = True
    max_output_mw: float = 50.0
    description: str = (
        "Uses thermoelectric elements embedded near the fleece layer to "
        "harvest small amounts of power from body + blanket heat gradient."
    )


@dataclass
class DualSolarPanelModule:
    """Dual solar harvesting module using nano-tech panels."""
    enabled: bool = True
    panel_count: int = 2
    panel_type: str = "flexible nano-structured thin-film solar"
    max_output_watts_per_panel: float = 2.0
    placement: str = "hood exterior and upper blanket quadrant"


@dataclass
class PowerModule:
    """Central power and regulation module."""
    input_voltage_min: float = 4.5
    input_voltage_max: float = 5.5
    connector_type: str = "USB-C or magnetic dock"
    battery_capacity_wh: float = 10.0
    thermal_harvesting: ThermalHarvestingModule = field(default_factory=ThermalHarvestingModule)
    dual_solar: DualSolarPanelModule = field(default_factory=DualSolarPanelModule)
    auto_shutdown_temp_c: float = 42.0  # hard safety cap
    eco_mode_max_temp_c: float = 35.0
    child_mode_max_temp_c: float = 37.0
    adult_mode_max_temp_c: float = 40.0


# =========================
#  ZIPPER & AUTHENTICATION
# =========================

@dataclass
class GoldBladeZipperBus:
    """Conductive zipper spine."""
    plating_material: str = "gold"
    segment_count: int = 60
    max_current_amps: float = 3.0
    corrosion_resistance_rating: str = "high"
    wash_cycle_rating: int = 200


@dataclass
class MagneticInlayStrip:
    """Magnetic inlay strip for soft docking and contact alignment."""
    magnet_type: str = "flexible ferrite with conductive pads"
    contact_pairs: int = 4
    retention_force_newtons: float = 5.0
    wash_safe_when_undocked: bool = True


@dataclass
class FingerprintProfile:
    """Stored fingerprint profile for adult supervision."""
    owner_label: str
    fingerprint_id_hash: str  # hashed representation
    is_primary: bool = False


@dataclass
class FingerprintZipperLock:
    """Fingerprint-gated zipper pull and dock."""
    sensor_type: str = "capacitive"
    max_stored_profiles: int = 3
    enrolled_profiles: List[FingerprintProfile] = field(default_factory=list)
    has_mechanical_override: bool = True
    override_tool_required: bool = True
    description: str = (
        "Zipper pull houses a fingerprint sensor and micro-latch. "
        "Docked at hood top; adult fingerprint required to unlock full disengagement."
    )


# =========================
#  CHILD SUPERVISED USAGE
# =========================

@dataclass
class ChildProfile:
    """Child usage profile for supervised mode."""
    age_band: AgeBand
    sensitivity: SensitivityLevel
    label: str = "Child Profile 1"


@dataclass
class UsageConfig:
    """Usage configuration and safety thresholds."""
    mode: UsageMode = UsageMode.CHILD_SUPERVISED
    child_profile: Optional[ChildProfile] = None
    allow_mode_change_without_fingerprint: bool = False

    def get_max_temp_c(self, power_module: PowerModule) -> float:
        """Return max allowed temperature based on mode and profile."""
        if self.mode == UsageMode.CHILD_SUPERVISED:
            return min(
                power_module.child_mode_max_temp_c,
                power_module.eco_mode_max_temp_c + 2.0
            )
        elif self.mode == UsageMode.ADULT_PERSONAL:
            return power_module.adult_mode_max_temp_c
        else:  # SERVICE_MODE
            return power_module.eco_mode_max_temp_c


# =========================
#  REPAIR CAPSULES & NANO-ANTS
# =========================

@dataclass
class NanoAntRepairCapsule:
    """
    Replaceable capsule containing 'nano-ants' – micro/nano repair agents
    (conceptual: conductive ink, micro-bots, or self-healing polymer).
    """
    capsule_id: str
    location: str  # e.g., "hood_left", "hood_right", "upper_quadrant"
    capacity_repair_events: int = 5
    remaining_repair_events: int = 5
    is_active: bool = True

    def register_repair_use(self) -> bool:
        """Consume one repair event if available."""
        if not self.is_active or self.remaining_repair_events <= 0:
            return False
        self.remaining_repair_events -= 1
        if self.remaining_repair_events == 0:
            self.is_active = False
        return True


@dataclass
class RepairProtocol:
    """Authenticated repair protocol for on-site repair events."""
    requires_adult_auth: bool = True
    logs_events: bool = True
    max_auto_repairs_per_day: int = 2

    def authorize_repair(
        self,
        usage_config: UsageConfig,
        fingerprint_verified: bool
    ) -> bool:
        """Check if a repair event is allowed."""
        if self.requires_adult_auth and not fingerprint_verified:
            return False
        if usage_config.mode == UsageMode.CHILD_SUPERVISED and not fingerprint_verified:
            return False
        return True


@dataclass
class RepairEvent:
    """Record of a repair event."""
    event_type: RepairEventType
    capsule_id: str
    fingerprint_verified: bool
    success: bool
    notes: str = ""


# =========================
#  MAIN BLANKET SYSTEM
# =========================

@dataclass
class NanoFiberThermalBlanketSystem:
    """Full system spec for the nano-fiber thermal blanket."""
    layers: BlanketLayerStack = field(default_factory=BlanketLayerStack)
    power: PowerModule = field(default_factory=PowerModule)
    zipper_bus: GoldBladeZipperBus = field(default_factory=GoldBladeZipperBus)
    magnetic_strip: MagneticInlayStrip = field(default_factory=MagneticInlayStrip)
    fingerprint_lock: FingerprintZipperLock = field(default_factory=FingerprintZipperLock)
    usage_config: UsageConfig = field(default_factory=UsageConfig)
    repair_capsules: List[NanoAntRepairCapsule] = field(default_factory=list)
    repair_protocol: RepairProtocol = field(default_factory=RepairProtocol)
    repair_log: List[RepairEvent] = field(default_factory=list)

    def register_child_profile(self, profile: ChildProfile):
        self.usage_config.mode = UsageMode.CHILD_SUPERVISED
        self.usage_config.child_profile = profile

    def authenticate_adult(self, fingerprint_hash: str) -> bool:
        """Simulate fingerprint verification."""
        for profile in self.fingerprint_lock.enrolled_profiles:
            if profile.fingerprint_id_hash == fingerprint_hash:
                return True
        return False

    def evaluate_thermal_event(self, sensed_temp_c: float) -> ThermalEventType:
        """Evaluate thermal state and decide if shutdown is needed."""
        max_temp = self.usage_config.get_max_temp_c(self.power)
        if sensed_temp_c <= max_temp:
            return ThermalEventType.NORMAL
        elif sensed_temp_c <= self.power.auto_shutdown_temp_c:
            return ThermalEventType.OVERHEAT_RISK
        else:
            return ThermalEventType.HARD_SHUTDOWN

    def perform_repair(
        self,
        event_type: RepairEventType,
        capsule_id: str,
        fingerprint_hash: Optional[str] = None
    ) -> RepairEvent:
        """Coordinate a repair event using nano-ant capsules."""
        fingerprint_verified = (
            fingerprint_hash is not None and self.authenticate_adult(fingerprint_hash)
        )
        authorized = self.repair_protocol.authorize_repair(
            self.usage_config,
            fingerprint_verified
        )

        capsule = next((c for c in self.repair_capsules if c.capsule_id == capsule_id), None)
        if not authorized or capsule is None:
            event = RepairEvent(
                event_type=event_type,
                capsule_id=capsule_id,
                fingerprint_verified=fingerprint_verified,
                success=False,
                notes="Unauthorized or capsule not found."
            )
            self.repair_log.append(event)
            return event

        success = capsule.register_repair_use()
        notes = "Repair executed." if success else "Capsule depleted or inactive."
        event = RepairEvent(
            event_type=event_type,
            capsule_id=capsule_id,
            fingerprint_verified=fingerprint_verified,
            success=success,
            notes=notes
        )
        self.repair_log.append(event)
        return event


# =========================
#  FACTORY / EXAMPLE SETUP
# =========================

def create_default_blanket_system() -> NanoFiberThermalBlanketSystem:
    """Create a default configured blanket system."""
    system = NanoFiberThermalBlanketSystem()

    # Example enrolled adult fingerprints (hashes are placeholders)
    system.fingerprint_lock.enrolled_profiles = [
        FingerprintProfile(owner_label="Parent A", fingerprint_id_hash="hash_parent_a", is_primary=True),
        FingerprintProfile(owner_label="Parent B", fingerprint_id_hash="hash_parent_b", is_primary=False),
    ]

    # Example child profile
    system.register_child_profile(
        ChildProfile(age_band=AgeBand.CHILD, sensitivity=SensitivityLevel.NORMAL)
    )

    # Example repair capsules in hood
    system.repair_capsules = [
        NanoAntRepairCapsule(capsule_id="hood_left", location="hood_left"),
        NanoAntRepairCapsule(capsule_id="hood_right", location="hood_right"),
    ]

    return system


if __name__ == "__main__":
    # Example: spin up system and simulate a repair event
    blanket = create_default_blanket_system()
    event = blanket.perform_repair(
        event_type=RepairEventType.GRID_BREAK,
        capsule_id="hood_left",
        fingerprint_hash="hash_parent_a"
    )
    print("Repair event:", event)
