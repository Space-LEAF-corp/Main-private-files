from dataclasses import dataclass, field
from typing import List, Dict


@dataclass
class Queen:
    id: str
    planet_id: str
    role: str
    abilities: List[str]


@dataclass
class Planet:
    id: str
    name: str
    role: str
    environment: str
    queen_id: str
    king_id: str
    lessons: List[str]
    rewards: Dict[str, str]


@dataclass
class Moonbase:
    id: str
    puzzles_completed: bool = False
    kodex_activated: bool = False

    def can_activate_kodex(self, environmental_damage_low: bool, data_collected: bool) -> bool:
        return self.puzzles_completed and environmental_damage_low and data_collected

    def activate_kodex(self) -> Dict[str, bool]:
        if not self.kodex_activated:
            self.kodex_activated = True
            return {
                "mars_unlocked": True,
                "neptune_unlocked": True,
                "jupiter_unlocked": True,
                "saturn_revealed": True
            }
        return {}


@dataclass
class SaturnSanctuary:
    id: str
    civilizations_present: List[str] = field(default_factory=list)
    forge_active: bool = False

    def register_civilization(self, civ_id: str):
        if civ_id not in self.civilizations_present:
            self.civilizations_present.append(civ_id)
        if len(self.civilizations_present) >= 3:
            self.forge_active = True

    def can_build_guardian(self) -> bool:
        return self.forge_active


@dataclass
class EarthCompanion:
    id: str = "queen_earth_black_widow"
    role: str = "digital_web_anchor"

    def guide_player(self, context: str) -> str:
        # Simple placeholder for guidance logic
        if context == "lost":
            return "A soft tug on your shoulder: a safe path reveals itself in the distance."
        if context == "danger":
            return "You feel her tense—something ahead is wrong. Move carefully."
        return "She rests quietly, watching, listening, weaving the web around your journey."


# Example: wiring the core pieces
moon = Moonbase(id="moon_kodex_sanctuary")
saturn = SaturnSanctuary(id="saturn_sanctuary_world")
earth_queen = EarthCompanion()

# Example flow (pseudo-runtime):
# 1) Player completes cave puzzles
moon.puzzles_completed = True

# 2) Check if Kodex can activate
if moon.can_activate_kodex(environmental_damage_low=True, data_collected=True):
    unlocks = moon.activate_kodex()
    # unlocks -> mars/neptune/jupiter/saturn

# 3) After visiting Mars, Neptune, Jupiter:
saturn.register_civilization("mars_spider_civ")
saturn.register_civilization("neptune_spider_civ")
saturn.register_civilization("jupiter_spider_civ")

if saturn.can_build_guardian():
    # Guardian assembly would be triggered here
    guardian_blueprint = "ultimate_spider_guardian_blueprint"