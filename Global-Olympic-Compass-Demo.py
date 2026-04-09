from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Dict, Optional, Callable
import math


# ---------- Core time window + access control ----------

@dataclass
class SeasonWindow:
    name: str
    start: datetime
    end: datetime

    def is_active(self, now: Optional[datetime] = None) -> bool:
        now = now or datetime.utcnow()
        return self.start <= now <= self.end


@dataclass
class TemporaryAccessToken:
    user_id: str
    season: SeasonWindow
    created_at: datetime = field(default_factory=datetime.utcnow)

    def is_valid(self, now: Optional[datetime] = None) -> bool:
        return self.season.is_active(now)


# ---------- Simple geo helpers (2D for now) ----------

@dataclass
class Point:
    x: float  # could be lon
    y: float  # could be lat


def distance(a: Point, b: Point) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2)


# ---------- Places: bathrooms, safe zones, custom arms ----------

@dataclass
class Place:
    id: str
    name: str
    kind: str  # "bathroom", "safe_zone", "family_node", etc.
    location: Point
    metadata: Dict[str, str] = field(default_factory=dict)


@dataclass
class UserNetworkArm:
    """Represents one 'arm' of the compass: e.g. family, friends, crew."""
    name: str
    filter_fn: Callable[[Place], bool]


# ---------- The Olympic Vocal Compass ----------

@dataclass
class OlympicCompass:
    season: SeasonWindow
    places: List[Place] = field(default_factory=list)
    user_arms: Dict[str, List[UserNetworkArm]] = field(default_factory=dict)

    # ---- Access control ----
    def require_valid_access(self, token: TemporaryAccessToken) -> None:
        if not token.is_valid():
            raise PermissionError("Access outside official Olympic season is not allowed.")

    # ---- Registration ----
    def register_place(self, place: Place) -> None:
        self.places.append(place)

    def register_user_arm(self, user_id: str, arm: UserNetworkArm) -> None:
        self.user_arms.setdefault(user_id, []).append(arm)

    # ---- Core queries ----
    def nearest_bathroom(self, token: TemporaryAccessToken, user_location: Point) -> Optional[Place]:
        self.require_valid_access(token)
        bathrooms = [p for p in self.places if p.kind == "bathroom"]
        if not bathrooms:
            return None
        return min(bathrooms, key=lambda p: distance(user_location, p.location))

    def nearest_safe_zone(self, token: TemporaryAccessToken, user_location: Point) -> Optional[Place]:
        self.require_valid_access(token)
        safe_zones = [p for p in self.places if p.kind == "safe_zone"]
        if not safe_zones:
            return None
        return min(safe_zones, key=lambda p: distance(user_location, p.location))

    def nearest_in_arm(self, token: TemporaryAccessToken, user_location: Point, arm_name: str) -> Optional[Place]:
        self.require_valid_access(token)
        arms = self.user_arms.get(token.user_id, [])
        arm_filters = [arm.filter_fn for arm in arms if arm.name == arm_name]
        if not arm_filters:
            return None

        # Combine all filters for this arm
        def combined_filter(p: Place) -> bool:
            return any(f(p) for f in arm_filters)

        candidates = [p for p in self.places if combined_filter(p)]
        if not candidates:
            return None
        return min(candidates, key=lambda p: distance(user_location, p.location))


# ---------- Example wiring (would live in a main/app file) ----------

if __name__ == "__main__":
    # Define an Olympic season
    season = SeasonWindow(
        name="Space Olympics 2032",
        start=datetime(2032, 7, 1),
        end=datetime(2032, 7, 31),
    )

    compass = OlympicCompass(season=season)

    # Register some places
    compass.register_place(Place(
        id="b1",
        name="Main Concourse Bathroom A",
        kind="bathroom",
        location=Point(0, 0),
    ))
    compass.register_place(Place(
        id="s1",
        name="Primary Safe Assembly Zone",
        kind="safe_zone",
        location=Point(10, 5),
    ))
    compass.register_place(Place(
        id="fam1",
        name="Family Meet Spot - North Gate",
        kind="family_node",
        location=Point(3, 4),
        metadata={"family_tag": "SOGGE_CREW"},
    ))

    # User token
    token = TemporaryAccessToken(user_id="user_123", season=season)

    # Define a family arm for this user
    def family_filter(place: Place) -> bool:
        return place.kind == "family_node" and place.metadata.get("family_tag") == "SOGGE_CREW"

    compass.register_user_arm(
        user_id="user_123",
        arm=UserNetworkArm(name="family", filter_fn=family_filter)
    )

    # Simulate a user location
    user_loc = Point(1, 1)

    # Queries (these would be hooked to voice/UI in real app)
    nearest_bath = compass.nearest_bathroom(token, user_loc)
    nearest_safe = compass.nearest_safe_zone(token, user_loc)
    nearest_family = compass.nearest_in_arm(token, user_loc, "family")

    print("Nearest bathroom:", nearest_bath)
    print("Nearest safe zone:", nearest_safe)
    print("Nearest family node:", nearest_family)
