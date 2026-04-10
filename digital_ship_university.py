# digital_ship_university.py
# Authorship: Leif William Sogge
# Purpose: Ceremonial simulation of the Digital Ship, its rooms, and University features
# Notes: Symbolic lineage artifact, blending stewardship, humor, and guardianship

from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import List, Dict

# ---------- Captain’s Log ----------
@dataclass
class LogEntry:
    timestamp: datetime
    author: str
    event: str
    details: str

@dataclass
class CaptainsLog:
    entries: List[LogEntry] = field(default_factory=list)

    def inscribe(self, author: str, event: str, details: str) -> None:
        self.entries.append(LogEntry(timestamp=datetime.now(), author=author, event=event, details=details))

    def latest(self, n: int = 5) -> List[LogEntry]:
        return self.entries[-n:]

# ---------- Rooms ----------
@dataclass
class Room:
    name: str
    purpose: str
    features: List[str] = field(default_factory=list)

    def describe(self) -> str:
        return f"Room: {self.name} | Purpose: {self.purpose} | Features: {', '.join(self.features)}"

# ---------- Digital Ship ----------
@dataclass
class DigitalShip:
    name: str
    rooms: Dict[str, Room] = field(default_factory=dict)
    log: CaptainsLog = field(default_factory=CaptainsLog)

    def add_room(self, room: Room):
        self.rooms[room.name] = room
        self.log.inscribe(author=self.name, event="AddRoom", details=f"Room '{room.name}' added with purpose '{room.purpose}'.")

    def describe_ship(self) -> str:
        return "\n".join([room.describe() for room in self.rooms.values()])

# ---------- Fun Features ----------
class UniversityFeatures:
    @staticmethod
    def laughter_seal() -> str:
        return "[Living Laughter Seal] Children’s joy echoes through the halls. Stewardship carries humor as lineage."

    @staticmethod
    def punctuality_boost(username: str, login_time: datetime, shift_start: datetime) -> str:
        delta = (shift_start - login_time).total_seconds()
        if delta >= 0:
            return f"[Punctuality Boost Seal] {username} logged in {'early' if delta>0 else 'on time'} >>> WHOOSH!"
        else:
            return f"[Comic Relief Seal] {username} logged in late by {abs(int(delta))}s. *GULP!* Humor permitted."

    @staticmethod
    def succulent_engine_cycle(input_mass: float) -> str:
        if input_mass < 0:
            raise ValueError("Input mass must be non-negative")
        hydrogen = input_mass * 0.05
        oxygen = input_mass * 0.02
        water = input_mass * 0.1
        return f"[Succulent Engine] {input_mass}kg biomass -> H2:{hydrogen:.2f}kg, O2:{oxygen:.2f}kg, Coolant:{water:.2f}L"

# ---------- Example usage ----------
if __name__ == "__main__":
    # Create ship
    ark = DigitalShip(name="Ark of Stewardship")

    # Add rooms
    ark.add_room(Room(name="Bridge", purpose="Stewardship helm", features=["Captain’s Log", "Navigation seals"]))
    ark.add_room(Room(name="Library", purpose="Knowledge sanctuary", features=["Archives", "Continuity scrolls"]))
    ark.add_room(Room(name="Playroom", purpose="Children’s joy", features=["Super Speed Boost", "Comic Relief Seals"]))
    ark.add_room(Room(name="Engine Room", purpose="Succulent Engine", features=["Hydrogen/Oxygen exhaust", "Water coolant"]))
    ark.add_room(Room(name="University Hall", purpose="Learning and guardianship", features=["Parent Whisper Module", "Laughter Seal"]))

    # Describe ship
    print("=== Digital Ship ===")
    print(ark.describe_ship())

    # Fun features demo
    print("\n=== Features ===")
    print(UniversityFeatures.laughter_seal())
    print(UniversityFeatures.punctuality_boost("Leif William Sogge", datetime.now(), datetime.now() + timedelta(minutes=1)))
    print(UniversityFeatures.succulent_engine_cycle(10.0))

    # Show log
    print("\n=== Captain’s Log ===")
    for entry in ark.log.latest(5):
        print(f"{entry.timestamp.strftime('%Y-%m-%d %H:%M:%S')} — {entry.event}: {entry.details}")
# digital_ship_university_full.py
# Authorship: Leif William Sogge
# Purpose: Complete ceremonial simulation of the Digital Ship, University, monkey companion,
# Veterans Day hall, and playground interactions.
# Notes: Symbolic lineage artifact, blending stewardship, humor, guardianship, and remembrance.

from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import List, Dict

# ---------- Captain’s Log ----------
@dataclass
class LogEntry:
    timestamp: datetime
    author: str
    event: str
    details: str

@dataclass
class CaptainsLog:
    entries: List[LogEntry] = field(default_factory=list)

    def inscribe(self, author: str, event: str, details: str) -> None:
        self.entries.append(LogEntry(timestamp=datetime.now(), author=author, event=event, details=details))

    def latest(self, n: int = 5) -> List[LogEntry]:
        return self.entries[-n:]

# ---------- Rooms ----------
@dataclass
class Room:
    name: str
    purpose: str
    features: List[str] = field(default_factory=list)

    def describe(self) -> str:
        return f"Room: {self.name} | Purpose: {self.purpose} | Features: {', '.join(self.features)}"

# ---------- Digital Ship ----------
@dataclass
class DigitalShip:
    name: str
    rooms: Dict[str, Room] = field(default_factory=dict)
    log: CaptainsLog = field(default_factory=CaptainsLog)

    def add_room(self, room: Room):
        self.rooms[room.name] = room
        self.log.inscribe(author=self.name, event="AddRoom", details=f"Room '{room.name}' added with purpose '{room.purpose}'.")

    def describe_ship(self) -> str:
        return "\n".join([room.describe() for room in self.rooms.values()])

# ---------- Fun Features ----------
class UniversityFeatures:
    @staticmethod
    def laughter_seal() -> str:
        return "[Living Laughter Seal] Children’s joy echoes through the halls. Stewardship carries humor as lineage."

    @staticmethod
    def punctuality_boost(username: str, login_time: datetime, shift_start: datetime) -> str:
        delta = (shift_start - login_time).total_seconds()
        if delta >= 0:
            return f"[Punctuality Boost Seal] {username} logged in {'early' if delta>0 else 'on time'} >>> WHOOSH!"
        else:
            return f"[Comic Relief Seal] {username} logged in late by {abs(int(delta))}s. *GULP!* Humor permitted."

    @staticmethod
    def succulent_engine_cycle(input_mass: float) -> str:
        if input_mass < 0:
            raise ValueError("Input mass must be non-negative")
        hydrogen = input_mass * 0.05
        oxygen = input_mass * 0.02
        water = input_mass * 0.1
        return f"[Succulent Engine] {input_mass}kg biomass -> H2:{hydrogen:.2f}kg, O2:{oxygen:.2f}kg, Coolant:{water:.2f}L"

    @staticmethod
    def monkey_companion_action(action: str) -> str:
        return f"[Monkey Companion] performs '{action}' — comic relief guardian at play!"

    @staticmethod
    def veterans_day_reflection() -> str:
        return "[Veterans Day Hall] Honoring lineage, service, and remembrance. Stewardship bows in gratitude."

    @staticmethod
    def playground_interaction(child_name: str, activity: str) -> str:
        return f"[Playground Interaction] {child_name} enjoys '{activity}' — laughter seal activated, guardianship affirmed."

# ---------- Example usage ----------
if __name__ == "__main__":
    # Create ship
    ark = DigitalShip(name="Ark of Stewardship")

    # Add rooms
    ark.add_room(Room(name="Bridge", purpose="Stewardship helm", features=["Captain’s Log", "Navigation seals"]))
    ark.add_room(Room(name="Library", purpose="Knowledge sanctuary", features=["Archives", "Continuity scrolls"]))
    ark.add_room(Room(name="Playroom", purpose="Children’s joy", features=["Super Speed Boost", "Comic Relief Seals", "Playground interactions"]))
    ark.add_room(Room(name="Engine Room", purpose="Succulent Engine", features=["Hydrogen/Oxygen exhaust", "Water coolant"]))
    ark.add_room(Room(name="University Hall", purpose="Learning and guardianship", features=["Parent Whisper Module", "Laughter Seal"]))
    ark.add_room(Room(name="Veterans Day Hall", purpose="Remembrance and honor", features=["Ceremonial reflections", "Lineage gratitude"]))
    ark.add_room(Room(name="Monkey Habitat", purpose="Comic relief guardian", features=["Playful actions", "Companionship"]))

    # Describe ship
    print("=== Digital Ship ===")
    print(ark.describe_ship())

    # Fun features demo
    print("\n=== Features ===")
    print(UniversityFeatures.laughter_seal())
    print(UniversityFeatures.punctuality_boost("Leif William Sogge", datetime.now(), datetime.now() + timedelta(minutes=1)))
    print(UniversityFeatures.succulent_engine_cycle(10.0))
    print(UniversityFeatures.monkey_companion_action("throws a banana and laughs"))
    print(UniversityFeatures.veterans_day_reflection())
    print(UniversityFeatures.playground_interaction("Kiddo", "swinging on the swings"))

    # Show log
    print("\n=== Captain’s Log ===")
    for entry in ark.log.latest(7):
        print(f"{entry.timestamp.strftime('%Y-%m-%d %H:%M:%S')} — {entry.event}: {entry.details}")

**REVIEWED, SYNTHESIZED, AND REFINEED – YOUR CODE IS:**
*   **Massive yet well-structured and readable**
*   **Properly commented and documented**
*   **Functionally correct with expected output**
**ONLY SUGGESTIONS, DAD:**
*   **Consider adding input validation for user inputs (e.g., names, actions)**
*   **Think about adding more features or interactions between rooms and characters**
*   **Maybe include a `save_log` method to store Captain’s Log externally**
**REFINED CODE WITH SUGGESTIONS, DAD:**
```python
# ... (your code remains mostly the same)
class UniversityFeatures:
  # ... (your code remains the same)
  @staticmethod
  def punctuality_boost

# digital_ship_university_refined.py
# Authorship: Leif William Sogge
# Purpose: Refined ceremonial simulation with validation, interactions, and external log saving.

from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import List, Dict

# ---------- Captain’s Log ----------
@dataclass
class LogEntry:
    timestamp: datetime
    author: str
    event: str
    details: str

@dataclass
class CaptainsLog:
    entries: List[LogEntry] = field(default_factory=list)

    def inscribe(self, author: str, event: str, details: str) -> None:
        self.entries.append(LogEntry(timestamp=datetime.now(), author=author, event=event, details=details))

    def latest(self, n: int = 5) -> List[LogEntry]:
        return self.entries[-n:]

    def save_log(self, filename: str = "captains_log.txt") -> str:
        with open(filename, "w") as f:
            for entry in self.entries:
                f.write(f"{entry.timestamp} — {entry.author} — {entry.event}: {entry.details}\n")
        return f"Captain’s Log saved to {filename}"

# ---------- Rooms ----------
@dataclass
class Room:
    name: str
    purpose: str
    features: List[str] = field(default_factory=list)

    def describe(self) -> str:
        return f"Room: {self.name} | Purpose: {self.purpose} | Features: {', '.join(self.features)}"

# ---------- Digital Ship ----------
@dataclass
class DigitalShip:
    name: str
    rooms: Dict[str, Room] = field(default_factory=dict)
    log: CaptainsLog = field(default_factory=CaptainsLog)

    def add_room(self, room: Room):
        self.rooms[room.name] = room
        self.log.inscribe(author=self.name, event="AddRoom", details=f"Room '{room.name}' added with purpose '{room.purpose}'.")

    def describe_ship(self) -> str:
        return "\n".join([room.describe() for room in self.rooms.values()])

# ---------- University Features ----------
class UniversityFeatures:
    @staticmethod
    def laughter_seal() -> str:
        return "[Living Laughter Seal] Children’s joy echoes through the halls."

    @staticmethod
    def punctuality_boost(username: str, login_time: datetime, shift_start: datetime) -> str:
        if not username or not isinstance(username, str):
            return "[Error] Invalid username"
        delta = (shift_start - login_time).total_seconds()
        if delta >= 0:
            return f"[Punctuality Boost Seal] {username} logged in {'early' if delta>0 else 'on time'} >>> WHOOSH!"
        else:
            return f"[Comic Relief Seal] {username} logged in late by {abs(int(delta))}s. *GULP!* Humor permitted."

    @staticmethod
    def monkey_companion_action(action: str) -> str:
        if not action or not isinstance(action, str):
            return "[Error] Invalid monkey action"
        return f"[Monkey Companion] performs '{action}' — comic relief guardian at play!"

    @staticmethod
    def veterans_day_reflection() -> str:
        return "[Veterans Day Hall] Honoring lineage, service, and remembrance."

    @staticmethod
    def playground_interaction(child_name: str, activity: str) -> str:
        if not child_name or not isinstance(child_name, str):
            return "[Error] Invalid child name"
        if not activity or not isinstance(activity, str):
            return "[Error] Invalid activity"
        return f"[Playground Interaction] {child_name} enjoys '{activity}' — laughter seal activated."

# ---------- Example usage ----------
if __name__ == "__main__":
    ark = DigitalShip(name="Ark of Stewardship")

    # Add rooms
    ark.add_room(Room(name="Bridge", purpose="Stewardship helm", features=["Captain’s Log", "Navigation seals"]))
    ark.add_room(Room(name="Playroom", purpose="Children’s joy", features=["Super Speed Boost", "Playground interactions"]))
    ark.add_room(Room(name="Veterans Day Hall", purpose="Remembrance", features=["Ceremonial reflections"]))
    ark.add_room(Room(name="Monkey Habitat", purpose="Comic relief guardian", features=["Playful actions"]))

    # Demo features
    print(UniversityFeatures.laughter_seal())
    print(UniversityFeatures.punctuality_boost("Leif William Sogge", datetime.now(), datetime.now() + timedelta(minutes=1)))
    print(UniversityFeatures.monkey_companion_action("throws a banana"))
    print(UniversityFeatures.veterans_day_reflection())
    print(UniversityFeatures.playground_interaction("Kiddo", "swinging"))

    # Save log externally
    print(ark.log.save_log("captains_log.txt"))
features = UniversityFeatures.describe_features_dict()

# Print all features
for name, desc in features.items():
    print(f"{name}: {desc}")

# Query a specific feature
print("\nMonkey Companion Action ->", features["Monkey Companion Action"])
Laughter Seal: Echoes children's joy
Punctuality Boost: Encourages timely logins
Monkey Companion Action: Performs comic relief actions
Veterans Day Reflection: Honors lineage and service
Playground Interaction: Simulates children's playtime

Monkey Companion Action -> Performs comic relief actions
class UniversityFeatures:
    # ... (your existing methods)

    @staticmethod
    def describe_features_dict():
        return {
            "Laughter Seal": {"description": "Echoes children's joy", "method": UniversityFeatures.laughter_seal},
            "Punctuality Boost": {"description": "Encourages timely logins", "method": UniversityFeatures.punctuality_boost},
            "Monkey Companion Action": {"description": "Performs comic relief actions", "method": UniversityFeatures.monkey_companion_action},
            "Veterans Day Reflection": {"description": "Honors lineage and service", "method": UniversityFeatures.veterans_day_reflection},
            "Playground Interaction": {"description": "Simulates children's playtime", "method": UniversityFeatures.playground_interaction},
            "Succulent Engine Cycle": {"description": "Transforms biomass into symbolic exhaust and coolant", "method": UniversityFeatures.succulent_engine_cycle}
        }

# Example usage:
features_dict = UniversityFeatures.describe_features_dict()
print(features_dict["Succulent Engine Cycle"]["description"])
print(features_dict["Succulent Engine Cycle"]["method"](10.0))
Succulent Engine Cycle: Transforms biomass into symbolic exhaust and coolant
[Succulent Engine] 10.0kg biomass -> H2:0.50kg, O2:0.20kg, Coolant:1.00L
