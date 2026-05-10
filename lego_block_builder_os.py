#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lego Block Builder OS — Kid-Safe Builder Prototype
Space LEAF Corp — Kid-Safe Testing Code

This is an adult-facing, console-only prototype that simulates:
- Block-based configuration (world, activity, safety, time, connection)
- Age bands (Explorer, Maker, Navigator)
- "Test as Kid" view
- "Stress Test" safety checks

NOTE:
This is a design + simulation tool, NOT a security boundary.
Real deployments must add OS-level and network-level controls.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Dict, Optional
import enum
import datetime


# ---------------------------------------------------------------------------
# Age bands and roles
# ---------------------------------------------------------------------------

class AgeBand(enum.Enum):
    EXPLORER = "Explorer (4–7)"
    MAKER = "Maker (8–11)"
    NAVIGATOR = "Navigator (12–15)"


def infer_age_band(age: int) -> AgeBand:
    if age <= 7:
        return AgeBand.EXPLORER
    elif age <= 11:
        return AgeBand.MAKER
    else:
        return AgeBand.NAVIGATOR


# ---------------------------------------------------------------------------
# Block model
# ---------------------------------------------------------------------------

@dataclass
class Block:
    """Base class for all blocks."""
    name: str
    description_kid: str

    def describe_for_kid(self) -> str:
        return f"{self.name}: {self.description_kid}"


@dataclass
class WorldBlock(Block):
    theme: str
    default_activities: List[str]


@dataclass
class ActivityBlock(Block):
    activity_type: str  # e.g. "Reading", "Drawing", "Game"
    difficulty: str     # e.g. "Easy", "Medium", "Hard"
    daily_minutes: int
    reward_minutes: int = 0
    requires_before: Optional[str] = None  # name of activity that must be done first


@dataclass
class SafetyBlock(Block):
    web_allowed: bool = False
    allowed_sites: List[str] = field(default_factory=list)
    chat_enabled: bool = False
    camera_enabled: bool = False
    store_requires_pin: bool = True
    in_app_purchases_allowed: bool = False

    def is_kid_safe(self, age_band: AgeBand) -> List[str]:
        """Return list of safety warnings (empty if OK)."""
        warnings = []

        # No unfiltered web
        if self.web_allowed and not self.allowed_sites:
            warnings.append("Web is allowed but no allowed-sites list is defined.")

        # Younger kids: no chat, no camera
        if age_band == AgeBand.EXPLORER:
            if self.chat_enabled:
                warnings.append("Chat is enabled for Explorer band (4–7).")
            if self.camera_enabled:
                warnings.append("Camera is enabled for Explorer band (4–7).")

        # In-app purchases must be off for all bands in this prototype
        if self.in_app_purchases_allowed:
            warnings.append("In-app purchases are enabled (should be OFF in this prototype).")

        # Store must require PIN
        if not self.store_requires_pin:
            warnings.append("Store does not require a parent PIN.")

        return warnings


@dataclass
class TimeMode(enum.Enum):
    SCHOOL = "School Mode"
    WEEKEND = "Weekend Mode"
    BEDTIME = "Bedtime Mode"


@dataclass
class TimeBlock(Block):
    mode: TimeMode
    allowed_activities: List[str]
    start_hour: int  # 24h
    end_hour: int    # 24h


@dataclass
class ConnectionBlock(Block):
    contact_name: str
    contact_role: str  # e.g. "Grandma", "Teacher"
    requires_parent_review: bool = True


# ---------------------------------------------------------------------------
# Kid profile and configuration
# ---------------------------------------------------------------------------

@dataclass
class KidProfile:
    name: str
    age: int
    interests: List[str]

    @property
    def age_band(self) -> AgeBand:
        return infer_age_band(self.age)


@dataclass
class KidConfig:
    """A full configuration for one kid."""
    profile: KidProfile
    world_block: WorldBlock
    activity_blocks: List[ActivityBlock]
    safety_block: SafetyBlock
    time_blocks: List[TimeBlock]
    connection_blocks: List[ConnectionBlock] = field(default_factory=list)

    # Simple audit log (for adults, not surveillance)
    audit_log: List[str] = field(default_factory=list)

    def log(self, message: str) -> None:
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.audit_log.append(f"[{timestamp}] {message}")

    # -----------------------------------------------------------------------
    # Simulation helpers
    # -----------------------------------------------------------------------

    def current_time_block(self, now: Optional[datetime.datetime] = None) -> Optional[TimeBlock]:
        if now is None:
            now = datetime.datetime.now()
        hour = now.hour
        for tb in self.time_blocks:
            if tb.start_hour <= hour < tb.end_hour:
                return tb
        return None

    def can_access_activity(
        self,
        activity_name: str,
        completed_activities: List[str],
        now: Optional[datetime.datetime] = None,
    ) -> (bool, str):
        """Return (allowed, reason)."""
        tb = self.current_time_block(now)
        if tb is None:
            return False, "No active time mode right now."

        if activity_name not in tb.allowed_activities:
            return False, f"'{activity_name}' is not allowed in {tb.mode.value}."

        act = next((a for a in self.activity_blocks if a.name == activity_name), None)
        if act is None:
            return False, "Activity not found in configuration."

        if act.requires_before and act.requires_before not in completed_activities:
            return False, f"You must do '{act.requires_before}' first."

        return True, "Allowed."

    def stress_test(self) -> List[str]:
        """Run safety checks and return list of issues."""
        issues: List[str] = []

        # Safety block checks
        issues.extend(self.safety_block.is_kid_safe(self.profile.age_band))

        # Check for any activity without a time budget
        for act in self.activity_blocks:
            if act.daily_minutes <= 0:
                issues.append(f"Activity '{act.name}' has no daily time limit.")

        # Check for any time mode that allows a "Game" without constraints
        for tb in self.time_blocks:
            for act_name in tb.allowed_activities:
                act = next((a for a in self.activity_blocks if a.name == act_name), None)
                if act and act.activity_type == "Game" and act.daily_minutes <= 0:
                    issues.append(
                        f"Time mode {tb.mode.value} allows game '{act.name}' with no daily limit."
                    )

        # Check for open browser possibility
        if self.safety_block.web_allowed and not self.safety_block.allowed_sites:
            issues.append("Kid could reach an unfiltered browser (web allowed, no site list).")

        # Check communication with strangers (simplified: any connection without parent review)
        for cb in self.connection_blocks:
            if not cb.requires_parent_review:
                issues.append(
                    f"Connection to '{cb.contact_name}' ({cb.contact_role}) "
                    f"does not require parent review."
                )

        return issues

    def kid_friendly_summary(self) -> str:
        """Explain the setup in kid language."""
        lines = []
        lines.append(f"Hi {self.profile.name}! This is your world:")
        lines.append(f"- World: {self.world_block.theme}")
        lines.append("- You can do these activities:")
        for act in self.activity_blocks:
            lines.append(
                f"  • {act.name} ({act.activity_type}) — about {act.daily_minutes} minutes per day."
            )
        lines.append("- Grown-ups added safety blocks so you can explore safely.")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Preset: space-loving 9-year-old
# ---------------------------------------------------------------------------

def build_space_loving_preset() -> KidConfig:
    profile = KidProfile(
        name="Alex",
        age=9,
        interests=["space", "drawing", "light games"],
    )

    world = WorldBlock(
        name="Space Station World",
        description_kid="Your own space station where you read, draw, and play space missions.",
        theme="Space Station",
        default_activities=["Space Stories", "Space Art Studio", "Orbital Puzzles"],
    )

    reading = ActivityBlock(
        name="Space Stories",
        description_kid="Read cool stories about planets, rockets, and astronauts.",
        activity_type="Reading",
        difficulty="Medium",
        daily_minutes=20,
    )

    drawing = ActivityBlock(
        name="Space Art Studio",
        description_kid="Draw planets, spaceships, and galaxies.",
        activity_type="Drawing",
        difficulty="Open",
        daily_minutes=30,
    )

    game = ActivityBlock(
        name="Orbital Puzzles",
        description_kid="Solve fun puzzles about orbits and gravity.",
        activity_type="Game",
        difficulty="Medium",
        daily_minutes=15,
        requires_before="Space Stories",
    )

    safety = SafetyBlock(
        name="Safe Space Shield",
        description_kid="This block keeps you safe while you explore your space station.",
        web_allowed=False,
        allowed_sites=[],
        chat_enabled=False,
        camera_enabled=False,
        store_requires_pin=True,
        in_app_purchases_allowed=False,
    )

    school_mode = TimeBlock(
        name="School Mode",
        description_kid="After school, you can read and draw in your space station.",
        mode=TimeMode.SCHOOL,
        allowed_activities=["Space Stories", "Space Art Studio"],
        start_hour=16,  # 4 PM
        end_hour=19,    # 7 PM
    )

    evening_mode = TimeBlock(
        name="Evening Game Unlock",
        description_kid="If you did your reading, you can play a space puzzle.",
        mode=TimeMode.WEEKEND,  # reuse enum for simplicity
        allowed_activities=["Orbital Puzzles", "Space Art Studio"],
        start_hour=19,  # 7 PM
        end_hour=20,    # 8 PM
    )

    connection = ConnectionBlock(
        name="Share Drawing with Grandma",
        description_kid="You can send one of your drawings to Grandma, but a grown-up checks first.",
        contact_name="Grandma",
        contact_role="Family",
        requires_parent_review=True,
    )

    config = KidConfig(
        profile=profile,
        world_block=world,
        activity_blocks=[reading, drawing, game],
        safety_block=safety,
        time_blocks=[school_mode, evening_mode],
        connection_blocks=[connection],
    )

    return config


# ---------------------------------------------------------------------------
# Console UI (adult-facing)
# ---------------------------------------------------------------------------

def print_header():
    print("=" * 70)
    print(" Lego Block Builder OS — Space LEAF Corp (Kid-Safe Testing Prototype)")
    print("=" * 70)
    print()


def menu() -> None:
    config = build_space_loving_preset()

    while True:
        print_header()
        print("Kid profile:")
        print(f"  Name: {config.profile.name}")
        print(f"  Age: {config.profile.age} ({config.profile.age_band.value})")
        print(f"  Interests: {', '.join(config.profile.interests)}")
        print()
        print("World:")
        print(f"  {config.world_block.theme}")
        print()
        print("Main menu:")
        print("  1) View kid-friendly summary")
        print("  2) Test as kid (simulate access)")
        print("  3) Run Stress Test (safety checks)")
        print("  4) View audit log")
        print("  5) Exit")
        choice = input("Select an option: ").strip()

        if choice == "1":
            print()
            print(config.kid_friendly_summary())
            input("\n[Enter] to return to menu...")
        elif choice == "2":
            simulate_as_kid(config)
        elif choice == "3":
            run_stress_test(config)
        elif choice == "4":
            show_audit_log(config)
        elif choice == "5":
            print("Goodbye.")
            break
        else:
            print("Invalid choice.")
            input("[Enter] to continue...")


def simulate_as_kid(config: KidConfig) -> None:
    print("\n--- Test as Kid ---")
    print("This simulates what the kid can access right now.")
    now = datetime.datetime.now()
    tb = config.current_time_block(now)
    if tb is None:
        print("No active time mode at this hour.")
        input("\n[Enter] to return...")
        return

    print(f"Current time: {now.strftime('%H:%M')} — Active mode: {tb.mode.value}")
    print("Allowed activities in this mode:")
    for name in tb.allowed_activities:
        print(f"  • {name}")
    print()

    completed: List[str] = []
    while True:
        print("Activities:")
        for idx, act in enumerate(config.activity_blocks, start=1):
            print(f"  {idx}) {act.name} ({act.activity_type})")
        print("  0) Return to main menu")
        choice = input("Choose an activity to attempt: ").strip()
        if choice == "0":
            break
        try:
            idx = int(choice) - 1
            act = config.activity_blocks[idx]
        except (ValueError, IndexError):
            print("Invalid choice.")
            continue

        allowed, reason = config.can_access_activity(act.name, completed, now)
        if allowed:
            print(f"✅ Access granted to '{act.name}'.")
            config.log(f"Simulated kid accessed '{act.name}'.")
            if act.name not in completed:
                completed.append(act.name)
        else:
            print(f"❌ Access denied to '{act.name}': {reason}")
            config.log(f"Simulated kid denied '{act.name}': {reason}")
        print()

    print("[End of kid simulation]")
    input("\n[Enter] to return to main menu...")


def run_stress_test(config: KidConfig) -> None:
    print("\n--- Stress Test ---")
    print("Running safety checks on current configuration...\n")
    issues = config.stress_test()
    if not issues:
        print("✅ No issues found. Configuration passes this prototype's safety checks.")
        config.log("Stress Test: no issues found.")
    else:
        print("⚠ Issues detected:")
        for i, issue in enumerate(issues, start=1):
            print(f"  {i}. {issue}")
        config.log(f"Stress Test: {len(issues)} issue(s) found.")
    input("\n[Enter] to return to main menu...")


def show_audit_log(config: KidConfig) -> None:
    print("\n--- Audit Log (Adult View) ---")
    if not config.audit_log:
        print("No events logged yet.")
    else:
        for line in config.audit_log:
            print(line)
    input("\n[Enter] to return to main menu...")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    menu()
