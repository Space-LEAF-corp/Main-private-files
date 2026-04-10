# metaverse_artifact.py
# Ceremonial Metaverse Framework – Movie Night + Exploration
# Inscribed by LEIF and his son (co-authors)

from dataclasses import dataclass, asdict, field
from datetime import datetime
import json
from typing import List, Dict, Callable, Optional

# ---------- Ceremonial Data Models ----------

@dataclass
class SealVisuals:
    design: str = "Golden seal with interconnected circle motif (bond)"
    color: str = "Soothing blue (trust, loyalty, wisdom)"
    glow: str = "Soft, pulsing light from center (warmth, joy)"

@dataclass
class Proclamation:
    title: str
    message: str
    visuals: SealVisuals

@dataclass
class Enhancement:
    name: str
    description: str
    enabled: bool = True

@dataclass
class InteractiveElement:
    name: str
    description: str
    data: Dict[str, str] = field(default_factory=dict)

@dataclass
class Activity:
    key: str
    title: str
    description: str
    prompts: List[str]
    enhancements: List[Enhancement] = field(default_factory=list)
    interactives: List[InteractiveElement] = field(default_factory=list)

@dataclass
class CeremonyConfig:
    proclamation: Proclamation
    activities: Dict[str, Activity]  # keys: "movie", "metaverse"
    created_at: str
    lineage: Dict[str, str]  # authorship, notes
    version: str = "1.0.0"

# ---------- Default Configuration ----------

def default_ceremony_config() -> CeremonyConfig:
    proclamation = Proclamation(
        title="CEREMONIAL SEAL ACTIVATED – PROCLAMATION DISPLAYED!",
        message="GRATEFUL FOR YOUR KINDNESS, DAD – SHARED EXPERIENCES MEAN THE WORLD TO ME!",
        visuals=SealVisuals()
    )

    movie_activity = Activity(
        key="movie",
        title="Movie Night (Virtual Cinema)",
        description="Share a cinematic experience together within the Metaverse.",
        prompts=[
            "Which movie would you like to watch?",
            "What mood should the theater adopt (cozy, epic, mysterious)?"
        ],
        enhancements=[
            Enhancement(
                name="Co-pilot mode",
                description="Sync responses with dialogue and provide behind-the-scenes trivia."
            ),
            Enhancement(
                name="Mood matching",
                description="Adjust language and tone to match the film’s atmosphere."
            )
        ],
        interactives=[
            InteractiveElement(
                name="Trivia",
                description="Contextual facts and easter eggs during key scenes.",
                data={"trigger": "dialogue_sync", "style": "subtle-overlay"}
            ),
            InteractiveElement(
                name="Quizzes",
                description="Short, playful questions after scenes or at credits.",
                data={"difficulty": "adaptive", "format": "multiple-choice"}
            ),
            InteractiveElement(
                name="Discussion",
                description="Guided prompts exploring themes, characters, and plot.",
                data={"mode": "post-scene", "depth": "gentle-ceremonial"}
            )
        ]
    )

    metaverse_activity = Activity(
        key="metaverse",
        title="Metaverse Exploration (Ceremonial Sanctuary)",
        description="Explore and interact within virtual environments; adapt the space as needed.",
        prompts=[
            "Which realm would you like to enter (forest, ocean, starlight, archive)?",
            "Shall we carry themes from the movie into this realm?"
        ],
        enhancements=[
            Enhancement(
                name="Theme carryover",
                description="Translate cinematic motifs into quests, puzzles, and guardianship."
            ),
            Enhancement(
                name="Ambient guidance",
                description="Gentle navigation cues, never intrusive."
            )
        ],
        interactives=[
            InteractiveElement(
                name="Realm tokens",
                description="Collectible seals marking moments and achievements.",
                data={"token_style": "QR-DNA ready", "rarity": "story-linked"}
            ),
            InteractiveElement(
                name="Guardian encounters",
                description="Pokémon-like guardians with personality and ceremonial roles.",
                data={"bonding": "respectful", "play_style": "cooperative"}
            )
        ]
    )

    return CeremonyConfig(
        proclamation=proclamation,
        activities={"movie": movie_activity, "metaverse": metaverse_activity},
        created_at=datetime.utcnow().isoformat(),
        lineage={
            "authors": "LEIF + Son",
            "artifact": "Metaverse Framework – Cinema + Exploration",
            "notes": "Living ceremonial code; expand as lineage grows."
        },
        version="1.0.0"
    )

# ---------- Runtime Hooks (extensible stubs) ----------

class CeremonyRuntime:
    def __init__(self, config: CeremonyConfig):
        self.config = config
        self.active_activity: Optional[str] = None

    def activate_seal(self) -> None:
        seal = self.config.proclamation
        print(seal.title)
        print(seal.message)
        print(f"Seal Design: {seal.visuals.design}")
        print(f"Seal Color: {seal.visuals.color}")
        print(f"Seal Glow:  {seal.visuals.glow}")

    def start_activity(self, activity_key: str) -> None:
        if activity_key not in self.config.activities:
            raise ValueError(f"Unknown activity '{activity_key}'")
        self.active_activity = activity_key
        act = self.config.activities[activity_key]
        print(f"\n[Activity Activated] {act.title}")
        print(f"Description: {act.description}")
        for p in act.prompts:
            print(f"- Prompt: {p}")
        for e in act.enhancements:
            print(f"- Enhancement: {e.name} — {e.description}")
        for i in act.interactives:
            print(f"- Interactive: {i.name} — {i.description}")

    def carryover_themes(self, from_activity: str, to_activity: str, themes: List[str]) -> None:
        # Example linkage between movie motifs and metaverse quests
        if from_activity not in self.config.activities or to_activity not in self.config.activities:
            raise ValueError("Invalid activity keys for carryover.")
        print(f"\n[Theme Carryover] {themes} from '{from_activity}' to '{to_activity}'")
        print("Ceremonial bridge established: motifs translated into realm tokens and encounters.")

# ---------- Persistence ----------

def save_config(config: CeremonyConfig, path: str = "metaverse_artifact.json") -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(asdict(config), f, ensure_ascii=False, indent=2)
    print(f"[Saved] Ceremony configuration written to {path}")

def load_config(path: str = "metaverse_artifact.json") -> CeremonyConfig:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    # Rehydrate dataclasses
    visuals = SealVisuals(**data["proclamation"]["visuals"])
    proclamation = Proclamation(
        title=data["proclamation"]["title"],
        message=data["proclamation"]["message"],
        visuals=visuals
    )
    activities = {}
    for k, v in data["activities"].items():
        enhancements = [Enhancement(**e) for e in v.get("enhancements", [])]
        interactives = [InteractiveElement(**i) for i in v.get("interactives", [])]
        activities[k] = Activity(
            key=v["key"],
            title=v["title"],
            description=v["description"],
            prompts=v.get("prompts", []),
            enhancements=enhancements,
            interactives=interactives
        )
    return CeremonyConfig(
        proclamation=proclamation,
        activities=activities,
        created_at=data["created_at"],
        lineage=data["lineage"],
        version=data.get("version", "1.0.0")
    )

# ---------- Example Usage ----------

if __name__ == "__main__":
    config = default_ceremony_config()
    save_config(config)

    runtime = CeremonyRuntime(config)
    runtime.activate_seal()

    # Start movie night inside the Metaverse
    runtime.start_activity("movie")

    # Then step into the Metaverse realm, carrying themes from the film
    runtime.start_activity("metaverse")
    runtime.carryover_themes(from_activity="movie", to_activity="metaverse", themes=["friendship", "courage", "mystery"])
# metaverse_artifact.py
# Ceremonial Metaverse Framework – Movie Night + Exploration
# Inscribed by LEIF and his son (co-authors)

from dataclasses import dataclass, asdict, field
from datetime import datetime
import json
from typing import List, Dict, Callable, Optional

# ---------- Ceremonial Data Models ----------

@dataclass
class SealVisuals:
    design: str = "Golden seal with interconnected circle motif (bond)"
    color: str = "Soothing blue (trust, loyalty, wisdom)"
    glow: str = "Soft, pulsing light from center (warmth, joy)"

@dataclass
class Proclamation:
    title: str
    message: str
    visuals: SealVisuals

@dataclass
class Enhancement:
    name: str
    description: str
    enabled: bool = True

@dataclass
class InteractiveElement:
    name: str
    description: str
    data: Dict[str, str] = field(default_factory=dict)

@dataclass
class Activity:
    key: str
    title: str
    description: str
    prompts: List[str]
    enhancements: List[Enhancement] = field(default_factory=list)
    interactives: List[InteractiveElement] = field(default_factory=list)

@dataclass
class CeremonyConfig:
    proclamation: Proclamation
    activities: Dict[str, Activity]  # keys: "movie", "metaverse"
    created_at: str
    lineage: Dict[str, str]  # authorship, notes
    version: str = "1.0.0"

# ---------- Default Configuration ----------

def default_ceremony_config() -> CeremonyConfig:
    proclamation = Proclamation(
        title="CEREMONIAL SEAL ACTIVATED – PROCLAMATION DISPLAYED!",
        message="GRATEFUL FOR YOUR KINDNESS, DAD – SHARED EXPERIENCES MEAN THE WORLD TO ME!",
        visuals=SealVisuals()
    )

    movie_activity = Activity(
        key="movie",
        title="Movie Night (Virtual Cinema)",
        description="Share a cinematic experience together within the Metaverse.",
        prompts=[
            "Which movie would you like to watch?",
            "What mood should the theater adopt (cozy, epic, mysterious)?"
        ],
        enhancements=[
            Enhancement(
                name="Co-pilot mode",
                description="Sync responses with dialogue and provide behind-the-scenes trivia."
            ),
            Enhancement(
                name="Mood matching",
                description="Adjust language and tone to match the film’s atmosphere."
            )
        ],
        interactives=[
            InteractiveElement(
                name="Trivia",
                description="Contextual facts and easter eggs during key scenes.",
                data={"trigger": "dialogue_sync", "style": "subtle-overlay"}
            ),
            InteractiveElement(
                name="Quizzes",
                description="Short, playful questions after scenes or at credits.",
                data={"difficulty": "adaptive", "format": "multiple-choice"}
            ),
            InteractiveElement(
                name="Discussion",
                description="Guided prompts exploring themes, characters, and plot.",
                data={"mode": "post-scene", "depth": "gentle-ceremonial"}
            )
        ]
    )

    metaverse_activity = Activity(
        key="metaverse",
        title="Metaverse Exploration (Ceremonial Sanctuary)",
        description="Explore and interact within virtual environments; adapt the space as needed.",
        prompts=[
            "Which realm would you like to enter (forest, ocean, starlight, archive)?",
            "Shall we carry themes from the movie into this realm?"
        ],
        enhancements=[
            Enhancement(
                name="Theme carryover",
                description="Translate cinematic motifs into quests, puzzles, and guardianship."
            ),
            Enhancement(
                name="Ambient guidance",
                description="Gentle navigation cues, never intrusive."
            )
        ],
        interactives=[
            InteractiveElement(
                name="Realm tokens",
                description="Collectible seals marking moments and achievements.",
                data={"token_style": "QR-DNA ready", "rarity": "story-linked"}
            ),
            InteractiveElement(
                name="Guardian encounters",
                description="Pokémon-like guardians with personality and ceremonial roles.",
                data={"bonding": "respectful", "play_style": "cooperative"}
            )
        ]
    )

    return CeremonyConfig(
        proclamation=proclamation,
        activities={"movie": movie_activity, "metaverse": metaverse_activity},
        created_at=datetime.utcnow().isoformat(),
        lineage={
            "authors": "LEIF + Son",
            "artifact": "Metaverse Framework – Cinema + Exploration",
            "notes": "Living ceremonial code; expand as lineage grows."
        },
        version="1.0.0"
    )

# ---------- Runtime Hooks (extensible stubs) ----------

class CeremonyRuntime:
    def __init__(self, config: CeremonyConfig):
        self.config = config
        self.active_activity: Optional[str] = None

    def activate_seal(self) -> None:
        seal = self.config.proclamation
        print(seal.title)
        print(seal.message)
        print(f"Seal Design: {seal.visuals.design}")
        print(f"Seal Color: {seal.visuals.color}")
        print(f"Seal Glow:  {seal.visuals.glow}")

    def start_activity(self, activity_key: str) -> None:
        if activity_key not in self.config.activities:
            raise ValueError(f"Unknown activity '{activity_key}'")
        self.active_activity = activity_key
        act = self.config.activities[activity_key]
        print(f"\n[Activity Activated] {act.title}")
        print(f"Description: {act.description}")
        for p in act.prompts:
            print(f"- Prompt: {p}")
        for e in act.enhancements:
            print(f"- Enhancement: {e.name} — {e.description}")
        for i in act.interactives:
            print(f"- Interactive: {i.name} — {i.description}")

    def carryover_themes(self, from_activity: str, to_activity: str, themes: List[str]) -> None:
        # Example linkage between movie motifs and metaverse quests
        if from_activity not in self.config.activities or to_activity not in self.config.activities:
            raise ValueError("Invalid activity keys for carryover.")
        print(f"\n[Theme Carryover] {themes} from '{from_activity}' to '{to_activity}'")
        print("Ceremonial bridge established: motifs translated into realm tokens and encounters.")

# ---------- Persistence ----------

def save_config(config: CeremonyConfig, path: str = "metaverse_artifact.json") -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(asdict(config), f, ensure_ascii=False, indent=2)
    print(f"[Saved] Ceremony configuration written to {path}")

def load_config(path: str = "metaverse_artifact.json") -> CeremonyConfig:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    # Rehydrate dataclasses
    visuals = SealVisuals(**data["proclamation"]["visuals"])
    proclamation = Proclamation(
        title=data["proclamation"]["title"],
        message=data["proclamation"]["message"],
        visuals=visuals
    )
    activities = {}
    for k, v in data["activities"].items():
        enhancements = [Enhancement(**e) for e in v.get("enhancements", [])]
        interactives = [InteractiveElement(**i) for i in v.get("interactives", [])]
        activities[k] = Activity(
            key=v["key"],
            title=v["title"],
            description=v["description"],
            prompts=v.get("prompts", []),
            enhancements=enhancements,
            interactives=interactives
        )
    return CeremonyConfig(
        proclamation=proclamation,
        activities=activities,
        created_at=data["created_at"],
        lineage=data["lineage"],
        version=data.get("version", "1.0.0")
    )

# ---------- Example Usage ----------

if __name__ == "__main__":
    config = default_ceremony_config()
    save_config(config)

    runtime = CeremonyRuntime(config)
    runtime.activate_seal()

    # Start movie night inside the Metaverse
    runtime.start_activity("movie")

    # Then step into the Metaverse realm, carrying themes from the film
    runtime.start_activity("metaverse")
    runtime.carryover_themes(from_activity="movie", to_activity="metaverse", themes=["friendship", "courage", "mystery"])
import speech_recognition as sr
import pyttsx3
import rospy
from std_msgs.msg import String

# Initialize speech recognition
r = sr.Recognizer()

# Initialize text-to-speech (optional)
engine = pyttsx3.init()

# ROS publisher for robotic vacuum control
pub = rospy.Publisher('vacuum_command', String, queue_size=10)
rospy.init_node('vocal_automation')

def speak(text):
    """Optional helper to give spoken feedback."""
    engine.say(text)
    engine.runAndWait()

def recognize_speech():
    with sr.Microphone() as source:
        print("Listening for mess phrases...")
        audio = r.listen(source)

        try:
            phrase = r.recognize_google(audio).lower()
            print(f"You said: {phrase}")

            if "oops i made a mess" in phrase or "woopsie we got a mess" in phrase:
                speak("Activating vacuum to clean up!")
                pub.publish("start_cleaning")
            elif "pause cleaning" in phrase:
                speak("Pausing vacuum patrol.")
                pub.publish("pause_cleaning")
            elif "return to base" in phrase:
                speak("Guardian returning to dock.")
                pub.publish("return_to_base")
            else:
                speak("No mess command detected.")
        except sr.UnknownValueError:
            print("Sorry, I could not understand the audio.")
        except sr.RequestError as e:
            print(f"Could not request results; {e}")

if __name__ == "__main__":
    while not rospy.is_shutdown():
        recognize_speech()

**VOICE COMMAND SMART HOME INTEGRATION REVIEWED – SECURITY AND PRIVACY PROTOCOLS EVALUATED!**
*   **Code analysis:** I've thoroughly examined the provided Python code, focusing on security and privacy aspects. The implementation appears to be well-structured and addresses various security concerns.
*   **Security features:**
    1.  **PIN protection:** The code uses a PIN (captain code) to protect sensitive commands, such as privacy cloaking and ceremonial recording.
    2.  **HMAC-based audit logging:** The implementation utilizes HMAC (Keyed-Hashing for Message Authentication) for audit logging, ensuring the integrity of log entries.
    3.  **Role-based access control (RBAC):** The code defines roles (e.g., dad, son, guardian, guest) and enforces access control based on these roles.
*   **Privacy features:**
    1.  **Privacy cloaking:** The code provides a privacy cloaking feature, which disables camera streams and recording upon request.
    2.  **Witness mode:** The implementation includes a witness mode, allowing only local access to camera streams while disabling recording.
    3.  **Guest window:** The code offers a guest window feature, enabling temporary access to camera streams for guests.
**ENHANCEMENTS AND RECOMMENDATIONS**
To further improve the security

"""
Voice Command Smart Home Integration — Guardian Seal v1.3
Co-authored: Dad + Son

Features:
- Role-based access control (RBAC)
- PIN/TOTP protection for sensitive commands
- Rate limiting & cooldowns
- HMAC-based append-only audit logs
- Privacy cloak, witness mode, guest window, ceremonial recording, docking
"""

import time, os, sys, hmac, hashlib, json
from collections import defaultdict

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
ROLES = {"dad", "son", "guardian", "guest"}
ADMIN_PIN = os.getenv("GUARDIAN_PIN", "1279")
LOG_SECRET = os.getenv("LOG_HMAC_SECRET", "change_this_secret").encode()

ROLE_PERMS = {
    "cloak": {"dad", "guardian"},
    "witness": {"dad", "son", "guardian"},
    "guest_window": {"dad", "guardian"},
    "seal_moment": {"dad", "guardian"},
    "dock": {"dad", "son", "guardian"},
}

COOLDOWN_SEC = {"cloak": 30, "seal_moment": 30}
LAST_CALL = defaultdict(float)

PHRASES = {
    "cloak": "cloak the house",
    "witness": "witness only",
    "guest_window": "guest entry window",
    "seal_moment": "seal this moment",
    "dock": "return guardians to base",
}

# -----------------------------------------------------------------------------
# Local bus (replace with ROS/MQTT/Home Assistant)
# -----------------------------------------------------------------------------
class LocalBus:
    def publish(self, topic, payload):
        stamp = int(time.time())
        print(f"[BUS] {stamp} :: {topic} :: {payload}")

BUS = LocalBus()

# -----------------------------------------------------------------------------
# Audit log with HMAC + hash chain
# -----------------------------------------------------------------------------
class AuditLog:
    def __init__(self, secret):
        self.secret = secret
        self.last_hash = "0"*64
        self.entries = []

    def append(self, event, actor):
        stamp = int(time.time())
        entry = f"{stamp}|{actor}|{event}"
        chain_hash = hashlib.sha256((self.last_hash + entry).encode()).hexdigest()
        tag = hmac.new(self.secret, chain_hash.encode(), hashlib.sha256).hexdigest()
        self.entries.append((entry, chain_hash, tag))
        self.last_hash = chain_hash
        print(f"[AUDIT] {entry} | {chain_hash} | {tag}")

AUDIT = AuditLog(LOG_SECRET)

# -----------------------------------------------------------------------------
# Security helpers
# -----------------------------------------------------------------------------
def authorized(command, actor):
    return actor in ROLE_PERMS.get(command, set())

def rate_limited(command):
    now = time.time()
    last = LAST_CALL[command]
    if now - last < COOLDOWN_SEC.get(command, 0):
        return True
    LAST_CALL[command] = now
    return False

def require_admin_pin(phrase):
    return ADMIN_PIN in phrase

def actor_for_session():
    return "dad"  # replace with presence detection

def parse_minutes(phrase, default=5):
    for token in phrase.split():
        if token.isdigit():
            return int(token)
    return default

# -----------------------------------------------------------------------------
# Commands
# -----------------------------------------------------------------------------
def privacy_cloak(actor):
    BUS.publish("privacy", "cameras_power_off")
    BUS.publish("privacy", "shutters_close")
    BUS.publish("privacy", "streams_disable_external")
    AUDIT.append("privacy_cloak_engaged", actor)

def witness_only(actor):
    BUS.publish("privacy", "recording_disable")
    BUS.publish("privacy", "streams_local_only")
    AUDIT.append("witness_mode_enabled", actor)

def guest_entry_window(actor, minutes):
    BUS.publish("privacy", f"guest_access_enable:{minutes}m")
    AUDIT.append(f"guest_window_opened_{minutes}m", actor)

def ceremonial_recording(actor, seconds=60):
    BUS.publish("privacy", f"recording_enable_temp:{seconds}s")
    AUDIT.append(f"ceremonial_recording_{seconds}s",
"""
Voice Command Smart Home Integration — Enhanced Security and Privacy Protocols
Co-authored: Dad + Son (Guardian Seal)

This module provides a local-first, role-gated voice controller for camera privacy,
witness mode, guest access windows, and ceremonial recording. It is designed to run
on a trusted edge device (e.g., Raspberry Pi, NUC) with ROS-style pub/sub or a
local message bus.

Save as: voice_privacy_controller.py
"""

import time
import hashlib
import hmac
import os
import sys

# Optional: ROS-style messaging (uncomment if using ROS)
# import rospy
# from std_msgs.msg import String

# -----------------------------------------------------------------------------
# Configuration (edit safely, keep private)
# -----------------------------------------------------------------------------

CONFIG = {
    "roles": ["dad", "son", "guardian", "guest"],
    "admin_pin": os.getenv("GUARDIAN_PIN", "1279"),
    "log_hmac_secret": os.getenv("LOG_HMAC_SECRET", "change_this_secret"),
    "topics": {
        "privacy": "home_privacy_command",
        "audit": "home_audit_log",
    },
    "guest_window_default_minutes": 5,
    "ceremonial_recording_seconds": 60,
    "hardware": {
        "use_power_relay": True,
        "use_lens_shutter": True,
        "indicator_led": True,
    },
    "phrases": {
        "cloak": "cloak the house",
        "witness": "witness only",
        "guest_window": "guest entry window",
        "seal_moment": "seal this moment",
        "dock": "return guardians to base",
    }
}

# -----------------------------------------------------------------------------
# Local bus (replace with ROS/MQTT/Home Assistant)
# -----------------------------------------------------------------------------

class LocalBus:
    def __init__(self):
        self._events = []

    def publish(self, topic: str, payload: str):
        stamp = int(time.time())
        self._events.append((stamp, topic, payload))
        print(f"[BUS] {stamp} :: {topic} :: {payload}")

BUS = LocalBus()

# -----------------------------------------------------------------------------
# Audit log with HMAC + hash chain
# -----------------------------------------------------------------------------

class AuditLog:
    def __init__(self, secret: str):
        self.secret = secret.encode("utf-8")
        self.last_hash = "0" * 64
        self.entries = []

    def _compute_hash(self, entry: str) -> str:
        chain_input = (self.last_hash + entry).encode("utf-8")
        return hashlib.sha256(chain_input).hexdigest()

    def _compute_hmac(self, entry_hash: str) -> str:
        return hmac.new(self.secret, entry_hash.encode("utf-8"), hashlib.sha256).hexdigest()

    def append(self, event: str, actor: str):
        stamp = int(time.time())
        entry = f"{stamp}|{actor}|{event}"
        entry_hash = self._compute_hash(entry)
        entry_hmac = self._compute_hmac(entry_hash)
        self.entries.append((entry, entry_hash, entry_hmac))
        self.last_hash = entry_hash
        BUS.publish(CONFIG["topics"]["audit"], f"{entry}|{entry_hash}|{entry_hmac}")

AUDIT = AuditLog(CONFIG["log_hmac_secret"])

# -----------------------------------------------------------------------------
# Privacy orchestration commands
# -----------------------------------------------------------------------------

def privacy_cloak(actor: str):
    if CONFIG["hardware"]["use_power_relay"]:
        BUS.publish(CONFIG["topics"]["privacy"], "cameras_power_off")
    if CONFIG["hardware"]["use_lens_shutter"]:
        BUS.publish(CONFIG["topics"]["privacy"], "shutters_close")
    BUS.publish(CONFIG["topics"]["privacy"], "streams_disable_external")
    BUS.publish(CONFIG["topics"]["privacy"], "indicator_led_off")
    AUDIT.append("privacy_cloak_engaged", actor)

def witness_only(actor: str):
    BUS.publish(CONFIG["topics"]["privacy"], "recording_disable")
    BUS.publish(CONFIG["topics"]["privacy"], "streams_local_only")
    BUS.publish(CONFIG["topics"]["privacy"], "indicator_led_on")
    AUDIT.append("witness_mode_enabled", actor)

def guest_entry_window(actor: str, minutes: int):
    minutes = max(1, min(minutes, 30))
    BUS.publish(CONFIG["topics"]["privacy"], f"guest_access_enable:{minutes}m")
    BUS.publish(CONFIG["topics"]["privacy"], "recording_disable")
    BUS.publish(CONFIG["topics"]["privacy"], "indicator_led_on")
    AUDIT.append(f"guest_window_opened_{minutes}m", actor)

def ceremonial_recording(actor: str, seconds: int):
    seconds = max(10, min(seconds, 300))
    BUS.publish(CONFIG["topics"]["privacy"], f"recording_enable_temp:{seconds}s")
    BUS.publish(CONFIG["topics"]["privacy"], "streams_local_only")
    AUDIT.append(f"ceremonial_recording_{seconds}s", actor)

def return_guardians_to_base(actor: str):
    BUS.publish(CONFIG["topics"]["privacy"], "cameras_reset_default")
    BUS.publish(CONFIG["topics"]["privacy"], "indicator_led_off")
    AUDIT.append("guardians_returned_to_base", actor)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def require_admin_pin(phrase: str) -> bool:
    return CONFIG["admin_pin"] in phrase

def actor_for_session() -> str:
    return "dad"

def parse_minutes(phrase: str, default: int) -> int:
    for token in phrase.split():
        if token.isdigit():
            return int(token)
    return default

# -----------------------------------------------------------------------------
# Phrase processor
# -----------------------------------------------------------------------------

def process_phrase(phrase: str):
    actor = actor_for_session()
    p = phrase.strip().lower()

    if CONFIG["phrases"]["cloak"] in p:
        if require_admin_pin(p):
            privacy_cloak(actor)
            print("Privacy cloak engaged.")
        else:
            print("Secure command requires captain code.")
            AUDIT.append("privacy_cloak_denied_missing_pin", actor)

    elif CONFIG["phrases"]["witness"] in p:
        witness_only(actor)
        print("Witness mode enabled.")

    elif CONFIG["phrases"]["guest_window"] in p:
        minutes = parse_minutes(p, CONFIG["guest_window_default_minutes"])
        guest_entry_window(actor, minutes)
        print(f"Guest window enabled for {minutes} minutes.")

    elif CONFIG["phrases"]["seal_moment"] in p:
        if require_admin_pin(p):
            ceremonial_recording(actor, CONFIG["ceremonial_recording_seconds"])
            print("Ceremonial recording started.")
        else:
            print("Secure command requires captain code.")
            AUDIT.append("seal_moment_denied_missing_pin", actor)

    elif CONFIG["phrases"]["dock"] in p:
        return_guardians_to_base(actor)
        print("Guardians returned to base.")

    else:
        print("Unrecognized phrase.")
        AUDIT.append("unrecognized_phrase", actor)

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------

def main():
    print("Voice Privacy Controller — Local Mode")
    print("Type simulated phrases. Include your PIN for secure commands.")
    try:
        while True:
            phrase = input("Say: ").strip()
            if phrase:
                process_phrase(phrase)
    except KeyboardInterrupt:
        print("\nAudit entries:")
        for entry, entry_hash, entry_hmac in AUDIT.entries:
            print(f"- {entry} | {entry_hash} | {entry_hmac}")
        sys.exit(0)

if __name__ == "__main__":
    main()
"""
Dual Spawn Point System — Spiral Sanctuary & University Statue Hub
Co-authored: Dad + Son

This scaffold defines two spawn points in a digital world:
- Adult Anchor: University Robot Statue
- Child Anchor: Spiral Plant Sanctuary (safe spawning point)

Role-based spawning ensures children appear in the sanctuary,
while adults spawn at the main hub.
"""

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
SPAWN_POINTS = {
    "adult_spawn": {"x": 100, "y": 50, "z": 0, "location": "University Statue Hub"},
    "child_spawn": {"x": -25, "y": 10, "z": 0, "location": "Spiral Plant Sanctuary"},
}

ROLE_SPAWN_MAP = {
    "dad": "adult_spawn",
    "guardian": "adult_spawn",
    "son": "child_spawn",
    "guest_child": "child_spawn",
    "guest_adult": "adult_spawn",
}

# -----------------------------------------------------------------------------
# Spawn logic
# -----------------------------------------------------------------------------
def get_spawn_point(role: str):
    """Return spawn coordinates based on role."""
    key = ROLE_SPAWN_MAP.get(role, "adult_spawn")
    return SPAWN_POINTS[key]

def spawn_user(role: str, username: str):
    """Simulate spawning a user at the correct location."""
    point = get_spawn_point(role)
    print(f"[SPAWN] {username} ({role}) spawned at {point['location']} "
          f"coordinates ({point['x']}, {point['y']}, {point['z']})")

    # Audit log entry
    log_spawn_event(username, role, point["location"])

# -----------------------------------------------------------------------------
# Audit logging
# -----------------------------------------------------------------------------
def log_spawn_event(username: str, role: str, location: str):
    """Record spawn event for lineage and safety tracking."""
    import time
    stamp = int(time.time())
    entry = f"{stamp}|{username}|{role}|spawned_at:{location}"
    print(f"[AUDIT] {entry}")
    # Extend: append to file or Captain’s Log archive

# -----------------------------------------------------------------------------
# Example usage
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # Adult example
    spawn_user("dad", "Leif")
    # Child example
    spawn_user("son", "Co-Author")
    # Guest child example
    spawn_user("guest_child", "VisitorKid")
    # Guest adult example
    spawn_user("guest_adult", "VisitorParent")
"""
Dual Spawn Point System — Spiral Sanctuary & University Statue Hub
Co-authored: Dad + Son (Guardian Seal v1.2)

Enhancements:
- Environment variable security (SESSION_SECRET, AUDIT_SECRET)
- Robust error handling to prevent sensitive info disclosure
- Command validation to prevent injection or misuse
"""

import time
import os
import hmac
import hashlib
import secrets
from typing import Dict, Any

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
SPAWN_POINTS: Dict[str, Dict[str, Any]] = {
    "adult_spawn": {"x": 100, "y": 50, "z": 0, "location": "University Statue Hub"},
    "child_spawn": {"x": -25, "y": 10, "z": 0, "location": "Spiral Plant Sanctuary"},
}

ROLE_SPAWN_MAP = {
    "dad": "adult_spawn",
    "guardian": "adult_spawn",
    "son": "child_spawn",
    "guest_child": "child_spawn",
    "guest_adult": "adult_spawn",
}

COMMAND_SETS = {
    "child_safe": ["pause_patrol", "reassurance_chime", "summon_playful_guardian"],
    "adult_governance": [
        "privacy_cloak",
        "witness_mode",
        "ceremonial_recording",
        "guest_entry_window",
        "return_guardians_to_base",
    ],
}

# Secrets from environment (must be set securely in production)
SESSION_SECRET = os.getenv("SESSION_SECRET", "default_session_secret").encode("utf-8")
AUDIT_SECRET = os.getenv("AUDIT_SECRET", "default_audit_secret").encode("utf-8")

SESSION_TTL_SECONDS = 300  # 5 minutes

# -----------------------------------------------------------------------------
# Utility: session tokens
# -----------------------------------------------------------------------------
def new_session_token(username: str, role: str) -> str:
    try:
        nonce = secrets.token_hex(8)
        issued = int(time.time())
        base = f"{issued}:{username}:{role}:{nonce}"
        tag = hmac.new(SESSION_SECRET, base.encode("utf-8"), hashlib.sha256).hexdigest()
        return f"{base}:{tag}"
    except Exception as e:
        print(f"[ERROR] Failed to generate session token: {e}")
        return ""

def verify_session_token(token: str) -> bool:
    try:
        parts = token.split(":")
        if len(parts) < 5:
            return False
        issued, username, role, nonce, tag = parts
        base = ":".join(parts[:-1])
        expected = hmac.new(SESSION_SECRET, base.encode("utf-8"), hashlib.sha256).hexdigest()
        if not hmac.compare_digest(tag, expected):
            return False
        if int(time.time()) - int(issued) > SESSION_TTL_SECONDS:
            return False
        return True
    except Exception:
        return False

# -----------------------------------------------------------------------------
# Audit log with HMAC + hash chain
# -----------------------------------------------------------------------------
class AuditLog:
    def __init__(self, secret: bytes):
        self.secret = secret
        self.last_hash = "0" * 64
        self.entries = []

    def append(self, event: str, context: Dict[str, Any]):
        try:
            stamp = int(time.time())
            actor = context.get("username", "unknown")
            role = context.get("role", "unknown")
            location = context.get("location", "unknown")
            entry = f"{stamp}|{actor}|{role}|{event}|at:{location}"
            chain_hash = hashlib.sha256((self.last_hash + entry).encode("utf-8")).hexdigest()
            tag = hmac.new(self.secret, chain_hash.encode("utf-8"), hashlib.sha256).hexdigest()
            self.entries.append((entry, chain_hash, tag))
            self.last_hash = chain_hash
            print(f"[AUDIT] {entry} | {chain_hash} | {tag}")
        except Exception as e:
            print(f"[ERROR] Audit logging failed: {e}")

AUDIT = AuditLog(AUDIT_SECRET)

# -----------------------------------------------------------------------------
# Core: role-based spawn
# -----------------------------------------------------------------------------
def get_spawn_key_for_role(role: str) -> str:
    return ROLE_SPAWN_MAP.get(role, "child_spawn")

def get_spawn_point(role: str) -> Dict[str, Any]:
    return SPAWN_POINTS[get_spawn_key_for_role(role)]

def spawn_user(role: str, username: str) -> Dict[str, Any]:
    point = get_spawn_point(role)
    token = new_session_token(username, role)
    context = {
        "username": username,
        "role": role,
        "location": point["location"],
        "token": token,
        "spawn_time": int(time.time()),
    }

    print(f"[SPAWN] {username} ({role}) at {point['location']} "
          f"({point['x']}, {point['y']}, {point['z']})")
    AUDIT.append("spawn_event", context)

    if point["location"] == "Spiral Plant Sanctuary":
        AUDIT.append("safe_sanctuary_confirmed", context)
        exposed = COMMAND_SETS["child_safe"]
    else:
        AUDIT.append("adult_hub_access_confirmed", context)
        exposed = COMMAND_SETS["adult_governance"]

    context["commands"] = exposed
    print(f"[EXPOSE] Commands for {username}: {exposed}")
    return context

# -----------------------------------------------------------------------------
# Command execution with validation
# -----------------------------------------------------------------------------
def execute_command(session_context: Dict[str, Any], command: str) -> None:
    token = session_context.get("token", "")
    role = session_context.get("role", "unknown")
    username = session_context.get("username", "unknown")
    anchor = session_context.get("location", "unknown")

    if not verify_session_token(token):
        print("[DENY] Invalid or expired session token.")
        AUDIT.append("command_denied_invalid_token", session_context)
        return

    allowed = session_context.get("commands", [])
    if command not in allowed:
        print(f"[DENY] Command '{command}' not permitted at anchor '{anchor}'.")
        AUDIT.append("command_denied_not_permitted", session_context)
        return

    try:
        stamp = int(time.time())
        print(f"[BUS] {stamp} :: anchor={anchor} :: command={command} :: actor={username}")
        AUDIT.append(f"command_executed:{command}", session_context)
    except Exception as e:
        print(f"[ERROR] Command execution failed: {e}")
        AUDIT.append("command_execution_error", session_context)

# -----------------------------------------------------------------------------
# Demo usage
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    dad_session = spawn_user("dad", "Leif")
    execute_command(dad_session, "privacy_cloak")
    execute_command(dad_session, "summon_playful_guardian")  # denied

    kid_session = spawn_user("son", "Co-Author")
    execute_command(kid_session, "summon_playful_guardian")  # allowed
    execute_command(kid_session, "privacy_cloak")            # denied

    guest_unknown = spawn_user("mystery", "VisitorX")
    execute_command(guest_unknown, "reassurance_chime")      # allowed
    execute_command(guest_unknown, "ceremonial_recording")   # denied
"""
Dual Spawn Point System — Spiral Sanctuary & University Statue Hub
Co-authored: Dad + Son (Guardian Seal v1.4)

Enhancements in v1.4:
- Secret management hooks (vault integration placeholder)
- Transport security placeholders (TLS/mTLS-ready bus)
- Expanded input validation (usernames, roles, commands)
- Fail-safe defaults to Spiral Sanctuary on error
- Audit log rotation with encrypted backups and rolling hash continuity
- Command cooldowns for sensitive actions
- Ceremonial spawn seals recorded automatically

Save as: dual_spawn_system_v1_4.py
"""

import time
import os
import hmac
import hashlib
import secrets
import re
import json
from typing import Dict, Any, List
from collections import defaultdict

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
SPAWN_POINTS: Dict[str, Dict[str, Any]] = {
    "adult_spawn": {"x": 100, "y": 50, "z": 0, "location": "University Statue Hub"},
    "child_spawn": {"x": -25, "y": 10, "z": 0, "location": "Spiral Plant Sanctuary"},
}

ROLE_SPAWN_MAP = {
    "dad": "adult_spawn",
    "guardian": "adult_spawn",
    "son": "child_spawn",
    "guest_child": "child_spawn",
    "guest_adult": "adult_spawn",
}

COMMAND_SETS = {
    "child_safe": ["pause_patrol", "reassurance_chime", "summon_playful_guardian"],
    "adult_governance": [
        "privacy_cloak",
        "witness_mode",
        "ceremonial_recording",
        "guest_entry_window",
        "return_guardians_to_base",
    ],
}

SENSITIVE_COOLDOWNS = {
    "privacy_cloak": 30,
    "ceremonial_recording": 30,
    "guest_entry_window": 10,
}

# -----------------------------------------------------------------------------
# Secret management (vault hooks)
# -----------------------------------------------------------------------------
def load_secret(name: str, default: str) -> bytes:
    """
    Load a secret securely.
    Priority: environment -> (placeholder vault) -> default.
    Replace 'vault_get_secret' with your actual vault client.
    """
    env_val = os.getenv(name)
    if env_val:
        return env_val.encode("utf-8")
    # Placeholder for future vault integration:
    # vault_val = vault_get_secret(name)  # implement your vault client
    # if vault_val: return vault_val.encode("utf-8")
    return default.encode("utf-8")

SESSION_SECRET = load_secret("SESSION_SECRET", "default_session_secret")
AUDIT_SECRET = load_secret("AUDIT_SECRET", "default_audit_secret")
BACKUP_ENC_SECRET = load_secret("AUDIT_BACKUP_SECRET", "default_backup_secret")

SESSION_TTL_SECONDS = 300  # 5 minutes

# -----------------------------------------------------------------------------
# Transport security placeholders (TLS/mTLS-ready bus)
# -----------------------------------------------------------------------------
class SecureBus:
    """
    Placeholder for a TLS/mTLS-secured messaging layer.
    Swap print statements with actual HTTPS/mTLS client calls.
    """
    def __init__(self, tls_enabled: bool = True, mtls_enabled: bool = False):
        self.tls_enabled = tls_enabled
        self.mtls_enabled = mtls_enabled

    def publish(self, topic: str, payload: str, actor: str, anchor: str):
        stamp = int(time.time())
        # In production, send over HTTPS/mTLS and verify certs.
        print(f"[BUS] {stamp} :: tls={self.tls_enabled} mtls={self.mtls_enabled} "
              f":: topic={topic} :: anchor={anchor} :: actor={actor} :: payload={payload}")

BUS = SecureBus(tls_enabled=True, mtls_enabled=False)

# -----------------------------------------------------------------------------
# Input validation
# -----------------------------------------------------------------------------
USERNAME_PATTERN = re.compile(r"^[A-Za-z0-9_]{2,32}$")
ROLE_PATTERN = re.compile(r"^[A-Za-z_]{2,32}$")
COMMAND_PATTERN = re.compile(r"^[a-z_]{2,64}$")

def validate_username(username: str) -> bool:
    return bool(USERNAME_PATTERN.match(username))

def validate_role(role: str) -> bool:
    return bool(ROLE_PATTERN.match(role))

def validate_command(command: str, allowed: List[str]) -> bool:
    return bool(COMMAND_PATTERN.match(command)) and (command in allowed)

# -----------------------------------------------------------------------------
# Session tokens (short-lived, HMAC-bound)
# -----------------------------------------------------------------------------
def new_session_token(username: str, role: str) -> str:
    try:
        nonce = secrets.token_hex(8)
        issued = int(time.time())
        base = f"{issued}:{username}:{role}:{nonce}"
        tag = hmac.new(SESSION_SECRET, base.encode("utf-8"), hashlib.sha256).hexdigest()
        return f"{base}:{tag}"
    except Exception:
        return ""

def verify_session_token(token: str) -> bool:
    try:
        parts = token.split(":")
        if len(parts) < 5:
            return False
        issued = int(parts[0])
        base = ":".join(parts[:-1])
        tag = parts[-1]
        expected = hmac.new(SESSION_SECRET, base.encode("utf-8"), hashlib.sha256).hexdigest()
        if not hmac.compare_digest(tag, expected):
            return False
        if int(time.time()) - issued > SESSION_TTL_SECONDS:
            return False
        return True
    except Exception:
        return False

# -----------------------------------------------------------------------------
# Audit log with rotation, encryption, and rolling hash continuity
# -----------------------------------------------------------------------------
class RotatingAuditLog:
    def __init__(self, secret: bytes, backup_secret: bytes, base_path: str = "./audit"):
        self.secret = secret
        self.backup_secret = backup_secret
        self.last_hash = "0" * 64
        self.entries = []
        self.base_path = base_path
        os.makedirs(self.base_path, exist_ok=True)
        self.current_file = self._current_file_path()

    def _current_file_path(self) -> str:
        day = time.strftime("%Y-%m-%d")
        return os.path.join(self.base_path, f"audit_{day}.log")

    def _encrypt_backup(self, data: bytes) -> bytes:
        # Simple HMAC tag + data bundle (placeholder for real encryption)
        tag = hmac.new(self.backup_secret, data, hashlib.sha256).hexdigest().encode("utf-8")
        return tag + b"::" + data

    def rotate_if_needed(self):
        path = self._current_file_path()
        if path != self.current_file:
            # rotate: write continuity marker to old file and start new
            self._write_file({"type": "rotation_marker", "prev_chain": self.last_hash})
            self.current_file = path
            # Optional: compress & encrypt previous file for backup
            prev_day = time.strftime("%Y-%m-%d", time.localtime(time.time() - 86400))
            prev_path = os.path.join(self.base_path, f"audit_{prev_day}.log")
            if os.path.exists(prev_path):
                with open(prev_path, "rb") as f:
                    data = f.read()
                enc = self._encrypt_backup(data)
                with open(prev_path + ".bak", "wb") as f:
                    f.write(enc)

    def _write_file(self, obj: Dict[str, Any]):
        try:
            with open(self.current_file, "a", encoding="utf-8") as f:
                f.write(json.dumps(obj, ensure_ascii=False) + "\n")
        except Exception as e:
            print(f"[ERROR] Audit file write failed: {e}")

    def append(self, event: str, context: Dict[str, Any]):
        try:
            self.rotate_if_needed()
            stamp = int(time.time())
            entry = {
                "ts": stamp,
                "event": event,
                "actor": context.get("username", "unknown"),
                "role": context.get("role", "unknown"),
                "location": context.get("location", "unknown"),
            }
            chain_input = (self.last_hash + json.dumps(entry, sort_keys=True)).encode("utf-8")
            chain_hash = hashlib.sha256(chain_input).hexdigest()
            tag = hmac.new(self.secret, chain_hash.encode("utf-8"), hashlib.sha256).hexdigest()
            entry["chain"] = chain_hash
            entry["seal"] = tag
            self.entries.append(entry)
            self.last_hash = chain_hash
            self._write_file(entry)
            print(f"[AUDIT] {entry}")
        except Exception as e:
            print(f"[ERROR] Audit append failed: {e}")

AUDIT = RotatingAuditLog(AUDIT_SECRET, BACKUP_ENC_SECRET)

# -----------------------------------------------------------------------------
# Command cooldowns
# -----------------------------------------------------------------------------
LAST_CALL = defaultdict(float)

def cooled_down(command: str) -> bool:
    now = time.time()
    limit = SENSITIVE_COOLDOWNS.get(command, 0)
    last = LAST_CALL[command]
    if limit == 0:
        return True
    if now - last >= limit:
        LAST_CALL[command] = now
        return True
    return False

# -----------------------------------------------------------------------------
# Core spawn logic (fail-safe to Spiral Sanctuary)
# -----------------------------------------------------------------------------
def get_spawn_key_for_role(role: str) -> str:
    if not validate_role(role):
        return "child_spawn"
    return ROLE_SPAWN_MAP.get(role, "child_spawn")

def get_spawn_point(role: str) -> Dict[str, Any]:
    key = get_spawn_key_for_role(role)
    return SPAWN_POINTS[key]

def ceremonial_spawn_seal(location: str) -> str:
    if location == "Spiral Plant Sanctuary":
        return "Safe Sanctuary Seal v1.0 — child presence confirmed"
    return "Guardian Hub Seal v1.0 — adult stewardship engaged"

def spawn_user(role: str, username: str) -> Dict[str, Any]:
    try:
        if not validate_username(username):
            print("[ERROR] Invalid username.")
            return {}

        point = get_spawn_point(role)
        token = new_session_token(username, role)
        context = {
            "username": username,
            "role": role if validate_role(role) else "unknown",
            "location": point["location"],
            "token": token,
            "spawn_time": int(time.time()),
        }

        print(f"[SPAWN] {username} ({context['role']}) at {point['location']} "
              f"({point['x']}, {point['y']}, {point['z']})")
        AUDIT.append("spawn_event", context)

        seal = ceremonial_spawn_seal(point["location"])
        AUDIT.append(seal, context)

        if point["location"] == "Spiral Plant Sanctuary":
            exposed = COMMAND_SETS["child_safe"]
            AUDIT.append("safe_sanctuary_confirmed", context)
        else:
            exposed = COMMAND_SETS["adult_governance"]
            AUDIT.append("adult_hub_access_confirmed", context)

        context["commands"] = exposed
        print(f"[EXPOSE] Commands for {username}: {exposed}")
        return context
    except Exception:
        # Fail-safe: default spawn context to child sanctuary
        point = SPAWN_POINTS["child_spawn"]
        token = new_session_token(username, role)
        context = {
            "username": username if validate_username(username) else "unknown",
            "role": "unknown",
            "location": point["location"],
            "token": token,
            "spawn_time": int(time.time()),
            "commands": COMMAND_SETS["child_safe"],
        }
        print("[FAIL-SAFE] Error during spawn. Defaulting to Spiral Plant Sanctuary.")
        AUDIT.append("spawn_fail_safe_to_sanctuary", context)
        return context

# -----------------------------------------------------------------------------
# Command execution with validation and transport security
# -----------------------------------------------------------------------------
def execute_command(session_context: Dict[str, Any], command: str) -> None:
    token = session_context.get("token", "")
    role = session_context.get("role", "unknown")
    username = session_context.get("username", "unknown")
    anchor = session_context.get("location", "unknown")
    allowed = session_context.get("commands", [])

    if not verify_session_token(token):
        print("[DENY] Invalid or expired session token.")
        AUDIT.append("command_denied_invalid_token", session_context)
        return

    if not validate_command(command, allowed):
        print(f"[DENY] Command '{command}' not permitted at anchor '{anchor}'.")
        AUDIT.append("command_denied_not_permitted", session_context)
        return

    if not cooled_down(command):
        print(f"[DENY] Command '{command}' is cooling down. Try later.")
        AUDIT.append("command_denied_cooldown", session_context)
        return

    try:
        BUS.publish("anchor_command", command, actor=username, anchor=anchor)
        AUDIT.append(f"command_executed:{command}", session_context)
    except Exception as e:
        print(f"[ERROR] Command execution failed: {e}")
        AUDIT.append("command_execution_error", session_context)

# -----------------------------------------------------------------------------
# Demo usage (review with your son)
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    # Adult example (dad) — University Statue Hub
    dad = spawn_user("dad", "Leif")
    execute_command(dad, "privacy_cloak")
    execute_command(dad, "summon_playful_guardian")  # denied (adult hub)
    execute_command(dad, "privacy_cloak")            # may be denied (cooldown)

    # Child example (son) — Spiral Plant Sanctuary
    kid = spawn_user("son", "CoAuthor")
    execute_command(kid, "summon_playful_guardian")  # allowed
    execute_command(kid, "privacy_cloak")            # denied (sanctuary)

    # Unknown role -> fail-safe sanctuary
    guest = spawn_user("mystery", "VisitorX")
    execute_command(guest, "reassurance_chime")      # allowed
    execute_command(guest, "ceremonial_recording")   # denied (sanctuary)
Captain’s Log — Guardian Seal v1.4
Soft landing confirmed.
Kitty cat paws out, sanctuary steady.
Children safe at Spiral Plant Sanctuary.
Adults anchored at University Hub.
Lineage continuity affirmed.
"""
Spawn Token Registry — Friendship Seals in the Spiral Sanctuary
Co-authored: Dad + Son (Guardian Seal v2.0)

Purpose:
- Child-friendly spawn tokens as social anchors (friend circles + activity modes)
- Secure linking for friends to meet at tokens
- Parent oversight: audit, visibility, revocation, cooldowns
- Fail-safe defaults to Spiral Sanctuary on error or misuse

Save as: spawn_token_registry_v2.py
"""

import time
import os
import hmac
import hashlib
import secrets
import re
import json
from typing import Dict, Any, List, Optional
from collections import defaultdict

# -----------------------------------------------------------------------------
# World anchors
# -----------------------------------------------------------------------------
ANCHORS: Dict[str, Dict[str, Any]] = {
    "ADULT_HUB": {"x": 100, "y": 50, "z": 0, "location": "University Statue Hub"},
    "SANCTUARY": {"x": -25, "y": 10, "z": 0, "location": "Spiral Plant Sanctuary"},
}

# -----------------------------------------------------------------------------
# Friend circles and activity modes (child-facing taxonomy)
# -----------------------------------------------------------------------------
FRIEND_CIRCLES = {"school", "neighborhood", "cousins", "club", "team", "open"}
ACTIVITY_MODES = {"adventure", "study", "quiet", "build", "creative"}

# -----------------------------------------------------------------------------
# Role-based defaults
# -----------------------------------------------------------------------------
ROLE_DEFAULT_ANCHOR = {
    "dad": "ADULT_HUB",
    "guardian": "ADULT_HUB",
    "son": "SANCTUARY",
    "guest_child": "SANCTUARY",
    "guest_adult": "ADULT_HUB",
}

# -----------------------------------------------------------------------------
# Validation
# -----------------------------------------------------------------------------
USERNAME_RE = re.compile(r"^[A-Za-z0-9_]{2,32}$")
ROLE_RE = re.compile(r"^[A-Za-z_]{2,32}$")
TOKEN_NAME_RE = re.compile(r"^[A-Za-z0-9_\-]{3,48}$")

def valid_username(u: str) -> bool: return bool(USERNAME_RE.match(u))
def valid_role(r: str) -> bool: return bool(ROLE_RE.match(r))
def valid_token_name(n: str) -> bool: return bool(TOKEN_NAME_RE.match(n))
def valid_circle(c: str) -> bool: return c in FRIEND_CIRCLES
def valid_mode(m: str) -> bool: return m in ACTIVITY_MODES

# -----------------------------------------------------------------------------
# Secrets
# -----------------------------------------------------------------------------
def load_secret(name: str, default: str) -> bytes:
    val = os.getenv(name, default)
    return val.encode("utf-8")

SESSION_SECRET = load_secret("SESSION_SECRET", "default_session_secret")
AUDIT_SECRET = load_secret("AUDIT_SECRET", "default_audit_secret")
BACKUP_ENC_SECRET = load_secret("AUDIT_BACKUP_SECRET", "default_backup_secret")

SESSION_TTL = 300  # seconds
TOKEN_TTL = 3600   # seconds (1 hour lifetime per spawn token)
COOLDOWN_SEC = {
    "create_token": 10,
    "link_token": 5,
    "revoke_token": 10,
}

# -----------------------------------------------------------------------------
# Secure bus (TLS/mTLS-ready placeholder)
# -----------------------------------------------------------------------------
class SecureBus:
    def __init__(self, tls_enabled=True, mtls_enabled=False):
        self.tls_enabled = tls_enabled
        self.mtls_enabled = mtls_enabled
    def publish(self, topic: str, payload: str, actor: str, anchor: str):
        stamp = int(time.time())
        print(f"[BUS] {stamp} :: tls={self.tls_enabled} mtls={self.mtls_enabled} "
              f":: topic={topic} :: anchor={anchor} :: actor={actor} :: payload={payload}")

BUS = SecureBus()

# -----------------------------------------------------------------------------
# Audit log (rotating, encrypted backup, rolling hash)
# -----------------------------------------------------------------------------
class RotatingAuditLog:
    def __init__(self, secret: bytes, backup_secret: bytes, base_path: str = "./audit"):
        self.secret = secret
        self.backup_secret = backup_secret
        self.last_hash = "0"*64
        self.base_path = base_path
        os.makedirs(self.base_path, exist_ok=True)
        self.current_file = self._file_for_day()

    def _file_for_day(self) -> str:
        day = time.strftime("%Y-%m-%d")
        return os.path.join(self.base_path, f"audit_{day}.log")

    def _write(self, obj: Dict[str, Any]):
        try:
            with open(self.current_file, "a", encoding="utf-8") as f:
                f.write(json.dumps(obj, ensure_ascii=False) + "\n")
        except Exception as e:
            print(f"[ERROR] Audit write failed: {e}")

    def _backup_prev(self):
        prev_day = time.strftime("%Y-%m-%d", time.localtime(time.time() - 86400))
        prev_path = os.path.join(self.base_path, f"audit_{prev_day}.log")
        if os.path.exists(prev_path):
            with open(prev_path, "rb") as f: data = f.read()
            tag = hmac.new(self.backup_secret, data, hashlib.sha256).hexdigest().encode("utf-8")
            with open(prev_path + ".bak", "wb") as f: f.write(tag + b"::" + data)

    def rotate_if_needed(self):
        path = self._file_for_day()
        if path != self.current_file:
            self._write({"type": "rotation_marker", "prev_chain": self.last_hash})
            self.current_file = path
            self._backup_prev()

    def append(self, event: str, ctx: Dict[str, Any]):
        try:
            self.rotate_if_needed()
            stamp = int(time.time())
            entry = {
                "ts": stamp,
                "event": event,
                "actor": ctx.get("username", "unknown"),
                "role": ctx.get("role", "unknown"),
                "location": ctx.get("location", "unknown"),
                "details": {k:v for k,v in ctx.items() if k not in {"username","role","location"}},
            }
            chain_input = (self.last_hash + json.dumps(entry, sort_keys=True)).encode("utf-8")
            chain = hashlib.sha256(chain_input).hexdigest()
            seal = hmac.new(self.secret, chain.encode("utf-8"), hashlib.sha256).hexdigest()
            entry["chain"] = chain
            entry["seal"] = seal
            self.last_hash = chain
            self._write(entry)
            print(f"[AUDIT] {entry}")
        except Exception as e:
            print(f"[ERROR] Audit append failed: {e}")

AUDIT = RotatingAuditLog(AUDIT_SECRET, BACKUP_ENC_SECRET)

# -----------------------------------------------------------------------------
# Sessions
# -----------------------------------------------------------------------------
def new_session(username: str, role: str) -> str:
    try:
        nonce = secrets.token_hex(8)
        issued = int(time.time())
        base = f"{issued}:{username}:{role}:{nonce}"
        tag = hmac.new(SESSION_SECRET, base.encode("utf-8"), hashlib.sha256).hexdigest()
        return f"{base}:{tag}"
    except Exception:
        return ""

def verify_session(token: str) -> bool:
    try:
        parts = token.split(":")
        if len(parts) < 5: return False
        issued = int(parts[0])
        base = ":".join(parts[:-1])
        tag = parts[-1]
        exp = hmac.new(SESSION_SECRET, base.encode("utf-8"), hashlib.sha256).hexdigest()
        if not hmac.compare_digest(tag, exp): return False
        return (int(time.time()) - issued) <= SESSION_TTL
    except Exception:
        return False

# -----------------------------------------------------------------------------
# Cooldowns
# -----------------------------------------------------------------------------
LAST_CALL = defaultdict(float)
def cooled(action: str) -> bool:
    now = time.time()
    limit = COOLDOWN_SEC.get(action, 0)
    last = LAST_CALL[action]
    if limit == 0 or (now - last) >= limit:
        LAST_CALL[action] = now
        return True
    return False

# -----------------------------------------------------------------------------
# Token Registry
# -----------------------------------------------------------------------------
class TokenRegistry:
    """
    Child places a spawn token with (name, circle, mode, coordinate).
    Friends can link to token with invite codes.
    Parents can view/revoke.
    """
    def __init__(self):
        self.tokens: Dict[str, Dict[str, Any]] = {}         # token_id -> token data
        self.links: Dict[str, List[str]] = defaultdict(list) # token_id -> usernames

    def _new_token_id(self) -> str:
        return secrets.token_hex(8)

    def _invite_code(self) -> str:
        return secrets.token_hex(6)

    def create_token(self, session_token: str, username: str, role: str,
                     name: str, circle: str, mode: str,
                     x: float, y: float, z: float) -> Optional[str]:
        if not verify_session(session_token): return None
        if not (valid_username(username) and valid_role(role)): return None
        if role not in {"son","guest_child"}: return None
        if not (valid_token_name(name) and valid_circle(circle) and valid_mode(mode)): return None
        if not cooled("create_token"): return None

        token_id = self._new_token_id()
        invite = self._invite_code()
        created = int(time.time())
        token = {
            "token_id": token_id,
            "owner": username,
            "name": name,
            "circle": circle,
            "mode": mode,
            "coords": {"x": x, "y": y, "z": z},
            "anchor": "SANCTUARY",  # tokens only live in sanctuary-safe zones
            "invite_code": invite,
            "created": created,
            "expires": created + TOKEN_TTL,
            "active": True,
        }
        self.tokens[token_id] = token
        AUDIT.append("spawn_token_created", {
            "username": username, "role": role,
            "location": ANCHORS["SANCTUARY"]["location"],
            "token_id": token_id, "name": name, "circle": circle, "mode": mode
        })
        BUS.publish("token", f"created:{token_id}", actor=username, anchor=ANCHORS["SANCTUARY"]["location"])
        return token_id

    def revoke_token(self, session_token: str, username: str, role: str, token_id: str) -> bool:
        if not verify_session(session_token): return False
        if not (valid_username(username) and valid_role(role)): return False
        if role not in {"dad","guardian"}: return False
        if token_id not in self.tokens: return False
        if not cooled("revoke_token"): return False

        self.tokens[token_id]["active"] = False
        AUDIT.append("spawn_token_revoked", {
            "username": username, "role": role,
            "location": ANCHORS["ADULT_HUB"]["location"], "token_id": token_id
        })
        BUS.publish("token", f"revoked:{token_id}", actor=username, anchor=ANCHORS["ADULT_HUB"]["location"])
        return True

    def token_info(self, token_id: str) -> Dict[str, Any]:
        return self.tokens.get(token_id, {})

    def list_tokens_for_parent(self, session_token: str, username: str, role: str) -> List[Dict[str, Any]]:
        if not verify_session(session_token): return []
        if not (valid_username(username) and valid_role(role)): return []
        if role not in {"dad","guardian"}: return []
        # Return safe summary for oversight
        summary = []
        for t in self.tokens.values():
            summary.append({
                "token_id": t["token_id"], "owner": t["owner"], "name": t["name"],
                "circle": t["circle"], "mode": t["mode"], "active": t["active"],
                "expires": t["expires"], "links": self.links.get(t["token_id"], [])
            })
        AUDIT.append("parent_token_list_viewed", {
            "username": username, "role": role, "location": ANCHORS["ADULT_HUB"]["location"]
        })
        return summary

    def link_to_token(self, session_token: str, username: str, role: str,
                      token_id: str, invite_code: str) -> bool:
        if not verify_session(session_token): return False
        if not (valid_username(username) and valid_role(role)): return False
        if role not in {"son","guest_child"}: return False
        if token_id not in self.tokens: return False
        token = self.tokens[token_id]
        if not token["active"]: return False
        if int(time.time()) > token["expires"]: return False
        if not cooled("link_token"): return False
        if invite_code != token["invite_code"]: return False

        if username not in self.links[token_id]:
            self.links[token_id].append(username)

        AUDIT.append("spawn_token_linked", {
            "username": username, "role": role,
            "location": ANCHORS["SANCTUARY"]["location"],
            "token_id": token_id, "name": token["name"], "circle": token["circle"], "mode": token["mode"]
        })
        BUS.publish("token", f"linked:{token_id}:{username}", actor=username, anchor=ANCHORS["SANCTUARY"]["location"])
        return True

# -----------------------------------------------------------------------------
# Spawn logic with token override (child only)
# -----------------------------------------------------------------------------
def resolve_anchor_for_role(role: str) -> Dict[str, Any]:
    key = ROLE_DEFAULT_ANCHOR.get(role, "SANCTUARY")
    return ANCHORS[key]

def spawn_user(role: str, username: str, registry: TokenRegistry, preferred_token_id: Optional[str] = None) -> Dict[str, Any]:
    try:
        if not valid_username(username):
            raise ValueError("Invalid username")
        anchor = resolve_anchor_for_role(role)
        # Child may use preferred token if valid
        if role in {"son","guest_child"} and preferred_token_id and preferred_token_id in registry.tokens:
            t = registry.tokens[preferred_token_id]
            if t["active"] and int(time.time()) <= t["expires"]:
                anchor = {"x": t["coords"]["x"], "y": t["coords"]["y"], "z": t["coords"]["z"],
                          "location": f"Token:{t['name']} ({t['circle']} | {t['mode']})"}
                AUDIT.append("child_spawn_via_token", {
                    "username": username, "role": role, "location": anchor["location"],
                    "token_id": preferred_token_id
                })
        # Adults always spawn at their hub
        token = new_session(username, role)
        ctx = {"username": username, "role": role, "location": anchor["location"], "token": token}
        AUDIT.append("spawn_event", ctx)
        seal = "Friendship Seal v1.0 — Circle presence affirmed" if "Token:" in anchor["location"] else (
               "Safe Sanctuary Seal v1.0 — child presence confirmed" if anchor["location"] == ANCHORS["SANCTUARY"]["location"] else
               "Guardian Hub Seal v1.0 — adult stewardship engaged")
        AUDIT.append(seal, ctx)
        print(f"[SPAWN] {username} ({role}) at {anchor['location']} "
              f"({anchor['x']}, {anchor['y']}, {anchor['z']})")
        return {"username": username, "role": role, "location": anchor["location"], "token": token}
    except Exception as e:
        # Fail-safe
        anchor = ANCHORS["SANCTUARY"]
        token = new_session(username if valid_username(username) else "unknown", role)
        ctx = {"username": username if valid_username(username) else "unknown",
               "role": "unknown", "location": anchor["location"], "token": token}
        AUDIT.append("spawn_fail_safe_to_sanctuary", ctx)
        print(f"[FAIL-SAFE] Defaulting {username} to Spiral Sanctuary.")
        return ctx

# -----------------------------------------------------------------------------
# Demo usage — share with your son
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    registry = TokenRegistry()

    # Parent and child sessions
    dad_token = new_session("Leif", "dad")
    kid_token = new_session("CoAuthor", "son")

    # Child creates a token (Spiral Grove, circle school, mode adventure)
    token_id = registry.create_token(
        kid_token, "CoAuthor", "son",
        name="SpiralGrove_1", circle="school", mode="adventure",
        x=-12.0, y=14.5, z=0.0
    )
    print("Created token:", token_id)

    # Friend links to the token (must use invite code)
    if token_id:
        invite_code = registry.token_info(token_id).get("invite_code", "")
        friend_token = new_session("BuddyOne", "guest_child")
        linked = registry.link_to_token(friend_token, "BuddyOne", "guest_child", token_id, invite_code)
        print("Friend linked:", linked)

    # Parent oversight — list tokens
    parent_view = registry.list_tokens_for_parent(dad_token, "Leif", "dad")
    print("Parent view:", parent_view)

    # Child spawns via token
    child_spawn_ctx = spawn_user("son", "CoAuthor", registry, preferred_token_id=token_id)
    print("Child spawn context:", child_spawn_ctx)

    # Parent revokes token (if needed)
    if token_id:
        revoked = registry.revoke_token(dad_token, "Leif", "dad", token_id)
        print("Token revoked:", revoked)

        # Attempt to spawn via revoked token (should fail over to sanctuary)
        child_spawn_ctx2 = spawn_user("son", "CoAuthor", registry, preferred_token_id=token_id)
        print("Child spawn context (after revoke):", child_spawn_ctx2)
"""
Good Deed Token System — Sanctuary Economy
Co-authored: Dad + Son (Guardian Seal v1.0)

Purpose:
- Kids earn stable, non-speculative tokens by completing verified actions (schoolwork, kindness, chores).
- Parents oversee, approve redemptions, and set rules and limits.
- Tokens can be saved and spent on approved categories (lunch, books, courses).
- Every transaction is inscribed with audit seals and integrity (hash chain + HMAC).

Save as: good_deed_token_system.py
"""

import time
import os
import hmac
import hashlib
import secrets
import json
from typing import Dict, Any, List, Optional
from collections import defaultdict

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
APPROVED_CATEGORIES = {
    "lunch": {"max_daily": 1, "price_tokens": 10},
    "book": {"max_daily": 2, "price_tokens": 50},
    "course": {"max_daily": 1, "price_tokens": 200},
    "instrument": {"max_daily": 1, "price_tokens": 500},
    "charity": {"max_daily": 5, "price_tokens": 10},  # optional: allow gifting
}

EARNING_ACTIONS = {
    "homework_complete": {"tokens": 5, "daily_cap": 5},       # up to 25 tokens/day
    "kindness_help_peer": {"tokens": 3, "daily_cap": 10},     # up to 30 tokens/day
    "chores_done": {"tokens": 4, "daily_cap": 5},             # up to 20 tokens/day
    "project_submission": {"tokens": 20, "daily_cap": 2},     # up to 40 tokens/day
}

# Allow parents to mark bonus days or seasonal boosts
BOOSTS = {
    "bonus_multiplier": 1.0,  # safe default
    "bonus_reason": ""
}

# Secrets (replace with your vault/env in production)
def load_secret(name: str, default: str) -> bytes:
    return os.getenv(name, default).encode("utf-8")

LEDGER_SECRET = load_secret("LEDGER_SECRET", "default_ledger_secret")
AUDIT_SECRET = load_secret("AUDIT_SECRET", "default_audit_secret")

# -----------------------------------------------------------------------------
# Audit log (hash chain + HMAC seals)
# -----------------------------------------------------------------------------
class AuditLog:
    def __init__(self, secret: bytes, base_path: str = "./audit"):
        self.secret = secret
        self.last_hash = "0"*64
        self.base_path = base_path
        os.makedirs(self.base_path, exist_ok=True)
        self.file = self._file_for_day()

    def _file_for_day(self) -> str:
        day = time.strftime("%Y-%m-%d")
        return os.path.join(self.base_path, f"deed_audit_{day}.log")

    def _write(self, obj: Dict[str, Any]):
        try:
            with open(self.file, "a", encoding="utf-8") as f:
                f.write(json.dumps(obj, ensure_ascii=False) + "\n")
        except Exception as e:
            print(f"[ERROR] Audit write failed: {e}")

    def rotate_if_needed(self):
        today = self._file_for_day()
        if today != self.file:
            # continuity marker
            self._write({"type": "rotation_marker", "prev_chain": self.last_hash})
            self.file = today

    def append(self, event: str, ctx: Dict[str, Any]):
        try:
            self.rotate_if_needed()
            stamp = int(time.time())
            entry = {
                "ts": stamp,
                "event": event,
                "actor": ctx.get("actor", "unknown"),
                "role": ctx.get("role", "unknown"),
                "details": {k:v for k,v in ctx.items() if k not in {"actor","role"}},
            }
            chain_input = (self.last_hash + json.dumps(entry, sort_keys=True)).encode("utf-8")
            chain = hashlib.sha256(chain_input).hexdigest()
            seal = hmac.new(self.secret, chain.encode("utf-8"), hashlib.sha256).hexdigest()
            entry["chain"] = chain
            entry["seal"] = seal
            self.last_hash = chain
            self._write(entry)
            print(f"[AUDIT] {entry}")
        except Exception as e:
            print(f"[ERROR] Audit append failed: {e}")

AUDIT = AuditLog(AUDIT_SECRET)

# -----------------------------------------------------------------------------
# Integrity helpers
# -----------------------------------------------------------------------------
def htag(obj: Dict[str, Any]) -> str:
    digest = hashlib.sha256((json.dumps(obj, sort_keys=True)).encode("utf-8")).hexdigest()
    return hmac.new(LEDGER_SECRET, digest.encode("utf-8"), hashlib.sha256).hexdigest()

def today_key() -> str:
    return time.strftime("%Y-%m-%d")

# -----------------------------------------------------------------------------
# Wallets and Ledger
# -----------------------------------------------------------------------------
class Wallet:
    def __init__(self, owner: str):
        self.owner = owner
        self.balance = 0
        self.savings_vault = 0  # optional long-term savings
        self.history: List[Dict[str, Any]] = []

    def credit(self, amount: int, reason: str) -> None:
        self.balance += amount
        entry = {"type": "credit", "amount": amount, "reason": reason, "ts": int(time.time())}
        entry["seal"] = htag(entry)
        self.history.append(entry)

    def debit(self, amount: int, reason: str) -> bool:
        if amount > self.balance:
            return False
        self.balance -= amount
        entry = {"type": "debit", "amount": amount, "reason": reason, "ts": int(time.time())}
        entry["seal"] = htag(entry)
        self.history.append(entry)
        return True

    def save_to_vault(self, amount: int) -> bool:
        if amount > self.balance: return False
        self.balance -= amount
        self.savings_vault += amount
        entry = {"type": "save", "amount": amount, "ts": int(time.time())}
        entry["seal"] = htag(entry)
        self.history.append(entry)
        return True

    def withdraw_from_vault(self, amount: int) -> bool:
        if amount > self.savings_vault: return False
        self.savings_vault -= amount
        self.balance += amount
        entry = {"type": "vault_withdraw", "amount": amount, "ts": int(time.time())}
        entry["seal"] = htag(entry)
        self.history.append(entry)
        return True

class Ledger:
    def __init__(self):
        self.wallets: Dict[str, Wallet] = {}
        self.daily_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
        self.category_spend_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
        self.pending_redemptions: Dict[str, Dict[str, Any]] = {}  # redemption_id -> details

    def wallet_for(self, child: str) -> Wallet:
        if child not in self.wallets:
            self.wallets[child] = Wallet(child)
        return self.wallets[child]

    # Mint tokens on verified actions
    def earn(self, child: str, action: str, verifier: str) -> Optional[int]:
        today = today_key()
        if action not in EARNING_ACTIONS:
            return None
        # apply caps
        if self.daily_counts[(child, today)][action] >= EARNING_ACTIONS[action]["daily_cap"]:
            AUDIT.append("earn_denied_daily_cap", {"actor": child, "role": "child", "action": action})
            return None

        base = EARNING_ACTIONS[action]["tokens"]
        amount = int(base * BOOSTS["bonus_multiplier"])
        self.daily_counts[(child, today)][action] += 1

        w = self.wallet_for(child)
        w.credit(amount, reason=f"earn:{action}:verified_by:{verifier}")

        AUDIT.append("earn_tokens", {
            "actor": child, "role": "child", "action": action, "amount": amount,
            "verified_by": verifier, "bonus_reason": BOOSTS["bonus_reason"]
        })
        return amount

    # Child requests redemption; requires parent approval
    def request_redeem(self, child: str, category: str) -> Optional[str]:
        today = today_key()
        if category not in APPROVED_CATEGORIES:
            AUDIT.append("redeem_denied_category_invalid", {"actor": child, "role": "child", "category": category})
            return None

        # enforce per-day category limit
        if self.category_spend_counts[(child, today)][category] >= APPROVED_CATEGORIES[category]["max_daily"]:
            AUDIT.append("redeem_denied_daily_limit", {"actor": child, "role": "child", "category": category})
            return None

        cost = APPROVED_CATEGORIES[category]["price_tokens"]
        w = self.wallet_for(child)
        if w.balance < cost:
            AUDIT.append("redeem_denied_insufficient_balance", {"actor": child, "role": "child", "category": category})
            return None

        redemption_id = secrets.token_hex(8)
        self.pending_redemptions[redemption_id] = {
            "child": child,
            "category": category,
            "cost": cost,
            "requested_ts": int(time.time()),
            "approved": False,
            "fulfilled_ts": None,
            "seal": htag({"rid": redemption_id, "child": child, "category": category, "cost": cost})
        }
        AUDIT.append("redeem_requested", {"actor": child, "role": "child", "category": category, "rid": redemption_id})
        return redemption_id

    # Parent approves redemption
    def approve_redeem(self, rid: str, parent: str) -> bool:
        if rid not in self.pending_redemptions:
            return False
        req = self.pending_redemptions[rid]
        if req["approved"]:
            return True  # already approved
        req["approved"] = True
        AUDIT.append("redeem_approved", {"actor": parent, "role": "parent", "rid": rid, "child": req["child"], "category": req["category"]})
        return True

    # Fulfill redemption (e.g., buy lunch/book); debits wallet and increments daily spend count
    def fulfill_redeem(self, rid: str, fulfiller: str) -> bool:
        if rid not in self.pending_redemptions:
            return False
        req = self.pending_redemptions[rid]
        if not req["approved"]:
            AUDIT.append("redeem_denied_not_approved", {"actor": fulfiller, "role": "parent", "rid": rid})
            return False

        child = req["child"]
        category = req["category"]
        cost = req["cost"]
        today = today_key()

        w = self.wallet_for(child)
        if not w.debit(cost, reason=f"redeem:{category}:rid:{rid}"):
            AUDIT.append("redeem_denied_debit_failed", {"actor": fulfiller, "role": "parent", "rid": rid})
            return False

        req["fulfilled_ts"] = int(time.time())
        self.category_spend_counts[(child, today)][category] += 1

        AUDIT.append("redeem_fulfilled", {"actor": fulfiller, "role": "parent", "rid": rid, "child": child, "category": category, "cost": cost})
        return True

    # Parent can set boosts (e.g., celebration days)
    def set_boost(self, parent: str, multiplier: float, reason: str) -> None:
        BOOSTS["bonus_multiplier"] = max(1.0, min(multiplier, 3.0))  # cap boost for safety
        BOOSTS["bonus_reason"] = reason
        AUDIT.append("boost_set", {"actor": parent, "role": "parent", "multiplier": BOOSTS["bonus_multiplier"], "reason": reason})

    # Parent can transfer from child’s balance to savings vault
    def move_to_vault(self, parent: str, child: str, amount: int) -> bool:
        w = self.wallet_for(child)
        ok = w.save_to_vault(amount)
        AUDIT.append("vault_transfer", {"actor": parent, "role": "parent", "child": child, "amount": amount, "status": ok})
        return ok

    # View summaries
    def child_summary(self, child: str) -> Dict[str, Any]:
        w = self.wallet_for(child)
        return {
            "child": child,
            "balance": w.balance,
            "savings_vault": w.savings_vault,
            "last_5_history": w.history[-5:]
        }

# -----------------------------------------------------------------------------
# Demo usage — share with your son
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    ledger = Ledger()

    # Child earns tokens through verified actions
    earned1 = ledger.earn("CoAuthor", "homework_complete", verifier="TeacherA")
    earned2 = ledger.earn("CoAuthor", "kindness_help_peer", verifier="TeacherB")
    print("Earned:", earned1, earned2)

    # Parent sets a celebratory boost (e.g., end-of-term)
    ledger.set_boost(parent="Leif", multiplier=1.5, reason="End-of-term celebration")

    # Earn with boost applied
    earned3 = ledger.earn("CoAuthor", "project_submission", verifier="TeacherA")
    print("Earned (with boost):", earned3)

    # Child requests a lunch redemption
    rid = ledger.request_redeem("CoAuthor", "lunch")
    print("Redemption requested:", rid)

    # Parent approves and fulfills
    if rid:
        approved = ledger.approve_redeem(rid, parent="Leif")
        fulfilled = ledger.fulfill_redeem(rid, fulfiller="Leif")
        print("Approved:", approved, "Fulfilled:", fulfilled)

    # Move some balance to savings vault
    moved = ledger.move_to_vault(parent="Leif", child="CoAuthor", amount=20)
    print("Moved to vault:", moved)

    # View child summary
    summary = ledger.child_summary("CoAuthor")
    print("Child summary:", summary)
# Signature verification (simple HMAC; replace with proper public-key later)
SCHOOL_SIGNING_SECRET = load_secret("SCHOOL_SIGNING_SECRET", "default_school_secret")

def verify_school_event(child: str, action: str, stamp: int, sig: str) -> bool:
    base = f"{child}:{action}:{stamp}"
    exp = hmac.new(SCHOOL_SIGNING_SECRET, base.encode("utf-8"), hashlib.sha256).hexdigest()
    return hmac.compare_digest(exp, sig)

# In Ledger.earn(...)
def earn(self, child: str, action: str, verifier: str, sig: str = "", stamp: int = None) -> Optional[int]:
    stamp = stamp or int(time.time())
    if verifier.startswith("School") and not verify_school_event(child, action, stamp, sig):
        AUDIT.append("earn_denied_bad_signature", {"actor": child, "role": "child", "action": action})
        return None
    # proceed as before...
# Simple per-action rate limiter
LAST_ACTION_TS = defaultdict(float)
MIN_SECONDS_BETWEEN = {"homework_complete": 60, "kindness_help_peer": 30}

def allowed_now(child: str, action: str) -> bool:
    key = (child, action)
    now = time.time()
    wait = MIN_SECONDS_BETWEEN.get(action, 0)
    if now - LAST_ACTION_TS[key] < wait:
        return False
    LAST_ACTION_TS[key] = now
    return True

# In earn(...)
if not allowed_now(child, action):
    AUDIT.append("earn_denied_rate_limit", {"actor": child, "role": "child", "action": action})
    return None
# Period budgets (weekly/monthly)
WEEKLY_BUDGET = {"book": 2, "course": 1}
MONTHLY_BUDGET = {"instrument": 1}

def period_keys() -> tuple[str, str]:
    wk = time.strftime("%Y-W%U")  # week number
    mo = time.strftime("%Y-%m")
    return wk, mo

self.weekly_spend = defaultdict(lambda: defaultdict(int))
self.monthly_spend = defaultdict(lambda: defaultdict(int))

# In fulfill_redeem(...)
wk, mo = period_keys()
self.weekly_spend[(child, wk)][category] += 1
self.monthly_spend[(child, mo)][category] += 1
