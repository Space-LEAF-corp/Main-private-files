"""
Captain's Log Entry Saver - Python Version
Saves ceremonial entries to captains-log/log.jsonl in newline-delimited JSON format.
Includes SHA-256 checksum for integrity verification.
"""
import os
import json
import datetime
import hashlib


LOG_DIR = "captains-log"
LOG_FILE = os.path.join(LOG_DIR, "log.jsonl")  # newline-delimited JSON


def ensure_dir(d):
    """Create directory if it doesn't exist."""
    if not os.path.exists(d):
        os.makedirs(d)


def add_checksum(entry):
    """Add SHA-256 checksum to entry for integrity verification."""
    # Create a copy without the checksum field to compute hash
    entry_copy = {k: v for k, v in entry.items() if k != "checksum_sha256"}
    s = json.dumps(entry_copy, sort_keys=True)
    entry["checksum_sha256"] = hashlib.sha256(s.encode("utf-8")).hexdigest()
    return entry


def save_entry(entry):
    """Save entry to log.jsonl with checksum."""
    ensure_dir(LOG_DIR)
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(json.dumps(entry, ensure_ascii=False) + "\n")
    print(f"Saved entry to {LOG_FILE}")


if __name__ == "__main__":
    # Example entry: ark-001 - The university became a moving ship
    entry = {
        "artifactId": "ark-001",
        "title": "The university became a moving ship",
        "timestamp": datetime.datetime.now().isoformat(),
        "author": "Leif William Sogge",
        "ceremonialContext": {
            "vessel": "Ark",
            "covenant": "Protection and peace",
            "status": "Consecrated and complete"
        },
        "entry": "On this day, the university becomes a moving ship—an Ark. Everyone aboard is safe.",
        "artifacts": [
            {
                "type": "image",
                "filename": "9ZYTqXhtqaFpyp6av5bHi.jpeg",
                "label": "Water-port vessel",
                "notes": "Twin circular engines; guardians stand watch."
            },
            {
                "type": "image",
                "filename": "F4TG8oJMN7gomr2TbCeWS.jpeg",
                "label": "Floating city above ocean",
                "notes": "Waterfalls cascade; sanctuary architecture."
            },
            {
                "type": "image",
                "filename": "tcSw9D3UpZBHqcBubnMuo.jpeg",
                "label": "Hybrid Pokémon (name reserved)",
                "notes": "Joy and flame—a guardian of balance."
            }
        ],
        "seals": {
            "externalAuthority": "God",
            "localSealApplied": False
        }
    }
    
    # Add checksum for integrity
    entry = add_checksum(entry)
    
    # Save entry
    save_entry(entry)
