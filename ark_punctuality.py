# ark_punctuality.py
# Authorship: Leif William Sogge
# Purpose: Adult-friendly punctuality seals + kid-friendly speed boost (simulated)
# Notes: Non-punitive, playful, lineage-preserving

from datetime import datetime, timedelta
from dataclasses import dataclass, field
from typing import List, Dict

# ---------- ACL (gentle, role-based) ----------
ACL: Dict[str, Dict] = {
    "Leif William Sogge": {"role": "admin", "permissions": ["read", "write", "execute"]},
    # Extend as needed:
    # "Mom Name": {"role": "guardian", "permissions": ["read"]},
    # "Dad Name": {"role": "guardian", "permissions": ["read"]},
    # "Team Member": {"role": "witness", "permissions": ["read"]},
}

def has_permission(username: str, permission: str) -> bool:
    user = ACL.get(username, {"role": "none", "permissions": []})
    return permission in user.get("permissions", [])

# ---------- Captain's Log ----------
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

# ---------- Comic relief / seals ----------
class Seals:
    @staticmethod
    def punctuality_boost(username: str, login_time: datetime, start_time: datetime) -> str:
        delta = (start_time - login_time).total_seconds()
        if delta >= 0:
            # On-time or early
            return (
                f"[Punctuality Boost Seal] {username} logged in {'early' if delta>0 else 'on time'}.\n"
                f"Speed lines: >>> WHOOSH! >>>\n"
                f"Recognition: Stewardship acknowledged, calm vigilance affirmed."
            )
        else:
            # Late
            return (
                f"[Comic Relief Seal] {username} logged in late by {abs(int(delta))}s.\n"
                f"Sound Effect: *GULP!*\n"
                f"Message: Humor permitted, learning embraced. Let's begin with joy."
            )

    @staticmethod
    def kid_speed_boost(username: str, destination_time: datetime) -> str:
        now = datetime.now()
        diff = (destination_time - now).total_seconds()
        if diff <= 0:
            arrival = now
        else:
            arrival = now + timedelta(seconds=diff * 0.1)  # 90% reduction (simulated)
        return (
            f"[Super Speed Boost] {username} arrives at {arrival.strftime('%Y-%m-%d %H:%M:%S')}!\n"
            f"Animation: >>> zoom zoom >>>\n"
            f"Safety: Arrived safely and ready to learn."
        )

# ---------- Adults: punctuality service ----------
@dataclass
class PunctualityService:
    log: CaptainsLog

    def process_login(self, username: str, login_time: datetime, shift_start: datetime) -> str:
        seal = Seals.punctuality_boost(username, login_time, shift_start)
        # Inscribe to Captain's Log
        self.log.inscribe(
            author=username,
            event="PunctualityCheck",
            details=f"Login: {login_time.isoformat()}, ShiftStart: {shift_start.isoformat()}, Seal: {seal.splitlines()[0]}"
        )
        return seal

# ---------- Example usage ----------
if __name__ == "__main__":
    log = CaptainsLog()
    service = PunctualityService(log=log)

    # Adult example
    user = "Leif William Sogge"
    start = datetime.now().replace(microsecond=0) + timedelta(minutes=1)  # shift starts in 1 min
    login = datetime.now().replace(microsecond=0)
    if has_permission(user, "read"):
        print(service.process_login(username=user, login_time=login, shift_start=start))
    else:
        print("[Access Notice] Read permission required to view seals.")

    # Late example (for demo)
    late_login = datetime.now().replace(microsecond=0) + timedelta(minutes=2)
    print(service.process_login(username=user, login_time=late_login, shift_start=start))

    # Kid-friendly speed boost
    destination = datetime.now() + timedelta(minutes=10)
    print(Seals.kid_speed_boost(username="Kiddo", destination_time=destination))

    # Show last 3 log entries
    for entry in log.latest(3):
        print(f"[Captain’s Log] {entry.timestamp.strftime('%Y-%m-%d %H:%M:%S')} — {entry.author} — {entry.event}: {entry.details}")
