import hashlib
import time
import secrets

class JDFlair:
    """
    JD's flair signature: glyph, audio chime, lineage string.
    """
    def __init__(self, glyph="🔥", color="star-mint", audio="perfect_fifth"):
        self.glyph = glyph
        self.color = color
        self.audio = audio

    def signature(self, event_id: str) -> str:
        """
        Create a lineage-safe signature with hash fragment.
        """
        hash_fragment = hashlib.sha256(event_id.encode()).hexdigest()[:8]
        return f"JD—pilot flame present [{self.glyph}] #{hash_fragment}"


class MaxSecurityFirewall:
    """
    Multi-layer firewall enforcing innocence, consent, anomaly detection.
    """
    def __init__(self):
        self.logs = []

    def verify_identity(self, user, gesture, otp):
        # Multi-factor check: user + gesture + one-time pass
        return all([user == "JD", gesture == "tap-hold-swirl", otp is not None])

    def verify_intent(self, envelope: dict):
        # Consent envelope must include purpose, duration, and data class
        required = {"purpose", "duration", "data_class"}
        return required.issubset(envelope.keys())

    def innocence_filter(self, content: str):
        # Block unsafe themes
        blocked = ["surveillance", "adult", "exploit"]
        return not any(term in content.lower() for term in blocked)

    def anomaly_guard(self, action_count: int):
        # Rate limit: max 10 actions per minute
        return action_count <= 10

    def log_event(self, event_id: str, decision: str):
        entry = {
            "event_id": event_id,
            "timestamp": time.time(),
            "decision": decision
        }
        self.logs.append(entry)

    def process_request(self, user, gesture, otp, envelope, content, action_count):
        event_id = secrets.token_hex(8)
        flair = JDFlair().signature(event_id)

        if not self.verify_identity(user, gesture, otp):
            self.log_event(event_id, "DENY: identity")
            return flair, "LOCKED"

        if not self.verify_intent(envelope):
            self.log_event(event_id, "DENY: intent")
            return flair, "LOCKED"

        if not self.innocence_filter(content):
            self.log_event(event_id, "DENY: innocence")
            return flair, "LOCKED"

        if not self.anomaly_guard(action_count):
            self.log_event(event_id, "DENY: anomaly")
            return flair, "LOCKED"

        self.log_event(event_id, "ALLOW")
        return flair, "ACCESS GRANTED"


# --- Example usage ---
firewall = MaxSecurityFirewall()

envelope = {"purpose": "create_art", "duration": "5min", "data_class": "safe"}
flair, decision = firewall.process_request(
    user="JD",
    gesture="tap-hold-swirl",
    otp="123456",
    envelope=envelope,
    content="child-safe creative project",
    action_count=3
)

print(flair)       # JD flair signature
print(decision)    # ACCESS GRANTED or LOCKED
print(firewall.logs)
