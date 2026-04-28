"""
Safe-Zone Bubble Demo
Public Transparency Edition
Space LEAF Corp — Stewardship Systems

This demo simulates:
- bubble initialization
- object compatibility checks
- private correction alerts
- non-escalation logic
- fallback protocols

No real detection models are included.
This code is safe for public GitHub release.
"""

class SafeZoneBubble:
    def __init__(self):
        self.active = False
        self.baseline = None

    def initialize_bubble(self):
        self.active = True
        self.baseline = {
            "crowd_density": "normal",
            "object_flow": "normal",
            "environment_status": "stable"
        }
        print("[Bubble] Safe-Zone initialized. Nonviolent boundary active.")

    def detect_object(self, obj):
        """
        obj = {"shape": "...", "material": "...", "intent": "unknown"}
        Only shape/material are evaluated.
        """
        incompatible_shapes = ["weapon_shape", "blade_shape"]
        incompatible_materials = ["weapon_metal"]

        if obj["shape"] in incompatible_shapes:
            return "incompatible"
        if obj["material"] in incompatible_materials:
            return "incompatible"

        return "compatible"

    def private_alert(self):
        print(
            "[Bubble] Private Notice: An item you are carrying is not "
            "permitted in this safe zone. Would you like assistance correcting this?"
        )

    def notify_security(self):
        print("[Bubble] Security Notice: Boundary inconsistency detected. Verification requested.")

    def process_object(self, obj):
        status = self.detect_object(obj)

        if status == "compatible":
            print("[Bubble] Object compatible. Entry permitted.")
            return

        # Incompatible object detected
        self.private_alert()
        self.notify_security()

    def fallback_mode(self, reason):
        self.active = False
        print(f"[Bubble] Fallback Mode Activated: {reason}")
        print("[Bubble] Automated assistance offline. Human stewards now guiding boundary.")


# Demo run
if __name__ == "__main__":
    bubble = SafeZoneBubble()
    bubble.initialize_bubble()

    test_object = {
        "shape": "weapon_shape",
        "material": "steel",
        "intent": "unknown"
    }

    bubble.process_object(test_object)
