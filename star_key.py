# -----------------------------------------
# star_key.py
# -----------------------------------------

class StarKey:
    def __init__(self):
        self.allowed_fields = {
            "bronze": ["energy_level"],
            "silver": ["energy_level", "focus_level"],
            "gold": ["energy_level", "focus_level", "social_load"]
        }

    def filter_context(self, tier, raw_context):
        """Return only the fields this tier is allowed to see."""
        allowed = self.allowed_fields.get(tier, [])
        return {k: v for k, v in raw_context.items() if k in allowed}