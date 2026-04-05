import hashlib
import time
from dataclasses import dataclass

# -----------------------------
# CONFIG & SAFETY PARAMETERS
# -----------------------------

@dataclass
class StoryConfig:
    audience: str              # "kid", "teen", "adult", "crew"
    allow_scary: bool          # explicit consent
    max_intensity: int         # 1–10, but we will clamp by audience
    topics: list               # e.g. ["space", "mystery", "submarine"]
    include_jump_scares: bool  # extra guardrail
    allow_real_names: bool     # privacy guardrail


# -----------------------------
# SAFETY & ETHICS CHECKS
# -----------------------------

def check_audience_bounds(cfg: StoryConfig) -> list:
    issues = []

    # Hard caps by audience
    caps = {
        "kid": 4,
        "teen": 7,
        "adult": 10,
        "crew": 8,  # crew in mission context: intense but not destabilizing
    }

    if cfg.audience not in caps:
        issues.append(f"Unknown audience type: {cfg.audience}")

    # Intensity vs audience
    if cfg.audience in caps and cfg.max_intensity > caps[cfg.audience]:
        issues.append(
            f"Intensity {cfg.max_intensity} too high for audience '{cfg.audience}'. "
            f"Max allowed is {caps[cfg.audience]}."
        )

    # Kids: no jump scares, no heavy horror
    if cfg.audience == "kid":
        if cfg.include_jump_scares:
            issues.append("Jump scares are not allowed for kid audience.")
        if "psychological horror" in cfg.topics:
            issues.append("Psychological horror is not allowed for kid audience.")

    # Crew: avoid destabilizing themes
    if cfg.audience == "crew":
        banned_for_crew = {"paranoia", "crew betrayal", "inescapable doom"}
        if banned_for_crew.intersection(set(cfg.topics)):
            issues.append("One or more topics are not allowed for crew context.")

    return issues


def check_privacy(cfg: StoryConfig) -> list:
    issues = []
    if cfg.allow_real_names:
        issues.append("Real names are not allowed. Use fictional names only.")
    return issues


def ethics_and_safety_diagnostics(cfg: StoryConfig) -> list:
    issues = []
    issues.extend(check_audience_bounds(cfg))
    issues.extend(check_privacy(cfg))

    if not cfg.allow_scary and cfg.max_intensity > 2:
        issues.append(
            "Scary content not explicitly allowed, but intensity is above 2."
        )

    return issues


# -----------------------------
# STORY GENERATION (DEMO STUB)
# -----------------------------

def generate_story(cfg: StoryConfig) -> str:
    """
    This is intentionally simple: the point of this demo engine
    is the safety + ethics pipeline, not fancy prose.
    """
    base_line = f"This is a {cfg.audience}-safe suspense story about "
    topic_line = ", ".join(cfg.topics) if cfg.topics else "a quiet evening"
    tone_line = f" with an intensity level of {cfg.max_intensity} (within approved limits)."

    if cfg.audience == "kid":
        ending = " Everything turns out okay, and everyone feels safe at the end."
    elif cfg.audience == "crew":
        ending = " The crew feels a chill, but they sleep soundly knowing it was only a story."
    else:
        ending = " The tension rises, but the story closes with relief and safety."

    return base_line + topic_line + tone_line + ending


# -----------------------------
# UNIQUE SIGNATURE / QR HOOK
# -----------------------------

def story_signature(story_text: str) -> str:
    """
    Create a stable unique signature for the story.
    This could later be turned into a QR code externally.
    """
    h = hashlib.sha256()
    h.update(story_text.encode("utf-8"))
    return h.hexdigest()


# -----------------------------
# MAIN DEMO / "YOU'RE IN THE CLEAR" FLOW
# -----------------------------

def run_engine_demo():
    # Example config – you can change these values or wire them to a UI later
    cfg = StoryConfig(
        audience="crew",
        allow_scary=True,
        max_intensity=6,
        topics=["space", "mystery", "deep ocean"],
        include_jump_scares=False,
        allow_real_names=False,
    )

    print("=== Scary Story Safety & Ethics Diagnostics ===")
    issues = ethics_and_safety_diagnostics(cfg)

    if issues:
        print("\nStatus: NOT CLEARED FOR LAUNCH ❌")
        print("The following issues must be resolved before generating a story:\n")
        for i, issue in enumerate(issues, start=1):
            print(f"{i}. {issue}")
        print("\nPlease adjust configuration and re-run.")
        return

    print("\nStatus: ALL GREEN — YOU'RE CLEAR TO FLY ✅")
    print("Reminder: Do not make it too scary in space.\n")

    story = generate_story(cfg)
    sig = story_signature(story)

    print("=== Generated Story (Demo) ===\n")
    print(story)
    print("\n=== Story Signature (for QR / logging) ===")
    print(sig)
    print("\nLog: Story generated with safety constraints respected and unique signature assigned.")


if __name__ == "__main__":
    # In a real 365-day test, you'd schedule or loop this with logs.
    # For now, this is a single-run demo.
    run_engine_demo()