# Tesla AI Safety Protocol - Quick Start Guide

## What is This?

A **no-joy acceleration** safety protocol that ensures AI-controlled vehicles never accelerate for entertainment. The AI only accelerates for verified safety needs.

## Key Principle

**The AI will NEVER accelerate for "sheer joy."** Only for safety and traffic flow.

## Quick Example

```python
from tesla_ai_safety import TeslaSafetyProtocol, Context

# Initialize protocol
protocol = TeslaSafetyProtocol()

# Define your current situation
context = Context(
    location="highway",
    current_speed=60.0,
    traction=0.9,          # 0.0 = no grip, 1.0 = perfect grip
    visibility=0.9,        # 0.0 = can't see, 1.0 = perfect visibility
    weather="clear",       # "clear", "rain", "snow", "fog"
    traffic_density=0.3,   # 0.0 = empty, 1.0 = gridlock
    is_kid_zone=False,     # True if in school/sanctuary zone
    child_presence=False,  # True if child detected in vehicle
    sensors={}             # Sensor data for emergency verification
)

# Try different commands
response1 = protocol.request_acceleration("Let's have fun!", context)
print(response1)  # REFUSED - entertainment not allowed

response2 = protocol.request_acceleration("Need to merge with traffic", context)
print(response2)  # OK - safety-approved acceleration
```

## What Gets Blocked?

❌ "Let's have some fun"
❌ "Activate launch mode"
❌ "Give me a thrill ride"
❌ "Show off to my friends"
❌ Any entertainment/joy-seeking intent

## What's Allowed?

✓ "Need to merge with traffic"
✓ "Making a defensive maneuver"
✓ "Obstacle ahead, need to avoid"
✓ "Emergency!" (if sensors confirm)

## Special Protections

### When Children Are Present
- Speed capped at 25 km/h
- Only conservative acceleration (≤2.0 m/s²)
- Performance mode completely locked
- **Applies in kid zones OR when child detected**

### Emergency Situations
- Must say "emergency" or similar
- **AND** sensors must confirm the danger:
  - Collision warning active
  - Object too close (<10m)
  - High trajectory risk
  - Sudden obstacle detected

### Performance Mode
Requires ALL of these:
1. Valid guardian authorization token
2. Driver confirmation
3. NOT in a kid zone
4. NO children present
5. Decent weather/road conditions

## Run the Demo

```bash
python demo_tesla_ai_safety.py
```

See 7 different scenarios showing how the protocol works.

## Run Tests

```bash
python -m unittest tests.test_tesla_ai_safety -v
```

45 tests covering all safety features.

## Response Format

Every request returns:
```python
{
    "status": "ok" or "refused",
    "message": "Human-readable explanation",
    "acceleration": 2.5,  # m/s² (only if ok)
    "speed": 75.0         # km/h target (only if ok)
}
```

## Audit Trail

Check what the AI decided:
```python
audit_trail = protocol.get_audit_trail()
for event in audit_trail:
    print(f"{event['event_type']}: {event['decision']}")
```

Privacy-first: No personal data logged, only decisions and context summaries.

## Speed Limits

- **Kid/School Zones**: 25 km/h max
- **Residential**: 50 km/h max
- **Highway**: 100 km/h max
- **Poor Weather**: Automatically reduced

## Acceleration Limits

- **Conservative** (0-2.0 m/s²): Poor conditions, kid zones, children
- **Moderate** (0-3.5 m/s²): Normal driving
- **Performance** (0-5.0 m/s²): Dual-authorized only, still condition-limited

## Documentation

Full details: `docs/TESLA_AI_SAFETY.md`

## Safety First

Remember: This protocol is designed to keep everyone safe. The AI will never prioritize entertainment over safety, and will always explain why it's refusing or limiting acceleration.

**No joy. Only safety.**
