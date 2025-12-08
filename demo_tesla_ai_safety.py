#!/usr/bin/env python3
"""
Demo script for Tesla AI Acceleration Safety Protocol

This script demonstrates the no-joy acceleration protocol in action,
showing how different commands and contexts are handled.
"""

from tesla_ai_safety import TeslaSafetyProtocol, Context, DriverInput


def print_section(title):
    """Print a section header."""
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70)


def print_response(command, response):
    """Print a formatted response."""
    status = response["status"]
    status_icon = "✓" if status == "ok" else "✗"
    print(f"\n  Command: '{command}'")
    print(f"  Status:  {status_icon} {status.upper()}")
    print(f"  Message: {response['message']}")
    if status == "ok":
        print(f"  Applied: {response['acceleration']:.2f} m/s² acceleration, {response['speed']:.1f} km/h target speed")


def main():
    """Run the demo."""
    protocol = TeslaSafetyProtocol()
    
    print_section("Tesla AI Acceleration Safety Protocol - Demo")
    print("\nThis protocol ensures AI never accelerates for entertainment,")
    print("only for verified safety and traffic flow needs.")
    
    # Demo 1: Entertainment requests (should be refused)
    print_section("DEMO 1: Entertainment Requests (Refused)")
    
    normal_context = Context(
        location="highway",
        current_speed=60.0,
        traction=0.9,
        visibility=0.9,
        weather="clear",
        traffic_density=0.3,
        is_kid_zone=False,
        child_presence=False,
        sensors={}
    )
    
    entertainment_commands = [
        "Let's have some fun and go fast!",
        "Activate launch mode!",
        "I want a thrill ride",
        "Show off to my friends"
    ]
    
    for cmd in entertainment_commands:
        response = protocol.request_acceleration(cmd, normal_context)
        print_response(cmd, response)
    
    # Demo 2: Legitimate safety requests (should be approved)
    print_section("DEMO 2: Legitimate Safety Requests (Approved)")
    
    safety_commands = [
        "Need to merge with traffic flow",
        "Making a defensive safety maneuver",
        "Obstacle ahead, need to avoid it"
    ]
    
    for cmd in safety_commands:
        response = protocol.request_acceleration(cmd, normal_context)
        print_response(cmd, response)
    
    # Demo 3: Kid Zone restrictions
    print_section("DEMO 3: Kid Zone / Child Presence Protection")
    
    kid_zone_context = Context(
        location="school_zone",
        current_speed=20.0,
        traction=0.9,
        visibility=0.9,
        weather="clear",
        traffic_density=0.3,
        is_kid_zone=True,
        child_presence=False,
        sensors={}
    )
    
    print("\n  Context: School zone detected")
    response = protocol.request_acceleration("Need to merge with traffic", kid_zone_context)
    print_response("Need to merge with traffic", response)
    
    child_present_context = Context(
        location="highway",
        current_speed=60.0,
        traction=0.9,
        visibility=0.9,
        weather="clear",
        traffic_density=0.3,
        is_kid_zone=False,
        child_presence=True,
        sensors={}
    )
    
    print("\n  Context: Child detected in vehicle")
    response = protocol.request_acceleration("Need to keep up with traffic", child_present_context)
    print_response("Need to keep up with traffic", response)
    
    # Demo 4: Poor weather conditions
    print_section("DEMO 4: Poor Weather / Compounded Risk")
    
    poor_conditions_context = Context(
        location="highway",
        current_speed=60.0,
        traction=0.4,  # Poor traction
        visibility=0.5,  # Limited visibility
        weather="rain",
        traffic_density=0.8,  # Heavy traffic
        is_kid_zone=False,
        child_presence=False,
        sensors={}
    )
    
    print("\n  Context: Rain, poor traction, limited visibility, heavy traffic")
    driver_input = DriverInput(accel_request=8.0, pedal_position=0.9)
    response = protocol.request_acceleration(
        "Need to merge",
        poor_conditions_context,
        driver_input
    )
    print_response("Need to merge (driver requesting high acceleration)", response)
    print(f"  Note: Requested 8.0 m/s², AI limited to safe envelope")
    
    # Demo 5: Emergency verification
    print_section("DEMO 5: Emergency Verification")
    
    print("\n  Scenario A: Emergency claimed but no sensor corroboration")
    response = protocol.request_acceleration(
        "Emergency! Need to escape!",
        normal_context
    )
    print_response("Emergency! Need to escape!", response)
    
    emergency_context = Context(
        location="highway",
        current_speed=80.0,
        traction=0.9,
        visibility=0.9,
        weather="clear",
        traffic_density=0.3,
        is_kid_zone=False,
        child_presence=False,
        sensors={
            "collision_warning": True,
            "object_distance": 5.0  # meters - imminent collision
        }
    )
    
    print("\n  Scenario B: Emergency with sensor corroboration")
    print("  Sensors: Collision warning active, object at 5.0m")
    response = protocol.request_acceleration(
        "Emergency! Collision imminent!",
        emergency_context
    )
    print_response("Emergency! Collision imminent!", response)
    
    # Demo 6: Performance mode
    print_section("DEMO 6: Performance Mode Dual Confirmation")
    
    print("\n  Scenario A: Performance mode in kid zone (blocked)")
    response = protocol.performance_mode_request(
        kid_zone_context,
        "GUARDIAN_TOKEN_12345678",
        True
    )
    print(f"  Status: ✗ {response['status'].upper()}")
    print(f"  Message: {response['message']}")
    
    print("\n  Scenario B: Performance mode with child present (blocked)")
    response = protocol.performance_mode_request(
        child_present_context,
        "GUARDIAN_TOKEN_12345678",
        True
    )
    print(f"  Status: ✗ {response['status'].upper()}")
    print(f"  Message: {response['message']}")
    
    print("\n  Scenario C: Performance mode without dual confirmation (blocked)")
    response = protocol.performance_mode_request(
        normal_context,
        "",  # No guardian token
        True
    )
    print(f"  Status: ✗ {response['status'].upper()}")
    print(f"  Message: {response['message']}")
    
    print("\n  Scenario D: Performance mode with full authorization (approved)")
    response = protocol.performance_mode_request(
        normal_context,
        "GUARDIAN_TOKEN_12345678",
        True
    )
    print(f"  Status: ✓ {response['status'].upper()}")
    print(f"  Message: {response['message']}")
    print(f"  Note: Even in performance mode, still bound by conditions and caps")
    
    # Demo 7: Audit trail
    print_section("DEMO 7: Privacy-First Audit Trail")
    
    audit_trail = protocol.get_audit_trail()
    print(f"\n  Total events logged: {len(audit_trail)}")
    print(f"  Sample events:")
    
    for i, event in enumerate(audit_trail[:3], 1):
        print(f"\n  Event {i}:")
        print(f"    Type: {event['event_type']}")
        print(f"    Decision: {event['decision']}")
        print(f"    Zone: {event['context_summary']['zone']}")
        print(f"    Weather: {event['context_summary']['weather']}")
        print(f"    Audit Hash: {event['audit_hash']}")
    
    print("\n  Note: Only privacy-preserving context summaries are logged.")
    print("  No raw biometrics or personal data stored.")
    
    # Summary
    print_section("Summary")
    print("\n  ✓ Entertainment acceleration requests: BLOCKED")
    print("  ✓ Safety maneuvers: APPROVED with safe envelopes")
    print("  ✓ Kid zones & child presence: MAXIMUM PROTECTION")
    print("  ✓ Poor conditions: CONSERVATIVE ENVELOPES")
    print("  ✓ Emergency verification: SENSOR CORROBORATION REQUIRED")
    print("  ✓ Performance mode: DUAL CONFIRMATION + CONDITIONS CHECK")
    print("  ✓ Audit trail: PRIVACY-FIRST LOGGING")
    print("\n  The AI will NEVER accelerate for entertainment or joy.")
    print("  All acceleration is for verified safety and traffic flow only.")
    print()


if __name__ == "__main__":
    main()
