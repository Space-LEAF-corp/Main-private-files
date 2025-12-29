# Morphine Grid Prototype - Test Suite

from morphine_grid import SecurityFramework

def run_tests():
    framework = SecurityFramework(perimeter_km=250)

    # Test 1: Basic activation
    report1 = framework.activate_protocol(user_id="Captain001", environment="Orbital Station")
    print("Test 1 - Basic Activation", report1)

    # Test 2: ICU Glasses disabled
    framework.glasses.alerts_enabled = False
    report2 = framework.activate_protocol(user_id="Captain002", environment="Training Grounds")
    print("Test 2 - ICU Glasses Disabled", report2)

    # Test 3: Geofence inactive
    framework.geofence.active = False
    report3 = framework.activate_protocol(user_id="Captain003", environment="Deep Space")
    print("Test 3 - Geofence Inactive", report3)

    # Test 4: Multiple users sequentially
    users = ["Alpha", "Bravo", "Charlie"]
    for user in users:
        report = framework.activate_protocol(user_id=user, environment="Simulation Chamber")
        print(f"Test 4 - User {user}", report)

if __name__ == "__main__":
    run_tests()
