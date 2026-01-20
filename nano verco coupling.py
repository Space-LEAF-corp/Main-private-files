import time
import random  # For simulating randomness in checks
from threading import Thread  # For background monitoring

class NanoVelcroCouplingSystem:
    def __init__(self, user_type='adult'):
        self.engaged = False
        self.user_type = user_type  # 'adult' or 'kid' for safety adjustments
        self.gyro_firewall_active = False
        self.encryption_key = self._generate_encryption_key()
        self.monitor_thread = None
        self.support_logs = []
        print(f"System initialized for {user_type} user.")

    def _generate_encryption_key(self):
        # Simulate backward encryption (reversed random string)
        key = ''.join(random.choice('abcdefghijklmnopqrstuvwxyz0123456789') for _ in range(16))
        return key[::-1]  # Reverse for "backward" encryption

    def engage_coupling(self):
        if self.engaged:
            print("Coupling already engaged.")
            return
        print("Engaging nano-velcro coupling...")
        time.sleep(1)  # Simulate engagement time
        self.engaged = True
        self._start_gyro_firewall()
        self._start_monitoring()
        print("Coupling engaged. Stable in space.")

    def disengage_with_button(self):
        if not self.engaged:
            print("Coupling not engaged.")
            return
        print("Button pressed: Disengaging nano-velcro coupling...")
        if self._perform_safety_checks():
            self.engaged = False
            self._stop_gyro_firewall()
            self._stop_monitoring()
            print("Coupling safely disengaged.")
        else:
            print("Disengagement failed: Safety checks not passed.")

    def _start_gyro_firewall(self):
        self.gyro_firewall_active = True
        print("Rotating gyroscopic firewall activated with backward encryption.")
        # Simulate rotation/encryption validation
        encrypted_data = self._encrypt_data("system_status: stable")
        print(f"Encrypted sample: {encrypted_data}")

    def _stop_gyro_firewall(self):
        self.gyro_firewall_active = False
        print("Gyroscopic firewall deactivated.")

    def _encrypt_data(self, data):
        # Simple mock backward encryption with key
        return ''.join(chr(ord(c) ^ ord(self.encryption_key[i % len(self.encryption_key)])) for i, c in enumerate(data))[::-1]

    def _perform_safety_checks(self):
        print("Running dual cross-validation and safety checks...")
        check1 = self._gyro_check()
        check2 = self._encryption_validation()
        if check1 and check2:
            print("All checks passed. Go for disengagement.")
            return True
        else:
            print("Checks failed. No-go: System remains engaged for safety.")
            return False

    def _gyro_check(self):
        # Simulate gyro stability check (random for demo)
        stability = random.uniform(0, 1)
        if stability > 0.2:  # Arbitrary threshold, lower for kids
            if self.user_type == 'kid':
                threshold = 0.1  # Extra safety for kids
                return stability > threshold
            return True
        return False

    def _encryption_validation(self):
        # Dual cross-validation: Encrypt and decrypt to verify
        test_data = "test_payload"
        encrypted = self._encrypt_data(test_data)
        decrypted = self._encrypt_data(encrypted)  # XOR reverses itself
        return decrypted == test_data

    def _start_monitoring(self):
        self.monitor_thread = Thread(target=self._monitor_support, daemon=True)
        self.monitor_thread.start()
        print("Support monitoring started.")

    def _stop_monitoring(self):
        # Thread will stop when engaged=False, but for simplicity, just log
        print("Support monitoring stopped.")
        print("Support logs:", self.support_logs)

    def _monitor_support(self):
        while self.engaged:
            status = "Stable" if random.random() > 0.1 else "Alert: Instability detected!"
            self.support_logs.append(status)
            print(f"Monitor: {status}")
            time.sleep(5)  # Check every 5 seconds

# Example usage:
system = NanoVelcroCouplingSystem(user_type='kid')  # Or 'adult'
system.engage_coupling()
time.sleep(10)  # Simulate time in space
system.disengage_with_button()