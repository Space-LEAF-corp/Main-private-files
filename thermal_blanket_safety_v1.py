"""
Thermal Blanket Safety Logic — Public Demo Build v1.0
Author: Leif William Sogge
Description:
    A simple, modular safety system for heated blankets during extreme cold.
    This public version demonstrates:
        - temperature monitoring
        - auto‑shutoff logic
        - overheat protection
        - cooldown cycles
        - user‑safe fallback states

    This is NOT a production system. It is a conceptual prototype intended
    for public review, education, and early-stage engineering reference.
"""

import time
from dataclasses import dataclass

# -----------------------------
# Configuration Parameters
# -----------------------------

@dataclass
class ThermalConfig:
    MIN_SAFE_TEMP: float = -40.0     # extreme cold threshold (°C)
    MAX_SAFE_TEMP: float = 55.0      # blanket surface limit (°C)
    TARGET_TEMP: float = 42.0        # ideal warming temperature (°C)
    COOLDOWN_TEMP: float = 38.0      # cooldown threshold (°C)
    CHECK_INTERVAL: float = 1.0      # seconds between sensor checks


# -----------------------------
# Blanket Controller
# -----------------------------

class HeatedBlanketController:
    def __init__(self, config: ThermalConfig):
        self.config = config
        self.heating = False
        self.last_temp = None

    # Simulated sensor read
    def read_temperature(self) -> float:
        # Placeholder for real hardware sensor
        raise NotImplementedError("Connect to actual temperature sensor")

    # Simulated heating element control
    def set_heating(self, state: bool):
        # Placeholder for GPIO or hardware relay
        self.heating = state
        print(f"[SYSTEM] Heating set to: {state}")

    # Core safety logic
    def safety_check(self, temp: float):
        self.last_temp = temp

        # Overheat protection
        if temp >= self.config.MAX_SAFE_TEMP:
            print("[ALERT] Overheat detected — shutting down heating element")
            self.set_heating(False)
            return

        # Maintain target temperature
        if temp < self.config.TARGET_TEMP and not self.heating:
            print("[INFO] Temperature low — activating heating")
            self.set_heating(True)

        # Cooldown cycle
        if temp > self.config.COOLDOWN_TEMP and self.heating:
            print("[INFO] Temperature above cooldown threshold — reducing heat")
            self.set_heating(False)

    # Main loop
    def run(self):
        print("[SYSTEM] Thermal Blanket Safety Logic Active")
        while True:
            try:
                temp = self.read_temperature()
                print(f"[SENSOR] Current temperature: {temp}°C")
                self.safety_check(temp)
                time.sleep(self.config.CHECK_INTERVAL)

            except KeyboardInterrupt:
                print("\n[SYSTEM] Manual shutdown requested")
                self.set_heating(False)
                break

            except Exception as e:
                print(f"[ERROR] Sensor or system failure: {e}")
                self.set_heating(False)
                break


# -----------------------------
# Entry Point (for engineers)
# -----------------------------

if __name__ == "__main__":
    config = ThermalConfig()
    controller = HeatedBlanketController(config)
    controller.run()
