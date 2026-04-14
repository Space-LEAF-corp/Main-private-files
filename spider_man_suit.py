import time
from typing import Optional, Any

try:
    import board  # pyright: ignore[reportMissingImports]
except ImportError:
    board = None
    print("Warning: board module not found. Board-related functionality will be unavailable.")

try:
    import busio  # pyright: ignore[reportMissingImports]
except ImportError:
    busio = None
    print("Warning: busio module not found. I2C communication will be unavailable.")

try:
    import adafruit_max30102  # pyright: ignore[reportMissingImports]
except ImportError:
    adafruit_max30102 = None
    print("Warning: adafruit_max30102 module not found. Heart rate/SpO2 sensor will be unavailable.")

import adafruit_ds18b20  # pyright: ignore[reportMissingImports]

try:
    import adafruit_mpu6050  # pyright: ignore[reportMissingImports]
except ImportError:
    adafruit_mpu6050 = None
    print("Warning: adafruit_mpu6050 module not found. Friction/impact detection will be unavailable.")

try:
    import adafruit_ssd1306  # pyright: ignore[reportMissingImports]
except ImportError:
    adafruit_ssd1306 = None
    print("Warning: adafruit_ssd1306 module not found. OLED HUD display will be unavailable.")

from adafruit_onewire.bus import OneWireBus  # type: ignore
from PIL import Image, ImageDraw, ImageFont  # type: ignore

# === SUIT MODE CONFIG ===
PEACEFUL_NEGOTIATION_MODE: bool = True  # Friendly neighborhood only. No killing.

# Initialize I2C bus and sensors
if busio is not None and board is not None:
    i2c: busio.I2C = busio.I2C(board.SCL, board.SDA)  # type: ignore
else:
    i2c = None
    print("Warning: I2C not initialized. Sensors depending on I2C will be unavailable.")

# Vital signs: Heart rate/SpO2
if adafruit_max30102 is not None and i2c is not None:
    pulse_sensor: Optional[Any] = adafruit_max30102.MAX30102(i2c)  # type: ignore
else:
    pulse_sensor = None

# Body temperature (connect DS18B20 to GPIO4)
if board is not None:
    ow_pin = getattr(board, "D4", None) or getattr(board, "GPIO4", None) or 4
else:
    ow_pin = 4
ow_bus: OneWireBus = OneWireBus(ow_pin)  # type: ignore
ds18 = adafruit_ds18b20.DS18B20(ow_bus, ow_bus.scan()[0])  # type: ignore

# Friction/impact: MPU6050 accelerometer
if adafruit_mpu6050 is not None and i2c is not None:
    mpu = adafruit_mpu6050.MPU6050(i2c)  # type: ignore
else:
    mpu = None

# HUD: OLED display (128x64, connected via I2C)
if adafruit_ssd1306 is not None and i2c is not None:
    oled: adafruit_ssd1306.SSD1306_I2C = adafruit_ssd1306.SSD1306_I2C(128, 64, i2c)  # type: ignore
else:
    oled = None
font: ImageFont.ImageFont = ImageFont.load_default()  # type: ignore

# Thresholds (Spider-Man tuning)
HEAT_THRESHOLD = 39.0   # Slightly lower—Spidey runs hot when stressed
ACCEL_THRESHOLD = 2.5   # g-force for swing/impact detection
STRESS_HEART_RATE = 120.0  # bpm

def read_vitals() -> tuple[float, float, float]:
    """Read heart rate, SpO2, and body temperature."""
    if pulse_sensor is not None:
        hr = pulse_sensor.heart_rate
        spo2 = pulse_sensor.spo2
        hr_val: float = float(hr) if hr is not None else 0.0
        spo2_val: float = float(spo2) if spo2 is not None else 0.0
    else:
        hr_val = 0.0
        spo2_val = 0.0

    temp_val: float = float(ds18.temperature) if ds18 is not None else 0.0  # type: ignore
    return hr_val, spo2_val, temp_val

def detect_friction() -> bool:
    """Detect high acceleration as proxy for swing/impact."""
    if mpu is None:
        return False
    accel: tuple[float, float, float] = mpu.acceleration  # type: ignore
    magnitude = (accel[0]**2 + accel[1]**2 + accel[2]**2)**0.5 / 9.81
    return magnitude > ACCEL_THRESHOLD

def compute_spidey_sense(hr: float, friction_alert: bool, temp: float) -> str:
    """
    Nonlethal 'Spidey-Sense' style guidance.
    Focus: awareness, de-escalation, and civilian safety.
    """
    if not PEACEFUL_NEGOTIATION_MODE:
        return "MODE ERROR"

    danger_score = 0

    if friction_alert:
        danger_score += 2
    if hr > STRESS_HEART_RATE:
        danger_score += 2
    if temp > HEAT_THRESHOLD:
        danger_score += 1

    if danger_score >= 4:
        return "HIGH TENSION: PROTECT + TALK"
    elif danger_score == 3:
        return "RISKY: WATCH + DE-ESCALATE"
    elif danger_score == 2:
        return "ALERT: STAY AWARE"
    elif danger_score == 1:
        return "MILD STRESS"
    else:
        return "ALL GOOD, NEIGHBOR"

def update_hud(hr: float, spo2: float, temp: float, friction_alert: bool):
    """Display data on OLED HUD with Spider-Man theming."""
    if oled is None:
        return

    image: Image.Image = Image.new("1", (oled.width, oled.height))  # type: ignore
    draw = ImageDraw.Draw(image)

    # Header
    draw.text((0, 0), "SPIDER-MAN HUD", font=font, fill=255)

    # Vital signs
    draw.text((0, 12), f"HR: {hr:.1f} bpm", font=font, fill=255)
    draw.text((0, 22), f"SpO2: {spo2:.1f}%", font=font, fill=255)
    draw.text((0, 32), f"Temp: {temp:.1f} C", font=font, fill=255)

    # Alerts
    y_alert = 42
    if temp > HEAT_THRESHOLD:
        draw.text((0, y_alert), "OVERHEAT RISK", font=font, fill=255)
        y_alert += 10
    if friction_alert:
        draw.text((0, y_alert), "SWING/IMPACT SPIKE", font=font, fill=255)
        y_alert += 10

    # Peaceful negotiation / Spidey-Sense
    sense = compute_spidey_sense(hr, friction_alert, temp)
    mode_label = "PEACEFUL MODE" if PEACEFUL_NEGOTIATION_MODE else "UNKNOWN MODE"
    draw.text((0, 54), f"{mode_label}: {sense}", font=font, fill=255)

    oled.image(image)  # type: ignore
    oled.show()        # type: ignore

# Main loop
if __name__ == "__main__":
    while True:
        hr, spo2, temp = read_vitals()
        friction_alert = detect_friction()

        print(
            f"[SPIDER-MAN] HR: {hr:.1f}, SpO2: {spo2:.1f}, "
            f"Temp: {temp:.1f}, Friction: {friction_alert}, "
            f"PeacefulMode: {PEACEFUL_NEGOTIATION_MODE}"
        )

        update_hud(hr, spo2, temp, friction_alert)
        time.sleep(1)


---
