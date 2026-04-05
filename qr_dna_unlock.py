# qr_dna_unlock.py
import qrcode
import json
import sys

# Your original data (with unicode dashes normalized to regular ones for safety)
data = {
    "user_id": "anonymous",
    "protocol": "qr_dna_unlock_v1",
    "transitions": [
        {"id": "T0", "from": "S0 - Idle", "to": "S1 - Awaiting QR Seal", "trigger": "User selects capsule", "conditions": "Capsule exists & not locked down", "actions": "Show 'Present QR Seal'", "log_schema": None},
        # ... copy-paste all other transitions here and replace — with - ...
        # (I replaced em-dashes with regular hyphens for widest compatibility)
        {"id": "T13", "from": "S9 - Lockdown", "to": "S0 - Idle", "trigger": "Lockdown cleared", "conditions": "Override or timer", "actions": "Reset counters; remain high-tension if desired", "log_schema": None}
    ]
}

# Minified JSON – smallest possible
json_payload = json.dumps(data, separators=(',', ':'), ensure_ascii=True)

print(f"Payload length: {len(json_payload)} characters")

# Try to create QR with highest error correction that still fits
qr = qrcode.QRCode(
    version=None,               # auto
    error_correction=qrcode.constants.ERROR_CORRECT_H,
    box_size=10,
    border=4,
)

qr.add_data(json_payload)

try:
    qr.make(fit=True)
    print(f"Created QR version {qr.version} with error correction H")
except Exception as e:
    print("High correction failed → trying medium correction...")
    qr = qrcode.QRCode(
        version=None,
        error_correction=qrcode.constants.ERROR_CORRECT_M,
        box_size=10,
        border=4,
    )
    qr.add_data(json_payload)
    qr.make(fit=True)
    print(f"Fallback: version {qr.version} with error correction M")

# Generate image
img = qr.make_image(fill_color="black", back_color="white")
img.save("qr_dna_unlock_v1.png")
print("Saved → qr_dna_unlock_v1.png")

# Optional: show in terminal if you have pillow + terminal image support
# img.show()