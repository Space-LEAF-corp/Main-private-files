# qr_full_protocol_low_ecc.py
import qrcode
import json

# Your exact structure (em-dashes and curly quotes preserved)
data = {
    "user_id": "anonymous",
    "protocol": "qr_dna_unlock_v1",
    "transitions": [
        {"id":"T0","from":"S0 — Idle","to":"S1 — Awaiting QR Seal","trigger":"User selects capsule","conditions":"Capsule exists & not locked down","actions":"Show “Present QR Seal”","log_schema":None},
        {"id":"T1","from":"S1 — Awaiting QR Seal","to":"S2 — Verifying QR Seal","trigger":"QR scanned","conditions":"QR payload decodable","actions":"Decode, parse, start verification","log_schema":None},
        {"id":"T2","from":"S2 — Verifying QR Seal","to":"S4 — Awaiting DNA Seal","trigger":"QR verification success","conditions":"sig, time, purpose, nonce, cid all valid","actions":"Mark nonce pending; show “Present DNA Seal”","log_schema":None},
        {"id":"T3","from":"S2 — Verifying QR Seal","to":"S3 — QR Seal Failed","trigger":"QR verification failure","conditions":"Any validation check fails","actions":"Increment failure counter; Frett tension++","log_schema":None},
        {"id":"T4","from":"S3 — QR Seal Failed","to":"S1 or S9","trigger":"Retry or threshold reached","conditions":"failure_count < threshold → S1; else → S9","actions":"Either allow retry or enter Lockdown","log_schema":None},
        {"id":"T5","from":"S4 — Awaiting DNA Seal","to":"S5 — Verifying DNA Seal","trigger":"DNA input started","conditions":"DNA mode configured","actions":"Capture biometric / gesture / passphrase","log_schema":None},
        {"id":"T6","from":"S5 — Verifying DNA Seal","to":"S7 — Capsule Unlocked","trigger":"DNA verification success","conditions":"DNA matches configured mode","actions":"Mark nonce used; decrypt; open read‐only; start timer","log_schema":None},
        {"id":"T7","from":"S5 — Verifying DNA Seal","to":"S6 — DNA Seal Failed","trigger":"DNA verification failure","conditions":"DNA mismatch","actions":"Increment DNA failure counter; Frett tension++","log_schema":None},
        {"id":"T8","from":"S6 — DNA Seal Failed","to":"S4 or S9","trigger":"Retry or threshold reached","conditions":"failure_count < threshold → S4; else → S9","actions":"Either allow retry or enter Lockdown","log_schema":None},
        {"id":"T9","from":"S7 — Capsule Unlocked","to":"S8 — Auto‐Relock","trigger":"Unlock timer expires","conditions":"unlock_window elapsed","actions":"Wipe temp key; close view","log_schema":None},
        {"id":"T10","from":"S8 — Auto‐Relock","to":"S0 — Idle","trigger":"Relock complete","conditions":"All cleanup done","actions":"Return to idle","log_schema":None},
        {"id":"T11","from":"S7 — Capsule Unlocked","to":"S8 — Auto‐Relock","trigger":"Manual relock","conditions":"User triggers close","actions":"Same as auto‐relock, user‐initiated","log_schema":None},
        {"id":"T12","from":"S3/S6/Any","to":"S9 — Lockdown","trigger":"Lockdown trigger","conditions":"Thresholds or Frett global alert","actions":"Freeze unlocks; require override","log_schema":None},
        {"id":"T13","from":"S9 — Lockdown","to":"S0 — Idle","trigger":"Lockdown cleared","conditions":"Override or timer","actions":"Reset counters; remain high‐tension if desired","log_schema":None}
    ]
}

payload = json.dumps(data, separators=(',', ':'), ensure_ascii=False)
print(f"Length: {len(payload)} characters")

qr = qrcode.QRCode(
    version=40,                     # forced maximum size
    error_correction=qrcode.constants.ERROR_CORRECT_L,  # lowest = biggest capacity
    box_size=10,
    border=4,
)

qr.add_data(payload)
qr.make(fit=False)  # we already know it fits version 40-L

img = qr.make_image(fill_color="black", back_color="white")
img.save("qr_dna_unlock_v1_FULL_LOW_ECC.png")
print("Saved as: qr_dna_unlock_v1_FULL_LOW_ECC.png")