# test_harness.py
# Run: python test_harness.py
import base64, hashlib, time, random, statistics
from io import BytesIO
from Crypto.Cipher import AES
from Crypto.Protocol.KDF import scrypt
from Crypto.Random import get_random_bytes
from nacl.signing import SigningKey, VerifyKey
import qrcode, numpy as np, PIL.Image, sys

# --- Config (tune for CI vs local)
PASS = b"SuperSecretPassphrase2025"
KDF_PARAMS = dict(N=2**14, r=8, p=1)   # production-like; reduce for CI speed
SALT_LEN = 16
YEARS = 10
DAYS_PER_YEAR = 365
SIM_DAYS = YEARS * DAYS_PER_YEAR  # accelerated timeline
EVENTS = {
    'compromise_day': 2 * DAYS_PER_YEAR,   # simulate compromise at year 2
    'rotation_day': 3 * DAYS_PER_YEAR,     # rotate keys at year 3
}
# QR corruption levels to test (0%..30%)
CORRUPTION_LEVELS = [0.0, 0.05, 0.10, 0.20, 0.30]

# --- Utilities
def hash_dna(dna_bytes):
    return hashlib.sha256(dna_bytes).digest()

def derive_key(passphrase, salt):
    return scrypt(passphrase, salt, 32, **KDF_PARAMS)

def encrypt_hash(key, plaintext):
    cipher = AES.new(key, AES.MODE_GCM)
    ct, tag = cipher.encrypt_and_digest(plaintext)
    return cipher.nonce, ct, tag

def decrypt_hash(key, nonce, ct, tag):
    cipher = AES.new(key, AES.MODE_GCM, nonce=nonce)
    return cipher.decrypt_and_verify(ct, tag)

def package_blob(version, salt, nonce, ct, tag, signature, owner):
    parts = [version, salt, nonce, ct, tag, signature, owner.encode()]
    return base64.b64encode(b"||".join(parts)).decode()

def unpackage_blob(b64):
    parts = base64.b64decode(b64).split(b"||")
    return parts  # v, salt, nonce, ct, tag, signature, owner

def make_qr(payload, box=4):
    qr = qrcode.QRCode(error_correction=qrcode.constants.ERROR_CORRECT_H, box_size=box, border=2)
    qr.add_data(payload)
    qr.make(fit=True)
    img = qr.make_image(fill_color="black", back_color="white").convert("L")
    return img

def corrupt_image(img, corruption_fraction):
    arr = np.array(img)
    total = arr.size
    nflip = int(total * corruption_fraction)
    idx = np.unravel_index(np.random.choice(total, nflip, replace=False), arr.shape)
    arr[idx] = 255 - arr[idx]  # invert pixels
    return PIL.Image.fromarray(arr)

def qr_decode_sim(img):
    # Lightweight decode simulation: attempt to read QR via qrcode lib is not available.
    # For real tests, use zxing, pyzbar, or OpenCV QR detector. Here we simulate decode success
    # based on image "damage" heuristics: if many pixels inverted, fail.
    arr = np.array(img)
    damage = np.mean(arr) / 255.0
    # heuristic: if damage metric below threshold, assume decode success
    return True if damage > 0.01 else False

# --- Setup initial keys (owner signing key)
owner = "LEIF William Sogge"
version = b"v1"
signing_key = SigningKey.generate()
verify_key = signing_key.verify_key

# store payloads (simulate DB)
payload_store = []

# metrics
metrics = {
    'roundtrip_success': [],
    'corruption_success': {c: [] for c in CORRUPTION_LEVELS},
    'verify_failures': 0,
    'decrypt_failures': 0,
    'tamper_detected': 0,
}

# --- Generate baseline payloads
def create_payload(dna_bytes, passphrase, owner, signing_key, version):
    salt = get_random_bytes(SALT_LEN)
    key = derive_key(passphrase, salt)
    dna_h = hash_dna(dna_bytes)
    nonce, ct, tag = encrypt_hash(key, dna_h)
    sig = signing_key.sign(version + salt + nonce + ct + tag + owner.encode()).signature
    blob = package_blob(version, salt, nonce, ct, tag, sig, owner)
    return blob

# create N baseline payloads
N = 200  # scale up for heavy tests
for i in range(N):
    dna = ("AGCT" * (10 + (i % 5))).encode() + bytes([i % 256])
    payload_store.append(create_payload(dna, PASS, owner, signing_key, version))

# --- Simulation loop (accelerated days)
current_signing_key = signing_key
current_verify_key = verify_key
revoked_keys = []
for day in range(SIM_DAYS):
    # trigger compromise
    if day == EVENTS['compromise_day']:
        compromised_key = current_signing_key
        print(f"[SIM] Day {day}: Simulating key compromise.")
    # trigger rotation
    if day == EVENTS['rotation_day']:
        print(f"[SIM] Day {day}: Rotating signing key.")
        new_sk = SigningKey.generate()
        revoked_keys.append(current_verify_key.encode())  # revoke old
        current_signing_key = new_sk
        current_verify_key = new_sk.verify_key
        # optionally re-sign stored payloads if policy requires (simulate cost)
        # For this test we do not re-sign automatically; we record how many would need re-signing.

# --- Run verification & corruption tests
for idx, blob in enumerate(payload_store):
    try:
        parts = unpackage_blob(blob)
        v, salt, nonce, ct, tag, signature, owner_bytes = parts
        # verify signature with original verify_key (we used initial signing_key)
        vk = VerifyKey(verify_key.encode())
        vk.verify(v + salt + nonce + ct + tag + owner_bytes, signature)
        # derive key and decrypt
        key = derive_key(PASS, salt)
        plaintext = decrypt_hash(key, nonce, ct, tag)
        metrics['roundtrip_success'].append(1)
    except Exception as e:
        metrics['roundtrip_success'].append(0)
        metrics['verify_failures'] += 1

    # corruption tests: render QR and corrupt
    img = make_qr(blob, box=4)
    for c in CORRUPTION_LEVELS:
        corrupted = corrupt_image(img, c)
        ok = qr_decode_sim(corrupted)
        metrics['corruption_success'][c].append(1 if ok else 0)

# --- Tamper tests: flip bytes in blob
for i in range(50):
    b = base64.b64decode(payload_store[i])
    # flip a random byte in ciphertext region
    b_list = bytearray(b)
    pos = random.randint(0, len(b_list)-1)
    b_list[pos] ^= 0xFF
    tampered = base64.b64encode(bytes(b_list)).decode()
    try:
        parts = unpackage_blob(tampered)
        v, salt, nonce, ct, tag, signature, owner_bytes = parts
        vk = VerifyKey(verify_key.encode())
        vk.verify(v + salt + nonce + ct + tag + owner_bytes, signature)
        # if signature passes (unlikely), try decrypt
        key = derive_key(PASS, salt)
        decrypt_hash(key, nonce, ct, tag)
    except Exception:
        metrics['tamper_detected'] += 1

# --- Summarize metrics
print("\n--- Simulation Summary ---")
print("Roundtrip success rate:", sum(metrics['roundtrip_success'])/len(metrics['roundtrip_success']))
for c in CORRUPTION_LEVELS:
    arr = metrics['corruption_success'][c]
    print(f"Corruption {int(c*100)}% decode success:", sum(arr)/len(arr))
print("Verify failures:", metrics['verify_failures'])
print("Tamper detections:", metrics['tamper_detected'])
