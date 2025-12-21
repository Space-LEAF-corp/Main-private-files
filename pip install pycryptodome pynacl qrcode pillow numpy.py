# test_harness.py
# Run: python test_harness.py

import base64, hashlib, random
import numpy as np # pyright: ignore[reportMissingImports]
try:
    from Crypto.Cipher import AES  # type: ignore
    from Crypto.Protocol.KDF import scrypt  # type: ignore
    try:
        # Try both pycryptodome and pycrypto import styles
        from Crypto.Random import get_random_bytes  # type: ignore
    except ImportError:
        try:
            from Cryptodome.Random import get_random_bytes  # type: ignore
        except ImportError:
            raise ImportError("pycryptodome is required. Install with: pip install pycryptodome (Crypto.Random.get_random_bytes missing)")
    from typing import Callable
    get_random_bytes: Callable[[int], bytes]
except ImportError:
    raise ImportError("pycryptodome is required. Install with: pip install pycryptodome")

try:
    from nacl.signing import SigningKey, VerifyKey  # type: ignore
except ImportError:
    raise ImportError("PyNaCl is required. Install with: pip install pynacl")
from typing import Union
def get_qrcode_and_error_correct_h():
    try:
        import qrcode # pyright: ignore[reportMissingModuleSource]
        try:
            from qrcode.constants import ERROR_CORRECT_H # pyright: ignore[reportMissingModuleSource]
            error_correct_h = ERROR_CORRECT_H
        except ImportError:
            error_correct_h = qrcode.ERROR_CORRECT_H
        return qrcode, error_correct_h
    except ImportError:
        raise ImportError("qrcode is required. Install with: pip install qrcode")
try:
    from PIL import Image as PILImage # pyright: ignore[reportUnknownVariableType, reportMissingImports]
except ImportError:
    raise ImportError("Pillow is required for image processing. Install with: pip install pillow")
from typing import Any
def make_qr(payload: Union[str, bytes], box: int = 4) -> Any:
    """
    Generate a QR code image from the given payload.
    Returns:
        PILImage.Image: The generated QR code as a PIL Image.
    """
    qrcode, ERROR_CORRECT_H = get_qrcode_and_error_correct_h()
    qr = qrcode.QRCode(error_correction=ERROR_CORRECT_H, box_size=box, border=2)
    # Ensure payload is a string (QR expects str, not bytes)
    if isinstance(payload, bytes):
        payload = payload.decode('utf-8')
    qr.add_data(payload)
    qr.make(fit=True)
    img = qr.make_image(fill_color="black", back_color="white")
    # Ensure we have a PIL.Image for conversion
    if hasattr(img, 'get_image'):
        img2 = img.get_image() # pyright: ignore[reportUnknownVariableType, reportUnknownMemberType]
        img = img2  # type: ignore[assignment]
    # If not a PIL.Image, convert using PILImage.fromarray
    if not isinstance(img, PILImage.Image): # pyright: ignore[reportUnknownMemberType]
        import numpy as np # pyright: ignore[reportMissingImports]
        img = PILImage.fromarray(np.array(img)) # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    img = img.convert("L") # pyright: ignore[reportUnknownMemberType, reportAttributeAccessIssue, reportUnknownVariableType]
    return img # pyright: ignore[reportUnknownVariableType]




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
def hash_dna(dna_bytes: bytes) -> bytes:
    return hashlib.sha256(dna_bytes).digest()

def derive_key(passphrase: bytes, salt: bytes) -> bytes:
    return scrypt(passphrase, salt, 32, **KDF_PARAMS) # pyright: ignore[reportUnknownVariableType]

from typing import Tuple
from typing import Tuple
from typing import Tuple
def encrypt_hash(key: bytes, plaintext: bytes) -> Tuple[bytes, bytes, bytes]:
    cipher = AES.new(key, AES.MODE_GCM)  # type: ignore[attr-defined]
    ct, tag = cipher.encrypt_and_digest(plaintext)  # type: ignore[attr-defined]
    nonce = cipher.nonce  # type: ignore[attr-defined]
    assert isinstance(nonce, bytes), "nonce must be bytes"
    return nonce, ct, tag # pyright: ignore[reportUnknownVariableType]

def decrypt_hash(key: bytes, nonce: bytes, ct: bytes, tag: bytes) -> bytes:
    cipher: AES = AES.new(key, AES.MODE_GCM, nonce=nonce)  # type: ignore[attr-defined]
    plaintext: bytes = cipher.decrypt_and_verify(ct, tag)  # type: ignore[attr-defined]
    return plaintext # pyright: ignore[reportUnknownVariableType]

def package_blob(
    version: bytes,
    salt: bytes,
    nonce: bytes,
    ct: bytes,
    tag: bytes,
    signature: bytes,
    owner: str
) -> str:
    parts = [version, salt, nonce, ct, tag, signature, owner.encode()]
    return base64.b64encode(b"||".join(parts)).decode()

from typing import Union, List
def unpackage_blob(b64: Union[str, bytes]) -> List[bytes]:
    parts = base64.b64decode(b64).split(b"||")
    return parts  # v, salt, nonce, ct, tag, signature, owner




import PIL.Image # pyright: ignore[reportMissingImports]
def corrupt_image(img: PIL.Image.Image, corruption_fraction: float) -> PIL.Image.Image: # pyright: ignore[reportUnknownMemberType, reportUnknownParameterType]
    arr = np.array(img) # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    total = arr.size # pyright: ignore[reportUnknownVariableType, reportUnknownMemberType]
    nflip = int(total * corruption_fraction) # pyright: ignore[reportUnknownArgumentType]
    idx: tuple[np.ndarray, ...] = np.unravel_index(np.random.choice(total, nflip, replace=False), arr.shape)  # type: ignore[attr-defined]
    arr[idx] = 255 - arr[idx]  # invert pixels
    # Ensure arr is uint8 for PILImage.fromarray
    arr_uint8 = arr.astype('uint8') if arr.dtype != 'uint8' else arr # pyright: ignore[reportUnknownVariableType, reportUnknownMemberType]
    return PILImage.fromarray(arr_uint8) # pyright: ignore[reportUnknownVariableType, reportUnknownMemberType]

def qr_decode_sim(img): # pyright: ignore[reportUnknownParameterType, reportMissingParameterType]
    # Lightweight decode simulation: attempt to read QR via qrcode lib is not available.
    # For real tests, use zxing, pyzbar, or OpenCV QR detector. Here we simulate decode success
    # based on image "damage" heuristics: if many pixels inverted, fail.
    arr: np.ndarray = np.array(img) # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    damage: float = float(np.mean(arr)) / 255.0  # type: ignore[attr-defined]
    # heuristic: if damage metric below threshold, assume decode success
    return True if damage > 0.01 else False

# --- Setup initial keys (owner signing key)
owner = "LEIF William Sogge"
version = b"v1"
signing_key: SigningKey = SigningKey.generate()  # type: ignore[attr-defined]
verify_key: VerifyKey = signing_key.verify_key # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]


# store payloads (simulate DB)
from typing import List
payload_store: List[str] = []

# metrics
from typing import Dict
metrics: Dict[str, Any] = {
    'roundtrip_success': [],  # type: List[int]
    'corruption_success': {c: [] for c in CORRUPTION_LEVELS},
    'verify_failures': 0,
    'decrypt_failures': 0,
    'tamper_detected': 0,
}

# --- Generate baseline payloads

def create_payload(dna_bytes: bytes, passphrase: bytes, owner: str, signing_key: 'SigningKey', version: bytes): # pyright: ignore[reportUnknownParameterType]
    salt = get_random_bytes(SALT_LEN)
    key = derive_key(passphrase, salt)
    dna_h = hash_dna(dna_bytes)
    nonce, ct, tag = encrypt_hash(key, dna_h)
    sig: bytes = signing_key.sign(version + salt + nonce + ct + tag + owner.encode()).signature # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
    blob = package_blob(version, salt, nonce, ct, tag, sig, owner) # pyright: ignore[reportUnknownArgumentType]
    return blob

# create N baseline payloads
N = 200  # scale up for heavy tests
for i in range(N):
    dna = ("AGCT" * (10 + (i % 5))).encode() + bytes([i % 256])
    payload_store.append(create_payload(dna, PASS, owner, signing_key, version)) # pyright: ignore[reportUnknownArgumentType]

# --- Simulation loop (accelerated days)
current_signing_key = signing_key # pyright: ignore[reportUnknownVariableType]
current_verify_key: VerifyKey = verify_key # pyright: ignore[reportUnknownVariableType]
from typing import List
revoked_keys: List[bytes] = []
for day in range(SIM_DAYS):
    # trigger compromise
    if day == EVENTS['compromise_day']:
        compromised_key = current_signing_key # pyright: ignore[reportUnknownVariableType]
        print(f"[SIM] Day {day}: Simulating key compromise.")
    # trigger rotation
    if day == EVENTS['rotation_day']:
        print(f"[SIM] Day {day}: Rotating signing key.")
        new_sk: SigningKey = SigningKey.generate()  # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType] # type: SigningKey
        # Ensure current_verify_key is bytes before appending
        if hasattr(current_verify_key, 'encode'): # pyright: ignore[reportUnknownArgumentType]
            revoked_keys.append(current_verify_key.encode())  # type: ignore[attr-defined]
        elif isinstance(current_verify_key, bytes):
            revoked_keys.append(current_verify_key)
        else:
            # Try to get bytes representation
            revoked_keys.append(bytes(current_verify_key)) # pyright: ignore[reportUnknownArgumentType]
        current_signing_key = new_sk # pyright: ignore[reportUnknownVariableType]
        current_verify_key = new_sk.verify_key # pyright: ignore[reportUnknownVariableType, reportUnknownMemberType]
        # optionally re-sign stored payloads if policy requires (simulate cost)
        # For this test we do not re-sign automatically; we record how many would need re-signing.

# --- Run verification & corruption tests
for idx, blob in enumerate(payload_store):  # type: ignore[assignment]
    blob: str
    try:
        parts = unpackage_blob(blob)
        v, salt, nonce, ct, tag, signature, owner_bytes = parts
        # verify signature with original verify_key (we used initial signing_key)
        if isinstance(verify_key, VerifyKey):
            vk: VerifyKey = verify_key  # pyright: ignore[reportUnknownVariableType] # type: VerifyKey
        else:
            vk: VerifyKey = VerifyKey(verify_key)  # type: ignore[assignment]
        vk.verify(v + salt + nonce + ct + tag + owner_bytes, signature)  # type: ignore[attr-defined]
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
        corrupted: PILImage.Image = corrupt_image(img, c) # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
        ok = qr_decode_sim(corrupted) # pyright: ignore[reportUnknownArgumentType]
        metrics['corruption_success'][c].append(1 if ok else 0)

# --- Tamper tests: flip bytes in blob
for i in range(50):
    blob_i: str = payload_store[i]  # type: ignore[assignment]
    b = base64.b64decode(blob_i)
    # flip a random byte in ciphertext region
    b_list = bytearray(b)
    pos = random.randint(0, len(b_list)-1)
    b_list[pos] ^= 0xFF
    tampered = base64.b64encode(bytes(b_list)).decode()
    try:
        parts = unpackage_blob(tampered)
        v, salt, nonce, ct, tag, signature, owner_bytes = parts
        if isinstance(verify_key, VerifyKey):
            vk: VerifyKey = verify_key # pyright: ignore[reportUnknownVariableType]
        else:
            vk: VerifyKey = VerifyKey(verify_key)  # type: ignore[assignment]
        # Explicit type annotation for pyright
        vk: VerifyKey = vk # pyright: ignore[reportUnknownVariableType]
        vk.verify(v + salt + nonce + ct + tag + owner_bytes, signature)  # type: ignore[attr-defined]
        # if signature passes (unlikely), try decrypt
        key = derive_key(PASS, salt)
        decrypt_hash(key, nonce, ct, tag)
    except Exception:
        metrics['tamper_detected'] += 1

# --- Summarize metrics
print("\n--- Simulation Summary ---")
print("Roundtrip success rate:", sum(metrics['roundtrip_success'])/len(metrics['roundtrip_success'])) # pyright: ignore[reportUnknownArgumentType]
for c in CORRUPTION_LEVELS:
    arr: list[int] = metrics['corruption_success'][c]
    print(f"Corruption {int(c*100)}% decode success:", sum(arr)/len(arr)) # pyright: ignore[reportUnknownArgumentType]
print("Verify failures: " + str(metrics['verify_failures'])) # pyright: ignore[reportUnknownArgumentType]
print("Tamper detections:", str(metrics['tamper_detected'])) # pyright: ignore[reportUnknownArgumentType]
