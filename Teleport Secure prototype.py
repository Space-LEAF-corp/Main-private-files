# Requires: pip install pycryptodome pynacl qrcode pillow
import base64, hashlib
from Crypto.Cipher import AES
from Crypto.Protocol.KDF import scrypt
from Crypto.Random import get_random_bytes
from nacl.signing import SigningKey, VerifyKey
from io import BytesIO
import qrcode

# --- Inputs (example)
fake_dna = ("AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT" * 10).encode()  # synthetic example only
passphrase = b"SuperSecretPassphrase2025"  # protect this; use secure vault in prod
owner = "LEIF William Sogge"
version = b"v1"

# --- 1) Hash DNA (non-reversible fingerprint)
dna_hash = hashlib.sha256(fake_dna).digest()

# --- 2) Derive symmetric key with scrypt (salted)
salt = get_random_bytes(16)
key = scrypt(passphrase, salt, 32, N=2**14, r=8, p=1)

# --- 3) AES-GCM encrypt the hash
cipher = AES.new(key, AES.MODE_GCM)
ciphertext, tag = cipher.encrypt_and_digest(dna_hash)
nonce = cipher.nonce

# --- 4) Sign the package with Ed25519 (author key)
# Generate or load a persistent signing key; here we generate for demo.
signing_key = SigningKey.generate()
verify_key = signing_key.verify_key
signature = signing_key.sign(version + salt + nonce + ciphertext + tag + owner.encode()).signature

# --- 5) Package: version||salt||nonce||ciphertext||tag||signature||owner  (base64)
blob = b"||".join([version, salt, nonce, ciphertext, tag, signature, owner.encode()])
qr_payload = base64.b64encode(blob).decode()

# --- 6) Make QR
qr = qrcode.QRCode(error_correction=qrcode.constants.ERROR_CORRECT_H, box_size=6, border=2)
qr.add_data(qr_payload)
qr.make(fit=True)
img = qr.make_image(fill_color="black", back_color="white")

# Save or display as needed (demo: base64 PNG)
buf = BytesIO()
img.save(buf, format="PNG")
png_b64 = base64.b64encode(buf.getvalue()).decode()
print("QR PNG base64 (short):", png_b64[:120], "...")

# --- 7) Verification snippet (any verifier)
def verify_payload(b64_payload, verify_key_bytes, passphrase):
    blob = base64.b64decode(b64_payload)
    parts = blob.split(b"||")
    v, salt, nonce, ciphertext, tag, signature, owner_bytes = parts
    # verify signature
    vk = VerifyKey(verify_key_bytes)
    vk.verify(v + salt + nonce + ciphertext + tag + owner_bytes, signature)
    # derive key and decrypt
    key = scrypt(passphrase, salt, 32, N=2**14, r=8, p=1)
    cipher = AES.new(key, AES.MODE_GCM, nonce=nonce)
    plaintext = cipher.decrypt_and_verify(ciphertext, tag)
    return plaintext, owner_bytes.decode()
