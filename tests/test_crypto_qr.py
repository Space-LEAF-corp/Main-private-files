# tests/test_crypto_qr.py
import base64, hashlib, os
from io import BytesIO
from Crypto.Cipher import AES
from Crypto.Protocol.KDF import scrypt
from Crypto.Random import get_random_bytes
from nacl.signing import SigningKey, VerifyKey
import qrcode, numpy as np, PIL.Image
import pytest
from pyzbar.pyzbar import decode as zbar_decode

# --- Config (CI-friendly defaults; increase KDF cost in prod)
PASS = b"CI_Test_Passphrase"
KDF_PARAMS = dict(N=2**12, r=8, p=1)   # reduce for CI speed; raise for prod
SALT_LEN = 16
OWNER = "LEIF William Sogge"
VERSION = b"v1"

# --- Helpers (same logic as production code)
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

def qr_decode_image(img):
    buf = BytesIO()
    img.save(buf, format="PNG")
    buf.seek(0)
    decoded = zbar_decode(PIL.Image.open(buf))
    if not decoded:
        return None
    return decoded[0].data.decode()

# --- Fixtures
@pytest.fixture(scope="module")
def signing_keypair():
    # In CI, generate ephemeral keys; in production, load persistent signing key from secure vault
    sk = SigningKey.generate()
    vk = sk.verify_key
    return sk, vk

# --- Tests
def test_roundtrip_encrypt_sign_verify(signing_keypair):
    sk, vk = signing_keypair
    dna = b"AGCT" * 16
    salt = get_random_bytes(SALT_LEN)
    key = derive_key(PASS, salt)
    dna_h = hash_dna(dna)
    nonce, ct, tag = encrypt_hash(key, dna_h)
    sig = sk.sign(VERSION + salt + nonce + ct + tag + OWNER.encode()).signature
    blob = package_blob(VERSION, salt, nonce, ct, tag, sig, OWNER)
    # QR encode/decode
    img = make_qr(blob)
    decoded = qr_decode_image(img)
    assert decoded is not None, "QR decode failed in CI environment"
    v, salt2, nonce2, ct2, tag2, sig2, owner2 = unpackage_blob(decoded)
    # verify signature
    vk.verify(v + salt2 + nonce2 + ct2 + tag2 + owner2, sig2)
    # decrypt
    key2 = derive_key(PASS, salt2)
    plaintext = decrypt_hash(key2, nonce2, ct2, tag2)
    assert plaintext == dna_h

def test_tamper_detection(signing_keypair):
    sk, vk = signing_keypair
    dna = b"AGCT" * 16
    salt = get_random_bytes(SALT_LEN)
    key = derive_key(PASS, salt)
    dna_h = hash_dna(dna)
    nonce, ct, tag = encrypt_hash(key, dna_h)
    sig = sk.sign(VERSION + salt + nonce + ct + tag + OWNER.encode()).signature
    blob = package_blob(VERSION, salt, nonce, ct, tag, sig, OWNER)
    raw = bytearray(base64.b64decode(blob))
    # flip a byte in ciphertext region
    pos = len(raw)//2
    raw[pos] ^= 0xFF
    tampered = base64.b64encode(bytes(raw)).decode()
    with pytest.raises(Exception):
        parts = unpackage_blob(tampered)
        v, salt2, nonce2, ct2, tag2, sig2, owner2 = parts
        vk.verify(v + salt2 + nonce2 + ct2 + tag2 + owner2, sig2)
        key2 = derive_key(PASS, salt2)
        decrypt_hash(key2, nonce2, ct2, tag2)

def test_wrong_passphrase_fails(signing_keypair):
    sk, vk = signing_keypair
    dna = b"AGCT" * 16
    salt = get_random_bytes(SALT_LEN)
    key = derive_key(PASS, salt)
    dna_h = hash_dna(dna)
    nonce, ct, tag = encrypt_hash(key, dna_h)
    sig = sk.sign(VERSION + salt + nonce + ct + tag + OWNER.encode()).signature
    blob = package_blob(VERSION, salt, nonce, ct, tag, sig, OWNER)
    # attempt decrypt with wrong passphrase
    parts = unpackage_blob(blob)
    v, salt2, nonce2, ct2, tag2, sig2, owner2 = parts
    wrong_key = derive_key(b"wrongpass", salt2)
    with pytest.raises(Exception):
        decrypt_hash(wrong_key, nonce2, ct2, tag2)