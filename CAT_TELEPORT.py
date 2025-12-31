try:
    import qrcode # pyright: ignore[reportMissingModuleSource]
    from qrcode.constants import ERROR_CORRECT_H # pyright: ignore[reportMissingModuleSource]
except ImportError:
    print("\n[ERROR] The 'qrcode' module is not installed or could not be resolved.\nPlease install it by running:\n    pip install qrcode[pil]\n")
    import sys
    sys.exit(1)
# NOTE: If you see an import error in your IDE (e.g., 'Import \"qrcode\" could not be resolved from source'), but the script runs fine, it means the package is installed for Python but not indexed by your IDE. As long as the script runs, you can ignore the warning.
import hashlib
from io import BytesIO
import base64

# Step 1: Fake cat DNA (in real life: sequence from Biopython or full genome)
fake_dna = "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT" * 100  # Unique to this cat

# Create secure DNA signature (SHA-256 hash—impossible to reverse)
dna_hash = hashlib.sha256(fake_dna.encode()).hexdigest()
print("Secure DNA Signature (Hash):", dna_hash)

# Step 2: ADD SECURITY - Encrypt the hash with a secret key
secret_key = "SuperSecretCatKey2025"  # Change this to something strong!

def simple_encrypt(data: str, key: str) -> str:
    # Simple XOR encryption (symmetric, for demo—real: use cryptography library + AES)
    return ''.join(chr(ord(d) ^ ord(k)) for d,k in zip(data, key * (len(data)//len(key) + 1)))

encrypted_sig = simple_encrypt(dna_hash, secret_key)

# Step 3: Pack into QR data (with identifier for teleport protocol)
qr_data = f"CAT_TELEPORT:{encrypted_sig}"

# Step 4: Generate secure QR code (high error correction = mutation-resistant!)
qr = qrcode.QRCode(
    error_correction=ERROR_CORRECT_H,  # Up to 30% damage fixable
    box_size=10,
    border=4,
)
qr.add_data(qr_data)
qr.make(fit=True)

img = qr.make_image(fill_color="black", back_color="white")

# Output as base64 so you can display/save it anywhere
buf = BytesIO()
img.save(buf)
buf.seek(0)
qr_base64 = base64.b64encode(buf.read()).decode('utf-8')

print("\nYour Secure Cat Teleport QR Ready! (Scan with key to reconstruct)")
# In a real app, display: <img src="data:image/png;base64,{qr_base64}" />

# Bonus: Decryption to verify (only works with correct key)
def simple_decrypt(encrypted: str, key: str) -> str:
    return simple_encrypt(encrypted, key)  # XOR is its own inverse

decrypted = simple_decrypt(encrypted_sig, secret_key)
print("\nDecrypted Hash Matches?", decrypted == dna_hash)  # Should be True
