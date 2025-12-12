import zlib
import hashlib
import os
import qrcode
import base64
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
from cryptography.hazmat.backends import default_backend
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Simulate larger human genome data (real: load full FASTA/FASTQ/variants VCF)
human_dna_chunk = "ACGT" * 1000000  # ~4MB raw demo — scales to real genomes
human_sequence = Seq(human_dna_chunk)

human_record = SeqRecord(
    human_sequence,
    id="HOMO_SAPIENS_001",
    name="You",
    description="Full secure multi-QR human genome teleport"
)

# Step 1: Prepare raw DNA bytes + compress
raw_dna_bytes = str(human_record.seq).encode('utf-8')
compressed_dna = zlib.compress(raw_dna_bytes, level=9)
print(f"Original: {len(raw_dna_bytes)/1e6:.2f} MB → Compressed: {len(compressed_dna)/1e6:.2f} MB")

# Step 2: AES-256-GCM Encryption (same as before)
password = b"UltraSecretHumanKey2025!"
salt = os.urandom(16)
kdf = PBKDF2HMAC(hashes.SHA256(), 32, salt, 600000, default_backend())
key = kdf.derive(password)
iv = os.urandom(12)

cipher = Cipher(algorithms.AES(key), modes.GCM(iv), backend=default_backend())
encryptor = cipher.encryptor()
ciphertext = encryptor.update(compressed_dna) + encryptor.finalize()
tag = encryptor.tag

encrypted_package = salt + iv + tag + ciphertext

# Step 3: Base64 encode the full encrypted blob
full_b64 = base64.b64encode(encrypted_package).decode('utf-8')

# Step 4: MULTI-QR SPLITTING
CHUNK_SIZE = 1800  # Safe bytes per QR (leaves room for header ~50 chars + overhead)
chunks = [full_b64[i:i + CHUNK_SIZE] for i in range(0, len(full_b64), CHUNK_SIZE)]
total = len(chunks)
print(f"Splitting into {total} secure QRs for teleport...")

qr_images = []
for idx, chunk in enumerate(chunks):
    part_payload = f"HUMAN_TELEPORT_V5_MULTI:{idx+1}/{total}:{chunk}"
    
    qr = qrcode.QRCode(
        error_correction=qrcode.constants.ERROR_CORRECT_H,  # Quantum noise resistant
        box_size=10,
        border=4,
    )
    qr.add_data(part_payload)
    qr.make(fit=True)
    img = qr.make_image(fill_color="black", back_color="white")
    qr_images.append(img)
    # In real script: img.save(f"human_teleport_part_{idx+1:04d}_of_{total}.png")
    print(f"Generated QR {idx+1}/{total}")

print("\n🔒 Multi-QR teleport pack complete! Scan all in order → reassemble → decrypt with password.")

# Bonus: Reassembly + Decryption (receiver side)
def reassemble_and_decrypt(scanned_payloads, password):
    # scanned_payloads = list of strings from scanning each QR
    full_b64 = ""
    for payload in sorted(scanned_payloads):  # Sort by part num if needed
        if not payload.startswith("HUMAN_TELEPORT_V5_MULTI:"):
            raise ValueError("Invalid QR")
        parts = payload.split(":", 2)
        chunk = parts[2]
        full_b64 += chunk
    
    data = base64.b64decode(full_b64)
    salt, iv, tag, ciphertext = data[:16], data[16:28], data[28:44], data[44:]
    
    kdf = PBKDF2HMAC(hashes.SHA256(), 32, salt, 600000, default_backend())
    key = kdf.derive(password)
    
    cipher = Cipher(algorithms.AES(key), modes.GCM(iv, tag), backend=default_backend())
    decryptor = cipher.decryptor()
    compressed = decryptor.update(ciphertext) + decryptor.finalize()
    
    original_dna_bytes = zlib.decompress(compressed)
    return original_dna_bytes.decode('utf-8')

# Test with fake scanned (in real: feed from scanner lib)
fake_scanned = [f"HUMAN_TELEPORT_V5_MULTI:{i+1}/{total}:{chunks[i]}" for i in range(total)]
recovered = reassemble_and_decrypt(fake_scanned, password)
print("Reconstruction SUCCESS!", recovered == str(human_record.seq))
