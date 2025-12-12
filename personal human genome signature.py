import zlib
import hashlib
import qrcode
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import base64
from io import BytesIO

# Simulate a chunk of human genome (real: use SeqIO.parse on a FASTA/FASTQ file)
human_dna_chunk = "ACGT" * 10000  # Scale up—real human chr1 alone is 248 million bp
human_sequence = Seq(human_dna_chunk)

human_record = SeqRecord(
    human_sequence,
    id="HOMO_SAPIENS_001",
    name="You",
    description="Human pre-teleport genome signature"
)

# Full sequence as bytes for compression
raw_dna_bytes = str(human_record.seq).encode('utf-8')

# Compact encoding: Compress heavily
compressed_dna = zlib.compress(raw_dna_bytes, level=9)  # gzip-level shrink
print(f"Original size: {len(raw_dna_bytes)} bytes")
print(f"Compressed size: {len(compressed_dna)} bytes")  # Way smaller!

# Secure hash of original
dna_hash = hashlib.sha256(raw_dna_bytes).hexdigest()

# Encrypt compressed data (demo XOR—use cryptography.fernet for real AES)
secret_key = "UltraSecretHumanKey2025"

def xor_encrypt(data, key):
    key_bytes = key.encode()
    return bytes(b ^ key_bytes[i % len(key_bytes)] for i, b in enumerate(data))

encrypted_compressed = xor_encrypt(compressed_dna, secret_key)

# QR payload
qr_data = f"HUMAN_TELEPORT_V3_BIO:{base64.b64encode(encrypted_compressed).decode()}"

# High-capacity QR
qr = qrcode.QRCode(
    version=40,  # Max size for big data
    error_correction=qrcode.constants.ERROR_CORRECT_H,  # Mutation-resistant
    box_size=10,
    border=4,
)
qr.add_data(qr_data)
qr.make(fit=True)
img = qr.make_image(fill_color="black", back_color="white")

# In real script: img.save("secure_human_teleport_qr.png")
print("Your Secure Human Genome QR is ready—scan with key to decrypt, decompress, reconstruct!")
