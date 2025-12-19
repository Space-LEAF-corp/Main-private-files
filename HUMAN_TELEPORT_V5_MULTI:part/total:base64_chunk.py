import zlib
import base64
try:
    from cryptography.hazmat.primitives import hashes  # type: ignore
    from cryptography.hazmat.primitives.ciphers import Cipher  # type: ignore
    from cryptography.hazmat.primitives.ciphers import modes  # type: ignore
    from cryptography.hazmat.primitives.ciphers import algorithms  # type: ignore
    from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC  # type: ignore
except ImportError as e:
    raise ImportError("cryptography module is required. Install with 'pip install cryptography'") from e
from typing import cast

from typing import Any



# Import Seq and SeqRecord from Bio.Seq with ImportError handling
try:
    try:

        # Import Seq from Bio.Seq with ImportError handling
        try:
            from Bio.Seq import Seq  # type: ignore[import]
        except ImportError as e:
            raise ImportError("Biopython is required. Install with 'pip install biopython'") from e
    except ImportError as e:
        raise ImportError("Biopython is required. Install with 'pip install biopython'") from e
    try:
        from Bio.SeqRecord import SeqRecord  # type: ignore
    except ImportError as e:
        raise ImportError("Biopython is required. Install with 'pip install biopython'") from e
    from Bio.Seq import Seq as _SeqType  # type: ignore
except ImportError as e:
    raise ImportError("Biopython is required. Install with 'pip install biopython'") from e

# Define password and salt before use
import os
password = b"UltraSecretHumanKey2025!"
salt = os.urandom(16)
from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC  # type: ignore
kdf: PBKDF2HMAC = PBKDF2HMAC(hashes.SHA256(), 32, salt, 600000)  # type: ignore
key = cast(bytes, kdf.derive(password))  # type: ignore[attr-defined]
iv = os.urandom(12)



from typing import cast



from cryptography.hazmat.primitives.ciphers import Cipher  # type: ignore
cipher = Cipher(algorithms.AES(key), modes.GCM(iv))  # type: ignore[attr-defined]
human_sequence: Seq = Seq("ACGT" * 1000)  # type: ignore[valid-type]
from Bio.SeqRecord import SeqRecord  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore
human_record: 'SeqRecord' = SeqRecord( # type: ignore
    human_sequence,
    id="HOMO_SAPIENS_001",
    name="You",
    description="Full secure multi-QR human genome teleport"
)

# Step 1: Prepare raw DNA bytes + compress
raw_dna_bytes = str(human_record.seq).encode('utf-8')  # type: ignore[arg-type]
compressed_dna = zlib.compress(raw_dna_bytes, level=9)
print(f"Original: {len(raw_dna_bytes)/1e6:.2f} MB → Compressed: {len(compressed_dna)/1e6:.2f} MB")



# Step 2: AES-256-GCM Encryption (same as before)
password = b"UltraSecretHumanKey2025!"
salt = os.urandom(16)
from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC  # type: ignore
from cryptography.hazmat.primitives.kdf import KeyDerivationFunction  # type: ignore
kdf: PBKDF2HMAC = PBKDF2HMAC(algorithm=hashes.SHA256(), length=32, salt=salt, iterations=600000)  # type: ignore[attr-defined]
key: bytes = kdf.derive(password)  # type: ignore
iv = os.urandom(12)



from cryptography.hazmat.primitives.ciphers.base import CipherContext  # type: ignore
encryptor: CipherContext = cipher.encryptor()  # type: ignore
ciphertext: bytes = encryptor.update(compressed_dna) + encryptor.finalize()  # type: ignore[attr-defined]
tag: bytes = encryptor.tag  # type: ignore[attr-defined]

encrypted_package: bytes = salt + iv + tag + ciphertext  # type: ignore[assignment]

# Step 3: Base64 encode the full encrypted blob
if not isinstance(encrypted_package, bytes):
    raise TypeError("encrypted_package must be bytes before base64 encoding.")
full_b64 = base64.b64encode(encrypted_package).decode('utf-8')

# Step 4: MULTI-QR SPLITTING


from typing import List, Any
CHUNK_SIZE = 1800  # Safe bytes per QR (leaves room for header ~50 chars + overhead)
assert isinstance(full_b64, str), "full_b64 must be a string before splitting into chunks"
chunks: List[str] = [full_b64[i:i + CHUNK_SIZE] for i in range(0, len(full_b64), CHUNK_SIZE)]
total = len(chunks)
print(f"Splitting into {total} secure QRs for teleport...")

qr_images: List[Any] = []










# Import qrcode and ERROR_CORRECT_H with proper error handling
try:
    import qrcode  # pyright: ignore[reportMissingModuleSource]
    try:
        ERROR_CORRECT_H = getattr(__import__("qrcode.constants", fromlist=["ERROR_CORRECT_H"]), "ERROR_CORRECT_H")
        _error_correct_h = ERROR_CORRECT_H
    except (ImportError, AttributeError, ModuleNotFoundError):
        _error_correct_h = None
        print("Warning: qrcode.constants could not be imported. QR generation will be skipped.")
except ImportError:
    qrcode = None
    _error_correct_h = None
    print("Warning: qrcode module not installed. Install with 'pip install qrcode[pil]'. QR generation will be skipped.")




for idx, chunk in enumerate(chunks):
    part_payload = f"HUMAN_TELEPORT_V5_MULTI:{idx+1}/{total}:{chunk}"
    if qrcode is not None and _error_correct_h is not None:
        qr = qrcode.QRCode(
            error_correction=_error_correct_h,  # Quantum noise resistant
            box_size=10,
            border=4,
        )
        qr.add_data(part_payload)
        qr.make(fit=True)
        img = qr.make_image(fill_color="black", back_color="white")
        qr_images.append(img)
        # In real script: img.save(f"human_teleport_part_{idx+1:04d}_of_{total}.png")
        print(f"Generated QR {idx+1}/{total}")
    else:
        print("qrcode module not installed, skipping QR generation for part", idx+1)

print("\n🔒 Multi-QR teleport pack complete! Scan all in order → reassemble → decrypt with password.")

# Bonus: Reassembly + Decryption (receiver side)
from typing import List

def reassemble_and_decrypt(scanned_payloads: List[str], password: bytes) -> str:
    # scanned_payloads = list of strings from scanning each QR
    # cryptography imports already handled at top
    full_b64 = ""
    for payload in sorted(scanned_payloads):  # Sort by part num if needed
        if not payload.startswith("HUMAN_TELEPORT_V5_MULTI:"):
            raise ValueError("Invalid QR")
        parts = payload.split(":", 2)
        chunk = parts[2]
        full_b64 += chunk
    # full_b64 is always a string here; no need for isinstance check
    data = base64.b64decode(full_b64.encode('utf-8'))
    salt, iv, tag, ciphertext = data[:16], data[16:28], data[28:44], data[44:]
    
    sha256_algorithm: hashes.HashAlgorithm = hashes.SHA256()  # type: ignore[valid-type]
    kdf: PBKDF2HMAC = PBKDF2HMAC( # pyright: ignore[reportUnknownVariableType]
        algorithm=sha256_algorithm,
        length=32,
        salt=salt,
        iterations=600000
    )
    key: bytes = kdf.derive(password)  # type: ignore
    
    cipher = Cipher(algorithms.AES(key), modes.GCM(iv, tag))  # type: ignore[attr-defined]
    from cryptography.hazmat.primitives.ciphers.base import CipherContext  # type: ignore
    decryptor: CipherContext = cipher.decryptor()  # type: ignore
    compressed: bytes = decryptor.update(ciphertext) + decryptor.finalize()  # type: ignore
    
    try:
        # Ensure compressed is bytes for type safety
        if not isinstance(compressed, bytes):
            raise TypeError("Decryption did not return bytes as expected.")
        original_dna_bytes = zlib.decompress(compressed)
    except zlib.error as e:
        raise RuntimeError(f"Decompression failed: {e}")
    return original_dna_bytes.decode('utf-8')

# Test with fake scanned (in real: feed from scanner lib)
fake_scanned = [f"HUMAN_TELEPORT_V5_MULTI:{i+1}/{total}:{chunks[i]}" for i in range(total)]
recovered = reassemble_and_decrypt(fake_scanned, password)
# Explicitly cast human_record.seq to str for type safety
try:
    from Bio.Seq import Seq  # type: ignore[import]
except ImportError as e:
    raise ImportError("Biopython is required. Install with 'pip install biopython'") from e
expected_seq: str = str(human_record.seq)  # type: ignore[arg-type]
print("Reconstruction SUCCESS!", recovered == expected_seq)
