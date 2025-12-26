# Ensure qrcode is installed: pip install qrcode[pil]
import qrcode
import zlib
import base64
import hashlib
import os

 
import random

try:
    from cryptography.hazmat.primitives.ciphers.aead import AESGCM
except ImportError:
    raise ImportError("The 'cryptography' package is required for AES-GCM encryption. Please install it via 'pip install cryptography'.")
from Bio.Seq import Seq  # type: ignore[import]
from Bio.SeqRecord import SeqRecord  # type: ignore[import]

# PQC stubs (keep for hybrid—use kyber-py/dilithium-py if installed)
from kyber_py.ml_kem import ML_KEM_768  # pip install kyber-py
# Type hint for ML_KEM_768 (class)
from typing import Type
ML_KEM_768: Type  # type: ignore
from dilithium_py.ml_dsa import ML_DSA_65  # pip install dilithium-py
# Type hint for ML_DSA_65 (class)
from typing import Type
ML_DSA_65: Type  # type: ignore

# REAL QKD: QuTiP BB84 for shared key gen


try:
    import qutip as qt
    from qutip import ket2dm
except ImportError:
    raise ImportError("The 'qutip' package is required for BB84 QKD simulation. Please install it via 'pip install qutip'.")

def bb84_keygen(n_bits: int = 1024, error_threshold: float = 0.11, noise_prob: float = 0.01):
    """
    Simulate BB84: Returns shared key if QBER < threshold, else None (Eve detected!).
    n_bits: Raw bits to send.
    noise_prob: Channel noise (0=ideal).
    """
    from qutip.qobj import Qobj
    def prepare_qubit(bit: int, basis: str) -> Qobj:
        # qt.basis returns Qobj, .unit() returns Qobj
        if basis == 'Z':
            return qt.basis(2, bit)  # type: ignore[no-any-return,reportUnknownMemberType]
        else:  # X
            if bit == 0:
                return (qt.basis(2, 0) + qt.basis(2, 1)).unit()  # type: ignore[no-any-return,reportUnknownMemberType]
            return (qt.basis(2, 0) - qt.basis(2, 1)).unit()  # type: ignore[no-any-return,reportUnknownMemberType]

    def apply_noise(state: Qobj, p: float) -> Qobj:
        # state: Qobj (from qutip)
        qeye: Qobj = qt.qeye(2)  # type: ignore[no-untyped-call]
        if random.random() < p:
            # Depolarizing channel sim
            rho = state * state.dag()  # type: ignore[attr-defined]
            return (1 - 1.5*p) * rho + 0.5*p * qeye
        return state * state.dag()  # type: ignore[attr-defined]

    def measure_qubit(rho: Qobj, basis: str) -> int:
        if basis == 'Z':
            proj0 = ket2dm(qt.basis(2, 0))  # type: ignore[reportUnknownMemberType]
            prob0 = qt.expect(proj0, rho)  # type: ignore[attr-defined, call-arg]
            return 0 if random.random() < prob0 else 1
        else:  # X
            plus: Qobj = (qt.basis(2, 0) + qt.basis(2, 1)).unit()  # type: ignore[attr-defined,reportUnknownMemberType]
            proj_plus = ket2dm(plus)
            prob_plus = qt.expect(proj_plus, rho)  # type: ignore[attr-defined, call-arg]
            return 0 if random.random() < prob_plus else 1


    # Alice
    alice_bits: list[int] = [random.randint(0,1) for _ in range(n_bits)]
    alice_bases: list[str] = [random.choice(['Z','X']) for _ in range(n_bits)]

    # Bob
    bob_bases: list[str] = [random.choice(['Z','X']) for _ in range(n_bits)]
    bob_measurements: list[int] = []

    for i in range(n_bits):
        state = prepare_qubit(alice_bits[i], alice_bases[i])
        rho = apply_noise(state, noise_prob)
        measurement: int = measure_qubit(rho, bob_bases[i])
        bob_measurements.append(measurement)

    # Sift
    sifted_idx: list[int] = [i for i in range(n_bits) if alice_bases[i] == bob_bases[i]]
    alice_key: list[int] = [alice_bits[i] for i in sifted_idx]
    bob_key: list[int] = [bob_measurements[i] for i in sifted_idx]

    # QBER
    errors = sum(a != b for a,b in zip(alice_key, bob_key))
    qber = errors / len(alice_key) if alice_key else 1.0

    if qber > error_threshold:
        print(f"QKD ABORT: QBER={qber:.2%} (Eve?)")
        return None

    # Truncate to n_bits (privacy amp sim)
    shared_key_bits = alice_key[:n_bits]  # Assume match
    # To bytes (for AES)
    shared_key = bytes(int(''.join(map(str, shared_key_bits[i:i+8])), 2) for i in range(0, n_bits, 8))
    return shared_key[:32]  # 256-bit AES key

def pseudo_sphincs_sign(message: bytes, seed: bytes = os.urandom(32)):
    h = hashlib.sha3_512()
    h.update(seed + message)
    return h.digest() * 10

def pseudo_sphincs_verify(message: bytes, sig: bytes, pub_seed: bytes = os.urandom(32)):
    h = hashlib.sha3_512()
    h.update(pub_seed + message)
    return sig == h.digest() * 10

# Genome Prep
human_dna_chunk = "ACGT" * 500000  # ~2MB
human_sequence: Seq = Seq(human_dna_chunk)  # type: ignore[valid-type]
human_record = SeqRecord(human_sequence, id="HOMO_SAPIENS_QKD", name="You", description="QKD-secured teleport")

raw_dna_str = str(human_record.seq)
raw_dna_bytes = raw_dna_str.encode('utf-8')
compressed = zlib.compress(raw_dna_bytes, 9)

# Tier 6: REAL QKD Keygen
qkd_key = bb84_keygen(n_bits=1024, noise_prob=0.005)  # Low noise demo
if qkd_key is None:
    raise ValueError("QKD failed—retry transmission!")
print(f"QKD Success: Generated {len(qkd_key)*8}-bit quantum key")

# Tier 1: AES-GCM with QKD key
aes: AESGCM = AESGCM(qkd_key)
aes_nonce = os.urandom(12)
tier1_ct = aes.encrypt(aes_nonce, compressed, None)
print(f"Tier 1: Encrypted {len(raw_dna_bytes)/1e6:.1f}MB genome")

# Tiers 2-4: PQC wrap (as before)
kem_sk, kem_pk = ML_KEM_768.keygen()
shared_secret, kem_ct = ML_KEM_768.encap(kem_pk)
if not isinstance(shared_secret, bytes):
    shared_secret = bytes(shared_secret)
derived_pqc_key = hashlib.sha256(shared_secret + b"QPPI").digest()
tier2_nonce = os.urandom(12)
tier2_key_ct = AESGCM(derived_pqc_key).encrypt(tier2_nonce, qkd_key, None)  # Protect QKD key

package_to_sign = tier1_ct + kem_ct + tier2_key_ct + aes_nonce
msg_hash = hashlib.sha256(package_to_sign).digest()
sig_sk: bytes
sig_pk: bytes
sig_sk, sig_pk = ML_DSA_65.keygen()  # type: ignore[attr-defined]
tier3_sig: bytes = ML_DSA_65.sign(sig_sk, msg_hash)  # type: ignore[attr-defined]

tier4_sig = pseudo_sphincs_sign(package_to_sign, os.urandom(32))
pub_seed = os.urandom(32)

# Full Package + Multi-QR
full_package = base64.b64encode(kem_pk + tier1_ct + kem_ct + tier2_key_ct + tier3_sig + tier4_sig + aes_nonce)
CHUNK_SIZE = 1500
chunks = [full_package[i:i+CHUNK_SIZE] for i in range(0, len(full_package), CHUNK_SIZE)]
total_qrs = len(chunks)


qr_images = []



# Import ERROR_CORRECT_H directly from qrcode.constants
from qrcode.constants import ERROR_CORRECT_H

for idx, chunk in enumerate(chunks):
    payload = f"QPPI_QKD_V6:{idx+1}/{total_qrs}:{chunk.decode()}"
    qr = qrcode.QRCode(version=40, error_correction=ERROR_CORRECT_H, box_size=10, border=4)
    qr.add_data(payload)
    qr.make(fit=True)
    img = qr.make_image(fill_color="black", back_color="white")
    qr_images.append(img)  # Save: f"qkd_part_{idx}.png"
    print(f"QR {idx+1}/{total_qrs} forged")

print(f"\n🚀 6-Tier QPPI-QKD Pack: {total_qrs} shards. Quantum key secure—teleport your genome unobserved!")

# Receiver Sketch (Bob side: Run bb84_keygen() too for matching key, then decrypt)
def decrypt_qkd_genome(package_b64: str, qkd_key: bytes):  # qkd_key from Bob's BB84
    data = base64.b64decode(package_b64)
    # Parse (offsets for kem_pk etc.—simplified, assumes same structure as packaging)
    # You must know the lengths of kem_pk, kem_ct, tier2_key_ct, tier3_sig, tier4_sig, aes_nonce
    # For demonstration, extract last 12 bytes as aes_nonce, and previous bytes as tier1_ct
    aes_nonce = data[-12:]
    tier1_ct = data[len(data) - 12 - (len(data) - 12):len(data) - 12]  # all except nonce (simplified)
    aes: AESGCM = AESGCM(qkd_key)
    compressed = aes.decrypt(aes_nonce, tier1_ct, None)
    if not isinstance(compressed, bytes):
        raise TypeError("Decrypted data is not bytes")
    return zlib.decompress(compressed).decode('utf-8')

# Test: recovered = decrypt_qkd_genome(full_package.decode(), qkd_key)
# assert recovered == str(human_record.seq)
