import hashlib
import os
from itertools import zip_longest

def rotate_left(data: bytes, n: int) -> bytes:
    n = n % (len(data) * 8)
    bitstring = int.from_bytes(data, "big")
    width = len(data) * 8
    bitstring = ((bitstring << n) & ((1 << width) - 1)) | (bitstring >> (width - n))
    return bitstring.to_bytes(len(data), "big")

def rotate_right(data: bytes, n: int) -> bytes:
    n = n % (len(data) * 8)
    bitstring = int.from_bytes(data, "big")
    width = len(data) * 8
    bitstring = (bitstring >> n) | ((bitstring << (width - n)) & ((1 << width) - 1))
    return bitstring.to_bytes(len(data), "big")

def interleave_bytes(a: bytes, b: bytes) -> bytes:
    out = bytearray()
    for x, y in zip_longest(a, b, fillvalue=0):
        out.append(x)
        out.append(y)
    return bytes(out)

def flip_bytes(data: bytes) -> bytes:
    # Bitwise NOT within byte range
    return bytes((~b) & 0xFF for b in data)

def anti_apple_dna_lock(user_id: str, device_secret: bytes = None) -> bytes:
    """
    Anti-Apple DNA Lock (toy demo):
    - user_id: stable identifier (string)
    - device_secret: per-device secret (bytes); if None, random for demo
    Returns: 32-byte key (256-bit)
    """
    if device_secret is None:
        device_secret = os.urandom(32)  # demo only

    # Step 1: derive two seeds (A, B)
    # A = stable seed (user + device)
    seed_a = hashlib.sha256(user_id.encode("utf-8") + device_secret).digest()
    # B = dynamic seed (nonce + user)
    nonce = os.urandom(16)
    seed_b = hashlib.sha256(nonce + user_id.encode("utf-8")).digest()

    # Step 2: split rotation
    r1 = seed_a[0] % 64  # small rotation factor from data
    r2 = seed_b[0] % 64
    a_rot = rotate_left(seed_a, r1)
    b_rot = rotate_right(seed_b, r2)

    # Step 3: DNA-style interleave
    helix = interleave_bytes(a_rot, b_rot)

    # Step 4: mix (simple hash-based "rounds")
    mixed = hashlib.sha256(helix).digest()
    mixed = hashlib.sha256(mixed + helix).digest()

    # Step 5: split, flip, fuse
    mid = len(mixed) // 2
    left = mixed[:mid]
    right = mixed[mid:]
    right_flipped = flip_bytes(right)
    fused = bytes(l ^ r for l, r in zip(left, right_flipped))

    # Step 6: final hardened key
    final_key = hashlib.sha256(fused + seed_a + seed_b).digest()
    return final_key

# Example usage
if __name__ == "__main__":
    user_id = "captain.leif"
    device_secret = os.urandom(32)  # in real use, this would be fixed per device
    key = anti_apple_dna_lock(user_id, device_secret)
    print("Anti-Apple DNA Lock key (hex):", key.hex())
