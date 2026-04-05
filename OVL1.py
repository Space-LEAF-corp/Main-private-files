import numpy as np
import hashlib
import base64
import os
from ecdsa import SigningKey, VerifyingKey, NIST256p
import ast


class SecureSatCom:
    """
    Demo-grade secure comms wrapper.
    - Shared secret for symmetric key derivation
    - ECDSA for authenticity
    - XOR + nonce as placeholder cipher (replace with AES in prod)
    """

    def __init__(self, shared_secret="leaf_starlink_secure_key_2026", sk: SigningKey | None = None):
        self.shared_secret = shared_secret.encode()

        # If a signing key is provided, use it; otherwise generate one (Starlink side)
        if sk is None:
            self.sk = SigningKey.generate(curve=NIST256p)
        else:
            self.sk = sk

        self.vk: VerifyingKey = self.sk.verifying_key

    def derive_key(self, nonce: bytes) -> bytes:
        """
        Derive a session key from shared secret + nonce (HKDF-like).
        """
        hash_obj = hashlib.sha256(self.shared_secret + nonce)
        return hash_obj.digest()[:16]  # 128-bit key

    def encrypt(self, message: str) -> str:
        """
        Simple XOR cipher with derived key + nonce.
        Returns base64-encoded payload: nonce || ciphertext
        """
        nonce = os.urandom(16)
        key = self.derive_key(nonce)

        msg_bytes = message.encode()
        key_stream = (key * (len(msg_bytes) // len(key) + 1))[:len(msg_bytes)]
        encrypted = bytes(a ^ b for a, b in zip(msg_bytes, key_stream))

        payload = nonce + encrypted
        return base64.b64encode(payload).decode()

    def decrypt(self, encrypted_message: str) -> str:
        """
        Decrypt counterpart: extract nonce, derive key, XOR.
        """
        payload = base64.b64decode(encrypted_message.encode())
        nonce, encrypted = payload[:16], payload[16:]
        key = self.derive_key(nonce)

        key_stream = (key * (len(encrypted) // len(key) + 1))[:len(encrypted)]
        decrypted_bytes = bytes(a ^ b for a, b in zip(encrypted, key_stream))
        return decrypted_bytes.decode()

    def sign_message(self, message: str) -> str:
        """Sign with ECDSA for authenticity."""
        return self.sk.sign(message.encode()).hex()

    @staticmethod
    def verify_signature(message: str, signature_hex: str, vk: VerifyingKey) -> bool:
        """Verify on receiver side."""
        sig = bytes.fromhex(signature_hex)
        return vk.verify(sig, message.encode())


class StarlinkNode:
    def __init__(self, position, velocity, com: SecureSatCom):
        self.position = np.array(position)  # [x, y, z] in km
        self.velocity = np.array(velocity)  # [vx, vy, vz] in km/s
        self.com = com

    def create_handover_packet(self, lane_id="L2A_S1"):
        packet = {
            "pos": self.position.tolist(),
            "vel": self.velocity.tolist(),
            "lane_id": lane_id,
            "timestamp": "2026-02-08T13:06:00Z",
        }
        return str(packet)

    def send_secure_packet(self):
        plain_packet = self.create_handover_packet()
        encrypted = self.com.encrypt(plain_packet)
        signature = self.com.sign_message(plain_packet)
        return encrypted, signature, self.com.vk  # share public key


class LeafNode:
    def __init__(self, position, velocity):
        self.position = np.array(position)
        self.velocity = np.array(velocity)

    def receive_and_process(self, encrypted_packet: str, signature: str, vk: VerifyingKey):
        # Decrypt
        com_dummy = SecureSatCom(shared_secret="leaf_starlink_secure_key_2026", sk=None)
        # We only use com_dummy for key derivation; vk is provided separately
        decrypted = com_dummy.decrypt(encrypted_packet)

        # Verify signature
        if not SecureSatCom.verify_signature(decrypted, signature, vk):
            raise ValueError("Signature invalid - transmission compromised!")

        # Parse packet
        packet = ast.literal_eval(decrypted)
        star_pos = np.array(packet["pos"])
        star_vel = np.array(packet["vel"])

        # Compute safety vectors
        rel_pos = self.position - star_pos
        rel_vel = self.velocity - star_vel
        dist = np.linalg.norm(rel_pos)

        dot = np.dot(rel_vel, rel_pos)
        rel_speed_sq = np.linalg.norm(rel_vel) ** 2

        safety_margin = 10.0  # km threshold for LEO lanes

        if rel_speed_sq > 0 and dot < 0:  # Approaching
            tca = -dot / rel_speed_sq
            closest_vector = rel_pos + rel_vel * tca
            closest_dist = np.linalg.norm(closest_vector)

            if closest_dist < safety_margin:
                evasion_dir = np.cross(rel_pos, rel_vel)
                if np.linalg.norm(evasion_dir) > 0:
                    evasion_dir /= np.linalg.norm(evasion_dir)
                evasion_thrust = evasion_dir * 0.5  # placeholder
            else:
                evasion_thrust = np.zeros(3)
        else:
            tca = float("inf")
            closest_dist = dist
            evasion_thrust = np.zeros(3)

        safety_vector = {
            "rel_pos": rel_pos.tolist(),
            "rel_vel": rel_vel.tolist(),
            "tca": float(tca),
            "closest_dist_km": float(closest_dist),
            "evasion_thrust_ms2": evasion_thrust.tolist(),
            "safe": closest_dist > safety_margin,
        }

        return packet, safety_vector


if __name__ == "__main__":
    # Shared comms instance (Starlink side)
    com = SecureSatCom(shared_secret="leaf_starlink_secure_key_2026")

    # Init nodes
    starlink = StarlinkNode(position=[7000, 0, 0], velocity=[0, 7.8, 0], com=com)
    leaf = LeafNode(position=[7050, 10, 5], velocity=[0, 7.7, 0.1])

    # Starlink sends
    encrypted, signature, vk = starlink.send_secure_packet()

    print("Secure Handover Demo:\n")
    print(f"Encrypted: {encrypted[:60]}...")
    print(f"Signature: {signature[:40]}...")

    # Leaf receives
    received_packet, safety = leaf.receive_and_process(encrypted, signature, vk)
    print(f"\nReceived Packet: {received_packet}")
    print(f"Safety Vectors: {safety}")