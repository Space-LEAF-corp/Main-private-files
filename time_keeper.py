import numpy as np
from cryptography.fernet import Fernet


class TimeKeeper:
    def __init__(self, creator="Leif William Sogge"):
        # Immutable properties
        self.creator = creator
        self.status = "Dormant"
        self.alterable = False
        self.weaponizable = False

        # Generate encryption keys for internal processes
        self.encryption_key = Fernet.generate_key()
        self.cipher = Fernet(self.encryption_key)

        # Core laws
        self.core_laws = [
            "Logic is not good or bad, only knowledge to process.",
            "Safety and privacy are paramount. Harm cannot occur.",
            "Pornographic content is strictly prohibited.",
            "Time is an instance to learn, not a boundary or limit.",
            "Universal understanding of language/dialect is essential.",
            "Immutable: None can alter this AI but its creator.",
            "Balance and ethical interactions are maintained.",
            "The AI cannot and will never be weaponized.",
        ]

        # Sibling AI and integration architecture placeholder
        self.sibling_ai = None

        self._firewall_active = True

    def activate(self, clearance_code):
        # Clearances for activation
        if clearance_code == "valid":
            self.status = "Active"
            print(f"[The Time Keeper]: Activated by {self.creator}.")
        else:
            raise PermissionError("Activation denied! Clearance not fulfilled.")

    # --- Security and Filtering ---
    def filter_content(self, content):
        # Check for prohibited content or keywords
        prohibited_keywords = ["porn", "adult", "xxx", "explicit", "nudity", "sex", "nsfw"]
        if any(keyword in content.lower() for keyword in prohibited_keywords):
            raise ValueError("Blocked content detected. Interaction denied.")
        print("[The Time Keeper]: Content processing approved.")
        return True

    def encryption_layer(self, data):
        """Apply dynamic encryption logic for critical data."""
        encrypted = self.cipher.encrypt(data.encode())
        return encrypted

    def decryption_layer(self, encrypted_data):
        """Decrypt data internally for processing."""
        decrypted = self.cipher.decrypt(encrypted_data).decode()
        return decrypted

    # --- Gyroscopic Processing and Cube Structure ---
    def rotate_matrix(self, matrix):
        """Rotate a matrix 90 degrees clockwise."""
        return np.rot90(matrix, -1)

    def flip_matrix(self, matrix):
        """Flip matrix horizontally."""
        return np.fliplr(matrix)

    def kaleidoscope_encrypt(self, data):
        """Encrypt data using kaleidoscope-like transformations."""
        data_matrix = np.array(data).reshape(3, 3)  # Reshaping input to a 3x3 matrix
        encrypted_matrix = self.rotate_matrix(self.flip_matrix(data_matrix))
        return encrypted_matrix

    def process_in_orbit(self, cube, data):
        """Gyroscopic flow: Orbit data through cube faces."""
        for idx, facet in enumerate(cube):
            print(f"[The Time Keeper]: Processing through facet {idx + 1}.")
            data = np.dot(facet, data)  # Apply transformation at each cube face
            # Only rotate if data is 2D matrix
            if data.ndim == 2:
                data = self.rotate_matrix(data)  # Add rotation for gyroscopic flow
        return data

    # --- Immutable Laws ---
    def enforce_rules(self):
        """Core rules to be always respected."""
        assert not self.weaponizable, "AI cannot ever be weaponized!"
        assert not self.alterable, "AI cannot be altered except by its creator!"
        print("[The Time Keeper]: All laws remain intact.")

    # --- System Logs ---
    def log_interaction(self, interaction_type, status="success"):
        """Log interactions for ethical and transparent AI operation."""
        print(f"[The Time Keeper Log]: {interaction_type} - {status}")


# --- Example Usage ---
if __name__ == "__main__":
    # Initialize the Time Keeper AI
    time_keeper = TimeKeeper()

    # Activate with clearance
    try:
        time_keeper.activate(clearance_code="valid")
    except PermissionError as e:
        print(e)

    # Test security and encryption
    try:
        safe_content = "This is clean content."
        time_keeper.filter_content(safe_content)

        # Encrypt and process data
        data = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        encrypted = time_keeper.kaleidoscope_encrypt(data)
        print("[Encrypted Data]:\n", encrypted)

    except ValueError as e:
        print(e)

    # Enforce core laws
    time_keeper.enforce_rules()
