import getpass
from dataclasses import dataclass, field
from typing import Dict, Any


# -----------------------------
# Data model: Leaf NFC Credential
# -----------------------------

@dataclass
class LeafCredential:
    identity_core: Dict[str, Any]
    nine_ports: Dict[str, Any]
    role: str  # "admin" or "user"

    def describe(self) -> None:
        print("\n=== LEAF NFC CREDENTIAL ===")
        print(f"Holder: {self.identity_core.get('display_name')}")
        print(f"User ID: {self.identity_core.get('user_id')}")
        print(f"Role:    {self.role}")
        print("\nNine Ports (Capabilities):")
        for port_name, port_value in self.nine_ports.items():
            print(f" - {port_name}: {port_value}")
        print("============================\n")


# -----------------------------
# Demo configuration
# -----------------------------

# In a real system, this would be securely stored / derived.
ADMIN_PASSWORD = "star-admin-demo"  # change this locally if you want

def build_admin_credential() -> LeafCredential:
    return LeafCredential(
        identity_core={
            "user_id": "space_leaf://user/leif",
            "display_name": "Leif Sogge",
            "issuer": "Space LEAF Corp",
            "lineage": ["docusign://demo-event", "local-device://seed"],
        },
        nine_ports={
            "port_1_auth_key": "admin-signing-key-ref",
            "port_2_encrypt_key": "admin-encryption-key-ref",
            "port_3_device_binding": "this-device-only",
            "port_4_guardian_profile": "PanTar-Guardian-Enabled",
            "port_5_permissions_profile": "full-access",
            "port_6_comm_channel": "local-only",
            "port_7_cloud_link": "disabled-in-demo",
            "port_8_recovery": "seed-phrase-stub",
            "port_9_personal_signature": "LWS-admin-style",
        },
        role="admin",
    )


def build_user_credential() -> LeafCredential:
    return LeafCredential(
        identity_core={
            "user_id": "space_leaf://user/demo-user",
            "display_name": "Standard User",
            "issuer": "Space LEAF Corp",
            "lineage": ["local-device://demo"],
        },
        nine_ports={
            "port_1_auth_key": "user-signing-key-ref",
            "port_2_encrypt_key": "user-encryption-key-ref",
            "port_3_device_binding": "this-device-only",
            "port_4_guardian_profile": "PanTar-Guardian-Required",
            "port_5_permissions_profile": "limited-access",
            "port_6_comm_channel": "local-only",
            "port_7_cloud_link": "disabled-in-demo",
            "port_8_recovery": "support-contact-stub",
            "port_9_personal_signature": "User-style",
        },
        role="user",
    )


# -----------------------------
# Authentication logic
# -----------------------------

def authenticate_as_admin() -> LeafCredential | None:
    print("\n[ADMIN PATH] Star key detected.")
    pwd = getpass.getpass("Enter admin password: ")
    if pwd == ADMIN_PASSWORD:
        print("Admin authentication successful.\n")
        return build_admin_credential()
    else:
        print("Admin authentication failed. Falling back to standard user.\n")
        return None


def authenticate_as_user() -> LeafCredential:
    print("\n[USER PATH] Standard user mode selected.")
    # In a real system, you might still ask for a PIN or local check.
    return build_user_credential()


# -----------------------------
# Main pre-demo flow
# -----------------------------

def main():
    print("======================================")
    print("  SPACE LEAF CORP - NFC eSIM PRE-DEMO ")
    print("======================================")
    print("This demo simulates a user-forged NFC credential.")
    print("Press '*' for ADMIN mode, or any other key for USER mode.")
    print("--------------------------------------")

    choice = input("Input: ").strip()

    # Default behavior: you usually log in as admin,
    # so we treat '*' as the privileged path.
    if choice == "*":
        cred = authenticate_as_admin()
        if cred is None:
            # Fallback to user if admin auth fails
            cred = authenticate_as_user()
    else:
        cred = authenticate_as_user()

    # Show the credential that would be "loaded"
    cred.describe()

    # Here is where you'd "boot" the rest of your system
    # using the credential's role and capabilities.
    if cred.role == "admin":
        print("Booting system in ADMIN mode (pre-demo)...")
        # Placeholder for your admin modules
        # e.g., load_admin_modules(cred)
    else:
        print("Booting system in USER mode (pre-demo)...")
        # Placeholder for your user modules
        # e.g., load_user_modules(cred)

    print("\nPre-demo complete. System ready.")


if __name__ == "__main__":
    main()