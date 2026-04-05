import time
import sys

# -----------------------------------------
# SPACE LEAF CORP — ACCESS PROTOCOL CHECK
# -----------------------------------------

ROLES = {
    "P": "Professional",
    "C": "Captain",
    "G": "Government Liaison",
    "S": "Student",
    "A": "Administrator",
    "M": "Maintenance",
    "E": "Emergency Bypass"
}

MAINTENANCE_MAX_TIME = 20  # seconds for demo purposes


def boot_sequence():
    print("\n🟢 SPACE LEAF CORP SYSTEM BOOTING...")
    time.sleep(0.5)
    print("🔧 Running diagnostics...")
    time.sleep(0.5)
    print("📡 Checking CARENY logic loop...")
    time.sleep(0.5)
    print("🌿 Core Law Verified: CARE")
    time.sleep(0.5)
    print("🚽 Loading DDD‐333‐X Access Protocol...\n")
    time.sleep(0.5)


def parse_access_code(code):
    try:
        ddd, triple, role = code.split("-")
    except ValueError:
        return None, None, None

    if ddd != "DDD" or triple != "333":
        return None, None, None

    return ddd, triple, role.upper()


def maintenance_session():
    print("\n🛠️ Maintenance session started.")
    print(f"⏱️ Auto‐logout in {MAINTENANCE_MAX_TIME} seconds...\n")

    for i in range(MAINTENANCE_MAX_TIME, 0, -1):
        sys.stdout.write(f"\r⏳ {i} seconds remaining...")
        sys.stdout.flush()
        time.sleep(1)

    print("\n\n🔒 Auto‐logout complete.")
    print("🧼 Maintenance session logged safely.\n")


def main():
    boot_sequence()

    access_code = input("Enter access code (format: DDD-333-X): ").strip()
    ddd, triple, role = parse_access_code(access_code)

    if not role or role not in ROLES:
        print("\n❌ Invalid access code.")
        return

    print(f"\n🔑 Access recognized: {ROLES[role]}")

    if role == "M":
        maintenance_session()
    elif role == "E":
        print("\n🚨 EMERGENCY BYPASS ACTIVATED")
        print("⚠️ Logging override event for audit.\n")
    else:
        print("\n🚪 Access granted.")
        print("📝 Session logged.\n")


if __name__ == "__main__":
    main()