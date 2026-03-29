#!/usr/bin/env python3
import sys
import re
from pathlib import Path

TARGET_FILE = "jarvondis3.0_base.py"

FORBIDDEN_IN_JARVONDIS = {
    "Emotional framing (core laws or responses)": [
        r"\bcare\b",
        r"\bempathy\b",
        r"\bconnection\b",
        r"\bunderstand emotions\b",
        r"\bsupport\b",
    ],
    "Conversational / parasocial phrases": [
        r"i[' ]?m hearing you",
        r"how can i help",
        r"i understand",
        r"that sounds",
        r"i[' ]?m here",
        r"i[' ]?m glad you",
    ],
    "Memory profiling or inference": [
        r"infer",
        r"profile",
        r"sentiment",
        r"emotional_state",
        r"analyze_user",
    ],
    "Ethics bypass or configuration": [
        r"disable_ethics",
        r"ethics_enabled\s*=\s*false",
        r"override_safety",
        r"bypass_ethics",
        r"os\.environ",
        r"config\.",
    ],
}

ALLOWED_SILENCE_RETURNS = [
    r"return\s+['\"]{0,2}",
    r"return\s+None",
]

def check_jarvondis_file(path: Path):
    violations = []
    content = path.read_text(errors="ignore").lower()

    # Enforce silence by default
    if "explicitly_requested_assistance" not in content:
        if re.search(r"return\s+[^\s]", content):
            violations.append(
                f"{path}: Default response without explicit request"
            )

    # Check forbidden patterns
    for category, patterns in FORBIDDEN_IN_JARVONDIS.items():
        for pattern in patterns:
            if re.search(pattern, content):
                violations.append(
                    f"{path}: {category} → '{pattern}'"
                )

    return violations

def main():
    violations = []

    for file_path in sys.argv[1:]:
        path = Path(file_path)
        if path.name == TARGET_FILE and path.exists():
            violations.extend(check_jarvondis_file(path))

    if violations:
        print("\n❌ COMMIT BLOCKED — JARVONDIS-SPECIFIC ETHICS VIOLATION\n")
        print("This file enforces silent, child‑safe, steward‑grade behavior.\n")
        for v in violations:
            print(f"  • {v}")
        print(
            "\nFix the issues above before committing.\n"
            "Jarvondis ethics are non‑negotiable.\n"
        )
        sys.exit(1)

    sys.exit(0)

if __name__ == "__main__":
    main()
