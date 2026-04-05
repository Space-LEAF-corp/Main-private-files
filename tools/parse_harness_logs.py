#!/usr/bin/env python3
"""
tools/parse_harness_logs.py

Parses harness or pytest output and prints a summary table and JSON report.
Expected to find lines like:
  Roundtrip success rate: 1.00
  Corruption 5% decode success: 0.98
  Verify failures: 0
  Tamper detections: 50
It also extracts timing and KDF parameter lines if present.
"""
import sys
import re
import json
from pathlib import Path
from collections import OrderedDict

PATTERNS = {
    "roundtrip": re.compile(r"Roundtrip success rate:\s*([0-9.]+)"),
    "corruption": re.compile(r"Corruption\s+(\d+)%\s+decode success:\s*([0-9.]+)"),
    "verify_failures": re.compile(r"Verify failures:\s*(\d+)"),
    "tamper": re.compile(r"Tamper detections:\s*(\d+)"),
    "kdf": re.compile(r"KDF_PARAMS.*N\s*=\s*2\*\*(\d+)"),
    "throughput": re.compile(r"throughput\s*`\(?ops/sec\)`?:\s*([0-9.]+)"),
    "rounds": re.compile(r"Generated\s+(\d+)\s+payloads"),
}

def parse_file(path):
    text = Path(path).read_text()
    summary = OrderedDict()
    m = PATTERNS["roundtrip"].search(text)
    summary["roundtrip_success_rate"] = float(m.group(1)) if m else None

    # corruption: collect multiple
    corruption = {}
    for m in PATTERNS["corruption"].finditer(text):
        pct = int(m.group(1))
        val = float(m.group(2))
        corruption[f"{pct}%"] = val
    summary["corruption_decode_success"] = corruption

    m = PATTERNS["verify_failures"].search(text)
    summary["verify_failures"] = int(m.group(1)) if m else 0
    m = PATTERNS["tamper"].search(text)
    summary["tamper_detections"] = int(m.group(1)) if m else 0
    m = PATTERNS["kdf"].search(text)
    summary["kdf_N_bits"] = int(m.group(1)) if m else None
    m = PATTERNS["throughput"].search(text)
    summary["throughput_ops_per_sec"] = float(m.group(1)) if m else None
    m = PATTERNS["rounds"].search(text)
    summary["generated_payloads"] = int(m.group(1)) if m else None
    return summary

def print_table(summary):
    print("\n--- Test Summary ---\n")
    print(f"**Roundtrip success rate**: {summary.get('roundtrip_success_rate')}")
    print(f"**Verify failures**: {summary.get('verify_failures')}")
    print(f"**Tamper detections**: {summary.get('tamper_detections')}")
    print(f"**Throughput ops/sec**: {summary.get('throughput_ops_per_sec')}")
    print(f"**KDF N exponent**: {summary.get('kdf_N_bits')}")
    print(f"**Generated payloads**: {summary.get('generated_payloads')}\n")
    print("**Corruption decode success by damage**:")
    for k, v in summary.get("corruption_decode_success", {}).items():
        print(f"  {k}: {v:.2f}")
    print("\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: python parse_harness_logs.py /path/to/harness_results.txt")
        sys.exit(2)
    path = sys.argv[1]
    summary = parse_file(path)
    print_table(summary)
    # write JSON report next to input file
    out = Path(path).with_suffix(".summary.json")
    out.write_text(json.dumps(summary, indent=2))
    print(f"JSON summary written to {out}")

if __name__ == "__main__":
    main()