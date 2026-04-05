#!/usr/bin/env python3
"""
tools/scan_qr_batch.py

Scans all PNG/JPG files in a folder and attempts to decode QR payloads using pyzbar and OpenCV fallback.
Outputs per-file result and an aggregate CSV summary suitable for the device matrix.
"""
import argparse
from pathlib import Path
from pyzbar.pyzbar import decode as zbar_decode
from PIL import Image
import cv2
import numpy as np
import csv
import base64
import sys

def decode_with_pyzbar(path):
    try:
        img = Image.open(path).convert("L")
        res = zbar_decode(img)
        if res:
            return res[0].data.decode()
    except Exception:
        return None
    return None

def decode_with_opencv(path):
    try:
        img = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
        detector = cv2.QRCodeDetector()
        data, points, _ = detector.detectAndDecode(img)
        if data:
            return data
    except Exception:
        return None
    return None

def analyze_payload(payload_b64):
    try:
        blob = base64.b64decode(payload_b64)
        parts = blob.split(b"||")
        return len(blob), len(parts)
    except Exception:
        return None, None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("folder", help="Folder with QR images")
    parser.add_argument("--out", default="scan_results.csv", help="CSV output")
    args = parser.parse_args()
    folder = Path(args.folder)
    files = list(folder.glob("*.png")) + list(folder.glob("*.jpg")) + list(folder.glob("*.jpeg"))
    rows = []
    success_count = 0
    for f in files:
        payload = decode_with_pyzbar(f)
        method = "pyzbar"
        if payload is None:
            payload = decode_with_opencv(f)
            method = "opencv" if payload else None
        success = bool(payload)
        if success:
            success_count += 1
            size_bytes, parts = analyze_payload(payload)
        else:
            size_bytes, parts = None, None
        rows.append({
            "file": f.name,
            "decoded": success,
            "method": method,
            "payload_size_bytes": size_bytes,
            "payload_parts": parts
        })
        print(f"{f.name}: decoded={success} method={method} size={size_bytes}")
    # write CSV
    with open(args.out, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["file","decoded","method","payload_size_bytes","payload_parts"])
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
    total = len(files)
    print(f"\nScanned {total} images. Success rate: {success_count}/{total} = {success_count/total if total else 0:.2f}")
    print(f"Results written to {args.out}")

if __name__ == "__main__":
    main()