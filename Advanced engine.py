#!/usr/bin/env python3
"""
advanced_random_security_qr_code_generator_engine.py

Modular QR generator with layered visual overlays.
Designed for mobile execution (Pythonista, Termux, etc.)
"""

import random
from datetime import datetime
import qrcode
from PIL import Image, ImageDraw


def generate_layered_qr(
    data: str,
    out_file: str = None,
    size: int = 800,
    overlay_intensity: int = 40,
    seed: int | None = None,
) -> str:
    """
    Generate a QR code with layered visual overlays.

    :param data: Text/URL to encode.
    :param out_file: Output PNG path. If None, auto-named.
    :param size: Final image size (square).
    :param overlay_intensity: 0–255, how strong overlays appear.
    :param seed: Optional random seed for reproducibility.
    :return: Path to saved PNG.
    """
    if seed is not None:
        random.seed(seed)

    # --- Base QR generation ---
    qr = qrcode.QRCode(
        version=None,
        error_correction=qrcode.constants.ERROR_CORRECT_H,
        box_size=10,
        border=4,
    )
    qr.add_data(data)
    qr.make(fit=True)

    base_img = qr.make_image(fill_color="black", back_color="white").convert("RGBA")
    base_img = base_img.resize((size, size), Image.NEAREST)

    # --- Overlay canvas ---
    overlay = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)

    # --- Layer 1: radial shading ---
    center = (size // 2, size // 2)
    max_radius = size // 2
    steps = 12

    for i in range(steps):
        radius = int(max_radius * (i + 1) / steps)
        alpha = int(overlay_intensity * (1 - i / steps))
        color = (0, 0, 0, alpha)
        bbox = [
            center[0] - radius,
            center[1] - radius,
            center[0] + radius,
            center[1] + radius,
        ]
        draw.ellipse(bbox, outline=color, width=2)

    # --- Layer 2: random glyphs ---
    for _ in range(25):
        shape = random.choice(["circle", "rect", "line"])
        x1, y1 = random.randint(0, size), random.randint(0, size)
        x2, y2 = random.randint(0, size), random.randint(0, size)
        alpha = random.randint(10, overlay_intensity)
        color = (0, 0, 0, alpha)

        if shape == "circle":
            r = random.randint(size // 40, size // 15)
            draw.ellipse([x1 - r, y1 - r, x1 + r, y1 + r], outline=color, width=2)
        elif shape == "rect":
            draw.rectangle([x1, y1, x2, y2], outline=color, width=2)
        else:
            draw.line([x1, y1, x2, y2], fill=color, width=2)

    # --- Composite final image ---
    final_img = Image.alpha_composite(base_img, overlay)

    # --- Save output ---
    if out_file is None:
        timestamp = datetime.utcnow().strftime("%Y%m%d_%H%M%S")
        out_file = f"advanced_qr_{timestamp}.png"

    final_img.save(out_file)
    return out_file


# --- Direct run block for Pythonista ---
if __name__ == "__main__":
    # Customize this line with your test input
    test_data = "https://example.com"
    print("[*] Generating advanced layered QR...")
    saved_path = generate_layered_qr(test_data)
    print(f"[+] Saved QR to: {saved_path}")
    print("[*] Scan it with your phone’s camera or QR app.")