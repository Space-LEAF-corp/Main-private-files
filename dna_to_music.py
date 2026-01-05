# dna_to_music.py
# Encodes frequencies into a DNA string, maps to MIDI notes and chords, writes MIDI, and optionally renders a simple organ-like WAV.
# Run locally. Requires: numpy, pretty_midi, mido, soundfile (optional for WAV)
# pip install numpy pretty_midi mido soundfile

import numpy as np
import pretty_midi
import math
import hashlib
import json
import soundfile as sf
import os

OUTPUT_DIR = "music_outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Safety note
with open(os.path.join(OUTPUT_DIR, "SAFETY_NOTE.txt"), "w") as fh:
    fh.write("This sandbox generates only MIDI, simple WAV renders, and textual artifacts. It does not model physiological or cognitive effects. Do not use these sounds on people without professional oversight and informed consent.\n")

# --- Encoding helpers ---
def freq_to_midi(f):
    return 69 + 12 * math.log2(f / 440.0)

def quantize_freq(f, step=0.01):
    return round(f / step) * step

def to_base4_bytes(value, bits=16):
    s = ''
    for i in range(bits//2):
        s = 'ACGT'[value & 0b11] + s
        value >>= 2
    return s

def freq_to_dna(f, bits=16, lo=20.0, hi=20000.0):
    v = int(((min(max(f, lo), hi) - lo) / (hi - lo)) * (2**bits - 1))
    return to_base4_bytes(v, bits=bits)

def dna_to_word(dna_sub):
    bank = ["whisper","pause","spark","riddle","breathe","orbit","echo","thread","glow","knot","lumen","tide"]
    h = int(hashlib.sha256(dna_sub.encode()).hexdigest(), 16)
    return bank[h % len(bank)]

# --- Music mapping ---
def build_midi_from_freqs(freqs, filename_midi, tempo_bpm=60, instrument_program=19):
    # instrument_program: General MIDI program number (0-based). 19 ~ Church Organ in many mappings.
    pm = pretty_midi.PrettyMIDI(initial_tempo=tempo_bpm)
    inst = pretty_midi.Instrument(program=instrument_program)
    base_time = 0.0
    total_dur = 8.0
    note_dur = total_dur / max(1, len(freqs))
    for f in freqs:
        midi_note = int(round(freq_to_midi(f)))
        midi_note = max(0, min(127, midi_note))
        note = pretty_midi.Note(velocity=80, pitch=midi_note, start=base_time, end=base_time + note_dur)
        inst.notes.append(note)
        base_time += note_dur
    pm.instruments.append(inst)
    pm.write(filename_midi)
    return filename_midi

# Simple additive organ-like render (mono) for preview only
def render_simple_organ_mono(freqs, filename_wav, sr=48000, dur=8.0):
    t = np.linspace(0, dur, int(sr*dur), endpoint=False)
    out = np.zeros_like(t)
    # additive partials typical of organ: fundamental + several harmonics with decreasing amplitude
    for f in freqs:
        for h, amp in [(1,1.0),(2,0.6),(3,0.35),(4,0.2),(6,0.12)]:
            out += amp * np.sin(2*np.pi*(f*h)*t)
    # normalize
    peak = np.max(np.abs(out))
    if peak > 0:
        out = out * (0.9 / peak)
    sf.write(filename_wav, out, sr)
    return filename_wav

# --- Main pipeline ---
def pipeline(freqs, user_name="listener"):
    # Encode DNA
    dna_pieces = [freq_to_dna(f) for f in freqs]
    dna = ''.join(dna_pieces)
    # Build MIDI
    midi_file = os.path.join(OUTPUT_DIR, f"{user_name}_generated_song.mid")
    build_midi_from_freqs(freqs, midi_file)
    # Optional WAV render with organ undertones
    wav_file = os.path.join(OUTPUT_DIR, f"{user_name}_generated_song.wav")
    render_simple_organ_mono(freqs, wav_file)
    # Build riddle
    tokens = [dna_to_word(dna[i:i+6]) for i in range(0, min(len(dna), 36), 6)]
    while len(tokens) < 3:
        tokens.append("echo")
    riddle = f"I carry a {tokens[0]} and a {tokens[1]}; between them hums a secret {tokens[2]}."
    # Manifest
    manifest = {
        "user": user_name,
        "frequencies": freqs,
        "dna": dna,
        "midi": midi_file,
        "wav_preview": wav_file,
        "riddle": riddle,
        "timestamp": datetime.utcnow().isoformat() + "Z"
    }
    manifest_path = os.path.join(OUTPUT_DIR, f"{user_name}_manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print("Generated:", manifest_path)
    print("Riddle:", riddle)
    return manifest_path

if __name__ == "__main__":
    from datetime import datetime
    # Example usage
    freqs = [45.0, 48.0, 42.0, 58.256]
    pipeline(freqs, user_name="leif_example")
    # Print consent template
    consent = """
Consent for voluntary listening trial
- This listening session is voluntary. You may stop at any time.
- The audio is experimental and intended for creative use only. No claims are made about altering brain function.
- Do not participate if you have epilepsy, are pregnant, have a pacemaker, or other medical concerns without medical advice.
- If you feel unwell, stop immediately.

Please initial to consent: ______
"""
    print(consent)
