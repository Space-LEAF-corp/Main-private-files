# signal_visualization.py
# Generates time-domain plots, spectrograms, FFTs, and a CSV summary for given frequencies.
# Run locally. Requires: numpy, scipy, matplotlib, pandas
# pip install numpy scipy matplotlib pandas

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import spectrogram, get_window
from scipy.fft import rfft, rfftfreq
import pandas as pd
from datetime import datetime

# --- Parameters ---
SR = 48000
FFT_LEN = 16384
SPEC_NPERSEG = 4096
SPEC_NOOVERLAP = SPEC_NPERSEG - 1024
DUR_SHORT = 1.0
DUR_LONG = 10.0
PEAK_NORMAL = 0.9
OUTPUT_DIR = "signal_outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Frequencies and binaural pairs
base_freqs = [45.0, 48.0, 42.0, 58.256]
pairs = []
for f in base_freqs:
    pairs.append(("single", f, f))
    pairs.append(("binaural", f, f + 3.0))

# Safety note file
with open(os.path.join(OUTPUT_DIR, "SAFETY_NOTE.txt"), "w") as fh:
    fh.write("This sandbox generates only raw audio signals and visualizations. It does not model physiological or cognitive effects. Do not use these signals on people without professional oversight and informed consent.\n")

# Utility functions
def make_tone(freq_left, freq_right, duration, sr=SR):
    t = np.arange(int(duration * sr)) / sr
    left = np.sin(2 * np.pi * freq_left * t)
    right = np.sin(2 * np.pi * freq_right * t)
    stereo = np.vstack([left, right]).T
    # normalize to peak
    peak = np.max(np.abs(stereo))
    if peak > 0:
        stereo = stereo * (PEAK_NORMAL / peak)
    return t, stereo

def plot_waveform(t, stereo, title, fname):
    plt.figure(figsize=(10,4))
    plt.plot(t, stereo[:,0], label='Left', alpha=0.8)
    plt.plot(t, stereo[:,1], label='Right', alpha=0.6)
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.close()

def plot_spectrogram(stereo, sr, title, fname):
    # compute spectrogram of the mean of channels for visualization
    mono = stereo.mean(axis=1)
    f, t, Sxx = spectrogram(mono, fs=sr, window='hann', nperseg=SPEC_NPERSEG, noverlap=SPEC_NOOVERLAP, scaling='density', mode='magnitude')
    plt.figure(figsize=(10,4))
    plt.pcolormesh(t, f, 20*np.log10(Sxx+1e-12), shading='gouraud', cmap='viridis')
    plt.ylim(0, sr/2)
    plt.xlabel("Time [s]")
    plt.ylabel("Frequency [Hz]")
    plt.title(title)
    plt.colorbar(label='Magnitude [dB]')
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.close()

def plot_fft(stereo, sr, title, fname):
    # FFT per channel
    left = stereo[:,0]
    right = stereo[:,1]
    N = len(left)
    L = rfft(left * np.hanning(N), n=FFT_LEN)
    R = rfft(right * np.hanning(N), n=FFT_LEN)
    freqs = rfftfreq(FFT_LEN, 1/sr)
    Lmag = 20*np.log10(np.abs(L)+1e-12)
    Rmag = 20*np.log10(np.abs(R)+1e-12)
    plt.figure(figsize=(10,4))
    plt.plot(freqs, Lmag, label='Left (dB)')
    plt.plot(freqs, Rmag, label='Right (dB)', alpha=0.8)
    plt.xlim(0, sr/2)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude [dB]")
    plt.title(title)
    # find top peaks
    # annotate top 5 peaks per channel
    def top_peaks(mag, freqs, topn=5):
        idx = np.argsort(mag)[-topn:]
        return sorted([(freqs[i], mag[i]) for i in idx], key=lambda x: x[0])
    Lpeaks = top_peaks(Lmag, freqs)
    Rpeaks = top_peaks(Rmag, freqs)
    for f,p in Lpeaks:
        plt.plot(f, p, 'ro', markersize=4)
        plt.text(f, p+3, f"{f:.1f}Hz", color='red', fontsize=8)
    for f,p in Rpeaks:
        plt.plot(f, p, 'go', markersize=4)
        plt.text(f, p+3, f"{f:.1f}Hz", color='green', fontsize=8)
    plt.legend()
    plt.tight_layout()
    plt.savefig(fname, dpi=150)
    plt.close()
    return Lpeaks, Rpeaks

# CSV metrics collection
rows = []

for typ, lf, rf in pairs:
    for dur in [DUR_SHORT, DUR_LONG]:
        t, stereo = make_tone(lf, rf, dur)
        title = f"{typ.upper()} {lf:.3f}Hz / {rf:.3f}Hz — {dur:.0f}s"
        base = f"{typ}_{lf:.3f}_{rf:.3f}_{int(dur)}s".replace('.', 'p')
        wf_name = os.path.join(OUTPUT_DIR, f"{base}_waveform.png")
        plot_waveform(t, stereo, title, wf_name)
        if dur == DUR_LONG:
            spec_name = os.path.join(OUTPUT_DIR, f"{base}_spectrogram.png")
            plot_spectrogram(stereo, SR, title + " Spectrogram", spec_name)
            fft_name = os.path.join(OUTPUT_DIR, f"{base}_fft.png")
            Lpeaks, Rpeaks = plot_fft(stereo, SR, title + " FFT", fft_name)
        else:
            # compute FFT metrics for short signal too
            fft_name = os.path.join(OUTPUT_DIR, f"{base}_fft.png")
            Lpeaks, Rpeaks = plot_fft(stereo, SR, title + " FFT", fft_name)

        # numeric metrics
        rms_left = np.sqrt(np.mean(stereo[:,0]**2))
        rms_right = np.sqrt(np.mean(stereo[:,1]**2))
        peak_left = np.max(np.abs(stereo[:,0]))
        peak_right = np.max(np.abs(stereo[:,1]))
        # format peaks as freq:mag_dB; take top 5
        def fmt_peaks(peaks):
            return ";".join([f"{f:.2f}:{m:.2f}" for f,m in peaks])
        rows.append({
            "signal_id": base,
            "type": typ,
            "base_frequency": lf,
            "left_freq": lf,
            "right_freq": rf,
            "duration_s": dur,
            "sample_rate": SR,
            "rms_left": rms_left,
            "rms_right": rms_right,
            "peak_left": peak_left,
            "peak_right": peak_right,
            "dominant_peaks_left": fmt_peaks(Lpeaks),
            "dominant_peaks_right": fmt_peaks(Rpeaks)
        })

# Save CSV
df = pd.DataFrame(rows)
csv_path = os.path.join(OUTPUT_DIR, "metrics_summary.csv")
df.to_csv(csv_path, index=False)
print("Outputs written to", OUTPUT_DIR)
print("CSV summary:", csv_path)
