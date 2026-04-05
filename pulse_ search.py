import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import datetime


class PulseSearch:
    def __init__(self):
        self.sample_rate = 44100  # Hz
        self.duration = 0.02      # seconds per pulse
        self.volume = 0.8         # [0.0, 1.0]

        # Multi-frequency sweep settings (Hz)
        self.sweep_bands = [
            (500, 1500),
            (1500, 3000),
            (3000, 6000),
        ]

        # Simulated environment
        self.echo_delay_sec = 0.004   # 4 ms delay
        self.echo_attenuation = 0.5   # echo strength
        self.noise_level = 0.1        # noise amplitude

    def _time_axis(self, duration=None):
        if duration is None:
            duration = self.duration
        n_samples = int(self.sample_rate * duration)
        t = np.arange(n_samples) / self.sample_rate
        return t

    def generate_tone_pulse(self, frequency):
        """Simple single-frequency sine pulse."""
        t = self._time_axis()
        pulse = self.volume * np.sin(2 * np.pi * frequency * t).astype(np.float32)
        return pulse

    def generate_chirp_pulse(self, f_start, f_end):
        """Frequency sweep (chirp) pulse."""
        t = self._time_axis()
        pulse = signal.chirp(t, f0=f_start, f1=f_end, t1=self.duration, method='linear')
        pulse = (self.volume * pulse).astype(np.float32)
        return pulse

    def simulate_channel(self, pulse):
        """
        Simulate a reflective environment:
        - delayed, attenuated echo
        - additive noise
        """
        n = len(pulse)
        delay_samples = int(self.echo_delay_sec * self.sample_rate)

        # Base received signal: just noise
        received = np.random.normal(0, self.noise_level, n)

        # Direct path (optional: comment out if you want only echo)
        received[:n] += pulse

        # Echo path (delayed, attenuated)
        if delay_samples < n:
            received[delay_samples:] += self.echo_attenuation * pulse[:n - delay_samples]

        return received.astype(np.float32)

    def detect_pulse(self, pulse, received):
        """
        Matched filter via correlation.
        Returns:
            corr: correlation array
            peak_index: index of max correlation
            peak_time: time (seconds) of detected peak
        """
        corr = signal.correlate(received, pulse, mode='full')
        lags = signal.correlation_lags(len(received), len(pulse), mode='full')
        peak_index = np.argmax(corr)
        peak_lag = lags[peak_index]
        peak_time = peak_lag / self.sample_rate
        return corr, peak_index, peak_time

    def send_and_receive(self, pulse):
        print("Pulse sent at", datetime.datetime.now())
        received = self.simulate_channel(pulse)
        return received

    def run_sweep_search(self, plot=True):
        """
        For each sweep band:
        - generate chirp pulse
        - send through simulated channel
        - detect echo via correlation
        - report detected delay
        """
        results = []

        for i, (f_start, f_end) in enumerate(self.sweep_bands, start=1):
            print(f"\n=== Sweep {i}: {f_start} Hz → {f_end} Hz ===")
            pulse = self.generate_chirp_pulse(f_start, f_end)
            received = self.send_and_receive(pulse)
            corr, peak_idx, peak_time = self.detect_pulse(pulse, received)

            print(f"Detected peak at ~{peak_time * 1000:.2f} ms")

            results.append({
                "band": (f_start, f_end),
                "peak_time_sec": peak_time,
                "pulse": pulse,
                "received": received,
                "corr": corr
            })

            if plot:
                self._plot_detection(pulse, received, corr, f_start, f_end)

        print("\nMulti-frequency pulse sweep search complete.")
        return results

    def _plot_detection(self, pulse, received, corr, f_start, f_end):
        t = self._time_axis()
        t_received = np.arange(len(received)) / self.sample_rate
        t_corr = np.arange(len(corr)) / self.sample_rate  # relative, just for visualization

        plt.figure(figsize=(12, 8))

        plt.subplot(3, 1, 1)
        plt.title(f"Pulse (Chirp {f_start}–{f_end} Hz)")
        plt.plot(t, pulse)
        plt.xlabel("Time [s]")
        plt.ylabel("Amplitude")

        plt.subplot(3, 1, 2)
        plt.title("Received Signal (with echo + noise)")
        plt.plot(t_received, received)
        plt.xlabel("Time [s]")
        plt.ylabel("Amplitude")

        plt.subplot(3, 1, 3)
        plt.title("Correlation (Matched Filter Output)")
        plt.plot(t_corr, corr)
        plt.xlabel("Relative Time [s]")
        plt.ylabel("Correlation")

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    ps = PulseSearch()
    # for now, just test a single chirp band:
    f_start, f_end = ps.sweep_bands[0]
    pulse = ps.generate_chirp_pulse(f_start, f_end)
    received = ps.simulate_channel(pulse)
    # later you can wire in plotting or detection here