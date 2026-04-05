"""
from SCIPY_COMPAT import signal
Lightweight, NumPy-based replacements for a few scipy.signal functions:
- chirp
- correlate
- correlation_lags

Designed for iOS / constrained environments where SciPy wheels can't be installed.
"""

import numpy as np


class _SignalCompat:
    # -----------------------------
    # CHIRP
    # -----------------------------
    @staticmethod
    def chirp(t, f0, f1, t1, method="linear", phi=0.0):
        """
        Generate a frequency sweep (chirp) signal.

        Parameters:
            t      : time array (seconds)
            f0     : starting frequency (Hz)
            f1     : ending frequency (Hz)
            t1     : time at which f1 is reached (seconds)
            method : only 'linear' supported here
            phi    : initial phase in degrees

        Returns:
            y      : chirp signal, same shape as t
        """
        t = np.asarray(t, dtype=float)
        phi_rad = np.deg2rad(phi)

        if method != "linear":
            raise ValueError("Only 'linear' chirp is implemented in SCIPY_COMPAT.")

        # Linear chirp: f(t) = f0 + (f1 - f0) * t / t1
        k = (f1 - f0) / float(t1)
        phase = 2.0 * np.pi * (f0 * t + 0.5 * k * t**2) + phi_rad
        y = np.cos(phase)
        return y.astype(np.float32)

    # -----------------------------
    # CORRELATE
    # -----------------------------
    @staticmethod
    def correlate(in1, in2, mode="full"):
        """
        Cross-correlation of two 1D sequences.

        Parameters:
            in1, in2 : input arrays
            mode     : 'full', 'same', or 'valid'

        Returns:
            correlation array
        """
        x = np.asarray(in1, dtype=float)
        y = np.asarray(in2, dtype=float)

        # Cross-correlation is convolution with one sequence reversed
        y_rev = y[::-1]
        corr_full = np.convolve(x, y_rev, mode="full")

        if mode == "full":
            return corr_full.astype(np.float32)
        elif mode == "same":
            # Center to match length of x
            start = (len(corr_full) - len(x)) // 2
            end = start + len(x)
            return corr_full[start:end].astype(np.float32)
        elif mode == "valid":
            # Valid region where signals fully overlap
            valid_len = max(len(x), len(y)) - min(len(x), len(y)) + 1
            start = (len(corr_full) - valid_len) // 2
            end = start + valid_len
            return corr_full[start:end].astype(np.float32)
        else:
            raise ValueError("mode must be 'full', 'same', or 'valid'")

    # -----------------------------
    # CORRELATION_LAGS
    # -----------------------------
    @staticmethod
    def correlation_lags(in1_len, in2_len, mode="full"):
        """
        Return lag indices for 1D cross-correlation.

        Parameters:
            in1_len : length of first sequence
            in2_len : length of second sequence
            mode    : 'full', 'same', or 'valid'

        Returns:
            lags : array of integer lags
        """
        if mode == "full":
            # From -(in2_len-1) to +(in1_len-1)
            lags = np.arange(-(in2_len - 1), in1_len)
        elif mode == "same":
            # Centered around zero, length = max(in1_len, in2_len)
            out_len = max(in1_len, in2_len)
            mid = (out_len - 1) // 2
            lags = np.arange(-mid, -mid + out_len)
        elif mode == "valid":
            # Only positions where signals fully overlap
            out_len = max(in1_len, in2_len) - min(in1_len, in2_len) + 1
            lags = np.arange(0, out_len)
        else:
            raise ValueError("mode must be 'full', 'same', or 'valid'")

        return lags.astype(int)


# Public interface: mimic scipy.signal namespace
signal = _SignalCompat()