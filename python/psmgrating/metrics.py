"""psmgrating.metrics -- the paper's three figures of merit, eqs. (172)-(174)."""
import numpy as np

__all__ = ["rms", "peak_error", "fwhm", "width_error"]


def rms(eta_model, eta_ref):
    """Root-mean-square departure over the band, eq. (172)."""
    m = np.asarray(eta_model, dtype=float)
    r = np.asarray(eta_ref, dtype=float)
    return float(np.sqrt(np.mean((m - r) ** 2)))


def peak_error(eta_model, eta_ref):
    """Error in the peak diffraction efficiency, eq. (173)."""
    return float(abs(np.max(eta_model) - np.max(eta_ref)))


def fwhm(lam, eta):
    """Full width at half maximum about the global maximum, located by linear
    interpolation between the bracketing grid points (text after eq. 174).
    Returns nan if the response does not fall to half on one side."""
    lam = np.asarray(lam, dtype=float)
    eta = np.asarray(eta, dtype=float)
    j = int(np.argmax(eta))
    half = eta[j] / 2.0
    left = np.nan
    for i in range(j, 0, -1):
        if eta[i - 1] <= half < eta[i] or eta[i - 1] < half <= eta[i]:
            frac = (eta[i] - half) / (eta[i] - eta[i - 1])
            left = lam[i] - frac * (lam[i] - lam[i - 1])
            break
    right = np.nan
    for i in range(j, len(eta) - 1):
        if eta[i + 1] <= half < eta[i] or eta[i + 1] < half <= eta[i]:
            frac = (eta[i] - half) / (eta[i] - eta[i + 1])
            right = lam[i] + frac * (lam[i + 1] - lam[i])
            break
    return float(right - left)


def width_error(lam, eta_model, eta_ref):
    """Fractional error in the band width, eq. (174)."""
    wr = fwhm(lam, eta_ref)
    wm = fwhm(lam, eta_model)
    return float(abs(wm - wr) / wr)
