"""psmgrating.chain -- the exactly unslanted grating as a stratified medium.

At K_x' = 0 every Floquet order shares one tangential wavevector, so the
orders are not independent exterior channels and |r_-1|^2 is not an
observable (paper Appendix A.6).  The plate is then a one-dimensional
stratified medium: one fringe period is built from 1024 sublayers sampled at
their midpoints with the exact layer propagator (eq. 188), the plate is
formed from that period matrix by binary exponentiation with rescaling, and
a separately constructed partial cell covers the remainder d - floor(d /
Lambda) Lambda.  The reflection coefficient is homogeneous in the total
matrix, so the rescaling cancels from it exactly.

Propagator conventions follow validation/rcw_tm/kx0_stratified/
kx0_stratified.py (validated to 7e-14 against an exact three-layer Abeles
slab and to 5e-11 against a DOP853 integration of the continuous profile):
  TE (sigma): (h, v) = (E, Z0 H_x), s = -1, q_j = k_y,j / k0
  TM (pi)   : (h, v) = (H, E_x / Z0), s = +1, q_j = k_y,j / (k0 eps_j)
  layer matrix [[cos D, -s i sin D / q_j], [-s i q_j sin D, cos D]], D = k_y,j t
"""
import numpy as np

from . import conventions as cv

__all__ = ["solve", "sweep"]

SUBLAYERS_PER_PERIOD = 1024


def _layer_products(k0, kx2, n0g, n1g, ks, y_mid, t, s):
    """Multiply the propagators of the sublayers centred at y_mid (ascending
    y), thickness t, for all wavelengths at once.  Returns (A, B, C, D)."""
    A = np.ones_like(k0, dtype=complex)
    B = np.zeros_like(A)
    C = np.zeros_like(A)
    D = np.ones_like(A)
    k02 = k0 * k0
    for y in y_mid:
        eps = (n0g + n1g * np.cos(ks * y)) ** 2
        ky = np.sqrt(k02 * eps - kx2 + 0j)
        q = ky / k0 if s < 0 else ky / (k0 * eps)
        dl = ky * t
        c = np.cos(dl)
        sn = np.sin(dl)
        m12 = -s * 1j * sn / q
        m21 = -s * 1j * q * sn
        A, B, C, D = (c * A + m12 * C, c * B + m12 * D,
                      m21 * A + c * C, m21 * B + c * D)
    return A, B, C, D


def _matmul(P, Qm):
    """2x2 product P @ Qm of per-wavelength matrices given as (A, B, C, D)."""
    A1, B1, C1, D1 = P
    A2, B2, C2, D2 = Qm
    return (A1 * A2 + B1 * C2, A1 * B2 + B1 * D2,
            C1 * A2 + D1 * C2, C1 * B2 + D1 * D2)


def _rescale(P):
    A, B, C, D = P
    scale = np.maximum.reduce([np.abs(A), np.abs(B), np.abs(C), np.abs(D)])
    scale = np.where(scale == 0.0, 1.0, scale)
    return (A / scale, B / scale, C / scale, D / scale), np.log(scale)


def solve(*, lam_c, lam_r, n0, n1, d, theta_c, theta_r, pol="sigma",
          D0=0.0, D1=0.0, chi0=None, chi1=None, sublayers=SUBLAYERS_PER_PERIOD):
    """Reflectance and transmittance of the exactly unslanted plate (Psi = 0)
    at the replay wavelengths lam_c (array).  Returns dict(eta_R, eta_T, r, t)."""
    lam = np.atleast_1d(np.asarray(lam_c, dtype=float))
    kk0 = 2.0 * np.pi / lam
    b = kk0 * n0
    kx = b * np.sin(theta_c)
    c0 = cv.chi_from_density(D0, lam, theta_c, 0.0, d) if chi0 is None else chi0 + 0.0 * lam
    c1 = cv.chi_from_density(D1, lam, theta_c, 0.0, d) if chi1 is None else chi1 + 0.0 * lam
    n0g = n0 - 1j * c0
    n1g = n1 - 1j * c1
    ks = cv.Ks(lam, lam_r, n0, theta_r)
    Lambda = cv.fringe_spacing(lam_r, n0, theta_r)
    s = 1.0 if pol == "pi" else -1.0
    kx2 = kx * kx

    t = Lambda / sublayers
    y_mid = (np.arange(sublayers) + 0.5) * t
    period = _layer_products(kk0, kx2, n0g, n1g, ks, y_mid, t, s)

    m = int(np.floor(d / Lambda))
    rem = d - m * Lambda
    result = (np.ones_like(kk0, dtype=complex), np.zeros_like(kk0, dtype=complex),
              np.zeros_like(kk0, dtype=complex), np.ones_like(kk0, dtype=complex))
    logscale = np.zeros_like(lam)
    base, base_log = _rescale(period)
    mm = m
    while mm:
        if mm & 1:
            result, lr = _rescale(_matmul(base, result))
            logscale = logscale + lr + base_log
        mm >>= 1
        if mm:
            base, lb2 = _rescale(_matmul(base, base))
            base_log = 2.0 * base_log + lb2
    if rem > 0.0:
        n_rem = max(1, int(np.ceil(rem / t)))
        t_rem = rem / n_rem
        y_rem = m * Lambda + (np.arange(n_rem) + 0.5) * t_rem
        Rm = _layer_products(kk0, kx2, n0g, n1g, ks, y_rem, t_rem, s)
        result = _matmul(Rm, result)

    A, B, C, D = result
    kz0 = b * np.cos(theta_c)
    q0 = kz0 / kk0 if s < 0 else kz0 / (kk0 * n0 * n0)
    sq = s * q0
    den = C - sq * A - sq * D + q0 * q0 * B
    r = (sq * A - C - sq * D + q0 * q0 * B) / den
    tamp = (A * (1.0 + r) + B * sq * (1.0 - r)) * np.exp(logscale)
    return dict(r=r, t=tamp, eta_R=np.abs(r) ** 2, eta_T=np.abs(tamp) ** 2)


def sweep(**kw):
    return solve(**kw)["eta_R"]
