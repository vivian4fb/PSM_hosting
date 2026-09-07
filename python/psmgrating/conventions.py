"""psmgrating.conventions -- single source of truth for the paper's conventions.

Paper: D. Brotherton-Ratcliffe and V. Sureshkumar, "Slanted planar volume
gratings: parallel-stacked-mirror models versus traditional coupled wave
descriptions" (2026).  Equation numbers below refer to that paper.

Units: metres and radians throughout.  Phasors exp(+i w t), exp(-i k.r); a
passive medium is n~ = n - i chi with chi > 0 (Section 2, footnote 1).

Geometry (paper Fig. 2): y is the DIRECTED fringe normal e_K and x lies along
the fringe planes; the slab frame (x', y') has y' along the boundary normal,
and the two are related by the slant Psi (eq. 20).  theta_c (replay) and
theta_r (recording) are measured from e_K; the reference wave makes the angle
theta_c - Psi with the boundary normal (eq. 21).  The scalar grating vector
K_s = 2 alpha beta cos(theta_r) (eq. 4) is SIGNED: it is negative wherever
theta_r is obtuse, which is the case over most of the transmission family.

Convention guards (each is a unit test in tests/test_conventions.py):
  * alpha = lambda_c / lambda_r                                   (eq. 5)
  * chi_0 = D_0 lambda_c cos(theta_c - Psi) / (2 pi d)            (eq. 143, 175)
    -- the direction cosine of the REFERENCE wave, not cos(theta_c) as in the
    August 2026 code, which gives gain for the obtuse theta_c of the
    transmission family.
  * geometry is decided by the sign of c_S^K (eq. 46), not of C_S/C_R.
"""
import numpy as np

__all__ = [
    "alpha", "beta", "k0", "Ks", "zeta", "cSK", "geometry", "chi_from_density",
    "density_from_chi", "fringe_spacing", "klein_cook", "internal_angle",
    "family_angles", "deg",
]

deg = np.pi / 180.0


def alpha(lam_c, lam_r):
    """Replay-to-recording wavelength ratio, eq. (5)."""
    return np.asarray(lam_c, dtype=float) / lam_r


def beta(lam_c, n0):
    """Propagation constant in the mean medium, eq. (3)."""
    return 2.0 * np.pi * n0 / np.asarray(lam_c, dtype=float)


def k0(lam_c):
    """Free-space wavenumber."""
    return 2.0 * np.pi / np.asarray(lam_c, dtype=float)


def Ks(lam_c, lam_r, n0, theta_r):
    """Signed scalar grating vector K_s = 2 alpha beta cos(theta_r), eq. (4)."""
    return 2.0 * alpha(lam_c, lam_r) * beta(lam_c, n0) * np.cos(theta_r)


def zeta(lam_c, lam_r, theta_c, theta_r):
    """The PSM* factor zeta = cos(theta_c) / (alpha cos(theta_r)), eq. (74).
    Unity at Bragg resonance and nowhere else."""
    return np.cos(theta_c) / (alpha(lam_c, lam_r) * np.cos(theta_r))


def cSK(lam_c, lam_r, theta_c, theta_r, psi):
    """y'-component of (k'_R - K') / beta, eq. (46).  Negative: reflection
    geometry (signal collected at the entry face).  Positive: transmission."""
    a = alpha(lam_c, lam_r)
    return np.cos(theta_c - psi) - 2.0 * a * np.cos(theta_r) * np.cos(psi)


def geometry(lam_c, lam_r, theta_c, theta_r, psi):
    """'reflection' or 'transmission' by the sign of c_S^K at the given replay
    wavelength (scalar).  For a sweep, classify at the Bragg wavelength and
    check sign agreement separately (the paper checked 3.3 million points)."""
    c = float(np.asarray(cSK(lam_c, lam_r, theta_c, theta_r, psi)).ravel()[0])
    return "reflection" if c < 0.0 else "transmission"


def chi_from_density(D, lam_c, theta_c, psi, d):
    """Imaginary index from Kogelnik's density parameter, eqs. (143)-(144),
    (175): chi = D lambda_c cos(theta_c - Psi) / (2 pi d).  The direction
    cosine is that of the reference wave, so D is the optical density
    accumulated along the reference ray."""
    return D * np.asarray(lam_c, dtype=float) * np.cos(theta_c - psi) / (2.0 * np.pi * d)


def density_from_chi(chi, lam_c, theta_c, psi, d):
    """Inverse of chi_from_density."""
    return 2.0 * np.pi * d * chi / (np.asarray(lam_c, dtype=float) * np.cos(theta_c - psi))


def fringe_spacing(lam_r, n0, theta_r):
    """Lambda = 2 pi / |K_s|, a positive length (text after eq. (5); eq. (170)
    at Bragg).  Fixed by the recording geometry."""
    return lam_r / (2.0 * n0 * abs(np.cos(theta_r)))


def klein_cook(lam_c, d, n0, Lambda, theta_int):
    """Klein-Cook parameter Q = 2 pi lambda_c d / (n0 Lambda^2 cos theta), eq.
    (171); theta is the internal replay angle from the boundary normal,
    theta_c - Psi.  The paper shows Q is NOT the right admissibility
    diagnostic; it is reported for context only."""
    return 2.0 * np.pi * np.asarray(lam_c, dtype=float) * d / (n0 * Lambda ** 2 * np.cos(theta_int))


def internal_angle(theta_air, n0):
    """Snell refraction of the external replay angle into the mean medium."""
    return np.arcsin(np.sin(theta_air) / n0)


def family_angles(theta_air, n0, psi):
    """The paper's Section 8 construction: one fixed reference beam at
    theta_air outside (theta_in inside), the slant swept; the angle the beams
    make with the directed fringe normal is theta_c = Psi + theta_in, and the
    grating is replayed by its own recording beam, theta_r = theta_c.
    Returns (theta_c, theta_r)."""
    th_in = internal_angle(theta_air, n0)
    theta_c = psi + th_in
    return theta_c, theta_c
