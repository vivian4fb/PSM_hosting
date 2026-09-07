"""psmgrating.twowave -- the four two-wave descriptions in closed form.

All four models solve the canonical constant-coefficient system (paper eqs.
32/34 for PSM and PSM*, eq. 99 for TCWT and Kogelnik):

    C_R dR/dy' + a R + i kappa S = 0
    C_S dS/dy' + (a + i theta) S + i kappa R = 0

with reflection boundary conditions R(0) = 1, S(d) = 0 or transmission
conditions R(0) = 1, S(0) = 0.  They differ only in the coefficient set.
The closed solutions are eqs. (48)-(51) and (56)-(57); the removable
singularity at Upsilon -> 0 is handled by the continuation of eq. (61) and the
paragraph after it.  Efficiencies carry the flux factor rho = |C_S / C_R|
(eqs. 50, 57).

Coefficient sets (sigma polarisation):
  PSM       eqs. (33), (35)          kappa_hat exact, including the chi terms
  PSM*      eq. (74): kappa_hat -> zeta kappa_hat, nothing else changes
  TCWT      eqs. (92), (100): a~ = a^ + i kappa^2 / beta
  Kogelnik  eq. (101): a~ = a^ (no kappa^2 term)
pi polarisation:
  PSM       kappa -> kappa cos(2 theta_c)                    eq. (41)
  PSM*      kappa -> zeta P_on kappa                          eqs. (85)-(86)
  TCWT      kappa -> kappa P_pi, loss a^ (not a~)             eqs. (125)-(126), (149)
  Kogelnik  kappa -> kappa cos(2 theta_c), loss a^            Section 6.4.2
"""
import numpy as np

from . import conventions as cv

__all__ = [
    "MODELS", "LABELS", "two_wave_amplitude", "two_wave_efficiency",
    "psm_coefficients", "cwt_coefficients", "on_shell_overlap", "coefficients",
    "efficiency", "sweep",
]

MODELS = ("psm", "psm_star", "kogelnik", "tcwt")
LABELS = {"psm": "PSM", "psm_star": "PSM*", "kogelnik": "Kogelnik", "tcwt": "TCWT"}

_SERIES_LIMIT = 1e-8   # |Upsilon d| below which the series of eq. (61) is used


def _continuation(Upsilon, d):
    """Return (Upsilon, x, E, F) with x = Upsilon d flipped so Re(x) >= 0,
    E = exp(-2x) and F = (1 - E)/Upsilon, the latter by its series
    d (2 - 2x + 4/3 x^2) below |x| = 1e-8 (paper, text after eq. 61).
    Every closed form is even in Upsilon, so the flip is free."""
    Upsilon = np.asarray(Upsilon, dtype=complex)
    x = Upsilon * d
    flip = x.real < 0.0
    Upsilon = np.where(flip, -Upsilon, Upsilon)
    x = np.where(flip, -x, x)
    E = np.exp(-2.0 * x)
    small = np.abs(x) < _SERIES_LIMIT
    safe_U = np.where(small, 1.0, Upsilon)
    F = np.where(small, d * (2.0 - 2.0 * x + (4.0 / 3.0) * x * x), (1.0 - E) / safe_U)
    return Upsilon, x, E, F


def two_wave_amplitude(CR, CS, a, theta, kappa, d, geometry):
    """Complex signal amplitude at the collection face.

    reflection   : S(0) of eq. (48)   [TCWT/Kogelnik: eq. (107)/(130)]
    transmission : S(d) of eq. (56)   [TCWT/Kogelnik: eq. (110)/(132)]
    """
    CR = np.asarray(CR, dtype=complex)
    CS = np.asarray(CS, dtype=complex)
    a = np.asarray(a, dtype=complex)
    theta = np.asarray(theta, dtype=complex)
    kappa = np.asarray(kappa, dtype=complex)
    Q = (CR - CS) * a + 1j * CR * theta                                # eq. (44)
    Ups2 = (Q * Q - 4.0 * CR * CS * kappa * kappa) / (4.0 * CR * CR * CS * CS)  # eq. (45)
    Ups = np.sqrt(Ups2)
    Ups, x, E, F = _continuation(Ups, d)
    if geometry == "reflection":
        # Upsilon coth(d Upsilon) = (1 + E) / F
        return -2j * CR * kappa / (Q - 2.0 * CS * CR * (1.0 + E) / F)
    if geometry == "transmission":
        # sinh(d Upsilon) / Upsilon = exp(x) F / 2; exponents combined before exp
        expo = -(CR + CS) * a * d / (2.0 * CR * CS) - 1j * theta * d / (2.0 * CS) + x
        return -(1j * kappa / CS) * np.exp(expo) * F / 2.0
    raise ValueError("geometry must be 'reflection' or 'transmission'")


def two_wave_efficiency(CR, CS, a, theta, kappa, d, geometry):
    """Diffraction efficiency rho |S|^2 with rho = |C_S / C_R| (eqs. 50, 57)."""
    S = two_wave_amplitude(CR, CS, a, theta, kappa, d, geometry)
    rho = np.abs(np.asarray(CS, dtype=complex) / np.asarray(CR, dtype=complex))
    return np.real(rho * np.abs(S) ** 2)


def on_shell_overlap(lam_c, lam_r, n0, theta_c, theta_r, psi):
    """P_on = k_c . k_S^on / beta^2, eq. (85): the face-parallel component of
    the signal is fixed by Floquet, k_S|| = k_c|| - K||, the modulus by the
    exterior dispersion relation |k_S^on| = beta, and the sign of the normal
    component by the geometry (negative for reflection).  Equals -cos(2
    theta_c) at Bragg resonance (eq. 148).  An evanescent signal (|k_S||| >
    beta) returns a complex value; that regime is outside the paper's scope
    and is flagged by the caller."""
    b = cv.beta(lam_c, n0)
    ks = cv.Ks(lam_c, lam_r, n0, theta_r)
    kcx = b * np.sin(theta_c - psi)
    kcy = b * np.cos(theta_c - psi)
    kSx = kcx + ks * np.sin(psi)           # K'_x' = -K_s sin(Psi), eq. (22)
    sgn = np.where(cv.cSK(lam_c, lam_r, theta_c, theta_r, psi) < 0.0, -1.0, 1.0)
    kSy = sgn * np.sqrt(b * b - kSx * kSx + 0j)
    return (kcx * kSx + kcy * kSy) / (b * b)


def psm_coefficients(lam_c, lam_r, n0, n1, chi0, chi1, theta_c, theta_r, psi,
                     pol="sigma", star=False):
    """PSM (eqs. 33, 35, 41) and PSM* (eqs. 74, 86) coefficient sets.
    Returns (C_R, C_S, a_bar, theta, kappa_hat)."""
    a_ = cv.alpha(lam_c, lam_r)
    b = cv.beta(lam_c, n0)
    cc = np.cos(theta_c)
    cr = np.cos(theta_r)
    CR = cc * np.cos(theta_c - psi) / (a_ * cr)
    CS = -cc * np.cos(theta_c + psi) / (a_ * cr)
    abar = chi0 * b * cc / (a_ * n0 * cr)
    theta = 2.0 * b * (cc / (a_ * cr) - 1.0) * cc * cc
    kap = (b / 2.0) * (n0 * n1 + chi0 * chi1 - 1j * (n0 * chi1 - n1 * chi0)) / (n0 * n0 + chi0 * chi0)
    if star:
        kap = cv.zeta(lam_c, lam_r, theta_c, theta_r) * kap
    if pol == "pi":
        if star:
            kap = kap * on_shell_overlap(lam_c, lam_r, n0, theta_c, theta_r, psi)
        else:
            kap = kap * np.cos(2.0 * theta_c)
    elif pol != "sigma":
        raise ValueError("pol must be 'sigma' or 'pi'")
    return CR, CS, abar, theta, kap


def cwt_coefficients(lam_c, lam_r, n0, n1, chi0, chi1, theta_c, theta_r, psi,
                     pol="sigma", model="tcwt"):
    """TCWT (eqs. 92, 100, 125-126) and Kogelnik (eq. 101, Section 6.4.2)
    coefficient sets.  Returns (c_R, c_S, a, theta, kappa)."""
    if model not in ("tcwt", "kogelnik"):
        raise ValueError("model must be 'tcwt' or 'kogelnik'")
    a_ = cv.alpha(lam_c, lam_r)
    b = cv.beta(lam_c, n0)
    kk0 = cv.k0(lam_c)
    cR = np.cos(theta_c - psi) + 0.0 * a_
    cS = np.cos(theta_c - psi) - 2.0 * a_ * np.cos(theta_r) * np.cos(psi)
    ahat = b * chi0 / n0
    kap = (kk0 / 2.0) * (n1 - 1j * chi1)
    thet = 2.0 * a_ * b * np.cos(theta_r) * (np.cos(theta_c) - a_ * np.cos(theta_r))
    if pol == "sigma":
        loss = ahat + 1j * kap * kap / b if model == "tcwt" else ahat + 0j
    elif pol == "pi":
        loss = ahat + 0j                                   # eq. (149)
        if model == "tcwt":
            kap = kap * (1.0 - 2.0 * a_ * np.cos(theta_r) * np.cos(theta_c))   # P_pi, eq. (126)
        else:
            kap = kap * np.cos(2.0 * theta_c)
    else:
        raise ValueError("pol must be 'sigma' or 'pi'")
    return cR, cS, loss, thet, kap


def coefficients(model, pol, lam_c, lam_r, n0, n1, chi0, chi1, theta_c, theta_r, psi):
    if model == "psm":
        return psm_coefficients(lam_c, lam_r, n0, n1, chi0, chi1, theta_c, theta_r, psi, pol, star=False)
    if model == "psm_star":
        return psm_coefficients(lam_c, lam_r, n0, n1, chi0, chi1, theta_c, theta_r, psi, pol, star=True)
    if model in ("tcwt", "kogelnik"):
        return cwt_coefficients(lam_c, lam_r, n0, n1, chi0, chi1, theta_c, theta_r, psi, pol, model)
    raise ValueError("model must be one of %s" % (MODELS,))


def _resolve_chi(D0, D1, chi0, chi1, lam_c, theta_c, psi, d):
    if chi0 is None:
        chi0 = cv.chi_from_density(D0, lam_c, theta_c, psi, d)
    if chi1 is None:
        chi1 = cv.chi_from_density(D1, lam_c, theta_c, psi, d)
    return chi0, chi1


def efficiency(model, pol, *, n0, n1, d, lam_r, theta_r, lam_c, theta_c, psi,
               D0=0.0, D1=0.0, chi0=None, chi1=None, geometry=None):
    """Diffraction efficiency of one two-wave model.  lam_c may be an array
    (a wavelength sweep); every other argument is a scalar.  The geometry is
    classified at lam_r (the recording wavelength) unless given."""
    lam_c = np.asarray(lam_c, dtype=float)
    chi0, chi1 = _resolve_chi(D0, D1, chi0, chi1, lam_c, theta_c, psi, d)
    if geometry is None:
        geometry = cv.geometry(lam_r, lam_r, theta_c, theta_r, psi)
    CR, CS, a, thet, kap = coefficients(model, pol, lam_c, lam_r, n0, n1, chi0, chi1,
                                        theta_c, theta_r, psi)
    return two_wave_efficiency(CR, CS, a, thet, kap, d, geometry)


def sweep(models, pol, *, lam_c, **params):
    """Wavelength sweep of several models: returns {model: eta array}."""
    return {m: efficiency(m, pol, lam_c=lam_c, **params) for m in models}
