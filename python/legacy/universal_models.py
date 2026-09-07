# VENDORED UNCHANGED from Code2/1_Optics/Research/DBR_volume_gratings/05 Aug files/code/universal_models.py
# Author: Vivian Amos Sureshkumar (reconstruction project, 2026). Kept for parity tests only;
# the package code in psmgrating/ follows the September 2026 paper conventions instead.
"""
universal_models.py  (Session 11, 2026-08-05)
=============================================
PSM and TCWT diffraction efficiency for slanted lossy volume gratings in the
UNIVERSAL slant convention of the published comparative paper
(Brotherton-Ratcliffe, Shi, Osanlou, Excell — PSM3_revised.pdf, Eqs. 2/3/9/12):

  * ONE set of coefficients covers reflection AND transmission:
        cR(PSM) =  cos(thetac) cos(thetac - psi) / (alpha cos(thetar))
        cS(PSM) = -cos(thetac) cos(thetac + psi) / (alpha cos(thetar))
        CT(PSM) =  2 beta (1 - cos(thetac)/(alpha cos(thetar))) cos^2(thetac)
    (and the TCWT/Kogelnik counterparts), with
        psi    = tilt of the grating vector K from the plate normal z
        thetar = recording angle measured from the FRINGE normal (K direction)
        thetac = replay    angle measured from the FRINGE normal
    Unslanted REFLECTION  grating:  psi = 0
    Unslanted TRANSMISSION grating: psi = 90 deg   (David B-R, email 2026-08-05;
    published paper: "For the transmission grating Snell's law limits the
    grating tilt angle to between 60 and 120 [deg]")

  * Only the BOUNDARY CONDITIONS distinguish the two regimes:
        reflection  : R(0) = 1, S(d) = 0   (signal exits the front face)
        transmission: R(0) = 1, S(0) = 0   (signal exits the back face)
    The applicable regime is dictated by the sign of cS/cR:
        cS/cR < 0  -> signal is backward-propagating  -> reflection BCs
        cS/cR > 0  -> signal is forward-propagating   -> transmission BCs

Method
------
Both PSM and TCWT are canonical 2x2 constant-coefficient linear ODE systems:

    cR dR/dz = -abar R                 - i kappa S
    cS dS/dz = -(abar + i CT) S        - i kappa R

The system is solved EXACTLY by the matrix exponential (scipy.linalg.expm);
no new closed forms are transcribed.  For reflection boundary conditions this
reproduces the validated closed-form (coth) results of psm_library.py to
machine precision (< 1e-15, verified in validate_05aug.py Check A); for
transmission boundary conditions it reproduces the published lossless closed
form Eq. (13) of PSM3_revised.pdf (Check B).

Pi polarisation: kappa -> kappa * cos(2 thetac) in BOTH theories (manuscript
Sec. 4.2/4.4 and Kogelnik ref 11).  thetac is measured from the fringe normal,
so 2*thetac is the R-S inter-beam angle in every geometry, reflection or
transmission — the double-angle factor carries over unchanged.

Loss parameters follow the current draft manuscript:
    chi0 = D0 * lamc * cos(thetac) / (2 pi d),   chi1 = chi1_factor * chi0
NOTE (flagged in notation_audit.md): this D0 definition contains cos(thetac)
with thetac measured from the FRINGE normal; for transmission geometries
(psi ~ 90 deg) it no longer represents an optical density along the plate
thickness.  Kept as-is pending David's decision.

Author: reconstruction project, Session 11. Validation: validate_05aug.py.
"""

import numpy as np
from scipy.linalg import expm

deg = np.pi / 180.0


def coefficients(lamc, lamr, n0, n1, Dparam, thetac, thetar, psi, d,
                 chi1_factor=0.0):
    """All PSM and TCWT canonical coefficients at one replay wavelength.

    Angles in radians; psi in the universal convention (reflection ~0,
    transmission ~90 deg).  Returns a dict.
    """
    alpha = lamc / lamr                      # legacy convention (replay/recording)
    beta  = 2.0 * np.pi * n0 / lamc
    chi0  = Dparam * lamc * np.cos(thetac) / (2.0 * np.pi * d)
    chi1  = chi1_factor * chi0

    # --- PSM (harmonic index) ------------------------------------------------
    abar  = (chi0 * beta * np.cos(thetac)) / (alpha * n0 * np.cos(thetar))
    khat  = beta * (n0*n1 + chi0*chi1 + 1j*(n0*chi1 - n1*chi0)) / \
            (2.0 * (n0**2 + chi0**2))
    CR_P  =  np.cos(thetac) * np.cos(thetac - psi) / (alpha * np.cos(thetar))
    CS_P  = -np.cos(thetac) * np.cos(thetac + psi) / (alpha * np.cos(thetar))
    CT_P  =  2.0 * beta * np.cos(thetac)**2 * \
             (1.0 - np.cos(thetac) / (alpha * np.cos(thetar)))

    # --- TCWT (harmonic index, Kogelnik-type) --------------------------------
    kappa = np.pi * (n1 + 1j * chi1) / lamc
    ahat  = -(beta * chi0 / n0) * (1.0 + 1j * chi0 / (2.0 * n0))
    atil  = ahat + 1j * kappa**2 / beta
    CR_K  = -np.cos(thetac - psi)
    CS_K  = -(np.cos(thetac - psi) - 2.0 * alpha * np.cos(thetar) * np.cos(psi))
    CT_K  =  2.0 * alpha * beta * np.cos(thetar) * \
             (np.cos(thetac) - alpha * np.cos(thetar))

    return dict(alpha=alpha, beta=beta, chi0=chi0, chi1=chi1,
                abar=abar, khat=khat, CR_P=CR_P, CS_P=CS_P, CT_P=CT_P,
                kappa=kappa, atil=atil, CR_K=CR_K, CS_K=CS_K, CT_K=CT_K)


def _eta_canonical(CR, CS, CT, aterm, kap, d, geometry):
    """Exact solution of the canonical 2x2 system for given BCs.

    reflection  : R(0)=1, S(d)=0 -> eta = -(CS/CR) |S(0)|^2
    transmission: R(0)=1, S(0)=0 -> eta = +(CS/CR) |S(d)|^2
    """
    A = np.array([[-aterm / CR,            -1j * kap / CR],
                  [-1j * kap / CS, -(aterm + 1j * CT) / CS]], dtype=complex)
    P = expm(A * d)
    if geometry == 'reflection':
        S0 = -P[1, 0] / P[1, 1]
        return float(np.real(-(CS / CR) * S0 * np.conj(S0)))
    elif geometry == 'transmission':
        Sd = P[1, 0]
        return float(np.real((CS / CR) * Sd * np.conj(Sd)))
    raise ValueError("geometry must be 'reflection' or 'transmission'")


def eta_two_wave(lamc, lamr, n0, n1, Dparam, thetac, thetar, psi, d,
                 theory='PSM', pol='sigma', geometry=None, chi1_factor=0.0):
    """PSM or TCWT diffraction efficiency, sigma or pi, either geometry.

    geometry=None auto-selects from the sign of cS/cR (the physical regime);
    passing it explicitly asserts the expectation (raises on mismatch).
    """
    c = coefficients(lamc, lamr, n0, n1, Dparam, thetac, thetar, psi, d,
                     chi1_factor)
    if theory.upper() == 'PSM':
        CR, CS, CT, aterm, kap = c['CR_P'], c['CS_P'], c['CT_P'], c['abar'], c['khat']
    elif theory.upper() == 'TCWT':
        CR, CS, CT, aterm, kap = c['CR_K'], c['CS_K'], c['CT_K'], c['atil'], c['kappa']
    else:
        raise ValueError("theory must be 'PSM' or 'TCWT'")

    if pol.lower() in ('pi', 'tm', 'p'):
        kap = kap * np.cos(2.0 * thetac)
    elif pol.lower() not in ('sigma', 'te', 's'):
        raise ValueError("pol must be 'sigma' or 'pi'")

    regime = 'reflection' if np.real(CS / CR) < 0 else 'transmission'
    if geometry is None:
        geometry = regime
    elif geometry != regime:
        raise ValueError(
            f"requested geometry '{geometry}' but cS/cR = {np.real(CS/CR):+.4f} "
            f"implies '{regime}' (psi = {psi/deg:.1f} deg)")
    return _eta_canonical(CR, CS, CT, aterm, kap, d, geometry)


def sweep_two_wave(lam_vec, lamr, n0, n1, Dparam, thetac, thetar, psi, d,
                   theory='PSM', pol='sigma', geometry=None, chi1_factor=0.0):
    """Wavelength sweep of eta_two_wave."""
    return np.array([eta_two_wave(l, lamr, n0, n1, Dparam, thetac, thetar,
                                  psi, d, theory, pol, geometry, chi1_factor)
                     for l in lam_vec])


# ----------------------------------------------------------------------------
# Published lossless closed forms (PSM3_revised.pdf Eqs. 10, 11, 13) — used
# only as independent cross-checks in validate_05aug.py, never for figures.
# ----------------------------------------------------------------------------
def eta_published_lossless(lamc, lamr, n0, n1, thetac, thetar, psi, d,
                           theory='PSM', geometry='transmission'):
    """Lossless closed forms: Eq.(10) reflection, Eq.(13) transmission,
    with Gamma^2 from Eq.(11); PSM or Kogelnik coefficients."""
    c = coefficients(lamc, lamr, n0, n1, 0.0, thetac, thetar, psi, d, 0.0)
    if theory.upper() == 'PSM':
        cR, cS, dth, kap = c['CR_P'], c['CS_P'], c['CT_P'], np.real(c['khat'])
    else:
        cR, cS, dth, kap = c['CR_K'], c['CS_K'], c['CT_K'], np.real(c['kappa'])
    Gam2 = -dth**2 / (4.0 * cS**2) - kap**2 / (cR * cS) + 0j
    Gam = np.sqrt(Gam2)
    if geometry == 'transmission':
        return float(np.real(kap**2 / (2.0 * Gam2 * cS * cR) *
                             (np.cosh(2.0 * d * Gam) - 1.0)))
    num = kap**2 * np.sinh(d * Gam)**2
    return float(np.real(num / (num - cR * cS * Gam2)))
