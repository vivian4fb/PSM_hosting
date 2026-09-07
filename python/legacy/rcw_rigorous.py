# VENDORED UNCHANGED from Code2/1_Optics/Research/DBR_volume_gratings/validation/rcw_tm/rcw_rigorous.py
# Author: Vivian Sureshkumar (reconstruction project, 2026). Kept for parity tests only;
# the package code in psmgrating/ follows the September 2026 paper conventions instead.
"""
rcw_rigorous.py
===============
Rigorous coupled-wave (RCW) solver for planar slanted volume gratings with
finite absorption, for BOTH polarisations:

  - TE (sigma):  field component Ey
  - TM (pi)   :  vector (Hy, Ex) first-order formulation  <-- the "Maxwell-exact"
                 TM treatment requested by D. Brotherton-Ratcliffe (email 2026-07-28):
                 naive scalar Hz/Hy second-order formulations of RCWT-pi are not
                 Maxwell-exact and converge slowly; the (hz, ex) vector form with a
                 full Toeplitz treatment of eps and 1/eps converges with few modes.

Physical model (identical geometry/conventions to the legacy DBR MATLAB codes
Lossy_RCW_Lam_with_slant*.m and psm_library.py):

  - grating occupies 0 <= z <= d, surrounded by index-matched LOSSLESS media n0
    (no Fresnel boundary reflections from the mean index step)
  - complex index inside the grating:
        n(r) = (n0 - i*chi0) + (n1 - i*chi1) * cos(K.r)
    with time convention exp(+i w t), spatial convention exp(-i k.r)
    (so Im(n) < 0 gives decay; chi0, chi1 >= 0)
  - chi0 = Dparam * lamc * cos(thetac) / (2*pi*d)   (Kogelnik density parameter)
  - permittivity eps(r) = n(r)^2 kept to ALL harmonics (m = 0, +-1, +-2) —
    the legacy codes drop the second harmonic of eps; here nothing is dropped.
  - 1/eps(r) expanded by FFT to machine precision (exact Toeplitz treatment)
  - incident plane wave k = beta*(sin(thetac-psi), cos(thetac-psi)),
    beta = 2*pi*n0/lamc
  - grating vector K = 2*beta_r*cos(thetar)*(-sin(psi), cos(psi)),
    beta_r = 2*pi*n0/lamr  (fixed by the recording geometry)
  - Bragg-reflected order is l = -1 (same labelling as the legacy codes)

Numerical method:
  - Floquet expansion F(x,z) = sum_l f_l(z) * exp(-i[(kx+l*Kx)x + l*Kz*z])
    (the slant phase l*Kz*z is folded into the basis so the ODE system has
    constant coefficients)
  - first-order system  psi' = M psi,  psi = [h; v] (TM: h=Hy, v=Ex/Z0;
    TE: h=Ey, v=Z0*Hx), solved by eigen-decomposition
  - boundary matching with DIRECTIONAL exponential anchoring (growing modes
    anchored at z=d, decaying at z=0) so evanescent orders never overflow —
    this permits large mode counts, unlike the legacy 7-mode codes.

Author: reconstruction project, 2026-07-28. Validation: run_tm_validation.py.
"""

import numpy as np


def _fourier_toeplitz(vals_on_grid, N):
    """Return the (2N-1)-long array of Fourier coefficients c_m, m=-(N-1)..(N-1),
    of a periodic function sampled on a uniform grid over [0, 2pi), using the
    convention  f(phi) = sum_m c_m exp(-i m phi).
    """
    Nfft = len(vals_on_grid)
    # c_m = (1/2pi) int f exp(+i m phi) dphi  ->  ifft convention
    coeffs = np.fft.ifft(vals_on_grid)
    ms = np.arange(-(N - 1), N)
    return coeffs[ms % Nfft]


def _toeplitz_from_coeffs(cm, N):
    """Build T[p,l] = c_{p-l} from coefficient array cm indexed m=-(N-1)..(N-1)."""
    idx = np.arange(N)
    return cm[(idx[:, None] - idx[None, :]) + (N - 1)]


def rcw_planar(lamc, lamr, n0, n1, Dparam, thetac, thetar, psi, d,
               pol='TM', chi1_factor=0.0, M=8, Nfft=4096):
    """
    Rigorous coupled-wave solution for a slanted lossy planar volume grating.

    Parameters
    ----------
    lamc, lamr : replay / recording wavelengths [m]
    n0, n1     : mean index, real index modulation
    Dparam     : Kogelnik photographic density D0 (chi0 = D*lamc*cos(thetac)/(2 pi d))
    thetac     : replay angle [rad]   (legacy convention, beam dir = thetac - psi)
    thetar     : recording half-angle [rad]
    psi        : grating slant angle [rad]
    d          : thickness [m]
    pol        : 'TM' (pi) or 'TE' (sigma)
    chi1_factor: chi1 = chi1_factor * chi0 (mixed/absorption grating)
    M          : Floquet orders l = -M..M  (2M+1 harmonics)
    Nfft       : sampling for the exact 1/eps Fourier expansion

    Returns
    -------
    dict with:
      ls      : order indices l = -M..M
      DE_R    : reflected diffraction efficiency per order
      DE_T    : transmitted diffraction efficiency per order
      eta_m1  : DE_R at l = -1 (Bragg-reflected order)
      sum_DE  : total (energy check; = 1 for lossless)
      R, T    : complex order amplitudes
    """
    lamc = float(lamc)
    k0 = 2.0 * np.pi / lamc
    beta = k0 * n0
    beta_r = 2.0 * np.pi * n0 / lamr

    chi0 = Dparam * lamc * np.cos(thetac) / (2.0 * np.pi * d)
    chi1 = chi1_factor * chi0

    # complex index (exp(+iwt), exp(-ik.r): absorbing => negative imaginary part)
    n0g = n0 - 1j * chi0
    n1g = n1 - 1j * chi1

    kx = beta * np.sin(thetac - psi)
    Kx = -2.0 * beta_r * np.cos(thetar) * np.sin(psi)
    Kz = 2.0 * beta_r * np.cos(thetar) * np.cos(psi)

    N = 2 * M + 1
    ls = np.arange(-M, M + 1)
    kxl = kx + ls * Kx

    # external z-wavenumbers in the index-matched lossless surround (index n0)
    kzl = np.sqrt(beta**2 - kxl**2 + 0j)
    kzl = np.where(kzl.imag > 0, -kzl, kzl)      # branch: Im(kzl) <= 0
    izero = M                                    # index of l = 0
    kz0 = kzl[izero].real
    eps_c = n0**2

    # ---- exact Fourier expansion of eps and 1/eps over one grating period ----
    phi = 2.0 * np.pi * np.arange(Nfft) / Nfft
    eps_prof = (n0g + n1g * np.cos(phi))**2
    cm = _fourier_toeplitz(eps_prof, N)          # eps harmonics (all kept)
    gm = _fourier_toeplitz(1.0 / eps_prof, N)    # exact 1/eps harmonics
    E = _toeplitz_from_coeffs(cm, N)
    G = _toeplitz_from_coeffs(gm, N)

    Lz = np.diag(1j * ls * Kz)
    X = np.diag(kxl)
    I = np.eye(N)

    if pol.upper() == 'TM':
        # h = Hy harmonics, v = Ex/Z0 harmonics
        #   h' = i l Kz h - i k0 (E v)
        #   v' = i l Kz v + (i/k0) X G X h - i k0 h
        Mmat = np.block([
            [Lz, -1j * k0 * E],
            [(1j / k0) * (X @ G @ X) - 1j * k0 * I, Lz]])
        q = kzl / (k0 * eps_c)
        sgn0, sgnd = +1.0, -1.0    # v(0) + Q h(0) = 2 q0 delta ; v(d) - Q h(d) = 0
    elif pol.upper() == 'TE':
        # h = Ey harmonics, v = Z0*Hx harmonics
        #   h' = i l Kz h + i k0 v
        #   v' = i l Kz v + i k0 (E h) - (i/k0) X^2 h
        Mmat = np.block([
            [Lz, 1j * k0 * I],
            [1j * k0 * E - (1j / k0) * (X @ X), Lz]])
        q = kzl / k0
        sgn0, sgnd = -1.0, +1.0    # v(0) - Q h(0) = -2 q0 delta ; v(d) + Q h(d) = 0
    else:
        raise ValueError("pol must be 'TM' or 'TE'")

    lam, W = np.linalg.eig(Mmat)
    Wh, Wv = W[:N, :], W[N:, :]

    # directional anchoring: growing modes (Re lam > 0) anchored at z=d
    up = lam.real > 0
    s0 = np.ones(2 * N, dtype=complex)               # column factors at z=0
    sd = np.ones(2 * N, dtype=complex)               # column factors at z=d
    s0[up] = np.exp(-lam[up] * d)
    sd[~up] = np.exp(lam[~up] * d)
    # (|s0|, |sd| <= 1 always: no overflow at any mode count)

    Q = np.diag(q)
    A = np.vstack([
        (Wv + sgn0 * Q @ Wh) * s0[None, :],
        (Wv + sgnd * Q @ Wh) * sd[None, :]])
    rhs = np.zeros(2 * N, dtype=complex)
    rhs[izero] = sgn0 * 2.0 * q[izero]

    c = np.linalg.solve(A, rhs)

    h0 = Wh @ (c * s0)
    hd = Wh @ (c * sd)

    delta = np.zeros(N, dtype=complex)
    delta[izero] = 1.0
    R = h0 - delta
    P = np.exp(-1j * ls * Kz * d)
    T = hd * P

    # diffraction efficiencies (index-matched real surround: eps_c cancels)
    DE_R = (kzl.real / kz0) * np.abs(R)**2
    DE_T = (kzl.real / kz0) * np.abs(T)**2

    # degenerate (unslanted, Kx = 0) case: all orders collinear -> coherent sums
    if abs(Kx) < 1e-9 * beta:
        R_coh = R.sum()
        T_coh = T.sum()                             # all kzl equal at Kx=0
        DE_R_coh = float(np.abs(R_coh)**2)
        DE_T_coh = float(np.abs(T_coh)**2)
    else:
        DE_R_coh = None
        DE_T_coh = None

    return dict(ls=ls, R=R, T=T, DE_R=DE_R, DE_T=DE_T,
                eta_m1=float(DE_R[izero - 1]),
                sum_DE=float(DE_R.sum() + DE_T.sum()),
                DE_R_coherent=DE_R_coh, DE_T_coherent=DE_T_coh,
                kzl=kzl, kxl=kxl)


def sweep_rcw(lam_c_vec, lamr, n0, n1, Dparam, thetac, thetar, psi, d,
              pol='TM', chi1_factor=0.0, M=8):
    """Wavelength sweep; returns eta(l=-1) array."""
    out = np.zeros(len(lam_c_vec))
    for ii, lamc in enumerate(lam_c_vec):
        out[ii] = rcw_planar(lamc, lamr, n0, n1, Dparam, thetac, thetar, psi, d,
                             pol=pol, chi1_factor=chi1_factor, M=M)['eta_m1']
    return out
