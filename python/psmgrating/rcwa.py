"""psmgrating.rcwa -- rigorous coupled-wave reference (paper Appendix A).

A direct numerical solution of the boundary-value problem the four two-wave
models approximate: Maxwell's equations in a slab 0 <= y' <= d of harmonic
complex index n~ = (n0 - i chi0) + (n1 - i chi1) cos(K.r), eq. (175),
embedded in an index-matched lossless exterior of index n0 on both sides
(Appendix A.1), expanded in Floquet orders |l| <= L (eq. 177), solved by
diagonalisation (eq. 180) with eigenmode anchoring (eqs. 189-190), and read
out through eqs. (185)-(186).

sigma (TE): (h, v) = (E_z, Z0 H_x); pi (TM): (h, v) = (H_z, E_x / Z0), eq.
(176).  The implementation follows Vivian's validated core
(validation/rcw_tm/rcw_rigorous.py, run_007/run_008/run_010) in the (h, -v)
gauge of eq. (179); efficiencies are gauge-independent.  Material Toeplitz
matrices are built from the closed-form harmonics of eps (eq. 181) and of
1/eps (validation/symbolic/gm_closed_form.py, exact to 1e-13 against the
4096-point FFT), so no FFT is needed and the code runs under Pyodide.

Only numpy is required.
"""
import numpy as np

from . import conventions as cv

__all__ = ["solve", "sweep", "eps_harmonics", "g_harmonic", "toeplitz"]


def eps_harmonics(a, b):
    """Fourier coefficients of eps = (a + b cos phi)^2, eq. (181)."""
    return {0: a * a + b * b / 2.0, 1: a * b, -1: a * b, 2: b * b / 4.0, -2: b * b / 4.0}


def _branch(a, b):
    w = np.sqrt(a * a - b * b + 0j)
    rho = (w - a) / b
    if abs(rho) >= 1.0:
        w = -w
        rho = (w - a) / b
    return w, rho


def g_harmonic(a, b, m):
    """Fourier coefficient g_m of 1/(a + b cos phi)^2 = rho^|m| (a + |m| w)/w^3
    with w = sqrt(a^2 - b^2) on the branch |rho| < 1, rho = (w - a)/b."""
    if b == 0:
        return 1.0 / (a * a) if m == 0 else 0.0
    w, rho = _branch(a, b)
    return rho ** abs(m) * (a + abs(m) * w) / w ** 3


def toeplitz(coeff, N):
    """T[p, l] = c_{p-l} from a callable coeff(m), m = -(N-1)..(N-1)."""
    ms = np.arange(-(N - 1), N)
    cm = np.array([coeff(int(m)) for m in ms], dtype=complex)
    idx = np.arange(N)
    return cm[(idx[:, None] - idx[None, :]) + (N - 1)]


def solve(*, lam_c, lam_r, n0, n1, chi0, chi1, theta_c, theta_r, psi, d,
          pol="sigma", L=10):
    """Rigorous solution at one replay wavelength.

    Returns a dict with the order indices `ls`, complex amplitudes `R`, `T`,
    per-order efficiencies `DE_R`, `DE_T` (eq. 186), `sum_DE` (eq. 187),
    the signal efficiency `eta` (order l = -1 at the face fixed by the sign
    of c_S^K, eq. 46), `geometry`, and the coherent sums at K_x = 0 (which
    are NOT the observable there: see chain.py and Appendix A.6).
    """
    lam_c = float(lam_c)
    kk0 = 2.0 * np.pi / lam_c
    b = kk0 * n0
    a_ = lam_c / lam_r
    ks = 2.0 * a_ * b * np.cos(theta_r)                    # signed, eq. (4)

    n0g = n0 - 1j * chi0
    n1g = n1 - 1j * chi1

    kx = b * np.sin(theta_c - psi)                          # eq. (88)
    Kx = -ks * np.sin(psi)                                  # eq. (89)
    Kz = ks * np.cos(psi)

    N = 2 * L + 1
    ls = np.arange(-L, L + 1)
    kxl = kx + ls * Kx
    kzl = np.sqrt(b * b - kxl * kxl + 0j)                   # eq. (183)
    kzl = np.where(kzl.imag > 0, -kzl, kzl)                 # Im(k_y,l) <= 0
    izero = L
    kz0 = kzl[izero].real
    eps_c = n0 * n0

    eh = eps_harmonics(n0g, n1g)
    E = toeplitz(lambda m: eh.get(m, 0.0), N)
    G = toeplitz(lambda m: g_harmonic(n0g, n1g, m), N)

    Lz = np.diag(1j * ls * Kz)
    X = np.diag(kxl)
    I = np.eye(N)

    if pol == "pi":
        Mmat = np.block([[Lz, -1j * kk0 * E],
                         [(1j / kk0) * (X @ G @ X) - 1j * kk0 * I, Lz]])
        q = kzl / (kk0 * eps_c)
        sgn0, sgnd = +1.0, -1.0
    elif pol == "sigma":
        Mmat = np.block([[Lz, 1j * kk0 * I],
                         [1j * kk0 * E - (1j / kk0) * (X @ X), Lz]])
        q = kzl / kk0
        sgn0, sgnd = -1.0, +1.0
    else:
        raise ValueError("pol must be 'sigma' or 'pi'")

    lam, W = np.linalg.eig(Mmat)
    Wh, Wv = W[:N, :], W[N:, :]

    up = lam.real > 0                                       # anchoring, eq. (189)
    s0 = np.ones(2 * N, dtype=complex)
    sd = np.ones(2 * N, dtype=complex)
    s0[up] = np.exp(-lam[up] * d)
    sd[~up] = np.exp(lam[~up] * d)

    Q = np.diag(q)
    A = np.vstack([(Wv + sgn0 * Q @ Wh) * s0[None, :],
                   (Wv + sgnd * Q @ Wh) * sd[None, :]])    # eq. (190)
    rhs = np.zeros(2 * N, dtype=complex)
    rhs[izero] = sgn0 * 2.0 * q[izero]
    c = np.linalg.solve(A, rhs)

    h0 = Wh @ (c * s0)
    hd = Wh @ (c * sd)
    delta = np.zeros(N, dtype=complex)
    delta[izero] = 1.0
    Ramp = h0 - delta
    Tamp = hd * np.exp(-1j * ls * Kz * d)

    DE_R = (kzl.real / kz0) * np.abs(Ramp) ** 2             # eq. (186)
    DE_T = (kzl.real / kz0) * np.abs(Tamp) ** 2

    geom = cv.geometry(lam_c, lam_r, theta_c, theta_r, psi)
    eta = float(DE_R[izero - 1]) if geom == "reflection" else float(DE_T[izero - 1])

    degenerate = abs(Kx) < 1e-9 * b
    out = dict(ls=ls, R=Ramp, T=Tamp, DE_R=DE_R, DE_T=DE_T,
               sum_DE=float(DE_R.sum() + DE_T.sum()), eta=eta, geometry=geom,
               kzl=kzl, kxl=kxl, degenerate=bool(degenerate))
    if degenerate:
        out["DE_R_coherent"] = float(np.abs(Ramp.sum()) ** 2)
        out["DE_T_coherent"] = float(np.abs(Tamp.sum()) ** 2)
    return out


def sweep(*, lam_c, lam_r, n0, n1, d, theta_c, theta_r, psi, pol="sigma", L=10,
          D0=0.0, D1=0.0, chi0=None, chi1=None):
    """Signal efficiency over an array of replay wavelengths.  chi0/chi1 follow
    the density parametrisation of eq. (175) unless given explicitly."""
    lam_c = np.atleast_1d(np.asarray(lam_c, dtype=float))
    out = np.empty(lam_c.shape)
    for i, lc in enumerate(lam_c):
        c0 = cv.chi_from_density(D0, lc, theta_c, psi, d) if chi0 is None else chi0
        c1 = cv.chi_from_density(D1, lc, theta_c, psi, d) if chi1 is None else chi1
        out[i] = solve(lam_c=lc, lam_r=lam_r, n0=n0, n1=n1, chi0=c0, chi1=c1,
                       theta_c=theta_c, theta_r=theta_r, psi=psi, d=d, pol=pol, L=L)["eta"]
    return out
