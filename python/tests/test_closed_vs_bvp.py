"""The closed forms against an independent matrix-exponential solution of the
same 2x2 boundary-value problem (scipy Pade expm), 200 random configurations
including lossy transmission gratings -- the case in which the legacy library
carried twice the correct attenuation exponent."""
import numpy as np
import pytest

from psmgrating import conventions as cv
from psmgrating import twowave as tw

scipy_linalg = pytest.importorskip("scipy.linalg")


def _bvp(CR, CS, a, theta, kap, d, geometry):
    A = np.array([[-a / CR, -1j * kap / CR],
                  [-1j * kap / CS, -(a + 1j * theta) / CS]], dtype=complex)
    P = scipy_linalg.expm(A * d)
    if geometry == "reflection":
        S0 = -P[1, 0] / P[1, 1]
        return abs(CS / CR) * abs(S0) ** 2
    Sd = P[1, 0]
    return abs(CS / CR) * abs(Sd) ** 2


def test_random_configurations():
    rng = np.random.default_rng(20260907)
    worst = 0.0
    n_ref = n_tr = 0
    for _ in range(200):
        n0 = rng.uniform(1.45, 1.7)
        n1 = rng.uniform(0.0, 0.08)
        d = rng.uniform(4e-6, 40e-6)
        lam_r = rng.uniform(450e-9, 650e-9)
        lam_c = lam_r * rng.uniform(0.85, 1.15)
        reflection = rng.random() < 0.5
        psi = (rng.uniform(-35, 10) if reflection else rng.uniform(60, 120)) * cv.deg
        theta_c = psi + rng.uniform(5, 40) * cv.deg
        theta_r = theta_c + rng.uniform(-3, 3) * cv.deg
        D0 = rng.uniform(0, 1.6)
        D1 = rng.uniform(0, D0) if rng.random() < 0.7 else 0.0
        chi0 = cv.chi_from_density(D0, lam_c, theta_c, psi, d)
        chi1 = cv.chi_from_density(D1, lam_c, theta_c, psi, d)
        geom = cv.geometry(lam_c, lam_r, theta_c, theta_r, psi)
        if geom == "reflection":
            n_ref += 1
        else:
            n_tr += 1
        for model in tw.MODELS:
            for pol in ("sigma", "pi"):
                CR, CS, a, th, kap = tw.coefficients(model, pol, lam_c, lam_r, n0, n1,
                                                     chi0, chi1, theta_c, theta_r, psi)
                closed = float(tw.two_wave_efficiency(CR, CS, a, th, kap, d, geom))
                ref = _bvp(complex(CR), complex(CS), complex(a), complex(th), complex(kap), d, geom)
                err = abs(closed - ref) / max(ref, 1e-12)
                worst = max(worst, err)
                assert err < 1e-9, (model, pol, geom, closed, ref)
    assert n_ref > 50 and n_tr > 50
    assert worst < 1e-9


def test_removable_singularity_is_smooth():
    # walk kappa through the locus Q^2 = 4 CR CS kappa^2 (Upsilon = 0) and
    # demand continuity of the closed form
    CR, CS, a, d = 0.9, 0.7, 0.0, 20e-6
    th = 2.0e4
    Q = 1j * CR * th
    k_star = np.sqrt(Q * Q / (4 * CR * CS))          # Upsilon = 0 exactly
    ks = k_star * (1 + np.linspace(-1e-6, 1e-6, 11))
    vals = np.array([float(tw.two_wave_efficiency(CR, CS, a, th, k, d, "transmission")) for k in ks])
    assert np.all(np.isfinite(vals))
    for k, v in zip(ks, vals):
        ref = _bvp(CR + 0j, CS + 0j, a + 0j, th + 0j, complex(k), d, "transmission")
        assert abs(v - ref) < 1e-10, (k, v, ref)
