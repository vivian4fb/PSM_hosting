"""Closed-form limits at exact Bragg resonance (paper Section 7)."""
import numpy as np

from psmgrating import conventions as cv
from psmgrating import twowave as tw

deg = cv.deg
LN3 = float(np.log(3.0))


def _bragg(psi_deg, theta_c_deg, **over):
    p = dict(n0=1.5, n1=0.045, d=20e-6, lam_r=532e-9, lam_c=532e-9,
             theta_c=theta_c_deg * deg, theta_r=theta_c_deg * deg, psi=psi_deg * deg)
    p.update(over)
    return p


def test_lossless_reflection_tanh2_eq141():
    p = _bragg(0.0, 10.0)
    b = cv.beta(p["lam_c"], p["n0"])
    kap = b * p["n1"] / (2 * p["n0"])
    expected = np.tanh(kap * p["d"] / np.cos(p["theta_c"])) ** 2
    for m in ("psm", "psm_star", "kogelnik"):
        eta = tw.efficiency(m, "sigma", **p)
        assert abs(eta - expected) < 1e-12, (m, eta, expected)
    # TCWT carries the second-order i kappa^2/beta term: close but not identical
    assert abs(tw.efficiency("tcwt", "sigma", **p) - expected) < 2e-4


def test_pi_reflection_tanh2_eq151_and_polarisation_null():
    p = _bragg(0.0, 10.0)
    b = cv.beta(p["lam_c"], p["n0"])
    kap = b * p["n1"] / (2 * p["n0"]) * abs(np.cos(2 * p["theta_c"]))
    expected = np.tanh(kap * p["d"] / np.cos(p["theta_c"])) ** 2
    for m in ("psm", "psm_star", "kogelnik", "tcwt"):
        eta = tw.efficiency(m, "pi", **p)
        tol = 1e-12 if m != "tcwt" else 2e-4
        assert abs(eta - expected) < tol, (m, eta, expected)
    # theta_c = 45 deg: the pi coupling vanishes identically (Section 7.2)
    p45 = _bragg(0.0, 45.0)
    for m in ("psm", "psm_star", "kogelnik", "tcwt"):
        assert tw.efficiency(m, "pi", **p45) < 1e-25


def test_psm_star_equals_psm_at_bragg_any_slant():
    for psi in (-25.0, -5.0, 7.5, 70.0, 85.0):
        for pol in ("sigma", "pi"):
            p = _bragg(psi, 33.86 + psi if psi > 30 else 20.0, n1=0.03, D0=0.5, D1=0.5)
            a = tw.efficiency("psm", pol, **p)
            s = tw.efficiency("psm_star", pol, **p)
            assert abs(a - s) < 1e-13, (psi, pol, a, s)


def test_absorption_reflection_limit_eq145():
    # 1/(7 + 4 sqrt 3) = 0.071797 as D0 = D1 -> infinity, unslanted sigma
    limit = 1.0 / (7.0 + 4.0 * np.sqrt(3.0))
    assert abs(limit - 0.071797) < 1e-6
    p = _bragg(0.0, 10.0, n1=0.0, D0=6.0, D1=6.0)
    for m in ("psm", "psm_star", "kogelnik", "tcwt"):
        eta = tw.efficiency(m, "sigma", **p)
        assert abs(eta - limit) < 2e-4, (m, eta)
    # finite density, eq. (145): D0 = D1 = 1
    D0 = 1.0
    x = np.sqrt(3.0) * D0 / 2.0
    finite = 1.0 / (4.0 + 4.0 * np.sqrt(3.0) / np.tanh(x) + 3.0 / np.tanh(x) ** 2)
    p1 = _bragg(0.0, 10.0, n1=0.0, D0=1.0, D1=1.0)
    for m in ("psm", "psm_star", "kogelnik"):
        eta = tw.efficiency(m, "sigma", **p1)
        assert abs(eta - finite) / finite < 1e-4, (m, eta, finite)


def test_pi_absorption_reflection_eq153():
    # theta_c = 10 deg: p = cos 20 deg gives 0.0622 rather than 0.0718 (Section 7.2)
    p = _bragg(0.0, 10.0, n1=0.0, D0=8.0, D1=8.0)
    pp = abs(np.cos(20.0 * deg))
    g = np.sqrt(1.0 - pp * pp / 4.0)
    limit = (1.0 - g) / (1.0 + g)
    assert abs(limit - 0.0622) < 5e-4
    for m in ("psm", "psm_star", "kogelnik"):
        eta = tw.efficiency(m, "pi", **p)
        assert abs(eta - limit) / limit < 2e-3, (m, eta, limit)


def test_transmission_absorption_optimum_1_over_27_eq162_164():
    # unslanted transmission grating: Psi = 90 deg, D0 = D1 = ln 3 -> 1/27
    p = _bragg(90.0, 20.0, n1=0.0, D0=LN3, D1=LN3)
    for m in ("psm", "psm_star", "kogelnik", "tcwt"):
        eta = tw.efficiency(m, "sigma", **p)
        assert abs(eta - 1.0 / 27.0) / (1.0 / 27.0) < 1e-4, (m, eta)
    # and it is a maximum: D0 = D1 = 0.8 and 1.4 give less
    for D in (0.8, 1.4):
        pD = _bragg(90.0, 20.0, n1=0.0, D0=D, D1=D)
        assert tw.efficiency("psm", "sigma", **pD) < 1.0 / 27.0
    # eq. (162): exp(-2 D0) sinh^2(D0/2)
    for D in (0.5, 1.0, 2.0):
        pD = _bragg(90.0, 20.0, n1=0.0, D0=D, D1=D)
        expected = np.exp(-2 * D) * np.sinh(D / 2) ** 2
        assert abs(tw.efficiency("psm", "sigma", **pD) - expected) / expected < 1e-4


def test_pi_transmission_absorption_eq167_169():
    theta_c = 20.0
    pp = abs(np.cos(2 * theta_c * deg))
    Dopt = np.log((2 + pp) / (2 - pp)) / pp                        # eq. (168)
    eta_max = pp ** 2 / (4 - pp ** 2) * ((2 - pp) / (2 + pp)) ** (2 / pp)   # eq. (169)
    p = _bragg(90.0, theta_c, n1=0.0, D0=Dopt, D1=Dopt)
    for m in ("psm", "psm_star", "kogelnik", "tcwt"):
        eta = tw.efficiency(m, "pi", **p)
        assert abs(eta - eta_max) / eta_max < 1e-4, (m, eta, eta_max)


def test_lossless_transmission_sin2():
    # eq. (161): eta = sin^2(kappa d / c) with c = cos(theta_c - Psi) = sin(theta_c) at Psi = 90
    p = _bragg(90.0, 20.0, n1=0.02)
    b = cv.beta(p["lam_c"], p["n0"])
    kap = b * p["n1"] / (2 * p["n0"])
    expected = np.sin(kap * p["d"] / np.sin(p["theta_c"])) ** 2
    for m in ("psm", "psm_star", "kogelnik"):
        assert abs(tw.efficiency(m, "sigma", **p) - expected) < 1e-12
