"""RCWA reference: energy conservation (eq. 187), convergence, weak-modulation
limit, absorption monotonicity, the K_x = 0 chain matrix (Appendix A.8)."""
import numpy as np

from psmgrating import chain, conventions as cv, presets, rcwa, twowave as tw

deg = cv.deg


def _cfg(psi_deg, n1=0.045, d=20e-6):
    return presets.configuration(presets.FAMILY1, psi_deg=psi_deg, d=d, n1=n1)


def test_energy_conservation_lossless():
    for psi, n1, d in ((7.5, 0.045, 20e-6), (-75.0, 0.045, 20e-6), (85.0, 0.08, 15e-6)):
        p = _cfg(psi, n1=n1, d=d)
        for pol in ("sigma", "pi"):
            for L in range(2, 8):
                for lc in (470e-9, 500e-9, 530e-9):
                    r = rcwa.solve(lam_c=lc, lam_r=p["lam_r"], n0=p["n0"], n1=p["n1"],
                                   chi0=0.0, chi1=0.0, theta_c=p["theta_c"], theta_r=p["theta_r"],
                                   psi=p["psi"], d=p["d"], pol=pol, L=L)
                    assert abs(r["sum_DE"] - 1.0) < 2e-12, (psi, pol, L, lc, r["sum_DE"])


def test_convergence_reflection_family():
    p = _cfg(-15.0, n1=0.065)
    for pol in ("sigma", "pi"):
        e3 = rcwa.sweep(lam_c=np.array([500e-9, 505e-9]), pol=pol, L=3, **p)
        e7 = rcwa.sweep(lam_c=np.array([500e-9, 505e-9]), pol=pol, L=7, **p)
        assert np.max(np.abs(e3 - e7)) < 1e-7


def test_weak_modulation_fourth_power():
    # |eta_RCWA - eta_Kogelnik| at Bragg falls as the fourth power of n1
    p = _cfg(-15.0, d=5e-6)
    diffs = []
    for n1 in (0.04, 0.02, 0.01, 0.005):
        q = dict(p, n1=n1)
        rig = rcwa.sweep(lam_c=np.array([500e-9]), pol="sigma", L=8, **q)[0]
        kog = float(tw.efficiency("kogelnik", "sigma", lam_c=500e-9, **q))
        diffs.append(abs(rig - kog))
    slopes = [np.log(diffs[i] / diffs[i + 1]) / np.log(2.0) for i in range(len(diffs) - 1)]
    assert all(d1 > d2 for d1, d2 in zip(diffs, diffs[1:])), diffs
    assert slopes[-1] > 3.0, (diffs, slopes)


def test_absorption_monotone():
    p = _cfg(-15.0, n1=0.045)
    prev = 2.0
    for D0 in (0.0, 0.3, 0.8, 1.5):
        chi0 = cv.chi_from_density(D0, 500e-9, p["theta_c"], p["psi"], p["d"])
        r = rcwa.solve(lam_c=500e-9, lam_r=p["lam_r"], n0=p["n0"], n1=p["n1"], chi0=chi0, chi1=0.0,
                       theta_c=p["theta_c"], theta_r=p["theta_r"], psi=p["psi"], d=p["d"], pol="sigma", L=6)
        assert r["sum_DE"] < prev + 1e-12
        prev = r["sum_DE"]
    assert prev < 0.9


def test_chain_matrix_energy_and_kogelnik_at_bragg():
    n0, n1, d, lam_r = 1.5, 0.045, 20e-6, 532e-9
    theta_c = theta_r = 10 * deg
    lam = np.linspace(520e-9, 545e-9, 7)
    for pol in ("sigma", "pi"):
        r = chain.solve(lam_c=lam, lam_r=lam_r, n0=n0, n1=n1, d=d, theta_c=theta_c,
                        theta_r=theta_r, pol=pol)
        assert np.max(np.abs(r["eta_R"] + r["eta_T"] - 1.0)) < 1e-10
    # sigma at Bragg vs Kogelnik tanh^2 (paper A.8 items 6-7; the project measured 6.9e-7)
    r = chain.solve(lam_c=np.array([lam_r]), lam_r=lam_r, n0=n0, n1=n1, d=d, theta_c=theta_c,
                    theta_r=theta_r, pol="sigma")
    kog = float(tw.efficiency("kogelnik", "sigma", n0=n0, n1=n1, d=d, lam_r=lam_r, theta_r=theta_r,
                              lam_c=lam_r, theta_c=theta_c, psi=0.0))
    assert abs(r["eta_R"][0] - kog) < 5e-4


def test_chain_matrix_matches_floquet_coherent_sum_at_kx0():
    n0, n1, d, lam_r = 1.5, 0.045, 20e-6, 532e-9
    theta_c = theta_r = 10 * deg
    lam = np.array([525e-9, 532e-9, 540e-9])
    st = chain.solve(lam_c=lam, lam_r=lam_r, n0=n0, n1=n1, d=d, theta_c=theta_c,
                     theta_r=theta_r, pol="sigma", D0=1.0)
    st2 = chain.solve(lam_c=lam, lam_r=lam_r, n0=n0, n1=n1, d=d, theta_c=theta_c,
                      theta_r=theta_r, pol="sigma", D0=1.0, sublayers=2048)
    for i, lc in enumerate(lam):
        chi0 = cv.chi_from_density(1.0, lc, theta_c, 0.0, d)
        r = rcwa.solve(lam_c=lc, lam_r=lam_r, n0=n0, n1=n1, chi0=chi0, chi1=0.0,
                       theta_c=theta_c, theta_r=theta_r, psi=0.0, d=d, pol="sigma", L=8)
        assert r["degenerate"]
        e1 = abs(r["DE_R_coherent"] - st["eta_R"][i])
        e2 = abs(r["DE_R_coherent"] - st2["eta_R"][i])
        # midpoint sublayers are second-order accurate: the residual is the
        # discretisation error of the plate, and doubling the sublayers should
        # cut it by about four (the project's Richardson study saw 1.4e-5 at 512)
        assert e1 < 5e-6, (lc, e1)
        assert e2 < e1 / 2.5, (lc, e1, e2)
