"""Convention guards: the paper's definitions, not the August 2026 code's."""
import numpy as np

from psmgrating import conventions as cv
from psmgrating import presets

deg = cv.deg


def test_alpha_is_replay_over_recording():
    assert cv.alpha(600e-9, 500e-9) == 1.2


def test_chi0_uses_reference_direction_cosine():
    fam = presets.FAMILY1
    p = presets.configuration(fam, psi_deg=85.0, d=20e-6)
    # theta_c = 85 + 33.86 = 118.86 deg is obtuse: cos(theta_c) < 0
    assert np.cos(p["theta_c"]) < 0
    chi0 = cv.chi_from_density(1.0, 500e-9, p["theta_c"], p["psi"], 20e-6)
    expected = 1.0 * 500e-9 * np.cos(33.86 * deg) / (2 * np.pi * 20e-6)
    assert chi0 > 0
    assert abs(chi0 - expected) / expected < 1e-3
    august = 1.0 * 500e-9 * np.cos(p["theta_c"]) / (2 * np.pi * 20e-6)
    assert august < 0                      # the August convention would give gain


def test_Ks_is_signed():
    fam = presets.FAMILY1
    p_t = presets.configuration(fam, psi_deg=85.0, d=20e-6)
    p_r = presets.configuration(fam, psi_deg=-29.0, d=20e-6)
    assert cv.Ks(500e-9, 500e-9, 1.5, p_t["theta_r"]) < 0
    assert cv.Ks(500e-9, 500e-9, 1.5, p_r["theta_r"]) > 0


def test_geometry_classification_by_cSK():
    fam = presets.FAMILY1
    for psi in fam["psi_transmission"]:
        p = presets.configuration(fam, psi_deg=psi, d=20e-6)
        for lc in (400e-9, 500e-9, 600e-9):
            assert cv.geometry(lc, p["lam_r"], p["theta_c"], p["theta_r"], p["psi"]) == "transmission"
    for psi in fam["psi_reflection"]:
        p = presets.configuration(fam, psi_deg=psi, d=20e-6)
        for lc in (400e-9, 500e-9, 600e-9):
            assert cv.geometry(lc, p["lam_r"], p["theta_c"], p["theta_r"], p["psi"]) == "reflection"


def test_zeta_is_one_at_bragg_and_departs_off_it():
    p = presets.configuration(presets.FAMILY1, psi_deg=-15.0, d=20e-6)
    assert abs(cv.zeta(500e-9, 500e-9, p["theta_c"], p["theta_r"]) - 1.0) < 1e-15
    assert abs(cv.zeta(550e-9, 500e-9, p["theta_c"], p["theta_r"]) - 1.0) > 0.05


def test_family_geometry_numbers_from_the_paper():
    # Section 8.1.1: 56.70 deg in air is 33.86 deg inside n0 = 1.5
    assert abs(cv.internal_angle(56.70 * deg, 1.5) / deg - 33.86) < 0.01
    # Table 1: at Psi = 90 the fringe spacing is 299 nm, at 55 it is 8.4 um
    p90 = presets.configuration(presets.FAMILY1, psi_deg=90.0, d=20e-6)
    p55 = presets.configuration(presets.FAMILY1, psi_deg=55.0, d=20e-6)
    assert abs(cv.fringe_spacing(500e-9, 1.5, p90["theta_r"]) * 1e9 - 299) < 1.0
    assert abs(cv.fringe_spacing(500e-9, 1.5, p55["theta_r"]) * 1e6 - 8.4) < 0.1
    # Table 1: Q_KC = 564 at Psi = 90, d = 20 um
    Lam = cv.fringe_spacing(500e-9, 1.5, p90["theta_r"])
    Q = cv.klein_cook(500e-9, 20e-6, 1.5, Lam, p90["theta_c"] - p90["psi"])
    assert abs(Q - 564) < 2
    # Section 8.2.1: the reflection fringe spacing stays between 167 and 194 nm
    for psi in (-31.0, -3.0):
        p = presets.configuration(presets.FAMILY1, psi_deg=psi, d=20e-6)
        assert 166e-9 < cv.fringe_spacing(500e-9, 1.5, p["theta_r"]) < 195e-9
    # Section 8.3: family 2, 40.00 deg in air is 23.38 deg inside n0 = 1.62
    assert abs(cv.internal_angle(40.00 * deg, 1.62) / deg - 23.38) < 0.01
