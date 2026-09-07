"""At Psi = 0 the August code's conventions coincide with the paper's for the
quantities that enter the efficiency; the sign conventions differ (conjugate
phasor) but eta must agree to machine precision for PSM (lossless and lossy)
and for TCWT (lossless only: the August a^ carries an extra chi0/(2 n0) term)."""
import numpy as np
import pytest

from psmgrating import conventions as cv
from psmgrating import twowave as tw

um = pytest.importorskip("legacy.universal_models")

deg = cv.deg


@pytest.mark.parametrize("D0,chi1f", [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)])
def test_psm_reflection_parity_at_psi0(D0, chi1f):
    n0, n1, d, lam_r = 1.5, 0.045, 20e-6, 532e-9
    theta_c, theta_r = 10 * deg, 20 * deg
    for lam_c in np.linspace(520e-9, 600e-9, 9):
        chi0 = cv.chi_from_density(D0, lam_c, theta_c, 0.0, d)
        chi1 = chi1f * chi0
        mine = tw.efficiency("psm", "sigma", n0=n0, n1=n1, d=d, lam_r=lam_r, theta_r=theta_r,
                             lam_c=lam_c, theta_c=theta_c, psi=0.0, chi0=chi0, chi1=chi1)
        legacy = um.eta_two_wave(lam_c, lam_r, n0, n1, D0, theta_c, theta_r, 0.0, d,
                                 theory="PSM", pol="sigma", chi1_factor=chi1f)
        assert abs(float(mine) - legacy) < 1e-12, (lam_c, mine, legacy)


def test_tcwt_lossless_parity_at_psi0():
    n0, n1, d, lam_r = 1.5, 0.045, 20e-6, 532e-9
    theta_c, theta_r = 10 * deg, 20 * deg
    for lam_c in np.linspace(520e-9, 600e-9, 9):
        mine = tw.efficiency("tcwt", "sigma", n0=n0, n1=n1, d=d, lam_r=lam_r, theta_r=theta_r,
                             lam_c=lam_c, theta_c=theta_c, psi=0.0)
        legacy = um.eta_two_wave(lam_c, lam_r, n0, n1, 0.0, theta_c, theta_r, 0.0, d,
                                 theory="TCWT", pol="sigma")
        assert abs(float(mine) - legacy) < 1e-12, (lam_c, mine, legacy)


def test_pi_psm_parity_at_psi0():
    n0, n1, d, lam_r = 1.5, 0.045, 20e-6, 532e-9
    theta_c, theta_r = 10 * deg, 20 * deg
    for lam_c in np.linspace(520e-9, 600e-9, 5):
        mine = tw.efficiency("psm", "pi", n0=n0, n1=n1, d=d, lam_r=lam_r, theta_r=theta_r,
                             lam_c=lam_c, theta_c=theta_c, psi=0.0)
        legacy = um.eta_two_wave(lam_c, lam_r, n0, n1, 0.0, theta_c, theta_r, 0.0, d,
                                 theory="PSM", pol="pi")
        assert abs(float(mine) - legacy) < 1e-12
