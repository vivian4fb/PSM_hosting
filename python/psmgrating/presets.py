"""psmgrating.presets -- the paper's two grating families and the figure
configurations whose printed figures of merit form the acceptance gate.

Family 1 (Section 8.1-8.2): n0 = 1.5, recorded at 500 nm with the reference
beam at 56.70 deg in air (33.86 deg inside), replayed with the same beam,
band 400-600 nm on 1201 points; thicknesses 5, 10, 20 um; phase n1 in
{0.03, 0.065, 0.08}; absorption D0 = D1 in {0.55, ln 3, 1.65}; mixed = both.
Transmission slants 90..55 deg by 2.5 (admissible >= 67.5); reflection slants
-31..-3 deg by 2.
Family 2 (Section 8.3): n0 = 1.62, 633 nm, 40.00 deg in air (23.38 inside),
band 533-733 nm, thicknesses 8, 30, 50 um; transmission 90..60 by 2.5
(admissible >= 77.5), reflection -24.5..+1.5 by 2.

Printed targets: the rms values quoted in the figure legends (PSM* and TCWT
against RCWA, |l| <= 10, over the full band).
"""
import numpy as np

from . import conventions as cv

__all__ = ["FAMILY1", "FAMILY2", "FAMILIES", "FIGURE_PRESETS", "configuration",
           "wavelength_grid", "preset_parameters", "LN3"]

LN3 = float(np.log(3.0))

FAMILY1 = dict(name="family1", n0=1.5, lam_r=500e-9, theta_air=56.70 * cv.deg,
               band=(400e-9, 600e-9), N=1201, thicknesses=(5e-6, 10e-6, 20e-6),
               n1_values=(0.03, 0.065, 0.08), D_values=(0.55, LN3, 1.65),
               psi_transmission=tuple(float(x) for x in np.arange(90.0, 54.9, -2.5)),
               psi_reflection=tuple(float(x) for x in np.arange(-31.0, -2.9, 2.0)),
               admissible_transmission_psi_min=67.5)

FAMILY2 = dict(name="family2", n0=1.62, lam_r=633e-9, theta_air=40.00 * cv.deg,
               band=(533e-9, 733e-9), N=1201, thicknesses=(8e-6, 30e-6, 50e-6),
               n1_values=(0.03, 0.065, 0.08), D_values=(0.55, LN3, 1.65),
               psi_transmission=tuple(float(x) for x in np.arange(90.0, 59.9, -2.5)),
               psi_reflection=tuple(float(x) for x in np.arange(-24.5, 1.6, 2.0)),
               admissible_transmission_psi_min=77.5)

FAMILIES = {"family1": FAMILY1, "family2": FAMILY2}


def configuration(family, *, psi_deg, d, n1=0.0, D0=0.0, D1=0.0):
    """Full parameter dict for one grating of a family (angles in radians)."""
    psi = psi_deg * cv.deg
    theta_c, theta_r = cv.family_angles(family["theta_air"], family["n0"], psi)
    return dict(n0=family["n0"], n1=n1, D0=D0, D1=D1, d=d, lam_r=family["lam_r"],
                theta_r=theta_r, theta_c=theta_c, psi=psi)


def wavelength_grid(family):
    lo, hi = family["band"]
    return np.linspace(lo, hi, family["N"])


FIGURE_PRESETS = [
    dict(id="fig06L", fig="Fig. 6 left", family="family1", psi_deg=85.0, d=20e-6, n1=0.065,
         D0=0.0, D1=0.0, pol="pi", rms_psm_star=0.0899, rms_tcwt=0.0285,
         caption="Phase transmission grating, pi, 2 propagating orders"),
    dict(id="fig06R", fig="Fig. 6 right", family="family1", psi_deg=65.0, d=20e-6, n1=0.065,
         D0=0.0, D1=0.0, pol="pi", rms_psm_star=0.3386, rms_tcwt=0.3454,
         caption="Phase transmission grating, pi, 7 propagating orders (inadmissible slant)"),
    dict(id="fig09L", fig="Fig. 9 left", family="family1", psi_deg=85.0, d=20e-6, n1=0.0,
         D0=LN3, D1=LN3, pol="sigma", rms_psm_star=2.6e-4, rms_tcwt=7.4e-5,
         caption="Absorption transmission grating, sigma"),
    dict(id="fig09R", fig="Fig. 9 right", family="family1", psi_deg=65.0, d=20e-6, n1=0.0,
         D0=LN3, D1=LN3, pol="sigma", rms_psm_star=0.0027, rms_tcwt=7.2e-5,
         caption="Absorption transmission grating, sigma (inadmissible slant)"),
    dict(id="fig12L", fig="Fig. 12 left", family="family1", psi_deg=85.0, d=20e-6, n1=0.065,
         D0=LN3, D1=LN3, pol="pi", rms_psm_star=0.0111, rms_tcwt=0.0037,
         caption="Mixed transmission grating, pi"),
    dict(id="fig12R", fig="Fig. 12 right", family="family1", psi_deg=65.0, d=20e-6, n1=0.065,
         D0=LN3, D1=LN3, pol="pi", rms_psm_star=0.0643, rms_tcwt=0.0677,
         caption="Mixed transmission grating, pi (inadmissible slant)"),
    dict(id="fig17L", fig="Fig. 17 left", family="family1", psi_deg=-29.0, d=20e-6, n1=0.065,
         D0=0.0, D1=0.0, pol="sigma", rms_psm_star=0.0341, rms_tcwt=0.0715,
         caption="Phase reflection grating, sigma, order leaves at -38 deg in air"),
    dict(id="fig17R", fig="Fig. 17 right", family="family1", psi_deg=-5.0, d=20e-6, n1=0.065,
         D0=0.0, D1=0.0, pol="sigma", rms_psm_star=0.0436, rms_tcwt=0.0788,
         caption="Phase reflection grating, sigma, order leaves at +37 deg in air"),
    dict(id="fig20L", fig="Fig. 20 left", family="family1", psi_deg=-29.0, d=20e-6, n1=0.0,
         D0=LN3, D1=LN3, pol="sigma", rms_psm_star=1.2e-5, rms_tcwt=3.5e-5,
         caption="Absorption reflection grating, sigma"),
    dict(id="fig20R", fig="Fig. 20 right", family="family1", psi_deg=-5.0, d=20e-6, n1=0.0,
         D0=LN3, D1=LN3, pol="sigma", rms_psm_star=1.6e-6, rms_tcwt=3.9e-5,
         caption="Absorption reflection grating, sigma"),
    dict(id="fig23L", fig="Fig. 23 left", family="family1", psi_deg=-29.0, d=20e-6, n1=0.065,
         D0=LN3, D1=LN3, pol="sigma", rms_psm_star=0.0042, rms_tcwt=0.0080,
         caption="Mixed reflection grating, sigma"),
    dict(id="fig23R", fig="Fig. 23 right", family="family1", psi_deg=-5.0, d=20e-6, n1=0.065,
         D0=LN3, D1=LN3, pol="sigma", rms_psm_star=0.0049, rms_tcwt=0.0087,
         caption="Mixed reflection grating, sigma"),
]


def preset_parameters(preset):
    """(family dict, parameter dict) for a figure preset."""
    fam = FAMILIES[preset["family"]]
    p = configuration(fam, psi_deg=preset["psi_deg"], d=preset["d"], n1=preset["n1"],
                      D0=preset["D0"], D1=preset["D1"])
    return fam, p
