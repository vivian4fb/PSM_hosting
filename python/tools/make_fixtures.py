"""Generate the JSON fixtures that tie the JavaScript port to the Python
reference, and the UI preset file.

    python tools/make_fixtures.py            (run from the python/ folder)

Outputs (repo root):
  fixtures/analytic_presets.json  12 figure presets x 25 wavelengths x 4 models
  fixtures/analytic_random.json   60 random configurations, both polarisations,
                                  with the full coefficient set of each model
  fixtures/rcwa_presets.json      12 presets x 5 wavelengths, L = 4 and 10
  presets/figures.json            UI presets (nm, um, degrees) with printed targets
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PY = os.path.dirname(HERE)
REPO = os.path.dirname(PY)
sys.path.insert(0, PY)

from psmgrating import conventions as cv, presets, rcwa, twowave as tw  # noqa: E402

FIX = os.path.join(REPO, "fixtures")
PRE = os.path.join(REPO, "presets")
os.makedirs(FIX, exist_ok=True)
os.makedirs(PRE, exist_ok=True)


def c2(z):
    z = complex(z)
    return [z.real, z.imag]


def analytic_presets():
    rows = []
    for pr in presets.FIGURE_PRESETS:
        fam, p = presets.preset_parameters(pr)
        lam = presets.wavelength_grid(fam)[::50]
        eta = {m: tw.efficiency(m, pr["pol"], lam_c=lam, **p).tolist() for m in tw.MODELS}
        rows.append(dict(id=pr["id"], pol=pr["pol"], params=p, lam=lam.tolist(), eta=eta))
    json.dump(rows, open(os.path.join(FIX, "analytic_presets.json"), "w"))
    print("analytic_presets.json:", len(rows), "presets")


def analytic_random():
    rng = np.random.default_rng(20260908)
    rows = []
    for _ in range(60):
        n0 = float(rng.uniform(1.45, 1.7))
        n1 = float(rng.uniform(0.0, 0.08))
        d = float(rng.uniform(4e-6, 40e-6))
        lam_r = float(rng.uniform(450e-9, 650e-9))
        lam_c = float(lam_r * rng.uniform(0.85, 1.15))
        reflection = rng.random() < 0.5
        psi = float((rng.uniform(-35, 10) if reflection else rng.uniform(60, 120)) * cv.deg)
        theta_c = float(psi + rng.uniform(5, 40) * cv.deg)
        theta_r = float(theta_c + rng.uniform(-3, 3) * cv.deg)
        D0 = float(rng.uniform(0, 1.6))
        D1 = float(rng.uniform(0, D0) if rng.random() < 0.7 else 0.0)
        chi0 = float(cv.chi_from_density(D0, lam_c, theta_c, psi, d))
        chi1 = float(cv.chi_from_density(D1, lam_c, theta_c, psi, d))
        geom = cv.geometry(lam_r, lam_r, theta_c, theta_r, psi)
        row = dict(params=dict(n0=n0, n1=n1, d=d, lam_r=lam_r, theta_r=theta_r, theta_c=theta_c,
                               psi=psi, D0=D0, D1=D1), lam_c=lam_c, chi0=chi0, chi1=chi1,
                   geometry=geom, coeff={}, eta={})
        for pol in ("sigma", "pi"):
            for m in tw.MODELS:
                CR, CS, a, th, kap = tw.coefficients(m, pol, lam_c, lam_r, n0, n1, chi0, chi1,
                                                     theta_c, theta_r, psi)
                row["coeff"]["%s:%s" % (m, pol)] = dict(CR=c2(CR), CS=c2(CS), a=c2(a), theta=c2(th), kappa=c2(kap))
                row["eta"]["%s:%s" % (m, pol)] = float(tw.efficiency(m, pol, lam_c=lam_c, **row["params"]))
        rows.append(row)
    json.dump(rows, open(os.path.join(FIX, "analytic_random.json"), "w"))
    print("analytic_random.json:", len(rows), "configurations")


def rcwa_presets():
    rows = []
    for pr in presets.FIGURE_PRESETS:
        fam, p = presets.preset_parameters(pr)
        lam = presets.wavelength_grid(fam)[::300]
        for L in (4, 10):
            eta = []
            sums = []
            orders = []
            for lc in lam:
                c0 = cv.chi_from_density(p["D0"], lc, p["theta_c"], p["psi"], p["d"])
                c1 = cv.chi_from_density(p["D1"], lc, p["theta_c"], p["psi"], p["d"])
                r = rcwa.solve(lam_c=lc, lam_r=p["lam_r"], n0=p["n0"], n1=p["n1"], chi0=c0, chi1=c1,
                               theta_c=p["theta_c"], theta_r=p["theta_r"], psi=p["psi"], d=p["d"],
                               pol=pr["pol"], L=L)
                eta.append(r["eta"])
                sums.append(r["sum_DE"])
                orders.append(dict(DE_R=r["DE_R"].tolist(), DE_T=r["DE_T"].tolist()))
            rows.append(dict(id=pr["id"], pol=pr["pol"], L=L, params=p, lam=lam.tolist(),
                             eta=eta, sum_DE=sums, orders=orders))
    json.dump(rows, open(os.path.join(FIX, "rcwa_presets.json"), "w"))
    print("rcwa_presets.json:", len(rows), "rows")


def ui_presets():
    out = dict(families={}, figures=[])
    for key, fam in presets.FAMILIES.items():
        out["families"][key] = dict(n0=fam["n0"], lam_r_nm=fam["lam_r"] * 1e9,
                                    theta_air_deg=fam["theta_air"] / cv.deg,
                                    theta_in_deg=cv.internal_angle(fam["theta_air"], fam["n0"]) / cv.deg,
                                    band_nm=[fam["band"][0] * 1e9, fam["band"][1] * 1e9], N=fam["N"],
                                    admissible_transmission_psi_min=fam["admissible_transmission_psi_min"])
    for pr in presets.FIGURE_PRESETS:
        out["figures"].append(dict(id=pr["id"], fig=pr["fig"], caption=pr["caption"], family=pr["family"],
                                   psi_deg=pr["psi_deg"], d_um=pr["d"] * 1e6, n1=pr["n1"], D0=pr["D0"],
                                   D1=pr["D1"], pol=pr["pol"], printed=dict(psm_star=pr["rms_psm_star"],
                                                                             tcwt=pr["rms_tcwt"])))
    json.dump(out, open(os.path.join(PRE, "figures.json"), "w"), indent=1)
    print("presets/figures.json:", len(out["figures"]), "figures")


if __name__ == "__main__":
    analytic_presets()
    analytic_random()
    rcwa_presets()
    ui_presets()
