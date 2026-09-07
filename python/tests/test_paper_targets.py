"""THE ACCEPTANCE GATE: the rms values printed in the paper's figure legends
(Figs. 6, 9, 12, 17, 20, 23), PSM* and TCWT against RCWA |l| <= 10 over the
full 400-600 nm band on 1201 points, reproduced to 2 per cent."""
import json
import os

import numpy as np
import pytest

from psmgrating import metrics, presets, rcwa, twowave as tw

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, "fixtures", "paper_targets_result.json")


def half_ulp(x):
    """Half a unit in the last printed significant figure of x (the paper
    prints 2-4 significant figures)."""
    s = "%.4g" % x
    mant = s.split("e")[0].replace("-", "").replace(".", "").lstrip("0")
    sig = len(mant.rstrip("0")) or 1
    exp10 = int(np.floor(np.log10(abs(x))))
    return 0.5 * 10.0 ** (exp10 - sig + 1)


@pytest.mark.parametrize("preset", presets.FIGURE_PRESETS, ids=[p["id"] for p in presets.FIGURE_PRESETS])
def test_figure_rms(preset):
    fam, p = presets.preset_parameters(preset)
    lam = presets.wavelength_grid(fam)
    pol = preset["pol"]
    rig = rcwa.sweep(lam_c=lam, pol=pol, L=10, **p)
    res = {}
    for m in ("psm_star", "tcwt", "psm", "kogelnik"):
        eta = tw.efficiency(m, pol, lam_c=lam, **p)
        res[m] = dict(rms=metrics.rms(eta, rig), peak=metrics.peak_error(eta, rig))
    row = dict(id=preset["id"], fig=preset["fig"], printed_psm_star=preset["rms_psm_star"],
               printed_tcwt=preset["rms_tcwt"], computed=res,
               rigorous_peak=float(np.max(rig)), rigorous_mean=float(np.mean(rig)))
    try:
        data = json.load(open(OUT)) if os.path.exists(OUT) else {}
    except Exception:
        data = {}
    data[preset["id"]] = row
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(data, open(OUT, "w"), indent=1)
    for m, key in (("psm_star", "rms_psm_star"), ("tcwt", "rms_tcwt")):
        printed = preset[key]
        got = res[m]["rms"]
        tol = max(0.02 * printed, half_ulp(printed))
        assert abs(got - printed) < tol, (preset["id"], m, got, printed, tol)
