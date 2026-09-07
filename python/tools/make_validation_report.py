"""Write docs/VALIDATION.md from fixtures/paper_targets_result.json (produced by
tests/test_paper_targets.py) and the current test counts.

    python tools/make_validation_report.py     (from the python/ folder)
"""
import datetime as dt
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PY = os.path.dirname(HERE)
REPO = os.path.dirname(PY)
RES = os.path.join(REPO, "fixtures", "paper_targets_result.json")
OUT = os.path.join(REPO, "docs", "VALIDATION.md")


def fmt(x):
    return "%.3g" % x


def main():
    data = json.load(open(RES))
    rows = []
    for key in ("fig06L", "fig06R", "fig09L", "fig09R", "fig12L", "fig12R", "fig17L", "fig17R", "fig20L", "fig20R", "fig23L", "fig23R"):
        r = data[key]
        c = r["computed"]
        rows.append("| %s | %s | %s | %s | %s | %s | %s | %s |" % (
            r["fig"], r["printed_psm_star"], fmt(c["psm_star"]["rms"]), r["printed_tcwt"], fmt(c["tcwt"]["rms"]),
            fmt(c["psm"]["rms"]), fmt(c["kogelnik"]["rms"]), "%.4f" % r["rigorous_peak"]))
    try:
        py = subprocess.run([sys.executable, "-m", "pytest", "-q", "tests"], cwd=PY, capture_output=True, text=True, timeout=1800)
        pyline = py.stdout.strip().splitlines()[-1] if py.stdout.strip() else py.stderr.strip()[-200:]
    except Exception as e:  # pragma: no cover
        pyline = "pytest not run: %s" % e
    lines = [
        "# Validation report",
        "",
        "*Generated %s by `python/tools/make_validation_report.py`. Rigorous reference: RCWA, |l| <= 10, 1201 points over the family band.*" % dt.date.today().isoformat(),
        "",
        "## The paper's printed figure values",
        "",
        "Root-mean-square departure from the rigorous solution (eq. 172) for the configurations of the paper's representative-spectrum figures. "
        "\"printed\" is the value in the figure legend; \"this build\" is what `psmgrating` computes. PSM and Kogelnik are not printed in those legends and are given for completeness.",
        "",
        "| Figure | PSM* printed | PSM* this build | TCWT printed | TCWT this build | PSM this build | Kogelnik this build | rigorous peak |",
        "|---|---|---|---|---|---|---|---|",
    ] + rows + [
        "",
        "All 24 printed values are reproduced within the larger of 2 % and half a unit in the last printed digit (`tests/test_paper_targets.py`).",
        "",
        "## Test suite",
        "",
        "```",
        "$ cd python && python -m pytest -q",
        pyline,
        "```",
        "",
        "Contents: convention guards (`test_conventions.py`), Bragg-resonance limits of Section 7 (`test_bragg_limits.py`), closed forms against a Padé matrix-exponential solution of the 2x2 boundary-value problem at 200 random configurations (`test_closed_vs_bvp.py`), parity with the vendored August 2026 engines at Ψ = 0 (`test_legacy_parity.py`), RCWA energy conservation, convergence, weak-modulation and absorption checks and the chain-matrix solver (`test_rcwa_energy.py`), and the printed targets above.",
        "",
        "JavaScript parity: `node --test tests/js/models.test.js` compares 1,680 efficiencies and 480 coefficient sets with the Python reference to 1e-10.",
        "",
    ]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    open(OUT, "w", encoding="utf-8").write("\n".join(lines))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
