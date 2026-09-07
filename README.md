# Volume Grating Explorer

An interactive, browser-native companion to

> D. Brotherton-Ratcliffe and V. Sureshkumar, *Slanted planar volume gratings:
> parallel-stacked-mirror models versus traditional coupled wave descriptions* (2026, manuscript).

Four two-wave descriptions of a slanted, absorbing planar volume grating — **PSM**, **PSM\***,
**Kogelnik** and a truncated coupled-wave theory (**TCWT**) — are plotted live against a rigorous
coupled-wave (RCWA) solution of the same boundary-value problem. Every curve is computed in the
reader's browser: the analytic models in JavaScript, the rigorous reference by the Python package
in this repository running unchanged under Pyodide. Every configuration is a permanent URL.

**Status (7 September 2026): public preview accompanying an unpublished manuscript.** The paper has
not yet been peer-reviewed; PSM\* and the comparisons shown are subject to review. `FEATURES.psmStar`
in `src/config.js` switches the PSM\* trace off should that ever be needed.

## What is validated, and how

This is an independent implementation written from the paper's equations. It is *not* the
authors' code deposit (release v1.1), which is a separate code base.

| Check | Where | Result (7 Sept 2026) |
|---|---|---|
| Closed forms vs a Padé matrix-exponential solution of the 2×2 boundary-value problem, 200 random configurations, all four models, both polarisations, both geometries | `python/tests/test_closed_vs_bvp.py` | relative error < 1e-9 |
| Bragg-resonance limits: tanh² (eq. 141), 1/(7+4√3) (eq. 145), π forms (151, 153), 1/27 at D₀ = ln 3 (162–164), π transmission optimum (167–169) | `python/tests/test_bragg_limits.py` | pass |
| Convention guards: α = λc/λr, χ₀ = D₀λc cos(θc−Ψ)/(2πd), signed Ks, geometry by c_S^K, the paper's Table 1 numbers | `python/tests/test_conventions.py` | pass |
| RCWA energy conservation (eq. 187), L = 2…7, both polarisations | `python/tests/test_rcwa_energy.py` | \|Ση − 1\| < 2e-12 |
| **The paper's printed figure values**: rms of PSM\* and TCWT against RCWA (\|l\| ≤ 10) for the configurations of Figs. 6, 9, 12, 17, 20 and 23 | `python/tests/test_paper_targets.py` | all 24 numbers reproduced to the printed precision (see `fixtures/paper_targets_result.json`) |
| JavaScript ↔ Python parity, 1,680 efficiencies and 480 coefficient sets | `tests/js/models.test.js` | agree to 1e-10 |

Run them:

```
cd python && python -m pytest -q          # needs numpy; scipy and pytest for the tests
node --test tests/js/models.test.js       # Node 18+
```

## Layout

```
index.html                 the explorer (no build step; open via any static server)
src/models/                complex.js geometry.js psm.js cwt.js twowave.js metrics.js   (JS port)
src/ui/                    controls, plots (uPlot), verdict rules, export, theme
src/rcwa/                  worker.js (Pyodide) and client.js
python/psmgrating/         conventions, twowave, rcwa, chain, metrics, presets   (the reference)
python/legacy/             the August 2026 engines, vendored unchanged for parity tests
python/tests/              the checks above
python/tools/make_fixtures.py   regenerates fixtures/ and presets/figures.json
fixtures/  presets/        parity fixtures and the twelve paper-figure presets
vendor/                    uPlot 1.6.31 (MIT)
```

Serve locally with `python -m http.server 8080` from the repository root and open
`http://localhost:8080/`. The rigorous solver loads Pyodide (about 12 MB, cached by the browser)
from the jsDelivr CDN on first use.

## Conventions

The single source of truth is `python/psmgrating/conventions.py` (mirrored by
`src/models/geometry.js`). Angles from the directed fringe normal; Ψ = 0 unslanted reflection,
Ψ = 90° unslanted transmission; geometry by the sign of c_S^K (eq. 46); densities along the
reference ray (eq. 175); index-matched exterior (Appendix A.1); at Ψ = 0 exactly the chain matrix
of Appendix A.6 replaces the Floquet solver.

## Licence

Code MIT (`LICENSE`); data files CC BY 4.0 (`LICENSE-DATA`). Cite the paper and `CITATION.cff`.
