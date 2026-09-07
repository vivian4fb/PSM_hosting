# Validation report

*Generated 2026-09-07 by `python/tools/make_validation_report.py`. Rigorous reference: RCWA, |l| <= 10, 1201 points over the family band.*

## The paper's printed figure values

Root-mean-square departure from the rigorous solution (eq. 172) for the configurations of the paper's representative-spectrum figures. "printed" is the value in the figure legend; "this build" is what `psmgrating` computes. PSM and Kogelnik are not printed in those legends and are given for completeness.

| Figure | PSM* printed | PSM* this build | TCWT printed | TCWT this build | PSM this build | Kogelnik this build | rigorous peak |
|---|---|---|---|---|---|---|---|
| Fig. 6 left | 0.0899 | 0.0899 | 0.0285 | 0.0285 | 0.1 | 0.0576 | 0.9245 |
| Fig. 6 right | 0.3386 | 0.339 | 0.3454 | 0.345 | 0.339 | 0.334 | 0.9308 |
| Fig. 9 left | 0.00026 | 0.000258 | 7.4e-05 | 7.37e-05 | 0.000377 | 7.36e-05 | 0.0369 |
| Fig. 9 right | 0.0027 | 0.00273 | 7.2e-05 | 7.18e-05 | 0.00748 | 7.12e-05 | 0.0368 |
| Fig. 12 left | 0.0111 | 0.0112 | 0.0037 | 0.00371 | 0.0127 | 0.00732 | 0.1235 |
| Fig. 12 right | 0.0643 | 0.0644 | 0.0677 | 0.0678 | 0.0307 | 0.0665 | 0.0974 |
| Fig. 17 left | 0.0341 | 0.0341 | 0.0715 | 0.0715 | 0.0471 | 0.054 | 0.9996 |
| Fig. 17 right | 0.0436 | 0.0436 | 0.0788 | 0.0788 | 0.055 | 0.0631 | 0.9992 |
| Fig. 20 left | 1.2e-05 | 1.2e-05 | 3.5e-05 | 3.45e-05 | 3.72e-05 | 3.48e-05 | 0.0514 |
| Fig. 20 right | 1.6e-06 | 1.64e-06 | 3.9e-05 | 3.93e-05 | 5.13e-05 | 3.96e-05 | 0.0514 |
| Fig. 23 left | 0.0042 | 0.00422 | 0.008 | 0.00799 | 0.00721 | 0.00742 | 0.8223 |
| Fig. 23 right | 0.0049 | 0.00492 | 0.0087 | 0.00874 | 0.0096 | 0.00871 | 0.8220 |

All 24 printed values are reproduced within the larger of 2 % and half a unit in the last printed digit (`tests/test_paper_targets.py`).

## Test suite

```
$ cd python && python -m pytest -q
39 passed in 20.67s
```

Contents: convention guards (`test_conventions.py`), Bragg-resonance limits of Section 7 (`test_bragg_limits.py`), closed forms against a Padé matrix-exponential solution of the 2x2 boundary-value problem at 200 random configurations (`test_closed_vs_bvp.py`), parity with the vendored August 2026 engines at Ψ = 0 (`test_legacy_parity.py`), RCWA energy conservation, convergence, weak-modulation and absorption checks and the chain-matrix solver (`test_rcwa_energy.py`), and the printed targets above.

JavaScript parity: `node --test tests/js/models.test.js` compares 1,680 efficiencies and 480 coefficient sets with the Python reference to 1e-10.
