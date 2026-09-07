# Conventions

Single source of truth: `python/psmgrating/conventions.py`, mirrored by `src/models/geometry.js`.
Equation numbers refer to Brotherton-Ratcliffe & Amos Sureshkumar (2026).

| Quantity | Definition | Eq. | Trap it replaces |
|---|---|---|---|
| Phasors | exp(+iωt), exp(−ik·r); passive medium ñ = n − iχ, χ > 0 | §2 fn. 1 | the August 2026 engines use the conjugate convention; efficiencies agree, amplitudes and the sign of ϑ do not |
| α | λc / λr | 5 | — |
| β | 2π n0 / λc | 3 | — |
| Ks | 2αβ cos θr, **signed**; negative wherever θr is obtuse | 4 | unsigned fringe spacing |
| θc, θr | measured from the **directed fringe normal** e_K; the reference wave makes θc − Ψ with the boundary normal | 21–22 | angles from the plate normal |
| Ψ | slant; Ψ = 0 unslanted reflection, Ψ = 90° unslanted transmission | 20 | — |
| Geometry | reflection if c_S^K = cos(θc − Ψ) − 2α cos θr cos Ψ < 0, else transmission | 46 | classification by sign(C_S/C_R) |
| χ0, χ1 | D λc cos(θc − Ψ)/(2πd): density along the **reference ray** | 143–144, 175 | D λc cos θc/(2πd), which gives gain for obtuse θc |
| ζ | cos θc/(α cos θr); PSM\* = PSM with κ̂ → ζκ̂ | 74 | — |
| π, PSM | κ̂π = κ̂σ cos 2θc (angle-only) | 41 | — |
| π, PSM\* | κ̂\*π = ζ P_on κ̂σ, P_on the on-shell projection | 85–86 | cos 2θc |
| π, TCWT | κπ = κ P_π, P_π = 1 − 2α cos θr cos θc; loss coefficient α̂, not α̃ | 125–126, 149 | Kogelnik's cos 2θc applied to TCWT |
| π, Kogelnik | κπ = κ cos 2θc; loss α̂ | §6.4.2 | — |
| σ, TCWT | α̃ = α̂ + iκ²/β | 100 | — |
| Flux factor | ρ = \|C_S / C_R\| | 50, 57 | — |
| Υ → 0 | (1 − e^{−2Υd})/Υ by its series below \|Υd\| = 10⁻⁸ | 61 | division by zero on the degenerate locus |
| Exterior | index-matched lossless n0 on both sides | A.1 | plate-surface Fresnel terms |
| Ψ = 0 | chain matrix, 1024 sublayers per period, binary exponentiation with rescaling | A.6 | the single Floquet order, which is not observable there |
| Truncation | production \|l\| ≤ 10; admissibility L_conv ≤ 4 for 10⁻⁴ agreement | §8 | the Klein–Cook parameter |

Family mode (the paper's Section 8 construction): one reference beam at θ_air outside, θ_in = asin(sin θ_air / n0) inside; θc = Ψ + θ_in; θr = θc; Bragg at λr.
