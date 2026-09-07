# VENDORED UNCHANGED from Code2/1_Optics/Research/DBR_volume_gratings/validation/symbolic/gm_closed_form.py
# Author: Vivian Amos Sureshkumar (reconstruction project, 2026). Kept for parity tests only;
# the package code in psmgrating/ follows the September 2026 paper conventions instead.
"""
TASK A4 - Closed-form Fourier coefficients of the reciprocal permittivity profile
================================================================================

Context (David's TM note, Eq 21-22): the RCW TM (pi) formulation needs the
Toeplitz matrix G built from the Fourier coefficients g_m of

    g(phi) = 1 / eps(phi) = 1 / (a + b*cos(phi))^2 ,

where a = n0g = n0 - i*chi0 and b = n1g = n1 - i*chi1 are the COMPLEX mean and
modulation index amplitudes (time convention exp(+i*w*t), Im n < 0 = loss).
The note remarks 1/eps has infinitely many harmonics; the existing codes
(matlab/rcw_pi/rcw_rigorous_planar.m lines 85-94) compute g_m by FFT:

    eps_prof = (n0g + n1g*cos(phi)).^2;   gf = ifft(1./eps_prof);   % Nfft=4096

This script DERIVES the exact closed form, valid for complex a, b, and
validates it against high-precision quadrature and the code's FFT route.

CLOSED FORMS (main result)
--------------------------
Let  w = sqrt(a^2 - b^2)  with the branch chosen so that

    rho = (w - a)/b        satisfies |rho| < 1.

Then, with the convention  f_m = (1/2pi) Int_0^{2pi} f(phi) e^{-i m phi} dphi :

    c_m = FourierCoeff[ 1/(a + b*cos(phi)) ]   = rho^{|m|} / w
    g_m = FourierCoeff[ 1/(a + b*cos(phi))^2 ] = rho^{|m|} * (a + |m|*w) / w^3

Derivation route:
  (1) z = e^{i phi} substitution: a + b*cos(phi) = (b/(2z))*(z-rho)(z-1/rho)
      with rho + 1/rho = -2a/b, rho*(1/rho) = 1; the residue at the interior
      pole z = rho gives c_m = rho^{|m|}/w   (rho - 1/rho = 2w/b).
  (2) g(phi) = -(d/da) 1/(a+b*cos(phi))  termwise (uniform convergence for
      |rho|<1), with  dw/da = a/w  and  drho/da = -rho/w :
      g_m = -(d/da) c_m = rho^{|m|} (a + |m| w)/w^3.

BRANCH RULE for complex a, b (b != 0, a^2 != b^2):
  The two square roots +-w give  rho(+w)*rho(-w) = (a^2-w^2)/b^2 = 1,
  i.e. the two sign choices give rho and 1/rho. Exactly one choice has
  |rho| < 1 unless |rho| = 1 (which happens iff a/b is real with |a/b|<=1,
  i.e. eps(phi) has a zero on the real phi axis - nonphysical/singular).
  Equivalent test:  |1/rho|^2 - |rho|^2 = (|w+a|^2 - |w-a|^2)/|b|^2
                    = 4*Re(w*conj(a))/|b|^2 ,
  so |rho| < 1  <=>  Re(w*conj(a)) > 0. For physical media (Re a > 0) the
  principal square root usually already satisfies this; the safe recipe used
  below is: take principal w, and flip its sign if |(w-a)/b| >= 1.

All prints are ASCII. Reproducible: fixed RNG seed 20260731.
Requires: sympy, mpmath, numpy.
"""

import json
import os
import random
import time

import numpy as np
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_MD = os.path.join(HERE, "gm_closed_form_results.md")

report = {}          # results container for the markdown note
t_start = time.time()


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def branch_w_rho_mp(a, b):
    """mpmath: principal sqrt, flipped if needed so |rho|<1. Returns (w, rho)."""
    w = mp.sqrt(a * a - b * b)
    rho = (w - a) / b
    if abs(rho) >= 1:
        w = -w
        rho = (w - a) / b
    return w, rho


def branch_w_rho_np(a, b):
    """numpy complex scalar version of the branch rule."""
    w = np.sqrt(complex(a) ** 2 - complex(b) ** 2)
    rho = (w - a) / b
    if abs(rho) >= 1:
        w = -w
        rho = (w - a) / b
    return w, rho


def c_closed(a, b, m, lib=mp):
    w, rho = (branch_w_rho_mp if lib is mp else branch_w_rho_np)(a, b)
    return rho ** abs(m) / w


def g_closed(a, b, m, lib=mp):
    w, rho = (branch_w_rho_mp if lib is mp else branch_w_rho_np)(a, b)
    return rho ** abs(m) * (a + abs(m) * w) / w ** 3


# ============================================================================
# A4-1  SYMBOLIC DERIVATION (sympy)
# ============================================================================
print("=" * 78)
print("A4-1  Symbolic derivation of c_m and g_m")
print("=" * 78)

a_s, b_s, z_s, w_s = sp.symbols("a b z w")
mod_poly = w_s ** 2 - (a_s ** 2 - b_s ** 2)     # relation defining w
rho_s = (w_s - a_s) / b_s                       # rho as function of (a,b,w)
inv_rho_s = -(w_s + a_s) / b_s                  # will prove this equals 1/rho


def is_zero_mod_w(expr):
    """True iff expr == 0 modulo w^2 = a^2 - b^2 (expr rational in a,b,z,w)."""
    num, _den = sp.fraction(sp.together(sp.expand(expr)))
    rem = sp.rem(sp.expand(num), mod_poly, w_s)
    return sp.simplify(rem) == 0


checks_A41 = {}

# --- (i) factorisation:  b*(z-rho)*(z-1/rho) = b*z^2 + 2*a*z + b  -----------
lhs = b_s * (z_s - rho_s) * (z_s - inv_rho_s)
rhs = b_s * z_s ** 2 + 2 * a_s * z_s + b_s
checks_A41["factorisation"] = is_zero_mod_w(lhs - rhs)
# and inv_rho_s really is 1/rho (i.e. rho * inv_rho = 1 mod w^2=a^2-b^2)
checks_A41["inv_rho_identity"] = is_zero_mod_w(rho_s * inv_rho_s - 1)
print("(i)   b(z-rho)(z-1/rho) = b z^2 + 2 a z + b        :",
      "PASS" if checks_A41["factorisation"] else "FAIL")
print("      rho * (-(w+a)/b) = 1  (so 1/rho = -(w+a)/b)  :",
      "PASS" if checks_A41["inv_rho_identity"] else "FAIL")

# --- (ii) pole-difference identities used by the residues -------------------
checks_A41["rho_minus_invrho"] = is_zero_mod_w((rho_s - inv_rho_s) - 2 * w_s / b_s)
checks_A41["rho_plus_invrho"] = sp.simplify((rho_s + inv_rho_s) + 2 * a_s / b_s) == 0
print("(ii)  rho - 1/rho = 2w/b ,  rho + 1/rho = -2a/b    :",
      "PASS" if (checks_A41["rho_minus_invrho"] and checks_A41["rho_plus_invrho"]) else "FAIL")

# --- (iii) residue derivation of c_m for concrete m = 0..4 ------------------
#  c_m = (1/2pi) Int  e^{-i m phi}/(a+b cos phi) dphi ;  z = e^{i phi},
#  evenness gives c_m = c_{-m}, so compute the m>=0 coefficient with e^{+im phi}:
#  c_m = (1/2pi i) Contour[ 2 z^m / (b (z-rho)(z-1/rho)) ] ; only pole inside
#  the unit circle (|rho|<1 branch) is the simple pole z = rho.
ok = True
for mm in range(5):
    integrand = 2 * z_s ** mm / (b_s * (z_s - rho_s) * (z_s - inv_rho_s))
    res = sp.residue(integrand, z_s, rho_s)
    ok = ok and is_zero_mod_w(res - rho_s ** mm / w_s)
checks_A41["residue_c_m_0to4"] = ok
print("(iii) Res_{z=rho} => c_m = rho^m / w   (m=0..4)    :",
      "PASS" if ok else "FAIL")

# --- (iv) closed-form geometric resummation (proves c_m for ALL m at once) --
#  Sum_{m=-inf}^{inf} (rho^{|m|}/w) z^m = (1/w)(1-rho^2)/((1-rho z)(1-rho/z))
#  must equal 1/(a + b(z+1/z)/2) as a rational identity mod w^2 = a^2-b^2.
S = (1 - rho_s ** 2) / ((1 - rho_s * z_s) * (1 - rho_s / z_s)) / w_s
T = 1 / (a_s + b_s * (z_s + 1 / z_s) / 2)
checks_A41["resummation_c"] = is_zero_mod_w(S - T)
print("(iv)  Sum rho^|m| z^m / w  == 1/(a+b cos)          :",
      "PASS" if checks_A41["resummation_c"] else "FAIL")

# --- (v) differentiation lemmas and chain rule for g_m ----------------------
W_expl = sp.sqrt(a_s ** 2 - b_s ** 2)
rho_expl = (W_expl - a_s) / b_s
lem1 = sp.simplify(sp.diff(W_expl, a_s) - a_s / W_expl) == 0          # w' = a/w
lem2 = sp.simplify(sp.diff(rho_expl, a_s) + rho_expl / W_expl) == 0   # rho' = -rho/w
checks_A41["lemma_dw_da"] = lem1
checks_A41["lemma_drho_da"] = lem2
print("(v)   dw/da = a/w ,  drho/da = -rho/w              :",
      "PASS" if (lem1 and lem2) else "FAIL")

#  chain rule with ABSTRACT functions r(a), W(a) obeying the two lemmas:
#  c_m = r^m/W  ==>  -dc_m/da = r^m (a + m W)/W^3   for symbolic m
aa = sp.symbols("alpha")                       # differentiation variable
m_s = sp.symbols("m", integer=True, nonnegative=True)
r_f = sp.Function("r")(aa)
W_f = sp.Function("W")(aa)
c_f = r_f ** m_s / W_f
minus_dc = -sp.diff(c_f, aa)
minus_dc = minus_dc.subs({sp.Derivative(r_f, aa): -r_f / W_f,
                          sp.Derivative(W_f, aa): aa / W_f})
target = r_f ** m_s * (aa + m_s * W_f) / W_f ** 3
checks_A41["chain_rule_gm_symbolic_m"] = sp.simplify(
    sp.powsimp(minus_dc - target, force=True)) == 0
print("(vi)  -d/da (r^m/W) = r^m (a+mW)/W^3  (symbolic m) :",
      "PASS" if checks_A41["chain_rule_gm_symbolic_m"] else "FAIL")

#  explicit-sqrt cross check for concrete m = 0..6
ok = True
for mm in range(7):
    gm_expl = -sp.diff(rho_expl ** mm / W_expl, a_s)
    tgt = rho_expl ** mm * (a_s + mm * W_expl) / W_expl ** 3
    ok = ok and (sp.simplify(sp.together(gm_expl - tgt)) == 0)
checks_A41["gm_explicit_m0to6"] = ok
print("(vii) explicit d/da check of g_m       (m=0..6)    :",
      "PASS" if ok else "FAIL")

# --- (viii) branch rule: the two roots +-w give rho and 1/rho ---------------
rho_plus = (w_s - a_s) / b_s
rho_minus = (-w_s - a_s) / b_s
checks_A41["branch_product"] = is_zero_mod_w(rho_plus * rho_minus - 1)
print("(viii) rho(+w)*rho(-w) = 1  (branch rule)          :",
      "PASS" if checks_A41["branch_product"] else "FAIL")

report["A41"] = checks_A41
assert all(checks_A41.values()), "A4-1 symbolic derivation failed: %s" % checks_A41


# ============================================================================
# A4-1b  HIGH-PRECISION NUMERIC VERIFICATION (mpmath, 50 digits, 24 draws)
# ============================================================================
print()
print("=" * 78)
print("A4-1b Numeric verification: 24 random complex draws, mpmath dps=50")
print("=" * 78)

mp.mp.dps = 50
rng = random.Random(20260731)
M_SET = (0, 1, 2, 3, 5, 8)
NPTS = 768                       # spectrally-exact trapezoid points (periodic)
NDRAWS = 24

max_err_c = mp.mpf(0)
max_err_g = mp.mpf(0)
branch_ok = True
draws = []

for idraw in range(NDRAWS):
    # draw a with Re(a)>0 (physical index) for 18 draws, arbitrary phase for 6
    while True:
        ra = 0.8 + 1.4 * rng.random()
        pha = (rng.random() - 0.5) * (2.0 if idraw < 18 else 6.0)
        a = mp.mpc(ra * mp.cos(pha), ra * mp.sin(pha))
        rb = (0.02 + 0.48 * rng.random()) * ra
        phb = 2 * mp.pi * rng.random()
        b = mp.mpc(rb * mp.cos(phb), rb * mp.sin(phb))
        wtry = mp.sqrt(a * a - b * b)
        r1 = abs((wtry - a) / b)
        if min(r1, 1 / r1) < 0.6:            # keep away from |rho| ~ 1
            break
    w, rho = branch_w_rho_mp(a, b)
    # --- branch rule checks
    rp = (mp.sqrt(a * a - b * b) - a) / b
    rm = (-mp.sqrt(a * a - b * b) - a) / b
    branch_ok &= abs(rp * rm - 1) < mp.mpf(10) ** (-45)
    branch_ok &= (abs(rho) < 1)
    branch_ok &= (mp.re(w * mp.conj(a)) > 0)     # equivalent criterion
    # --- spectrally exact trapezoid Fourier coefficients (error ~ rho^NPTS)
    phis = [2 * mp.pi * k / NPTS for k in range(NPTS)]
    fvals = [1 / (a + b * mp.cos(p)) for p in phis]
    gvals = [1 / (a + b * mp.cos(p)) ** 2 for p in phis]
    for m in M_SET:
        em = [mp.e ** (-1j * m * p) for p in phis]
        c_num = sum(f * e for f, e in zip(fvals, em)) / NPTS
        g_num = sum(g * e for g, e in zip(gvals, em)) / NPTS
        max_err_c = max(max_err_c, abs(c_num - rho ** m / w))
        max_err_g = max(max_err_g, abs(g_num - rho ** m * (a + m * w) / w ** 3))
    draws.append((complex(a), complex(b), float(abs(rho))))

print("draws                       : %d  (seed 20260731)" % NDRAWS)
print("|rho| range                 : %.4f .. %.4f"
      % (min(d[2] for d in draws), max(d[2] for d in draws)))
print("max |c_m^num - c_m^closed|  : %s" % mp.nstr(max_err_c, 3))
print("max |g_m^num - g_m^closed|  : %s" % mp.nstr(max_err_g, 3))
print("branch rule (product=1, unique |rho|<1, Re(w conj a)>0) :",
      "PASS" if branch_ok else "FAIL")

# independent adaptive-quadrature spot check (not DFT based) on 5 draws
max_err_quad = mp.mpf(0)
for (a_c, b_c, _r) in draws[:5]:
    a = mp.mpc(a_c)
    b = mp.mpc(b_c)
    w, rho = branch_w_rho_mp(a, b)
    for m in (0, 3, 8):
        cq = mp.quad(lambda p: mp.e ** (-1j * m * p) / (a + b * mp.cos(p)),
                     [0, mp.pi, 2 * mp.pi], maxdegree=9) / (2 * mp.pi)
        gq = mp.quad(lambda p: mp.e ** (-1j * m * p) / (a + b * mp.cos(p)) ** 2,
                     [0, mp.pi, 2 * mp.pi], maxdegree=9) / (2 * mp.pi)
        max_err_quad = max(max_err_quad, abs(cq - rho ** m / w),
                           abs(gq - rho ** m * (a + m * w) / w ** 3))
print("adaptive mp.quad spot check (5 draws x m=0,3,8), max err : %s"
      % mp.nstr(max_err_quad, 3))

report["A41b"] = dict(ndraws=NDRAWS,
                      max_err_c=mp.nstr(max_err_c, 3),
                      max_err_g=mp.nstr(max_err_g, 3),
                      max_err_quad=mp.nstr(max_err_quad, 3),
                      branch_ok=bool(branch_ok))
assert max_err_c < mp.mpf(10) ** (-40) and max_err_g < mp.mpf(10) ** (-40)
assert max_err_quad < mp.mpf(10) ** (-38)
assert branch_ok


# ============================================================================
# A4-2  VALIDATION vs FFT (numpy, N = 4096, m = -8..8)
# ============================================================================
print()
print("=" * 78)
print("A4-2  Closed form vs FFT (N=4096, m=-8..8), float64")
print("=" * 78)

NFFT = 4096
MS = np.arange(-8, 9)

cases = {
    "real   (n0=1.5,        n1=0.045)": (1.5 + 0j, 0.045 + 0j),
    "complex(n0g=1.5-2e-3i, n1g=0.045-2e-3i)": (1.5 - 0.002j, 0.045 - 0.002j),
}
report["A42"] = {}
phi = 2 * np.pi * np.arange(NFFT) / NFFT
for label, (a, b) in cases.items():
    gprof = 1.0 / (a + b * np.cos(phi)) ** 2
    # convention identical to the MATLAB code: gf = ifft(1./eps_prof)
    gf = np.fft.ifft(gprof)
    g_fft = gf[np.mod(MS, NFFT)]
    g_cf = np.array([g_closed(a, b, int(m), lib=np) for m in MS])
    err_abs = float(np.max(np.abs(g_fft - g_cf)))
    err_rel = float(np.max(np.abs(g_fft - g_cf) / np.abs(g_cf)))
    w, rho = branch_w_rho_np(a, b)
    print("%-42s max abs err = %.3e   max rel err = %.3e"
          % (label, err_abs, err_rel))
    report["A42"][label] = dict(err_abs=err_abs, err_rel=err_rel,
                                g0=complex(g_cf[8 + 0]),
                                rho=complex(rho), w=complex(w))
    assert err_abs < 1e-13, "FFT mismatch"

# note: ifft vs fft convention is immaterial because g(phi) is even (g_m=g_-m);
# verify that explicitly:
a, b = cases["complex(n0g=1.5-2e-3i, n1g=0.045-2e-3i)"]
gprof = 1.0 / (a + b * np.cos(phi)) ** 2
sym_err = float(np.max(np.abs(np.fft.ifft(gprof)[np.mod(MS, NFFT)]
                              - np.fft.fft(gprof)[np.mod(MS, NFFT)] / NFFT)))
print("evenness check |ifft - fft/N| over m=-8..8         : %.3e" % sym_err)


# ============================================================================
# A4-3  MATCH VIVIAN'S CONSTRUCTION (FFT of 1/eps as coded in rcw_rigorous_planar.m)
# ============================================================================
print()
print("=" * 78)
print("A4-3  Exact match to the code's construction (ifft(1./eps_prof))")
print("=" * 78)
# rcw_rigorous_planar.m lines 85-94:
#   phi      = 2*pi*(0:Nfft-1)/Nfft;         Nfft = 4096
#   eps_prof = (n0g + n1g*cos(phi)).^2;
#   gf = ifft(1./eps_prof);                  gm = gf(mod(ms,Nfft)+1)
#   G  = Toeplitz  G(p,l) = g_{p-l}
# Reproduce literally (eps first, then reciprocal, then ifft) and compare.
report["A43"] = {}
for label, (a, b) in cases.items():
    eps_prof = (a + b * np.cos(phi)) ** 2
    gf = np.fft.ifft(1.0 / eps_prof)               # the code's route
    # Toeplitz for M=8 (N=17 modes) from code route vs closed form
    Mtr = 8
    N = 2 * Mtr + 1
    ms_full = np.arange(-(N - 1), N)
    gm_code = gf[np.mod(ms_full, NFFT)]
    gm_exact = np.array([g_closed(a, b, int(m), lib=np) for m in ms_full])
    pp, ll = np.meshgrid(np.arange(N), np.arange(N), indexing="ij")
    G_code = gm_code[pp - ll + N - 1]
    G_exact = gm_exact[pp - ll + N - 1]
    errG = float(np.max(np.abs(G_code - G_exact)))
    # bonus: the code's eps harmonics cf = ifft(eps_prof) vs exact
    cf = np.fft.ifft(eps_prof)
    e_exact = np.zeros(2 * NFFT, dtype=complex)    # index shift trick
    e_m = {0: a * a + b * b / 2, 1: a * b, -1: a * b,
           2: b * b / 4, -2: b * b / 4}
    err_e = max(abs(cf[np.mod(m, NFFT)] - e_m.get(m, 0.0))
                for m in range(-8, 9))
    print("%-42s max |G_code - G_exact| = %.3e ; eps harm err = %.3e"
          % (label, errG, err_e))
    report["A43"][label] = dict(errG=errG, err_e=float(err_e))
    assert errG < 1e-13
print("=> the code's FFT-built Toeplitz G is EXACT to machine precision;")
print("   supports David's Section-10 row: reciprocal profile via FFT is the")
print("   correct finite-truncation construction (Li: Fourier-factorise 1/eps,")
print("   i.e. build G from the coefficients of 1/eps, never from inv(E)).")


# ============================================================================
# A4-4  THE WRONG CONSTRUCTION:  G_wrong = inv(E_truncated)
# ============================================================================
print()
print("=" * 78)
print("A4-4  Wrong construction G_wrong = inv(E_trunc) vs exact G, M = 1..6")
print("=" * 78)

def toeplitz_from(coeff_fn, Mtr):
    N = 2 * Mtr + 1
    ms_full = np.arange(-(N - 1), N)
    v = np.array([coeff_fn(int(m)) for m in ms_full])
    pp, ll = np.meshgrid(np.arange(N), np.arange(N), indexing="ij")
    return v[pp - ll + N - 1]

wrong_cases = {
    "weak   n0=1.5, n1=0.045 (baseline)": (1.5 + 0j, 0.045 + 0j),
    "strong n0=1.5, n1=0.3": (1.5 + 0j, 0.3 + 0j),
    "complex baseline n0g=1.5-2e-3i, n1g=0.045-2e-3i": (1.5 - 0.002j, 0.045 - 0.002j),
}
report["A44"] = {}
for label, (a, b) in wrong_cases.items():
    e_m = {0: a * a + b * b / 2, 1: a * b, -1: a * b, 2: b * b / 4, -2: b * b / 4}
    rows = []
    w, rho = branch_w_rho_np(a, b)
    for Mtr in range(1, 7):
        E = toeplitz_from(lambda m: e_m.get(m, 0.0), Mtr)
        G_wrong = np.linalg.inv(E)
        G_exact = toeplitz_from(lambda m: g_closed(a, b, m, lib=np), Mtr)
        err_max = float(np.max(np.abs(G_wrong - G_exact)))
        ctr = Mtr                                            # central index
        err_ctr = float(abs(G_wrong[ctr, ctr] - G_exact[ctr, ctr]))
        rows.append((Mtr, err_max, err_ctr))
    print("%s   (|rho| = %.5f)" % (label, abs(rho)))
    print("   M   max|G_wrong-G_exact|   central-elem err")
    for Mtr, e1, e2 in rows:
        print("   %d        %.3e            %.3e" % (Mtr, e1, e2))
    report["A44"][label] = dict(rho=float(abs(rho)), rows=rows)
print("=> error is tiny for weak modulation but orders of magnitude larger at")
print("   n1=0.3, and it never reaches machine precision at fixed M: inv() of")
print("   a truncated E is NOT the Toeplitz of 1/eps (Li's factorisation rules).")


# ============================================================================
# A4-5  DECAY RATE |rho| : why small M suffices in the TM eigenproblem
# ============================================================================
print()
print("=" * 78)
print("A4-5  Harmonic decay |g_m| ~ |rho|^m : truncation payoff")
print("=" * 78)
report["A45"] = {}
for label, (a, b) in cases.items():
    w, rho = branch_w_rho_np(a, b)
    g0 = g_closed(a, b, 0, lib=np)
    decs = [(m, abs(g_closed(a, b, m, lib=np) / g0)) for m in range(0, 9)]
    print("%s" % label)
    print("   |rho| = %.6f   (decay factor per harmonic ~ 1/%.1f)"
          % (abs(rho), 1 / abs(rho)))
    print("   m : " + "  ".join("%d" % m for m, _ in decs))
    print("   |g_m/g_0| : " + "  ".join("%.2e" % d for _, d in decs))
    report["A45"][label] = dict(rho=float(abs(rho)),
                                decay=[(m, float(d)) for m, d in decs])
print("=> baseline |rho| ~ 0.015: g_m falls ~66x per harmonic, so a Toeplitz")
print("   truncation M=3-4 already carries |g_m/g_0| ~ 1e-6..1e-7, consistent")
print("   with the observed M ~ 3-4 convergence of the TM RCW solver.")


# ============================================================================
# write the results note (ASCII markdown)
# ============================================================================
runtime = time.time() - t_start


def fmt_c(zc):
    return "%.15g%+.15gi" % (zc.real, zc.imag)


lines = []
A = lines.append
A("# A4 - Closed-form Fourier coefficients of 1/eps (reciprocal profile)")
A("")
A("Script: `validation/symbolic/gm_closed_form.py` (this note is auto-generated")
A("by that script; rerun it to regenerate). Runtime %.1f s. Seed 20260731." % runtime)
A("")
A("## Result (new derivation, valid for COMPLEX index / lossy gratings)")
A("")
A("With a = n0 - i chi0, b = n1 - i chi1, w = sqrt(a^2 - b^2) on the branch")
A("with |rho| < 1 where rho = (w - a)/b, and Fourier convention")
A("f_m = (1/2pi) Int_0^{2pi} f(phi) exp(-i m phi) dphi :")
A("")
A("    c_m = FourierCoeff[ 1/(a + b cos phi)   ] = rho^|m| / w")
A("    g_m = FourierCoeff[ 1/(a + b cos phi)^2 ] = rho^|m| (a + |m| w) / w^3")
A("")
A("Branch rule: the two roots +-w give rho(+w) rho(-w) = 1 exactly, so the sign")
A("choices give rho and 1/rho; pick the sign with |rho| < 1, equivalently")
A("Re(w conj(a)) > 0 (proved: |1/rho|^2 - |rho|^2 = 4 Re(w conj(a))/|b|^2).")
A("|rho| = 1 only if a/b is real in [-1,1] (eps has a real zero - excluded).")
A("Practical recipe: principal sqrt, flip sign if |(w-a)/b| >= 1.")
A("")
A("## A4-1 Derivation checks (sympy, all PASS)")
A("")
for k, v in report["A41"].items():
    A("- %s : %s" % (k, "PASS" if v else "FAIL"))
A("")
A("Route: z = exp(i phi) turns a + b cos phi into (b/2z)(z-rho)(z-1/rho);")
A("residue at the interior simple pole z = rho gives c_m (rho - 1/rho = 2w/b);")
A("then g(phi) = -(d/da) 1/(a+b cos phi) termwise with dw/da = a/w and")
A("drho/da = -rho/w gives g_m = -(d/da) c_m (chain rule proved for symbolic m).")
A("")
A("## A4-1b High-precision numeric verification (mpmath, dps=50)")
A("")
A("- %d random complex draws (18 physical Re a > 0, 6 arbitrary phase)" % NDRAWS)
A("- 768-point spectrally-exact periodic trapezoid, m in {0,1,2,3,5,8}")
A("- max |c_m num - closed| = %s" % report["A41b"]["max_err_c"])
A("- max |g_m num - closed| = %s" % report["A41b"]["max_err_g"])
A("- independent adaptive mp.quad spot check: max err = %s"
  % report["A41b"]["max_err_quad"])
A("- branch rule held on every draw (product = 1, unique |rho|<1, sign test)")
A("")
A("## A4-2 Closed form vs FFT (N = 4096, m = -8..8, float64)")
A("")
for label, d in report["A42"].items():
    A("- %s : max abs err %.3e, max rel err %.3e"
      % (label, d["err_abs"], d["err_rel"]))
A("")
A("## A4-3 Match to the repo's construction")
A("")
A("`matlab/rcw_pi/rcw_rigorous_planar.m` lines 85-94 build gf = ifft(1./eps_prof)")
A("with eps_prof = (n0g + n1g cos phi).^2, Nfft = 4096, then Toeplitz G(p,l) =")
A("g_{p-l}. Reproducing that literally (eps first, then reciprocal, then ifft):")
A("")
for label, d in report["A43"].items():
    A("- %s : max |G_code - G_exact| = %.3e (eps harmonics err %.3e)"
      % (label, d["errG"], d["err_e"]))
A("")
A("The code's FFT-built G is therefore EXACT to machine precision - it equals")
A("the Toeplitz of the true Fourier coefficients of 1/eps. This supports")
A("David's Section-10 table row: the reciprocal profile is handled by the")
A("correct finite-truncation construction (Li's rule: Fourier-factorise 1/eps).")
A("")
A("## A4-4 The wrong construction G_wrong = inv(E_truncated)")
A("")
for label, d in report["A44"].items():
    A("%s (|rho| = %.5f)" % (label, d["rho"]))
    A("")
    A("| M | max element err | central element err |")
    A("|---|-----------------|---------------------|")
    for Mtr, e1, e2 in d["rows"]:
        A("| %d | %.3e | %.3e |" % (Mtr, e1, e2))
    A("")
A("inv(E_trunc) != Toeplitz(1/eps): the max element error sits on the corner")
A("rows and PLATEAUS with M (baseline ~4e-4, n1=0.3 ~1.9e-2, never reaching")
A("machine precision) because the boundary rows of a truncated Toeplitz inverse")
A("solve a different (hard-truncated) problem - the numerical face of Li's")
A("Fourier factorisation rules. The CENTRAL element does converge, roughly like")
A("|rho|^(2M) with algebraic prefactors (baseline ratio per M ~ 2.5e3, close to")
A("1/|rho|^2 = 4.4e3; n1=0.3 ratio ~ 55-70 vs 1/|rho|^2 = 98), so the error is")
A("negligible for weak modulation but ~1e-3 at n1=0.3 and small M. Building G")
A("from the exact g_m (closed form or fine FFT) is both exact and cheap.")
A("")
A("## A4-5 Decay payoff")
A("")
for label, d in report["A45"].items():
    A("- %s : |rho| = %.6f, |g_m/g_0| at m = 4 is %.2e, at m = 8 is %.2e"
      % (label, d["rho"], d["decay"][4][1], d["decay"][8][1]))
A("")
A("Baseline |rho| ~ 0.0150: each harmonic falls by ~66x, so M = 3-4 keeps all")
A("neglected |g_m/g_0| below ~1e-6, matching the observed M ~ 3-4 convergence")
A("of the TM (pi) RCW solver.")
A("")

with open(RESULTS_MD, "w", encoding="ascii") as fh:
    fh.write("\n".join(lines))

print()
print("Results note written to: %s" % RESULTS_MD)
print("ALL A4 CHECKS PASSED  (runtime %.1f s)" % runtime)
