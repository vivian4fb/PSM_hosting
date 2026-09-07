// "Which model?" — the paper's practical rules applied to the current grating,
// plus live figures of merit against the rigorous curve when it is available.
import * as G from '../models/geometry.js';
import { rms, peakError, widthError } from '../models/metrics.js';
import { LABELS } from '../models/twowave.js';
import { RCWA } from '../config.js';

const fmt = (x, d = 4) => (Number.isFinite(x) ? x.toFixed(d) : '—');
const fmtE = (x) => (Number.isFinite(x) ? x.toExponential(2) : '—');

function propagatingOrders(p, lam_c, L) {
  const b = G.beta(lam_c, p.n0);
  const ks = G.Ks(lam_c, p.lam_r, p.n0, p.theta_r);
  const kx = b * Math.sin(p.theta_c - p.psi);
  const Kx = -ks * Math.sin(p.psi);
  let n = 0;
  for (let l = -L; l <= L; l++) if (Math.abs(kx + l * Kx) < b) n++;
  return n;
}

export function analyse(p, lam, curves, models, pol, extra) {
  const lam_r = p.lam_r;
  const cSK0 = G.cSK(lam_r, lam_r, p.theta_c, p.theta_r, p.psi);
  const geom = cSK0 < 0 ? 'reflection' : 'transmission';
  const phase = p.D0 === 0 && p.D1 === 0;
  const absorptionOnly = p.n1 === 0;
  const Lambda = G.fringeSpacing(p.lam_r, p.n0, p.theta_r);
  const Q = G.kleinCook(lam_r, p.d, p.n0, Lambda, p.theta_c - p.psi);
  const orders = propagatingOrders(p, lam_r, 40);
  const ref = curves.rcwa;
  const peakRef = ref ? Math.max(...ref) : (curves.psm_star ? Math.max(...curves.psm_star) : NaN);
  const zetaRange = [G.zeta(lam[0], lam_r, p.theta_c, p.theta_r), G.zeta(lam[lam.length - 1], lam_r, p.theta_c, p.theta_r)];
  const chi0 = G.chiFromDensity(p.D0, lam_r, p.theta_c, p.psi, p.d);
  const scopeOk = p.n1 / p.n0 < 1 - Math.abs(Math.sin(p.theta_c));

  // sign agreement between C_S (PSM mirror ray) and c_S^K across the band
  const CSsign = Math.sign(-Math.cos(p.theta_c + p.psi) * Math.cos(p.theta_c) / Math.cos(p.theta_r));
  let disagree = 0;
  for (let i = 0; i < lam.length; i++) if (Math.sign(G.cSK(lam[i], lam_r, p.theta_c, p.theta_r, p.psi)) !== CSsign) disagree++;

  const lines = [];
  const bullet = (cls, text) => lines.push({ cls, text });

  if (geom === 'reflection') {
    bullet('ok', `Reflection geometry (c_S^K = ${fmt(cSK0, 3)} at λr): PSM* is the model to prefer. In the paper it is closer to the rigorous solution than either coupled-wave model in 810/810 configurations of family 1 and 750/756 of family 2 (Tables 2–3).`);
    if (phase && !absorptionOnly && Number.isFinite(peakRef) && peakRef > 0.999) {
      bullet('warn', `Saturated phase reflection grating (peak η = ${fmt(peakRef, 4)} > 0.999): PSM and PSM* are then close and the published PSM is marginally ahead (Section 8.3). The whole broadband error sits in the sidelobes.`);
    }
    if (pol === 'sigma') {
      bullet(phase ? 'warn' : 'info', phase
        ? 'σ, reflection, phase: TCWT\'s iκ²/β term is harmful here (21 % worse than Kogelnik). Between the coupled-wave pair, use Kogelnik\'s coefficients (Section 8.2.2).'
        : 'σ, reflection with absorption: the κ² term helps slightly; TCWT is marginally the better coupled-wave model (Section 8.2.3–8.2.4).');
    } else {
      bullet('info', 'π, reflection: Kogelnik\'s angle-only cos 2θc is the better polarisation prescription for the coupled-wave pair (factors 1.7–2.9 over TCWT\'s band-varying P_π); PSM* uses the on-shell overlap P_on (Section 5.5).');
    }
  } else {
    bullet('ok', `Transmission geometry (c_S^K = ${fmt(cSK0, 3)} at λr): the coupled-wave lineshape is the one to prefer. In the paper it is closer than PSM* in every transmission cell, by median factors 1.5–10; PSM* still predicts the peak efficiency best (Section 8.1).`);
    bullet('info', pol === 'sigma'
      ? 'σ, transmission: TCWT (with the iκ²/β term) is closer than Kogelnik by 2 % (phase) to 19 % (mixed); for pure absorption they are indistinguishable.'
      : 'π, transmission: TCWT\'s band-varying projection P_π is the better prescription, by factors 1.4–3.1 over Kogelnik\'s angle-only cos 2θc (Section 8.1.3).');
  }

  bullet('info', `Fringe spacing Λ = ${(Lambda * 1e9).toFixed(0)} nm; Klein–Cook Q = ${Q.toFixed(0)}; ${orders} nominal orders propagate inside the medium at λr. The paper shows Q is not the diagnostic: the truncation the rigorous solution needs (L_conv ≤ 4) decides whether any two-wave model applies (Section 8, conclusion 4).`);
  bullet('info', `ζ = cos θc/(α cos θr) runs from ${fmt(zetaRange[0], 3)} to ${fmt(zetaRange[1], 3)} across the band (unity at Bragg); PSM* differs from PSM by exactly this factor on the coupling constant (eq. 74).`);
  bullet('info', `Mean loss χ0/n0 = ${fmtE(chi0 / p.n0)} at λr; the two-wave reductions are joint weak-loss, weak-modulation reductions (Section 6.1).`);

  if (!scopeOk) bullet('bad', `Outside the model's scope: n1/n0 = ${fmt(p.n1 / p.n0, 3)} exceeds 1 − |sin θc| = ${fmt(1 - Math.abs(Math.sin(p.theta_c)), 3)} (eq. 87): a ray no longer exists at every depth and the change of variables fails.`);
  if (pol === 'pi' && Math.abs(Math.abs(p.theta_c) - Math.PI / 4) < 2 * G.deg) bullet('warn', 'π polarisation near θc = 45°: the coupling vanishes (Brewster-like null, Section 7.2). Every model returns ≈ 0.');
  if (Math.abs(p.psi) < 0.01 * G.deg) bullet('warn', 'Ψ = 0 exactly: all nominal orders share one tangential wavevector; the rigorous curve shown is the stratified-medium (chain-matrix) reflectance of Appendix A.6, the physically observable quantity.');
  if (disagree > 0) bullet('bad', `sign(C_S) disagrees with sign(c_S^K) at ${disagree} of ${lam.length} wavelengths: the closed forms are to be read as restricted to configurations where the two agree (Section 4). Treat the curves there with caution.`);
  if (extra && extra.admissibility) {
    const a = extra.admissibility;
    bullet(a.maxDiff <= RCWA.admissibilityTol ? 'ok' : 'bad',
      `Admissibility: max|η(L=4) − η(L=8)| = ${fmtE(a.maxDiff)} over the band ${a.maxDiff <= RCWA.admissibilityTol ? '≤' : '>'} 1e-4. ${a.maxDiff <= RCWA.admissibilityTol ? 'A half-width of four suffices: the two-wave comparison is meaningful here.' : 'The rigorous solution needs more than four orders: no two-wave model is useful here (the paper excludes such slants).'}`);
  }

  // figures of merit
  const rows = [];
  if (ref) {
    for (const m of models) {
      const c = curves[m];
      if (!c) continue;
      rows.push({ model: LABELS[m], rms: rms(c, ref), peak: peakError(c, ref), width: widthError(Array.from(lam), c, ref) });
    }
  }
  return { geom, lines, rows, peakRef, Lambda, Q, orders };
}

export function render(el, tableEl, result, printed) {
  el.innerHTML = '';
  for (const l of result.lines) {
    const p = document.createElement('p');
    p.className = `v-${l.cls}`;
    p.textContent = l.text;
    el.append(p);
  }
  tableEl.innerHTML = '';
  if (!result.rows.length) {
    tableEl.innerHTML = '<p class="muted">Figures of merit appear when the rigorous curve is available.</p>';
    return;
  }
  const t = document.createElement('table');
  t.innerHTML = '<thead><tr><th>model</th><th>E_rms (eq. 172)</th><th>E_peak (eq. 173)</th><th>E_W (eq. 174)</th><th>paper</th></tr></thead>';
  const tb = document.createElement('tbody');
  for (const r of result.rows) {
    const tr = document.createElement('tr');
    const key = r.model === 'PSM*' ? 'psm_star' : r.model === 'TCWT' ? 'tcwt' : null;
    const pv = printed && key ? printed[key] : null;
    const ok = pv != null ? Math.abs(r.rms - pv) <= Math.max(0.02 * pv, 0.5 * Math.pow(10, Math.floor(Math.log10(pv)) - 1)) : null;
    tr.innerHTML = `<td>${r.model}</td><td>${fmtE(r.rms)}</td><td>${fmtE(r.peak)}</td><td>${fmtE(r.width)}</td>` +
      `<td>${pv != null ? `${pv} ${ok ? '✓' : '✗'}` : ''}</td>`;
    tb.append(tr);
  }
  t.append(tb);
  tableEl.append(t);
}
