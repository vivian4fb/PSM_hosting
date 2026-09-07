// Explorer state: UI units (nm, um, degrees) <-> URL hash <-> SI parameters.
import * as G from './models/geometry.js';

export const DEFAULTS = {
  v: 1,
  preset: 'fig17L',
  n0: 1.5, d_um: 20, lr_nm: 500, th_air: 56.7,
  adv: 0, tc: 0, tr: 0,
  psi: -29, n1: 0.065, D0: 0, D1: 0,
  pol: 'sigma', a_nm: 400, b_nm: 600, N: 1201,
  models: 'psm_star,tcwt,psm,kogelnik', rcwa: 1, L: 10,
};

const NUMERIC = new Set(['v', 'n0', 'd_um', 'lr_nm', 'th_air', 'adv', 'tc', 'tr', 'psi', 'n1', 'D0', 'D1', 'a_nm', 'b_nm', 'N', 'rcwa', 'L']);

export function fromHash(hash) {
  const s = { ...DEFAULTS };
  const q = (hash || '').replace(/^#/, '');
  if (!q) return s;
  for (const part of q.split('&')) {
    const [k, v] = part.split('=');
    if (!(k in DEFAULTS)) continue;
    const val = decodeURIComponent(v ?? '');
    s[k] = NUMERIC.has(k) ? Number(val) : val;
    if (NUMERIC.has(k) && !Number.isFinite(s[k])) s[k] = DEFAULTS[k];
  }
  return s;
}

export function toHash(s) {
  const parts = [];
  for (const k of Object.keys(DEFAULTS)) {
    if (k === 'v') { parts.push('v=1'); continue; }
    if (s[k] === DEFAULTS[k]) continue;
    parts.push(`${k}=${encodeURIComponent(s[k])}`);
  }
  return '#' + parts.join('&');
}

// SI parameter set used by the models and the worker.
export function toParams(s) {
  const n0 = s.n0;
  let theta_c;
  let theta_r;
  if (s.adv) {
    theta_c = s.tc * G.deg;
    theta_r = s.tr * G.deg;
  } else {
    const ang = G.familyAngles(s.th_air * G.deg, n0, s.psi * G.deg);
    theta_c = ang.theta_c;
    theta_r = ang.theta_r;
  }
  return {
    n0, n1: s.n1, d: s.d_um * 1e-6, lam_r: s.lr_nm * 1e-9,
    theta_r, theta_c, psi: s.psi * G.deg, D0: s.D0, D1: s.D1,
  };
}

export function wavelengths(s) {
  const N = Math.max(3, Math.round(s.N));
  const a = s.a_nm * 1e-9;
  const b = s.b_nm * 1e-9;
  const lam = new Float64Array(N);
  for (let i = 0; i < N; i++) lam[i] = a + ((b - a) * i) / (N - 1);
  return lam;
}

export function modelsOf(s) {
  return s.models.split(',').filter(Boolean);
}

export function applyPreset(s, fig, families) {
  const fam = families[fig.family];
  return {
    ...s,
    preset: fig.id,
    n0: fam.n0, lr_nm: fam.lam_r_nm, th_air: fam.theta_air_deg, adv: 0,
    psi: fig.psi_deg, d_um: fig.d_um, n1: fig.n1, D0: fig.D0, D1: fig.D1,
    pol: fig.pol, a_nm: fam.band_nm[0], b_nm: fam.band_nm[1], N: fam.N,
    rcwa: 1, L: 10,
  };
}
