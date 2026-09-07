// The canonical two-wave system solved in closed form (paper eqs. 44-61 and
// 105-112), the model dispatcher and the wavelength sweep. Mirrors
// python/psmgrating/twowave.py exactly; parity is enforced by tests/js.
import * as Z from './complex.js';
import * as G from './geometry.js';
import { psmCoefficients } from './psm.js';
import { cwtCoefficients } from './cwt.js';

export const MODELS = ['psm', 'psm_star', 'kogelnik', 'tcwt'];
export const LABELS = { psm: 'PSM', psm_star: 'PSM*', kogelnik: 'Kogelnik', tcwt: 'TCWT' };

const SERIES_LIMIT = 1e-8;

// Upsilon flipped so Re(Upsilon d) >= 0; E = exp(-2x); F = (1 - E)/Upsilon with
// the series d (2 - 2x + 4/3 x^2) below |x| = 1e-8 (paper, after eq. 61).
export function continuation(Ups, d) {
  let x = Z.scale(Ups, d);
  if (x.re < 0) {
    Ups = Z.neg(Ups);
    x = Z.neg(x);
  }
  const E = Z.exp(Z.scale(x, -2));
  let F;
  if (Z.abs(x) < SERIES_LIMIT) {
    const x2 = Z.mul(x, x);
    F = Z.scale(Z.add(Z.sub(Z.C(2, 0), Z.scale(x, 2)), Z.scale(x2, 4 / 3)), d);
  } else {
    F = Z.div(Z.sub(Z.C(1, 0), E), Ups);
  }
  return { Ups, x, E, F };
}

// Complex signal amplitude: S(0) (reflection, eq. 48) or S(d) (transmission, eq. 56).
export function amplitude(CR, CS, a, theta, kappa, d, geometry) {
  const Q = Z.add(Z.mul(Z.sub(CR, CS), a), Z.mulI(Z.mul(CR, theta)));            // eq. 44
  const CRCS = Z.mul(CR, CS);
  const num = Z.sub(Z.mul(Q, Q), Z.scale(Z.mul(CRCS, Z.mul(kappa, kappa)), 4));
  const den = Z.scale(Z.mul(CRCS, CRCS), 4);
  const { x, E, F } = continuation(Z.sqrt(Z.div(num, den)), d);                     // eq. 45
  if (geometry === 'reflection') {
    const t = Z.div(Z.scale(Z.mul(CRCS, Z.add(Z.C(1, 0), E)), 2), F);              // 2 CS CR Ups coth(d Ups)
    return Z.div(Z.scale(Z.mulI(Z.mul(CR, kappa)), -2), Z.sub(Q, t));
  }
  if (geometry === 'transmission') {
    const e1 = Z.scale(Z.div(Z.mul(Z.add(CR, CS), a), CRCS), -d / 2);
    const e2 = Z.scale(Z.mulI(Z.div(theta, CS)), -d / 2);
    const expo = Z.add(Z.add(e1, e2), x);
    const pre = Z.neg(Z.mulI(Z.div(kappa, CS)));
    return Z.scale(Z.mul(pre, Z.mul(Z.exp(expo), F)), 0.5);
  }
  throw new Error("geometry must be 'reflection' or 'transmission'");
}

// eta = |C_S / C_R| |S|^2 (eqs. 50, 57)
export function efficiencyFromCoefficients(c, d, geometry) {
  const S = amplitude(c.CR, c.CS, c.a, c.theta, c.kappa, d, geometry);
  return Z.abs(Z.div(c.CS, c.CR)) * Z.abs2(S);
}

export function coefficients(model, pol, lam_c, lam_r, n0, n1, chi0, chi1, theta_c, theta_r, psi) {
  if (model === 'psm') return psmCoefficients(lam_c, lam_r, n0, n1, chi0, chi1, theta_c, theta_r, psi, pol, false);
  if (model === 'psm_star') return psmCoefficients(lam_c, lam_r, n0, n1, chi0, chi1, theta_c, theta_r, psi, pol, true);
  if (model === 'tcwt' || model === 'kogelnik') return cwtCoefficients(lam_c, lam_r, n0, n1, chi0, chi1, theta_c, theta_r, psi, pol, model);
  throw new Error(`unknown model ${model}`);
}

// p: {n0, n1, d, lam_r, theta_r, theta_c, psi, D0, D1, chi0?, chi1?, geometry?}
export function efficiency(model, pol, lam_c, p) {
  const chi0 = p.chi0 != null ? p.chi0 : G.chiFromDensity(p.D0 || 0, lam_c, p.theta_c, p.psi, p.d);
  const chi1 = p.chi1 != null ? p.chi1 : G.chiFromDensity(p.D1 || 0, lam_c, p.theta_c, p.psi, p.d);
  const geometry = p.geometry || G.geometry(p.lam_r, p.lam_r, p.theta_c, p.theta_r, p.psi);
  const c = coefficients(model, pol, lam_c, p.lam_r, p.n0, p.n1, chi0, chi1, p.theta_c, p.theta_r, p.psi);
  return efficiencyFromCoefficients(c, p.d, geometry);
}

// Wavelength sweep: returns {model: Float64Array}
export function sweep(models, pol, lam, p) {
  const out = {};
  for (const m of models) {
    const arr = new Float64Array(lam.length);
    for (let i = 0; i < lam.length; i++) arr[i] = efficiency(m, pol, lam[i], p);
    out[m] = arr;
  }
  return out;
}
