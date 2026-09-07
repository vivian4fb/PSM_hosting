// TCWT (eqs. 92, 100, 125-126, 149) and Kogelnik (eq. 101, Section 6.4.2) coefficient sets.
import * as Z from './complex.js';
import * as G from './geometry.js';

export function cwtCoefficients(lam_c, lam_r, n0, n1, chi0, chi1, theta_c, theta_r, psi, pol, model) {
  if (model !== 'tcwt' && model !== 'kogelnik') throw new Error("model must be 'tcwt' or 'kogelnik'");
  const a = G.alpha(lam_c, lam_r);
  const b = G.beta(lam_c, n0);
  const kk0 = G.k0(lam_c);
  const cR = Math.cos(theta_c - psi);
  const cS = Math.cos(theta_c - psi) - 2 * a * Math.cos(theta_r) * Math.cos(psi);
  const ahat = (b * chi0) / n0;
  let kappa = Z.C((kk0 / 2) * n1, (-kk0 / 2) * chi1);
  const theta = 2 * a * b * Math.cos(theta_r) * (Math.cos(theta_c) - a * Math.cos(theta_r));
  let loss;
  if (pol === 'sigma') {
    loss = model === 'tcwt' ? Z.add(Z.C(ahat), Z.scale(Z.mulI(Z.mul(kappa, kappa)), 1 / b)) : Z.C(ahat);
  } else if (pol === 'pi') {
    loss = Z.C(ahat);                                                  // eq. 149
    kappa = model === 'tcwt'
      ? Z.scale(kappa, 1 - 2 * a * Math.cos(theta_r) * Math.cos(theta_c))   // P_pi, eq. 126
      : Z.scale(kappa, Math.cos(2 * theta_c));
  } else {
    throw new Error("pol must be 'sigma' or 'pi'");
  }
  return { CR: Z.C(cR), CS: Z.C(cS), a: loss, theta: Z.C(theta), kappa };
}
