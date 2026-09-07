// PSM (eqs. 33, 35, 41) and PSM* (eqs. 74, 85-86) coefficient sets.
import * as Z from './complex.js';
import * as G from './geometry.js';

// P_on = k_c . k_S^on / beta^2, eq. (85). Complex when the signal is evanescent.
export function onShellOverlap(lam_c, lam_r, n0, theta_c, theta_r, psi) {
  const b = G.beta(lam_c, n0);
  const ks = G.Ks(lam_c, lam_r, n0, theta_r);
  const kcx = b * Math.sin(theta_c - psi);
  const kcy = b * Math.cos(theta_c - psi);
  const kSx = kcx + ks * Math.sin(psi);                       // K'_x' = -K_s sin(Psi)
  const sgn = G.cSK(lam_c, lam_r, theta_c, theta_r, psi) < 0 ? -1 : 1;
  const kSy = Z.scale(Z.sqrt(Z.C(b * b - kSx * kSx, 0)), sgn);
  return Z.scale(Z.add(Z.C(kcx * kSx, 0), Z.scale(kSy, kcy)), 1 / (b * b));
}

export function psmCoefficients(lam_c, lam_r, n0, n1, chi0, chi1, theta_c, theta_r, psi, pol, star) {
  const a = G.alpha(lam_c, lam_r);
  const b = G.beta(lam_c, n0);
  const cc = Math.cos(theta_c);
  const cr = Math.cos(theta_r);
  const CR = (cc * Math.cos(theta_c - psi)) / (a * cr);
  const CS = (-cc * Math.cos(theta_c + psi)) / (a * cr);
  const abar = (chi0 * b * cc) / (a * n0 * cr);
  const theta = 2 * b * (cc / (a * cr) - 1) * cc * cc;
  const den = n0 * n0 + chi0 * chi0;
  let kappa = Z.C(((b / 2) * (n0 * n1 + chi0 * chi1)) / den, (-(b / 2) * (n0 * chi1 - n1 * chi0)) / den);
  if (star) kappa = Z.scale(kappa, G.zeta(lam_c, lam_r, theta_c, theta_r));
  if (pol === 'pi') {
    kappa = star
      ? Z.mul(kappa, onShellOverlap(lam_c, lam_r, n0, theta_c, theta_r, psi))
      : Z.scale(kappa, Math.cos(2 * theta_c));
  } else if (pol !== 'sigma') {
    throw new Error("pol must be 'sigma' or 'pi'");
  }
  return { CR: Z.C(CR), CS: Z.C(CS), a: Z.C(abar), theta: Z.C(theta), kappa };
}
