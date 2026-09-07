// Geometric and material conventions of Brotherton-Ratcliffe & Amos (2026).
// Mirrors python/psmgrating/conventions.py; equation numbers are the paper's.
// Units: metres and radians.

export const deg = Math.PI / 180;

export const alpha = (lam_c, lam_r) => lam_c / lam_r;                          // eq. 5
export const beta = (lam_c, n0) => (2 * Math.PI * n0) / lam_c;                 // eq. 3
export const k0 = (lam_c) => (2 * Math.PI) / lam_c;
export const Ks = (lam_c, lam_r, n0, theta_r) =>
  2 * alpha(lam_c, lam_r) * beta(lam_c, n0) * Math.cos(theta_r);              // eq. 4, signed
export const zeta = (lam_c, lam_r, theta_c, theta_r) =>
  Math.cos(theta_c) / (alpha(lam_c, lam_r) * Math.cos(theta_r));              // eq. 74
export const cSK = (lam_c, lam_r, theta_c, theta_r, psi) =>
  Math.cos(theta_c - psi) - 2 * alpha(lam_c, lam_r) * Math.cos(theta_r) * Math.cos(psi); // eq. 46
export const geometry = (lam_c, lam_r, theta_c, theta_r, psi) =>
  cSK(lam_c, lam_r, theta_c, theta_r, psi) < 0 ? 'reflection' : 'transmission';
export const chiFromDensity = (D, lam_c, theta_c, psi, d) =>
  (D * lam_c * Math.cos(theta_c - psi)) / (2 * Math.PI * d);                  // eqs. 143-144, 175
export const densityFromChi = (chi, lam_c, theta_c, psi, d) =>
  (2 * Math.PI * d * chi) / (lam_c * Math.cos(theta_c - psi));
export const fringeSpacing = (lam_r, n0, theta_r) =>
  lam_r / (2 * n0 * Math.abs(Math.cos(theta_r)));                            // eq. 170 at Bragg
export const kleinCook = (lam_c, d, n0, Lambda, theta_int) =>
  (2 * Math.PI * lam_c * d) / (n0 * Lambda * Lambda * Math.cos(theta_int));  // eq. 171
export const internalAngle = (theta_air, n0) => Math.asin(Math.sin(theta_air) / n0);
export function familyAngles(theta_air, n0, psi) {
  const th = psi + internalAngle(theta_air, n0);
  return { theta_c: th, theta_r: th };
}
