// Minimal complex arithmetic. Values are plain {re, im} objects; every
// function returns a new object. Branches (sqrt, exp) follow numpy so the
// JavaScript port reproduces the Python reference bit-for-bit where it matters.

export const C = (re, im = 0) => ({ re, im });
export const toC = (x) => (typeof x === 'number' ? C(x, 0) : x);
export const add = (a, b) => C(a.re + b.re, a.im + b.im);
export const sub = (a, b) => C(a.re - b.re, a.im - b.im);
export const mul = (a, b) => C(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
export const scale = (a, s) => C(a.re * s, a.im * s);
export const neg = (a) => C(-a.re, -a.im);
export const conj = (a) => C(a.re, -a.im);
export const mulI = (a) => C(-a.im, a.re);            // i * a
export const abs2 = (a) => a.re * a.re + a.im * a.im;
export const abs = (a) => Math.hypot(a.re, a.im);

export function div(a, b) {
  const d = b.re * b.re + b.im * b.im;
  return C((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
}

// Principal square root (numpy convention, including the sign of -0 imaginary parts).
export function sqrt(z) {
  const r = Math.hypot(z.re, z.im);
  const re = Math.sqrt(Math.max(0, (r + z.re) / 2));
  let im = Math.sqrt(Math.max(0, (r - z.re) / 2));
  if (z.im < 0 || Object.is(z.im, -0)) im = -im;
  return C(re, im);
}

export function exp(z) {
  const e = Math.exp(z.re);
  return C(e * Math.cos(z.im), e * Math.sin(z.im));
}
