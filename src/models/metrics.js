// The paper's figures of merit, eqs. (172)-(174). Mirrors python/psmgrating/metrics.py.

export function rms(model, ref) {
  let s = 0;
  for (let i = 0; i < ref.length; i++) {
    const e = model[i] - ref[i];
    s += e * e;
  }
  return Math.sqrt(s / ref.length);
}

export function peakError(model, ref) {
  let mm = -Infinity;
  let mr = -Infinity;
  for (let i = 0; i < ref.length; i++) {
    if (model[i] > mm) mm = model[i];
    if (ref[i] > mr) mr = ref[i];
  }
  return Math.abs(mm - mr);
}

// Full width at half maximum about the global maximum, linear interpolation
// between bracketing grid points; NaN if the response does not fall to half.
export function fwhm(lam, eta) {
  let j = 0;
  for (let i = 1; i < eta.length; i++) if (eta[i] > eta[j]) j = i;
  const half = eta[j] / 2;
  let left = NaN;
  for (let i = j; i > 0; i--) {
    if ((eta[i - 1] <= half && half < eta[i]) || (eta[i - 1] < half && half <= eta[i])) {
      const frac = (eta[i] - half) / (eta[i] - eta[i - 1]);
      left = lam[i] - frac * (lam[i] - lam[i - 1]);
      break;
    }
  }
  let right = NaN;
  for (let i = j; i < eta.length - 1; i++) {
    if ((eta[i + 1] <= half && half < eta[i]) || (eta[i + 1] < half && half <= eta[i])) {
      const frac = (eta[i] - half) / (eta[i] - eta[i + 1]);
      right = lam[i] + frac * (lam[i + 1] - lam[i]);
      break;
    }
  }
  return right - left;
}

export function widthError(lam, model, ref) {
  const wr = fwhm(lam, ref);
  const wm = fwhm(lam, model);
  return Math.abs(wm - wr) / wr;
}
