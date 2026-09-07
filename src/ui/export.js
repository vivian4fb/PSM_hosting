// CSV / PNG export and share link.
import { VERSION } from '../config.js';
import { LABELS } from '../models/twowave.js';

export function csvText(state, lamNm, curves, models, geom) {
  const cols = ['lambda_nm'];
  const keys = [];
  if (curves.rcwa) { cols.push('rcwa'); keys.push('rcwa'); }
  for (const m of models) if (curves[m]) { cols.push(m); keys.push(m); }
  for (const m of models) if (curves[m] && curves.rcwa) cols.push(`delta_${m}`);
  const head = [
    `# Volume Grating Explorer v${VERSION} — independent implementation of Brotherton-Ratcliffe & Amos (2026)`,
    `# geometry=${geom} ${Object.entries(state).map(([k, v]) => `${k}=${v}`).join(' ')}`,
    `# link=${location.href}`,
    cols.join(','),
  ];
  const lines = [];
  for (let i = 0; i < lamNm.length; i++) {
    const row = [lamNm[i].toFixed(4)];
    for (const k of keys) row.push(curves[k][i].toPrecision(10));
    for (const m of models) if (curves[m] && curves.rcwa) row.push((curves[m][i] - curves.rcwa[i]).toPrecision(6));
    lines.push(row.join(','));
  }
  return head.concat(lines).join('\n') + '\n';
}

export function download(name, blob) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = name;
  document.body.append(a);
  a.click();
  setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
}

export function pythonSnippet(state, params, models, pol) {
  const p = params;
  const m = models.map((x) => `"${x}"`).join(', ');
  return [
    '# pip install psmgrating  (or run from python/ in the repository)',
    'import numpy as np',
    'from psmgrating import twowave, rcwa',
    `lam = np.linspace(${state.a_nm}e-9, ${state.b_nm}e-9, ${state.N})`,
    `p = dict(n0=${p.n0}, n1=${p.n1}, d=${p.d.toExponential(6)}, lam_r=${p.lam_r.toExponential(6)},`,
    `         theta_r=${p.theta_r.toPrecision(10)}, theta_c=${p.theta_c.toPrecision(10)}, psi=${p.psi.toPrecision(10)}, D0=${p.D0}, D1=${p.D1})`,
    `eta = twowave.sweep([${m}], "${pol}", lam_c=lam, **p)   # ${models.map((x) => LABELS[x]).join(', ')}`,
    `ref = rcwa.sweep(lam_c=lam, pol="${pol}", L=${state.L}, **p)          # rigorous reference`,
  ].join('\n');
}
