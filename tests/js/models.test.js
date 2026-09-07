// Parity of the JavaScript models with the Python reference (fixtures made by
// python/tools/make_fixtures.py). Run: node --test tests/js/
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';

import * as Z from '../../src/models/complex.js';
import { MODELS, coefficients, efficiency, sweep } from '../../src/models/twowave.js';
import { rms, fwhm } from '../../src/models/metrics.js';

const load = (name) => JSON.parse(readFileSync(new URL(`../../fixtures/${name}`, import.meta.url)));

function close(a, b, absTol = 1e-10, relTol = 1e-9) {
  return Math.abs(a - b) <= absTol || Math.abs(a - b) <= relTol * Math.max(Math.abs(a), Math.abs(b));
}

test('complex sqrt and exp follow numpy branches', () => {
  const s = Z.sqrt(Z.C(-4, 0));
  assert.ok(close(s.re, 0) && close(s.im, 2));
  const s2 = Z.sqrt(Z.C(-4, -0));
  assert.ok(close(s2.im, -2));
  const e = Z.exp(Z.C(0, Math.PI));
  assert.ok(close(e.re, -1) && close(e.im, 0, 1e-15));
});

test('figure presets: four models on 25 wavelengths each', () => {
  const rows = load('analytic_presets.json');
  assert.equal(rows.length, 12);
  let n = 0;
  for (const row of rows) {
    const eta = sweep(MODELS, row.pol, row.lam, row.params);
    for (const m of MODELS) {
      for (let i = 0; i < row.lam.length; i++) {
        assert.ok(close(eta[m][i], row.eta[m][i]), `${row.id} ${m} @${row.lam[i]}: js ${eta[m][i]} py ${row.eta[m][i]}`);
        n++;
      }
    }
  }
  assert.ok(n >= 1200);
});

test('random configurations: coefficients and efficiencies', () => {
  const rows = load('analytic_random.json');
  assert.equal(rows.length, 60);
  for (const row of rows) {
    const p = row.params;
    for (const pol of ['sigma', 'pi']) {
      for (const m of MODELS) {
        const c = coefficients(m, pol, row.lam_c, p.lam_r, p.n0, p.n1, row.chi0, row.chi1, p.theta_c, p.theta_r, p.psi);
        const ref = row.coeff[`${m}:${pol}`];
        for (const k of ['CR', 'CS', 'a', 'theta', 'kappa']) {
          assert.ok(close(c[k].re, ref[k][0], 1e-12, 1e-11) && close(c[k].im, ref[k][1], 1e-12, 1e-11),
            `${m}:${pol} ${k}: js (${c[k].re},${c[k].im}) py (${ref[k][0]},${ref[k][1]})`);
        }
        const eta = efficiency(m, pol, row.lam_c, p);
        assert.ok(close(eta, row.eta[`${m}:${pol}`]), `${m}:${pol} eta js ${eta} py ${row.eta[`${m}:${pol}`]}`);
      }
    }
  }
});

test('metrics: rms and fwhm behave', () => {
  const lam = [0, 1, 2, 3, 4];
  const eta = [0, 0.5, 1, 0.5, 0];
  assert.ok(close(fwhm(lam, eta), 2));
  assert.ok(close(rms([1, 1], [0, 0]), 1));
});
