// Web Worker: the rigorous reference computed by the Python package itself,
// byte-identical to the CI-tested files, running under Pyodide. Chunked so
// that cancel messages are honoured between chunks.
let pyodide = null;
let ready = false;
let current = null;   // id of the job in flight

function post(msg) { self.postMessage(msg); }

async function init(cfg) {
  try {
    importScripts(cfg.indexURL + 'pyodide.js');
    pyodide = await loadPyodide({ indexURL: cfg.indexURL });
    await pyodide.loadPackage('numpy');
    pyodide.FS.mkdirTree('/py/psmgrating');
    for (const f of cfg.files) {
      const r = await fetch(cfg.base + 'python/psmgrating/' + f, { cache: 'no-cache' });
      if (!r.ok) throw new Error(`cannot fetch ${f}: ${r.status}`);
      pyodide.FS.writeFile('/py/psmgrating/' + f, await r.text());
    }
    pyodide.runPython(`
import sys, json
sys.path.insert(0, '/py')
import numpy as np
from psmgrating import rcwa, chain, conventions as cv

def sweep_chunk(params_json, pol, lam_json, L):
    p = json.loads(params_json)
    lam = np.array(json.loads(lam_json), dtype=float)
    if abs(p['psi']) < 1e-9:
        r = chain.solve(lam_c=lam, lam_r=p['lam_r'], n0=p['n0'], n1=p['n1'], d=p['d'],
                        theta_c=p['theta_c'], theta_r=p['theta_r'], pol=pol, D0=p['D0'], D1=p['D1'])
        return [float(x) for x in r['eta_R']]
    out = rcwa.sweep(lam_c=lam, lam_r=p['lam_r'], n0=p['n0'], n1=p['n1'], d=p['d'],
                     theta_c=p['theta_c'], theta_r=p['theta_r'], psi=p['psi'], pol=pol, L=int(L),
                     D0=p['D0'], D1=p['D1'])
    return [float(x) for x in out]
`);
    ready = true;
    console.log('rcwa worker ready, pyodide ' + pyodide.version);
    post({ type: 'ready', version: pyodide.version });
  } catch (err) {
    console.log('rcwa worker init error: ' + err);
    post({ type: 'error', message: String(err && err.message ? err.message : err) });
  }
}

async function sweep(job) {
  if (!ready) { post({ type: 'error', id: job.id, message: 'Pyodide not ready' }); return; }
  current = job.id;
  const fn = pyodide.globals.get('sweep_chunk');
  const paramsJson = JSON.stringify(job.params);
  const lam = job.lam;
  const chunk = job.chunk || 25;
  const eta = new Float64Array(lam.length);
  const t0 = performance.now();
  try {
    for (let s = 0; s < lam.length; s += chunk) {
      if (current !== job.id) { fn.destroy(); return; }         // cancelled
      const part = Array.from(lam.subarray(s, Math.min(lam.length, s + chunk)));
      const res = fn(paramsJson, job.pol, JSON.stringify(part), job.L);
      const arr = res.toJs();
      res.destroy();
      for (let i = 0; i < arr.length; i++) eta[s + i] = arr[i];
      post({ type: 'progress', id: job.id, done: Math.min(lam.length, s + chunk), total: lam.length,
             eta: eta.slice(0, Math.min(lam.length, s + chunk)) });
      await new Promise((r) => setTimeout(r, 0));               // let cancel messages arrive
    }
    fn.destroy();
    post({ type: 'done', id: job.id, eta, ms: performance.now() - t0, L: job.L, N: lam.length });
  } catch (err) {
    post({ type: 'error', id: job.id, message: String(err && err.message ? err.message : err) });
  }
}

self.onmessage = (ev) => {
  const m = ev.data;
  if (m.type === 'init') init(m);
  else if (m.type === 'sweep') sweep(m);
  else if (m.type === 'cancel') { if (current === m.id) current = null; }
};
