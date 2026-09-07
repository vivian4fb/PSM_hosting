// Main-thread manager for the rigorous-reference worker: one job at a time,
// stale jobs cancelled, results delivered progressively.
import { PYODIDE, RCWA } from '../config.js';

const PY_FILES = ['__init__.py', 'conventions.py', 'twowave.py', 'rcwa.py', 'chain.py', 'metrics.py', 'presets.py'];

export class RcwaClient {
  constructor(handlers) {
    this.h = handlers;              // {onStatus, onProgress, onDone, onError}
    this.worker = null;
    this.ready = false;
    this.nextId = 1;
    this.active = null;
    this.timer = null;
  }

  start(base) {
    if (this.worker) return;
    this.worker = new Worker(new URL('./worker.js', import.meta.url));
    this.worker.onmessage = (ev) => this._onMessage(ev.data);
    this.worker.onerror = (e) => this.h.onError(`worker: ${e.message || e}`);
    this.h.onStatus('loading the rigorous solver (Pyodide + numpy, ~12 MB once)…');
    this.worker.postMessage({ type: 'init', indexURL: PYODIDE.indexURL, base, files: PY_FILES });
  }

  _onMessage(m) {
    if (m.type === 'ready') {
      this.ready = true;
      this.h.onStatus(`rigorous solver ready (Pyodide ${m.version})`);
      if (this.pending) { const p = this.pending; this.pending = null; this.request(p.params, p.pol, p.lam, p.L, p.tag); }
      return;
    }
    if (m.type === 'error') { this.h.onError(m.message); return; }
    if (m.id !== (this.active && this.active.id)) return;     // stale
    if (m.type === 'progress') this.h.onProgress(this.active, m.done, m.total, m.eta);
    else if (m.type === 'done') { const a = this.active; this.active = null; this.h.onDone(a, m.eta, m.ms); }
  }

  // Debounced request; returns immediately.
  schedule(params, pol, lam, L, tag) {
    clearTimeout(this.timer);
    this.timer = setTimeout(() => this.request(params, pol, lam, L, tag), RCWA.debounceMs);
  }

  request(params, pol, lam, L, tag) {
    if (!this.ready) { this.pending = { params, pol, lam, L, tag }; return; }
    this.cancel();
    const id = this.nextId++;
    this.active = { id, params, pol, lam, L, tag };
    this.worker.postMessage({ type: 'sweep', id, params, pol, lam, L, chunk: RCWA.chunk });
  }

  cancel() {
    if (this.active) { this.worker.postMessage({ type: 'cancel', id: this.active.id }); this.active = null; }
  }
}
