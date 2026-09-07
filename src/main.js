// Wiring: state <-> controls <-> models <-> plots <-> rigorous worker.
import { FEATURES, MODEL_ORDER, RCWA as RC } from './config.js';
import { fromHash, toHash, toParams, wavelengths, modelsOf, applyPreset } from './state.js';
import { LABELS, sweep } from './models/twowave.js';
import { Plots } from './ui/plot.js';
import { buildControls, syncControls } from './ui/controls.js';
import { initTheme, currentTheme } from './ui/theme.js';
import { analyse, render } from './ui/verdict.js';
import { csvText, download, pythonSnippet } from './ui/export.js';
import { RcwaClient } from './rcwa/client.js';

const $ = (id) => document.getElementById(id);
const available = MODEL_ORDER.filter((m) => FEATURES.psmStar || m !== 'psm_star');

let state = fromHash(location.hash);
if (!FEATURES.psmStar) state.models = modelsOf(state).filter((m) => m !== 'psm_star').join(',');
let presets = { families: {}, figures: [] };
let lam = wavelengths(state);
let curves = {};
let rcwaCache = { key: null, eta: null, L: null, N: null, preview: false };
let admissibility = null;
let admState = null;

const plots = new Plots($('plot-top'), $('plot-bot'), currentTheme());

function rcKey(s) {
  const p = toParams(s);
  return JSON.stringify([p, s.pol, s.a_nm, s.b_nm, s.N]);
}

function setStatus(text, busy) {
  $('status').textContent = text;
  $('prog').hidden = !busy;
}

function lamNm() {
  return Array.from(lam, (x) => x * 1e9);
}

function interpolate(srcLam, srcEta, dstLam) {
  const out = new Float64Array(dstLam.length);
  let j = 0;
  for (let i = 0; i < dstLam.length; i++) {
    const x = dstLam[i];
    while (j < srcLam.length - 2 && srcLam[j + 1] < x) j++;
    const x0 = srcLam[j];
    const x1 = srcLam[j + 1];
    const t = x1 === x0 ? 0 : (x - x0) / (x1 - x0);
    out[i] = srcEta[j] * (1 - t) + srcEta[j + 1] * t;
  }
  return out;
}

function currentPrinted() {
  const fig = presets.figures.find((f) => f.id === state.preset);
  if (!fig) return null;
  // printed values apply only if the configuration is unchanged
  const s = applyPreset(state, fig, presets.families);
  const same = ['n0', 'lr_nm', 'th_air', 'adv', 'psi', 'd_um', 'n1', 'D0', 'D1', 'pol', 'a_nm', 'b_nm', 'N'].every((k) => Number(s[k]) === Number(state[k]) || s[k] === state[k]);
  return same && rcwaCache.L === 10 && !rcwaCache.preview ? fig.printed : null;
}

function updateVerdict() {
  const models = modelsOf(state).filter((m) => available.includes(m));
  const res = analyse(toParams(state), lam, curves, models, state.pol, { admissibility });
  render($('verdict'), $('metrics'), res, currentPrinted());
}

function redraw() {
  const models = modelsOf(state).filter((m) => available.includes(m));
  plots.update(lamNm(), curves, new Set([...models, 'rcwa']));
  updateVerdict();
}

function recompute(live) {
  lam = wavelengths(state);
  const p = toParams(state);
  const models = modelsOf(state).filter((m) => available.includes(m));
  const t0 = performance.now();
  curves = sweep(models, state.pol, lam, p);
  const dt = performance.now() - t0;
  const k = rcKey(state);
  if (rcwaCache.key === k && rcwaCache.eta) {
    curves.rcwa = rcwaCache.eta;
  } else {
    curves.rcwa = null;
    admissibility = null;
    if (state.rcwa) scheduleRcwa(live);
    else setStatus(`analytic models: ${models.length} × ${lam.length} points in ${dt.toFixed(1)} ms`, false);
  }
  redraw();
  history.replaceState(null, '', toHash(state));
}

function scheduleRcwa(live) {
  const p = toParams(state);
  const key = rcKey(state);
  if (lam.length > RC.previewPoints) {
    const preview = wavelengths({ ...state, N: RC.previewPoints });
    client.schedule(p, state.pol, preview, Math.min(state.L, 4), { stage: 'preview', key });
  } else {
    client.schedule(p, state.pol, lam, state.L, { stage: 'full', key });
  }
  setStatus(live ? 'rigorous reference queued…' : 'rigorous reference queued…', true);
}

const client = new RcwaClient({
  onStatus: (t) => setStatus(t, false),
  onError: (m) => setStatus(`rigorous solver error: ${m}`, false),
  onProgress: (job, done, total) => {
    $('prog').value = Math.round((100 * done) / total);
    setStatus(`${job.tag.stage === 'preview' ? 'preview' : job.tag.stage === 'full' ? 'rigorous' : job.tag.stage} L = ${job.L}: ${done}/${total}`, true);
  },
  onDone: (job, eta, ms) => {
    if (job.tag.stage === 'adm4' || job.tag.stage === 'adm8') { onAdmissibility(job, eta); return; }
    if (job.tag.key !== rcKey(state)) return;                       // stale
    if (job.tag.stage === 'preview') {
      rcwaCache = { key: job.tag.key, eta: interpolate(job.lam, eta, lam), L: job.L, N: job.lam.length, preview: true };
      curves.rcwa = rcwaCache.eta;
      redraw();
      setStatus(`preview L = ${job.L}, ${job.lam.length} points in ${(ms / 1000).toFixed(1)} s — computing full grid at L = ${state.L}…`, true);
      client.request(toParams(state), state.pol, lam, state.L, { stage: 'full', key: job.tag.key });
    } else {
      rcwaCache = { key: job.tag.key, eta, L: job.L, N: job.lam.length, preview: false };
      curves.rcwa = eta;
      redraw();
      setStatus(`rigorous reference: L = ${job.L}, ${job.lam.length} points in ${(ms / 1000).toFixed(1)} s`, false);
    }
  },
});

function onAdmissibility(job, eta) {
  if (!admState || admState.key !== job.tag.key) return;
  if (job.tag.stage === 'adm4') {
    admState.eta4 = eta;
    setStatus('admissibility: L = 8 pass…', true);
    client.request(toParams(state), state.pol, admState.lam, 8, { stage: 'adm8', key: admState.key });
  } else {
    let maxDiff = 0;
    for (let i = 0; i < eta.length; i++) maxDiff = Math.max(maxDiff, Math.abs(eta[i] - admState.eta4[i]));
    admissibility = { maxDiff };
    admState = null;
    setStatus(`admissibility check done: max|η(4) − η(8)| = ${maxDiff.toExponential(2)}`, false);
    updateVerdict();
  }
}

function onChange(key, value, live) {
  state = { ...state, [key]: value };
  if (key !== 'models' && key !== 'rcwa' && key !== 'L') state.preset = state.preset && key === 'preset' ? value : 'custom';
  if (key === 'adv') { const box = $('adv-box'); if (box) box.hidden = !value; }
  if (key === 'L' || key === 'rcwa') { rcwaCache = { key: null }; }
  if (live && lam.length > 601) {
    // cheap redraw while dragging: analytic only on the full grid is still fast
    recompute(true);
  } else {
    recompute(false);
  }
  syncPresetSelect();
}

function syncPresetSelect() {
  const sel = $('preset');
  sel.value = presets.figures.some((f) => f.id === state.preset) ? state.preset : 'custom';
  const fig = presets.figures.find((f) => f.id === state.preset);
  $('preset-caption').textContent = fig ? `${fig.fig}: ${fig.caption}. Printed rms — PSM* ${fig.printed.psm_star}, TCWT ${fig.printed.tcwt}.` : 'Custom configuration — every slider is live.';
}

async function loadPresets() {
  try {
    const r = await fetch(new URL('../presets/figures.json', import.meta.url));
    presets = await r.json();
  } catch (e) {
    presets = { families: {}, figures: [] };
  }
  const sel = $('preset');
  sel.innerHTML = '';
  const custom = document.createElement('option');
  custom.value = 'custom';
  custom.textContent = 'custom';
  sel.append(custom);
  for (const f of presets.figures) {
    const o = document.createElement('option');
    o.value = f.id;
    o.textContent = `${f.fig} — ${f.caption}`;
    sel.append(o);
  }
  sel.addEventListener('change', () => {
    const fig = presets.figures.find((f) => f.id === sel.value);
    if (!fig) return;
    state = applyPreset(state, fig, presets.families);
    if (!FEATURES.psmStar) state.models = modelsOf(state).filter((m) => m !== 'psm_star').join(',');
    rcwaCache = { key: null };
    syncControls(state);
    recompute(false);
    syncPresetSelect();
  });
  if (!location.hash) {
    const fig = presets.figures.find((f) => f.id === state.preset);
    if (fig) state = applyPreset(state, fig, presets.families);
  }
}

function wireButtons() {
  $('btn-csv').addEventListener('click', () => {
    const models = modelsOf(state).filter((m) => available.includes(m));
    const geom = analyse(toParams(state), lam, curves, models, state.pol, {}).geom;
    download('volume-grating-explorer.csv', new Blob([csvText(state, lamNm(), curves, models, geom)], { type: 'text/csv' }));
  });
  $('btn-png').addEventListener('click', async () => {
    const url = plots.toPNG();
    const blob = await (await fetch(url)).blob();
    download('volume-grating-explorer.png', blob);
  });
  $('btn-link').addEventListener('click', async () => {
    try { await navigator.clipboard.writeText(location.href); setStatus('link copied', false); }
    catch (e) { setStatus(location.href, false); }
  });
  $('btn-py').addEventListener('click', async () => {
    const models = modelsOf(state).filter((m) => available.includes(m));
    const txt = pythonSnippet(state, toParams(state), models, state.pol);
    try { await navigator.clipboard.writeText(txt); setStatus('Python snippet copied', false); }
    catch (e) { setStatus('clipboard unavailable', false); }
  });
  $('btn-verify').addEventListener('click', () => {
    state = { ...state, rcwa: 1, L: 10 };
    rcwaCache = { key: null };
    syncControls(state);
    client.cancel();
    client.request(toParams(state), state.pol, lam, 10, { stage: 'full', key: rcKey(state) });
    setStatus('rigorous reference at L = 10 on the full grid…', true);
    history.replaceState(null, '', toHash(state));
  });
  $('btn-adm').addEventListener('click', () => {
    const grid = wavelengths({ ...state, N: RC.previewPoints });
    admState = { key: rcKey(state), lam: grid, eta4: null };
    client.cancel();
    client.request(toParams(state), state.pol, grid, 4, { stage: 'adm4', key: admState.key });
    setStatus('admissibility: L = 4 pass…', true);
  });
}

async function main() {
  initTheme((theme) => { plots.build(theme); redraw(); });
  await loadPresets();
  buildControls($('ctl-root'), state, onChange, { models: available, labels: LABELS });
  plots.build(currentTheme());
  window.addEventListener('resize', () => plots.resize());
  wireButtons();
  syncPresetSelect();
  const base = new URL('../', import.meta.url).href;
  client.start(base);
  recompute(false);
}

main();
