// Builds the control panel from a spec and keeps it in sync with the state.

const SPEC = [
  { group: 'Plate', items: [
    { key: 'n0', label: 'n₀ mean index', min: 1.3, max: 2.0, step: 0.001 },
    { key: 'd_um', label: 'd thickness (µm)', min: 1, max: 60, step: 0.1 },
    { key: 'lr_nm', label: 'λr recording wavelength (nm)', min: 350, max: 800, step: 1 },
    { key: 'th_air', label: 'reference beam in air (°)', min: 0, max: 85, step: 0.01 },
  ] },
  { group: 'Grating', items: [
    { key: 'psi', label: 'Ψ slant (°)', min: -40, max: 120, step: 0.5, hint: '0 = unslanted reflection, 90 = unslanted transmission' },
    { key: 'n1', label: 'n₁ index modulation', min: 0, max: 0.1, step: 0.001 },
    { key: 'D0', label: 'D₀ mean optical density', min: 0, max: 2.5, step: 0.01 },
    { key: 'D1', label: 'D₁ density modulation', min: 0, max: 2.5, step: 0.01 },
  ] },
  { group: 'Replay', items: [
    { key: 'pol', label: 'polarisation', radio: [['sigma', 'σ (TE)'], ['pi', 'π (TM)']] },
    { key: 'a_nm', label: 'band start (nm)', min: 300, max: 900, step: 1 },
    { key: 'b_nm', label: 'band end (nm)', min: 300, max: 900, step: 1 },
    { key: 'N', label: 'points', select: [301, 601, 1201] },
  ] },
];

const ADVANCED = [
  { key: 'tc', label: 'θc replay angle from fringe normal (°)', min: -60, max: 150, step: 0.01 },
  { key: 'tr', label: 'θr recording angle from fringe normal (°)', min: -60, max: 150, step: 0.01 },
];

function numericRow(item, state, onChange) {
  const row = document.createElement('div');
  row.className = 'ctl';
  const lab = document.createElement('label');
  lab.htmlFor = `in-${item.key}`;
  lab.textContent = item.label;
  const range = document.createElement('input');
  range.type = 'range';
  range.id = `rg-${item.key}`;
  range.min = item.min; range.max = item.max; range.step = item.step;
  range.value = state[item.key];
  range.setAttribute('aria-label', item.label);
  const num = document.createElement('input');
  num.type = 'number';
  num.id = `in-${item.key}`;
  num.min = item.min; num.max = item.max; num.step = item.step;
  num.value = state[item.key];
  range.addEventListener('input', () => { num.value = range.value; onChange(item.key, Number(range.value), true); });
  range.addEventListener('change', () => onChange(item.key, Number(range.value), false));
  num.addEventListener('change', () => {
    const v = Number(num.value);
    if (Number.isFinite(v)) { range.value = v; onChange(item.key, v, false); }
  });
  row.append(lab, range, num);
  if (item.hint) {
    const h = document.createElement('div');
    h.className = 'hint';
    h.textContent = item.hint;
    row.append(h);
  }
  return row;
}

function radioRow(item, state, onChange) {
  const row = document.createElement('div');
  row.className = 'ctl';
  const lab = document.createElement('span');
  lab.textContent = item.label;
  row.append(lab);
  const wrap = document.createElement('div');
  wrap.className = 'radios';
  for (const [val, text] of item.radio) {
    const l = document.createElement('label');
    const r = document.createElement('input');
    r.type = 'radio'; r.name = `rd-${item.key}`; r.value = val; r.checked = state[item.key] === val;
    r.addEventListener('change', () => { if (r.checked) onChange(item.key, val, false); });
    l.append(r, document.createTextNode(' ' + text));
    wrap.append(l);
  }
  row.append(wrap);
  return row;
}

function selectRow(item, state, onChange) {
  const row = document.createElement('div');
  row.className = 'ctl';
  const lab = document.createElement('label');
  lab.htmlFor = `in-${item.key}`;
  lab.textContent = item.label;
  const sel = document.createElement('select');
  sel.id = `in-${item.key}`;
  for (const v of item.select) {
    const o = document.createElement('option');
    o.value = v; o.textContent = v; o.selected = Number(state[item.key]) === Number(v);
    sel.append(o);
  }
  sel.addEventListener('change', () => onChange(item.key, Number(sel.value), false));
  row.append(lab, sel);
  return row;
}

export function buildControls(root, state, onChange, opts) {
  root.innerHTML = '';
  for (const g of SPEC) {
    const card = document.createElement('section');
    card.className = 'card';
    const h = document.createElement('h3');
    h.textContent = g.group;
    card.append(h);
    for (const item of g.items) {
      if (item.radio) card.append(radioRow(item, state, onChange));
      else if (item.select) card.append(selectRow(item, state, onChange));
      else card.append(numericRow(item, state, onChange));
    }
    if (g.group === 'Plate') {
      const adv = document.createElement('div');
      adv.className = 'ctl';
      const l = document.createElement('label');
      const c = document.createElement('input');
      c.type = 'checkbox'; c.checked = !!state.adv; c.id = 'in-adv';
      c.addEventListener('change', () => onChange('adv', c.checked ? 1 : 0, false));
      l.append(c, document.createTextNode(' advanced: set θc and θr directly (family mode: θc = Ψ + θin, θr = θc)'));
      adv.append(l);
      card.append(adv);
      const advBox = document.createElement('div');
      advBox.id = 'adv-box';
      advBox.hidden = !state.adv;
      for (const item of ADVANCED) advBox.append(numericRow(item, state, onChange));
      card.append(advBox);
    }
    root.append(card);
  }
  // models and reference
  const card = document.createElement('section');
  card.className = 'card';
  const h = document.createElement('h3');
  h.textContent = 'Models';
  card.append(h);
  const chosen = new Set(state.models.split(','));
  const list = document.createElement('div');
  list.className = 'checks';
  for (const m of opts.models) {
    const l = document.createElement('label');
    const c = document.createElement('input');
    c.type = 'checkbox'; c.value = m; c.checked = chosen.has(m);
    c.addEventListener('change', () => {
      const cur = new Set(document.querySelectorAll('.checks input:checked'));
      onChange('models', Array.from(cur).map((x) => x.value).join(','), false);
    });
    l.append(c, document.createTextNode(' ' + opts.labels[m]));
    list.append(l);
  }
  card.append(list);
  const rc = document.createElement('div');
  rc.className = 'ctl';
  const rl = document.createElement('label');
  const rcb = document.createElement('input');
  rcb.type = 'checkbox'; rcb.checked = !!state.rcwa; rcb.id = 'in-rcwa';
  rcb.addEventListener('change', () => onChange('rcwa', rcb.checked ? 1 : 0, false));
  rl.append(rcb, document.createTextNode(' rigorous reference (RCWA, runs in your browser)'));
  rc.append(rl);
  card.append(rc);
  card.append(selectRow({ key: 'L', label: 'Floquet orders |l| ≤ L', select: [3, 4, 6, 8, 10] }, state, onChange));
  root.append(card);
}

export function syncControls(state) {
  for (const k of Object.keys(state)) {
    const rg = document.getElementById(`rg-${k}`);
    const num = document.getElementById(`in-${k}`);
    if (rg) rg.value = state[k];
    if (num && num.tagName === 'INPUT' && num.type === 'number') num.value = state[k];
    if (num && num.tagName === 'SELECT') num.value = state[k];
  }
  for (const r of document.querySelectorAll('input[name="rd-pol"]')) r.checked = r.value === state.pol;
  const adv = document.getElementById('in-adv');
  if (adv) adv.checked = !!state.adv;
  const box = document.getElementById('adv-box');
  if (box) box.hidden = !state.adv;
  const chosen = new Set(state.models.split(','));
  for (const c of document.querySelectorAll('.checks input')) c.checked = chosen.has(c.value);
  const rcb = document.getElementById('in-rcwa');
  if (rcb) rcb.checked = !!state.rcwa;
}
