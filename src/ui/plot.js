// Two synced uPlot panels: efficiency vs replay wavelength, and |Δη| against the rigorous curve.
import { SERIES } from '../config.js';

const FLOOR = 1e-5;   // log floor of the difference strip, as in the paper's Fig. 5

function palette(theme, key) {
  return theme === 'dark' ? SERIES[key].dark : SERIES[key].light;
}

function axisColor(theme) {
  return theme === 'dark' ? '#a0a6b0' : '#5c5c57';
}

function gridColor(theme) {
  return theme === 'dark' ? 'rgba(255,255,255,0.08)' : 'rgba(0,0,0,0.08)';
}

export class Plots {
  constructor(topEl, botEl, theme) {
    this.topEl = topEl;
    this.botEl = botEl;
    this.theme = theme;
    this.keys = ['rcwa', 'psm_star', 'tcwt', 'psm', 'kogelnik'];
    this.top = null;
    this.bot = null;
    this.sync = uPlot.sync('vge');
  }

  _series(theme, isDiff) {
    const out = [{ label: 'λ (nm)' }];
    for (const k of this.keys) {
      if (isDiff && k === 'rcwa') continue;
      const s = SERIES[k];
      out.push({
        label: isDiff ? `|Δη| ${s.label}` : s.label,
        stroke: palette(theme, k),
        width: isDiff ? 1.5 : s.width,
        dash: s.dash.length ? s.dash : undefined,
        spanGaps: false,
        value: (u, v) => (v == null ? '—' : isDiff ? v.toExponential(2) : v.toFixed(4)),
      });
    }
    return out;
  }

  _opts(theme, isDiff, width) {
    const ax = axisColor(theme);
    const gr = gridColor(theme);
    return {
      width,
      height: isDiff ? 170 : 340,
      cursor: { sync: { key: 'vge' }, drag: { x: true, y: false } },
      legend: { show: !isDiff },
      scales: isDiff ? { x: { time: false }, y: { distr: 3, range: (u, min, max) => [FLOOR, Math.max(1, max || 1)] } }
                     : { x: { time: false }, y: { range: (u, min, max) => [0, Math.min(1.0, Math.max(1e-3, (max || 0) * 1.08))] } },
      axes: [
        { stroke: ax, grid: { stroke: gr }, ticks: { stroke: gr }, label: isDiff ? 'replay wavelength λc (nm)' : '' },
        { stroke: ax, grid: { stroke: gr }, ticks: { stroke: gr },
          label: isDiff ? '|Δη| vs RCWA' : 'diffraction efficiency η', size: 64,
          values: isDiff ? (u, vals) => vals.map((v) => (v == null ? '' : v.toExponential(0))) : undefined },
      ],
      series: this._series(theme, isDiff),
    };
  }

  build(theme) {
    this.theme = theme;
    if (this.top) this.top.destroy();
    if (this.bot) this.bot.destroy();
    const w = Math.max(320, this.topEl.clientWidth);
    this.top = new uPlot(this._opts(theme, false, w), this.data || [[], [], [], [], [], []], this.topEl);
    this.bot = new uPlot(this._opts(theme, true, w), this.diff || [[], [], [], [], []], this.botEl);
    this.sync.sub(this.top);
    this.sync.sub(this.bot);
  }

  // lamNm: Float64Array; curves: {rcwa?: Float64Array|null, psm_star..: Float64Array}
  update(lamNm, curves, shown) {
    const x = Array.from(lamNm);
    const data = [x];
    for (const k of this.keys) {
      const c = curves[k];
      data.push(c && shown.has(k) ? Array.from(c).map((v) => (Number.isFinite(v) ? v : null)) : x.map(() => null));
    }
    const diff = [x];
    const ref = curves.rcwa;
    for (const k of this.keys) {
      if (k === 'rcwa') continue;
      const c = curves[k];
      if (!ref || !c || !shown.has(k)) { diff.push(x.map(() => null)); continue; }
      diff.push(Array.from(c).map((v, i) => (Number.isFinite(v) && Number.isFinite(ref[i]) ? Math.max(FLOOR, Math.abs(v - ref[i])) : null)));
    }
    this.data = data;
    this.diff = diff;
    if (!this.top) this.build(this.theme);
    this.top.setData(data);
    this.bot.setData(diff);
  }

  resize() {
    const w = Math.max(320, this.topEl.clientWidth);
    if (this.top) this.top.setSize({ width: w, height: 340 });
    if (this.bot) this.bot.setSize({ width: w, height: 170 });
  }

  // Compose both canvases into one PNG (2x device pixels already applied by uPlot)
  toPNG() {
    const a = this.top.ctx.canvas;
    const b = this.bot.ctx.canvas;
    const c = document.createElement('canvas');
    c.width = Math.max(a.width, b.width);
    c.height = a.height + b.height;
    const ctx = c.getContext('2d');
    ctx.fillStyle = this.theme === 'dark' ? '#1a1b1e' : '#ffffff';
    ctx.fillRect(0, 0, c.width, c.height);
    ctx.drawImage(a, 0, 0);
    ctx.drawImage(b, 0, a.height);
    return c.toDataURL('image/png');
  }
}
