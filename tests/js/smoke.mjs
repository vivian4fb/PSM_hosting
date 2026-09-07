// End-to-end smoke test with the local Chrome: loads the explorer, waits for the
// rigorous solver (Pyodide) to finish at L = 10 on the Fig. 17 left preset, and
// checks that the printed paper values are reproduced. Not part of `node --test`.
//   node tests/js/smoke.mjs [url] [screenshot.png]
import puppeteer from 'puppeteer-core';
import { existsSync } from 'node:fs';

const url = process.argv[2] || 'http://127.0.0.1:8080/';
const shot = process.argv[3] || 'smoke.png';
const candidates = [
  'C:/Program Files/Google/Chrome/Application/chrome.exe',
  'C:/Program Files (x86)/Google/Chrome/Application/chrome.exe',
  'C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe',
];
const exe = candidates.find((c) => existsSync(c));
if (!exe) { console.error('no Chrome/Edge found'); process.exit(2); }

const browser = await puppeteer.launch({ executablePath: exe, headless: true, args: ['--no-sandbox', '--disable-gpu'] });
const page = await browser.newPage();
await page.setViewport({ width: 1400, height: 1900 });
const errors = [];
page.on('pageerror', (e) => errors.push(`pageerror: ${e.message}`));
page.on('console', (m) => { if (m.type() === 'error') errors.push(`console.error: ${m.text()}`); });

const t0 = Date.now();
await page.goto(url, { waitUntil: 'load' });
await page.waitForFunction(() => /rigorous solver ready|rigorous reference:|solver error|preview L/.test(document.getElementById('status').textContent), { timeout: 180000 });
console.log('status after load:', await page.$eval('#status', (e) => e.textContent), `(${((Date.now() - t0) / 1000).toFixed(1)} s)`);
await page.waitForFunction(() => /rigorous reference: L = 10|solver error/.test(document.getElementById('status').textContent), { timeout: 300000 });
console.log('status:', await page.$eval('#status', (e) => e.textContent), `(${((Date.now() - t0) / 1000).toFixed(1)} s)`);
const rows = await page.$$eval('#metrics tbody tr', (trs) => trs.map((tr) => Array.from(tr.children).map((td) => td.textContent)));
console.table(rows);
await page.screenshot({ path: shot, fullPage: true });
if (errors.length) { console.log('ERRORS:'); for (const e of errors) console.log(' ', e); }
const ok = rows.some((r) => r[0] === 'PSM*' && r[4].includes('✓')) && rows.some((r) => r[0] === 'TCWT' && r[4].includes('✓'));
console.log(ok ? 'SMOKE PASS: printed values reproduced in the browser' : 'SMOKE FAIL');
await browser.close();
process.exit(ok && errors.length === 0 ? 0 : 1);
