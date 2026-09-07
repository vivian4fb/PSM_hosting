// Light / dark theme: prefers-color-scheme by default, data-theme override remembered per browser.
export function currentTheme() {
  const forced = document.documentElement.getAttribute('data-theme');
  if (forced) return forced;
  return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
}

export function initTheme(onChange) {
  let saved = null;
  try { saved = localStorage.getItem('vge-theme'); } catch (e) { saved = null; }
  if (saved === 'dark' || saved === 'light') document.documentElement.setAttribute('data-theme', saved);
  const btn = document.getElementById('theme-toggle');
  const label = () => { btn.textContent = currentTheme() === 'dark' ? 'Light mode' : 'Dark mode'; };
  label();
  btn.addEventListener('click', () => {
    const next = currentTheme() === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', next);
    try { localStorage.setItem('vge-theme', next); } catch (e) { /* ignore */ }
    label();
    onChange(next);
  });
  if (window.matchMedia) {
    window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => {
      if (!document.documentElement.getAttribute('data-theme')) { label(); onChange(currentTheme()); }
    });
  }
}
