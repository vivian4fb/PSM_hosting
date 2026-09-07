// Build-time switches and constants for the explorer.

export const FEATURES = {
  // PSM* (paper Section 5) is unpublished until the preprint is public.
  // The private preview for the authors carries it; flip to false for any
  // public deployment made before preprint day.
  psmStar: true,
};

export const PYODIDE = {
  version: '0.28.3',
  indexURL: 'https://cdn.jsdelivr.net/pyodide/v0.28.3/full/',
};

export const RCWA = {
  previewPoints: 301,   // first pass while dragging
  chunk: 25,            // wavelengths per worker message
  debounceMs: 400,
  admissibilityTol: 1e-4,   // paper, head of Section 8: |eta(L) - eta(L_conv)| <= 1e-4
};

// Series identities (fixed order, never cycled). Hue identities follow the paper's figures.
export const SERIES = {
  rcwa:     { label: 'RCWA (rigorous)', light: '#A63D2F', dark: '#CF4F55', width: 3,   dash: [] },
  psm_star: { label: 'PSM*',            light: '#178A5E', dark: '#1BA07F', width: 2,   dash: [] },
  tcwt:     { label: 'TCWT',            light: '#2F4FB8', dark: '#4F73E3', width: 2,   dash: [8, 5] },
  psm:      { label: 'PSM',             light: '#2F4FB8', dark: '#4F73E3', width: 2,   dash: [2, 4] },
  kogelnik: { label: 'Kogelnik',        light: '#D97A1A', dark: '#BD8C10', width: 2,   dash: [8, 4, 2, 4] },
};

export const MODEL_ORDER = ['psm_star', 'tcwt', 'psm', 'kogelnik'];
export const VERSION = '0.1.0';
