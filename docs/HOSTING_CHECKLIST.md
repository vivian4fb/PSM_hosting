# Hosting

The explorer is served by **GitHub Pages** from the public repository
`https://github.com/vivian4fb/PSM_hosting` (branch `main`, root, Jekyll disabled by `.nojekyll`):

    https://vivian4fb.github.io/PSM_hosting/

Decision 2026-09-07: published with PSM* included, before the preprint, so that David can review it.

## Updating the site

Edit, run the tests, commit, push. Pages rebuilds in about a minute.

    cd C:\Users\vivia\Repos\psm-explorer
    cd python && python -m pytest -q && cd ..
    node --test tests/js/models.test.js
    git add -A && git commit -m "..." && git push

## If the site ever has to go private

Make the repository private (Settings → General → Danger zone). GitHub Pages stops serving it on a
free account; copies already downloaded cannot be recalled. A private alternative is Cloudflare
Pages connected to the private repo (free tier, build command empty, output directory `/`).

## Local check

    python -m http.server 8080      (repository root) → http://localhost:8080/
    node tests/js/smoke.mjs http://127.0.0.1:8080/ out.png
