# Web page for public sharing

Static one-pager for external audiences (hiring panels, collaborators).

## Contents
- `index.html` — self-contained page with inline CSS, references figures/
- `vercel.json` — deployment config (clean URLs, cache headers)
- `figures/` — 4 selected figures at ~300 KB each

## Deploy to Vercel (fastest path, ~5 minutes)

1. Push the whole project to a GitHub repo (public or private both work)
2. Vercel dashboard → **New Project** → import the repo
3. In the setup screen, set **Root Directory** = `web`
4. Framework preset: **Other** (no build step needed)
5. Deploy → get a `*.vercel.app` URL immediately
6. Add custom domain if desired

## Iterating

Any edit to `index.html` + push to GitHub auto-triggers redeploy on the same URL.

## Local preview

```
cd web
python -m http.server 8000
```

Then open http://localhost:8000
