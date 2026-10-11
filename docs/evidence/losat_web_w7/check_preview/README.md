# The preview check (S16 item 7)

What to run after the Owner has deployed LOSAT Web to a Cloudflare Pages preview
(`docs/web/deploy_cloudflare.md`), and what it checks. It runs the same way against a local
site (a dry run, below). It reads from the site and writes only a record on this machine; it
sends no research data and no GA4 event (unless you ask for one, check 5).

## Run it

Needs `web/app` installed (`npm ci`) and Playwright's browsers (`npx playwright install
chromium firefox webkit`; on a host where WebKit's libraries are missing, set
`LOSAT_WEB_WEBKIT_EXECUTABLE` to a launcher, as the E2E tests do). The check has no
dependencies of its own: `node_modules` here is a link to `web/app/node_modules`.

```bash
cd web/app
LOSAT_WEB_PREVIEW_URL=https://<branch>.<project>.pages.dev \
LOSAT_WEB_PREVIEW_HEADERS=dist/_headers \
LOSAT_WEB_PREVIEW_OUT=$HOME/losat-preview-check/<release-id> \
npx playwright test -c ../../docs/evidence/losat_web_w7/check_preview/playwright.config.ts
```

| Variable | Meaning |
|---|---|
| `LOSAT_WEB_PREVIEW_URL` | required; the origin of the site (https, or http for localhost) |
| `LOSAT_WEB_PREVIEW_OUT` | required; the directory for the record (created) |
| `LOSAT_WEB_PREVIEW_HEADERS` | the `_headers` file the site was deployed with; default `web/app/public/_headers`. Use the `dist/_headers` that was uploaded |
| `LOSAT_WEB_PREVIEW_SEND_VISIT` | `1` sends one valid event to `/api/visit` (check 5). Default: not sent |
| `LOSAT_WEB_WEBKIT_EXECUTABLE` | a WebKit launcher (as in `web/app/playwright.config.ts`) |

The check runs one test at a time (a deployed site, and a search needs the machine). The whole
check takes a few minutes on a preview (the "Offline ready" step downloads the release, about
8 MB with the engine). The exit status is 0 only if every test passed; skipped tests say why.

**Dry run** against a local site: `npm run build`, `node build/site-server.ts dist --port 4665`
(an engine build: `LOSAT_WEB_REACTORS=<reactors> npm run build`, port 4765), then
`LOSAT_WEB_PREVIEW_URL=http://localhost:4665` and `LOSAT_WEB_PREVIEW_HEADERS=dist/_headers`.
The local server has no Functions, so check 5 is skipped unless something answers `/api/visit`.

## What it checks

Project `http` (Node's `fetch`, no browser; `http.spec.ts`):

- **Headers of each kind of response** (test `headers: <kind>`): `/`, `/boot.js`, `/boot.css`,
  `/releases.json`, `/sw.js` (skipped if the site has none), the release's `release-manifest.json`,
  its entry script and entry style sheet, its engine worker script, its two engine modules
  (skipped for a FakeEngine site, whose manifest has `engine: null`), a missing file at the site
  root and a missing file under `/r/<release-id>/`. For each 200: status, the headers of the
  `_headers` rules for that path (`headersFor` of `web/app/build/headers.ts`) all equal, the
  content type (`application/wasm` for the modules), an `ETag`. Then the conditional requests:
  `If-None-Match` (with `If-Modified-Since` too when there is a `Last-Modified`), and
  `If-Modified-Since` alone when there is one. They must give 304 with no body, and COOP, COEP,
  CORP and the CSP on the 304 equal the rules'. **A 304 without COEP/CORP is a failure**
  (WebKit refuses an Engine worker re-created after a cancel when the 304 lacks them: W3 gate
  record, decision 12). The other headers of a 304 or of a 404 (`Cache-Control`,
  `X-Content-Type-Options`, `Referrer-Policy`) are recorded as `notes`, not failures. A missing
  file must be 404 (not the boot page with 200: Pages' single-page mode without the top-level
  `404.html`). A 404 under `/r/` that carries the immutable `Cache-Control` is noted.
- **Bytes** (`the files of the latest release are served byte for byte`): every file of the
  release manifest, and the shell files listed in `/releases.json`, have the listed size and
  SHA-256 as served (what the service worker checks before "Offline ready"). A host that
  rewrites a file (minifies, injects a script) fails here whatever the headers say.
- **Boot page** (`the boot page loads nothing from another origin ...`): no `src`/`href` to another
  origin, no Cloudflare beacon (Web Analytics) or Rocket Loader, no inline script.
- **Releases** (`every release of the site index has its manifest`): each release in
  `/releases.json` answers its manifest with 200 and its own ID.
- **`/api/visit`** (only when `/releases.json` has `analytics`): a POST with a body that is not the
  event gives 400, a GET gives 405; the answers' own headers are recorded (`_headers` does not
  apply to a Function). With `LOSAT_WEB_PREVIEW_SEND_VISIT=1`, one valid event
  (release `preview-check`, a random client ID) must give 204. If the project has the GA4
  secrets in this environment, that event reaches GA4: the preview should not have them.

Projects `chromium`, `firefox`, `webkit` (`browser.spec.ts`):

- **Start** (`the app starts through the boot page, isolated, on the latest release`):
  `crossOriginIsolated === true`, `SharedArrayBuffer` exists, the entry script of the latest
  release was requested, no boot message, and the release shown (the `release-id` in the header,
  else the pin in `localStorage`) equals `/releases.json`'s `latest`.
- **Cancel** (`a search cancelled while it runs, then a new search completes`): a BLASTP search of
  a proteome is started, cancelled while running, and a new search completes. Skipped for a
  FakeEngine site (no search long enough to cancel). It uses the helpers of
  `web/app/tests/e2e/support/search.ts` and the proteome `LOSAT/tests/fasta/NZ_CP006932.faa`, so
  it needs the repository checkout.
- **Offline** (`"Offline ready" appears; in Chromium an offline reload starts the app`): skipped if
  the site has no `/sw.js`. Otherwise the header's `offline-status` reaches `data-state="ready"` ("Offline ready", up to 3 minutes); in Chromium only,
  the context goes offline and a reload starts the app with `crossOriginIsolated === true`.
  (Firefox and WebKit cannot go offline reliably under Playwright; the local V-OFF test goes
  offline by stopping the server.)

## The record

In `LOSAT_WEB_PREVIEW_OUT`: `summary.json` (counts per project and one line per test: status, skip
reason, first line of the error), `http.json`, `chromium.json`, `firefox.json`, `webkit.json` (every
test with its status, duration, error and the record it attached: each request's status and all
its headers, the conditional requests, the failures and the notes). `test-results/` holds
Playwright's traces of failed tests. The files hold URLs, statuses and response headers; no
research data.

**What to send back to the session:** the preview URL, the release ID and the branch that was
deployed, the wrangler version, and the five JSON files (`summary.json`, `http.json`,
`chromium.json`, `firefox.json`, `webkit.json`). If a test failed, its `failures` and `notes`
in `http.json` say which header differs.

## Files

`playwright.config.ts` (projects, output), `support.ts` (the site, the record, a fetch helper),
`http.spec.ts`, `browser.spec.ts`, `record-reporter.ts` (writes the JSON files), `globals.d.ts`
and `tsconfig.json` (typecheck: `cd web/app && npx tsc -p
../../docs/evidence/losat_web_w7/check_preview/tsconfig.json`), `package.json` (marks the folder
as ES modules), `node_modules` (a link).
