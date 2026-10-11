# Deploying LOSAT Web to Cloudflare Pages (for the Owner)

The Owner deploys; no agent logs in to Cloudflare, runs `wrangler` against an account or sends
anything to a service. This document is the procedure, and what to run on a preview before the
site goes to production. It follows the layout in [`release_manifest.md`](release_manifest.md)
(site, releases, site index, `_headers`) and the plan's §5.9 (offline use and explicit update) and
§5.10 (GA4). Cloudflare facts below were read from Cloudflare's documentation on 2026-10-11 and
are cited with their pages at the end; anything not found there is marked "not verified".

## 1. Choices

- **Direct upload, not Pages' Git build.** The site needs the two engine modules (the Rust/Wasm
  reactors, `web/adapter/tools/build_reactors.py`: Rust 1.92.0 with the targets `wasm32-wasip1`
  and `wasm32-wasip1-threads`). Pages' build image (v3) lists Node.js, Python, Ruby, Go and Bun
  and no Rust, and a build times out after 20 minutes. Installing Rust inside a Pages build is
  possible in principle but was not tried, and the engine builds already run in
  `.github/workflows/web.yml`. So the site is built on the Owner's machine (or in CI) and the
  finished directory is uploaded. Cloudflare says that a Direct Upload project cannot be switched
  to Git integration later.
- **Wrangler.** `npx wrangler@4.149.0 ...` (the latest on the npm registry on 2026-10-11; needs
  Node >= 22). Pin the version in every command and write it in the deployment notes. The app
  itself needs Node >= 22.18 (`web/app/package.json`); CI uses Node 24, use the same.
- **Projects and branches.** One Pages project (name the Owner's), created with a production
  branch. A deployment to the production branch is production (`<project>.pages.dev` and custom
  domains); a deployment to any other branch is a preview, at
  `https://<branch>.<project>.pages.dev` (branch name with non-alphanumeric characters turned into
  `-`; the alias follows the latest deployment of that branch) and at a per-deployment URL
  `https://<hash>.<project>.pages.dev`. Previews are public by default; Cloudflare Access can
  restrict them (project settings; this does not cover `*.pages.dev` production). Every preview
  answer carries `X-Robots-Tag: noindex`.
- **No `wrangler.toml`.** `web/app` has none, so the dashboard settings (variables, secrets,
  compatibility) are the configuration. Do not add one without deciding that it becomes the
  source of truth (Cloudflare then reads variables from the file).

## 2. One-time setup (Owner)

```bash
npx wrangler@4.149.0 login        # or CLOUDFLARE_API_TOKEN and CLOUDFLARE_ACCOUNT_ID in the environment
npx wrangler@4.149.0 pages project create <project> --production-branch main
```

For a Direct Upload project the production branch cannot be changed in the dashboard (only through
the API's Update Project, `production_branch`), so choose it now. Then, in the dashboard (Workers &
Pages > the project > Settings):

- **Variables and Secrets** (`web/app/functions/api/visit.ts`, the GA4 Function of
  [`analytics.md`](analytics.md)): `GA4_MEASUREMENT_ID` and `GA4_API_SECRET`, the secret with
  "Encrypt", **in the Production environment only**. The Function reads them as `context.env`;
  without them it still answers 204 and forwards nothing, so a preview never measures. Check the
  environment column after saving. (`wrangler pages secret put <KEY> --project-name <project>`
  exists too; which environment it writes was not verified, so use the dashboard.)
- Do **not** turn on Web Analytics for the project: Cloudflare then injects a beacon script into
  HTML, which loads from another origin and which the Content-Security-Policy refuses
  (web/AGENTS.md: no script from another origin in a document that handles research data). The
  preview check fails if it finds one.
- Optionally, Access for the preview environment.

## 3. Build a release

From a clean checkout of the commit to release: no uncommitted or untracked file in `web/app` or in
the other files the build reads (`web/app/build/inputs.ts`); the release ID ends in `-dirty`
otherwise, and a `-dirty` release must not be deployed. A tree without git cannot be built
without `LOSAT_WEB_RELEASE_ID`.

```bash
# 1. The reactors, as .github/workflows/web.yml's adapter job builds them (Rust 1.92.0, both wasm targets).
RUSTUP_TOOLCHAIN=1.92.0 python3 web/adapter/tools/build_reactors.py \
    --target-dir <scratch>/adapter-target --output-dir <scratch>/reactors

# 2. The site.
cd web/app
npm ci
git status --porcelain                      # must print nothing (the release ID checks build/inputs.ts's paths)
LOSAT_WEB_REACTORS=<scratch>/reactors \
LOSAT_WEB_ANALYTICS_ENDPOINT=/api/visit \
    npm run build                           # writes dist/
python3 -c "import json; print(json.load(open('dist/releases.json'))['latest'])"   # must not end in -dirty or contain -fake
```

- `LOSAT_WEB_REACTORS` is required: without it the build is a FakeEngine site (release ID
  `<commit>-fake`), which must not be deployed.
- Build with `NODE_ENV` unset (or `production`): any other value makes a development build with
  the build machine's paths in the bundle, so `vite.config.ts` refuses it (`unset NODE_ENV` first
  in a shell or a tool that exports it).
- `LOSAT_WEB_ANALYTICS_ENDPOINT=/api/visit` turns measurement on in the site index
  (`/releases.json`). It is a property of the site, not of the release. Set it for the preview and
  the production site alike (S16 decision 9: the same `dist/` is promoted, below); the preview has
  no GA4 secrets, so its Function forwards nothing, and the preview check can exercise
  `/api/visit`. [`analytics.md`](analytics.md) says the same. Leave it out to deploy a site that
  measures nothing.
- Release IDs: `<commit 12>-<engine digest 8>`, the first 12 hex digits of the commit in every
  clone (a CI build and a local build of one commit get the same ID and the same bytes).
  `LOSAT_WEB_RELEASE_ID=<name>` names one by hand
  (pattern in `release_manifest.md`); a release is never rebuilt under an ID that is already
  served with other bytes (`assemble` refuses it).

### Keep the older releases

A tab that runs an older release keeps loading that release's files, so a new deployment must
still serve them. Keep, outside the repository, the directory of every site that was uploaded
(for example `~/losat-web-sites/<release-id>/`), and assemble the new site on top of the last
production one:

```bash
npm run assemble -- --site dist --previous ~/losat-web-sites/<last production release-id>
cp -a dist ~/losat-web-sites/<new release-id>      # keep what you upload
```

`assemble` checks every older release against its manifest, adds it to `dist/r/`, and rewrites
`dist/releases.json` (latest = the new build). **How many releases to keep is the Owner's decision
(pending).** The command keeps every release it is given; to drop old ones, leave them out of
`--previous` (a tab pinned to a dropped release is told that the release is no longer served and
offered the latest). Cloudflare's limits that bear on this: 20,000 files per deployment on the
free plan (100,000 on paid plans) and 25 MiB per file. A release with the engine is about 15 files
and 8 MB (the larger module is 3.9 MB), so neither limit is near; the Pages limits page lists no
total size.

Check the result locally before uploading:

```bash
node build/site-server.ts dist --port 4665 &     # serves dist/ as Pages would (no Functions)
LOSAT_WEB_PREVIEW_URL=http://localhost:4665 LOSAT_WEB_PREVIEW_HEADERS=dist/_headers \
LOSAT_WEB_PREVIEW_OUT=$HOME/losat-preview-check/local \
  npx playwright test -c ../../docs/evidence/losat_web_w7/check_preview/playwright.config.ts
```

## 4. Deploy the preview, check it, promote it

```bash
cd web/app            # so that functions/ next to dist/ is uploaded (see below)
npx wrangler@4.149.0 pages deploy dist --project-name <project> --branch preview \
    --commit-hash "$(git rev-parse HEAD)" --commit-message "LOSAT Web <release-id>"
```

- Run it in `web/app`: Wrangler uploads a `functions/` folder that exists where the command is
  run (the dashboard's drag-and-drop upload does not compile Functions). Wrangler
  generates a `_routes.json` for a project with `functions/` (its content was not verified: look
  at the deploy output; Cloudflare says that without a narrower rule every request invokes the
  Functions, which would count static files against the free plan's 100,000 requests a day, while
  requests to static assets that do not invoke a Function are free and unlimited).
  The Function imports `../../src/domain/analytics`; Wrangler is expected to bundle it, which was
  not tried here (if the deploy cannot build the Function, it must be made self-contained).
- `dist/_headers` is read from the root of the uploaded directory and is not served.
- The commit message is stored with the deployment: no secrets in it.
- Wrangler uploads at most 20,000 files and 25 MiB per file.

Then run the preview check (section 5) against `https://preview.<project>.pages.dev` (the alias),
or against the per-deployment URL that Wrangler prints. If it passes and the Owner has looked at
the preview by hand, promote **the same directory** (no rebuild: the bytes the check saw are the
bytes of production):

```bash
npx wrangler@4.149.0 pages deploy dist --project-name <project> --branch main
```

and run the check once more against `https://<project>.pages.dev`. Production's secrets come into
play only now; do not run the check with `LOSAT_WEB_PREVIEW_SEND_VISIT=1` against production
unless one real event in GA4 is acceptable.

## 5. The preview check

`docs/evidence/losat_web_w7/check_preview/` ([README](../evidence/losat_web_w7/check_preview/README.md)):

```bash
cd web/app
LOSAT_WEB_PREVIEW_URL=https://preview.<project>.pages.dev \
LOSAT_WEB_PREVIEW_HEADERS=dist/_headers \
LOSAT_WEB_PREVIEW_OUT=$HOME/losat-preview-check/<release-id> \
npx playwright test -c ../../docs/evidence/losat_web_w7/check_preview/playwright.config.ts
```

(Add `LOSAT_WEB_WEBKIT_EXECUTABLE=<launcher>` on a host where Playwright's WebKit lacks libraries.)
It checks, in Chromium, Firefox and WebKit, that the app starts through the boot page with
`crossOriginIsolated === true` on the latest release; for each kind of response (site files,
manifest, entry and worker scripts, both engine modules, a missing file) the headers of the 200 or
404 and of conditional requests that should give 304 (COOP, COEP, CORP and the CSP must be on the
304 too: WebKit refuses a re-created Engine worker otherwise; the `/r/*` immutable and site-file
`no-cache` `Cache-Control`s on the 200s); that every file is served byte for byte as the manifest
lists it and the boot page has no foreign or injected script; a search cancelled while running
followed by a search that completes (engine builds); "Offline ready" and, in Chromium, an offline
reload; and `/api/visit` (400 for a refused body; a valid event only with
`LOSAT_WEB_PREVIEW_SEND_VISIT=1`).

**Send back to the session**: the preview URL, the release ID, the branch, the Wrangler version,
the date, and `summary.json`, `http.json`, `chromium.json`, `firefox.json`, `webkit.json` from
`LOSAT_WEB_PREVIEW_OUT` (they hold only URLs, statuses and response headers).

## 6. What to expect from Pages (and what the check settles)

- **`_headers`** ([`public/_headers`](../../web/app/public/_headers)): at most 100 rules (the file
  has 8) and 2,000 characters a line (the CSP line is about 200). Every matching rule applies and a
  header set by two rules is joined with `, `; the file avoids that (`Cache-Control` appears only
  in the specific rules). `_headers` **does not apply to responses of Pages Functions**: the
  `/api/visit` answers carry the headers the Function sets itself (`Cache-Control: no-store`,
  CORP, `nosniff`) and no COOP/COEP; that is fine for a same-origin `fetch` from the isolated
  page. Cloudflare's page does not say whether the rules reach `304` and `404` answers; the check
  tests it, and a failure there is the finding to send back (the local server and Vite servers
  always add the headers).
- **Top-level `404.html`.** Without it Pages assumes a single-page application and answers every
  unknown path with `/` and status 200; the build writes `404.html` (`src/boot/not-found.html`),
  so a missing file is 404. Pages looks for the nearest `404.html` up the tree. Pages redirects
  `/404.html` itself to `/404`, which is why it is not in the shell.
- **`.html` redirects.** `/name.html` redirects to `/name`, `/dir/index.html` to `/dir/`; nothing
  links to those.
- **Caching.** Without a `Cache-Control` rule Pages sends `public, max-age=0, must-revalidate`;
  `/r/*` is immutable for a year, the site files `no-cache` (site-level files are revalidated at
  every load; `/sw.js` too: the browser also revalidates it itself). Pages sends an `ETag` on every
  200 and answers a matching `If-None-Match` with 304. Cloudflare may compress text and make the
  ETag weak (`W/"..."`); the check sends the ETag as received.
- **A 404 under `/r/` and the immutable header.** The rule matches by path, so a missing file
  under `/r/` may get `Cache-Control: public, max-age=31536000, immutable`, and a browser or CDN
  may then remember "missing" for a year. Release IDs are not reused, so this matters only for a
  release that is removed and later added again, or requested before it is deployed (a tab
  pinned to an ID that a later deployment adds). The check records what Pages does in the
  `notes` of `headers: missing-release-file`; nothing in `_headers` can tell a 404 from a 200,
  so if Pages does it, the answer is to deploy a release before anything can ask for it.
- **Previews and production share `dist/`.** The site-level files carry no environment; what
  differs is the Function's secrets (Production only) and the URL.
- **Not verified**: the exact order and form of Pages' answers for `HEAD` or ranges, whether
  `_headers` reach 304/404 answers (the check will say), and how wrangler's `pages secret put`
  chooses an environment.

## Sources (read 2026-10-11)

- Headers: <https://developers.cloudflare.com/pages/configuration/headers/>
- Serving Pages (404.html, single-page mode, redirects, ETag): <https://developers.cloudflare.com/pages/configuration/serving-pages/>
- Limits: <https://developers.cloudflare.com/pages/platform/limits/>
- Direct Upload: <https://developers.cloudflare.com/pages/get-started/direct-upload/>
- Wrangler `pages` commands: <https://developers.cloudflare.com/workers/wrangler/commands/pages/>
- Preview deployments: <https://developers.cloudflare.com/pages/configuration/preview-deployments/>
- Branch deployment controls: <https://developers.cloudflare.com/pages/configuration/branch-build-controls/>
- Build image (Node.js 22.16.0 default; no Rust listed): <https://developers.cloudflare.com/pages/configuration/build-image/>
- Functions routing and `_routes.json`: <https://developers.cloudflare.com/pages/functions/routing/>
- Functions bindings (variables and secrets): <https://developers.cloudflare.com/pages/functions/bindings/>
- Functions pricing and limits: <https://developers.cloudflare.com/pages/functions/pricing/>
- Wrangler on the npm registry (4.149.0, Node >= 22): <https://registry.npmjs.org/wrangler/latest>
