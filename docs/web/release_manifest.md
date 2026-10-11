# LOSAT Web site layout, release manifest (format 1), site index and boot page

A LOSAT Web site serves every release in a directory of its own, `/r/<release-id>/`, and a small
boot page at `/` that starts the release this browser has pinned (plan §5.9, design §14, S16
item 1). A release never changes after it is built, and a browser never switches releases silently:
the pin changes only when the user chooses another release. This document is the contract of the
layout; `web/app/src/domain/release-manifest.ts` and `web/app/src/domain/site-index.ts` check the
two JSON files, `web/app/build/site.ts` and `web/app/build/site-files.ts` write them, and
`web/app/src/boot/` is the boot page.

None of these files holds research data: they describe the application's own files.

## The site

`npm run build` (in `web/app`) writes a site into `dist/` (or `vite build --outDir <dir>`), after
removing the directory (it refuses a directory that is not empty and has no `releases.json`,
`index.html` or `_headers`):

```
dist/
  _headers                       the response headers (copied from public/_headers, their only source)
  404.html                       the page of a missing file (src/boot/not-found.html)
  index.html                     the boot page (src/boot/index.html): no inline script or style
  boot.js                        the boot script (src/boot/boot.ts, one ES module)
  boot.css                       the style sheet of the boot page's messages and of 404.html
  releases.json                  the site index
  sw.js                          the service worker (src/sw/worker.ts; offline use and the explicit update)
  r/<release-id>/                one release: the Vite build with base /r/<release-id>/
    release-manifest.json
    assets/...                   JavaScript, CSS, workers, and every file the code imports
    engine/...                   the two engine modules (a build with LOSAT_WEB_REACTORS)
```

- **Site-level files** (everything outside `r/`) are not versioned: each deployment replaces them.
  The files of `web/app/public/` are copied to the site root, never into a release.
- **A release directory** holds everything the release loads, and nothing of a release lives
  outside it. Files that belong to a release are imported by the code (Vite `?url` or
  `new URL(..., import.meta.url)`), so they get content-hashed names under `assets/` and enter the
  manifest without a list kept by hand. The build writes no `index.html` into a release: the boot
  page loads the entry that the manifest names.
- **The release ID** (`<release-id>`) is `LOSAT_WEB_RELEASE_ID` if set, else
  `<commit 12>-<engine digest 8, or fake>[-dirty]` (`web/app/build/release.ts`; S16 decision 1);
  it matches `^[0-9a-z][0-9a-z.-]{0,63}$`. The application knows its own as `__LOSAT_RELEASE_ID__`.
  - `<commit 12>` is the first 12 hex digits of `git rev-parse HEAD`: a fixed length, never git's
    own abbreviation, whose length depends on the number of objects in the clone.
  - `-dirty` marks a build whose inputs differ from the commit (uncommitted or untracked files).
    The inputs are listed once, in `web/app/build/inputs.ts`: `web/app`, and the files of the
    repository that the build plugins read (`docs/web/verification_cells.tsv`, the certification
    records of the verification table, `LOSAT/tests/wasi_shared_memory.js`). The plugins take
    their paths from that list. A `-dirty` release must not be deployed.
  - A build without a commit (no git) stops unless `LOSAT_WEB_RELEASE_ID` names its release; only
    the development server and Vitest run without one (their ID is `unknown-<engine>`).
- **Determinism.** A release directory is immutable and cached for a year, so its ID must fix its
  bytes. No file holds a time; lists are sorted. Two builds of the same
  tree with the same engine modules write the same bytes (`tests/unit/release-bytes.test.ts` builds
  twice and compares every file).
- **Older releases.** A build writes one release. `build/assemble-site.ts` copies the releases of
  earlier sites into a new one (below), so that tabs that still run them keep their files.

The development server (`npm run dev`) is not a site: it serves `index.html` -> `src/main.ts` at
`/`, without the boot page and without release directories.

## The release manifest (`r/<release-id>/release-manifest.json`)

A JSON object, written with two-space indentation and a final line feed. Every field below is
required. A reader ignores fields it does not know, so a later release may add fields within
format 1; a change that an older reader would misread needs format 2.

| Field | Value |
|---|---|
| `format` | `1` |
| `release` | the release ID; it must be the ID of the directory the manifest is in |
| `app` | `{ "version", "build" }`: `package.json`'s version and the first 8 hex digits of the build's commit (`unknown` without git); the app writes the same `build` into session files (`__LOSAT_APP_BUILD__`) |
| `engine` | `{ "serial", "threads" }`: the engine build names of the two modules (`engineBuildName`, as Run details and session files record them); `null` in a FakeEngine build |
| `entry` | `{ "script", "styles" }`: the module script and the style sheets (in order) that the boot page adds; each path is one of `assets` |
| `assets` | every file of the release directory except the manifest: `{ "path", "size", "sha256" }`, sorted by `path` (code-unit order), each path once |

- `path` is relative to the release directory: segments of letters, digits, `_`, `.`, `~` and `-`
  separated by `/`, none empty and none starting with `.`; at most 512 characters.
- `size` is the file's length in bytes; `sha256` is the lower-case hex SHA-256 of its bytes, as the
  site serves them (the threaded engine module with its shared-memory guard).
- At most 100,000 assets and 100 style sheets; text fields at most 200 characters.

The service worker (WP-U) fetches every listed file, checks its size and SHA-256, and only then
calls the release stored ("Offline ready"). `build/assemble-site.ts` refuses a release directory
that has a file the manifest does not list, lacks a listed file, or has a listed file with other
bytes.

## The site index (`/releases.json`)

Written when the site is built or assembled, never by a release (S16 decision 2). Fields that a
reader does not know are ignored.

| Field | Value |
|---|---|
| `format` | `1` |
| `latest` | the release a first visit pins and that an update offers; one of `releases` |
| `releases` | 1 to 1,000 `{ "id", "version", "build" }`: every release the site serves (`version` and `build` are the manifest's `app`) |
| `analytics` | optional; `{ "endpoint" }`: the measurement endpoint, a same-origin path under `/api/` (`^/api/[a-z0-9][a-z0-9/-]{0,63}$`). Present only when `LOSAT_WEB_ANALYTICS_ENDPOINT` is set when the site is built or assembled; anything else stops the build. Without it the site measures nothing |
| `shell` | optional; the site-level files a browser needs to start offline: `{ "path", "size", "sha256" }` for `/` (the bytes of `index.html`), `/boot.css` and `/boot.js`, sorted by path. Paths are absolute site paths |

`404.html` is not in the shell: a browser that is offline has nothing to gain from it, and Cloudflare
Pages redirects `/404.html` to `/404` (below).

## The boot page (`/`, `/boot.js`)

The boot page is static HTML with a `<div id="boot" hidden>` for its messages and the
`<div id="app">` in which the application mounts. `boot.js` does, in order:

1. Read the pinned release ID from `localStorage['losat-web:release']` (S16 decision 3). A read
   that throws, or a value that is not a release ID, is no pin.
2. Unless `navigator.onLine === false`, fetch `/releases.json` (`cache: 'no-store'`, 5 s timeout).
   A failure is not an error. With a pin, this request runs alongside step 4 and its answer is used
   only if the pinned release cannot start.
3. Choose the release: the pin if there is one; else the index's `latest`, which is then pinned (a
   write that throws is not an error: this visit starts the release anyway). With neither:
   - offline: "LOSAT Web needs a network connection the first time it is opened."
   - online, the index not read: "LOSAT Web could not read this site's list of releases. Reload
     the page to try again."
4. Fetch `/r/<id>/release-manifest.json` (default cache mode, so that the service worker answers it
   when the release is stored; 20 s timeout) and check it (format 1, `release` = `<id>`, entry paths
   among the assets). When it cannot be used:
   - 404 or 410, an answer that is not JSON (a host that answers every path with a page), or a
     manifest that fails the check: "This browser uses LOSAT Web release `<id>`, which this site no
     longer serves." With a site index whose `latest` is another release, a button "Use release
     `<latest>`" writes the pin and reloads; if the write throws, the page says that the browser did
     not let LOSAT Web save the choice. Never a silent switch.
   - no answer while offline (the release is not stored): "This browser uses LOSAT Web release
     `<id>`, which is not stored for offline use. Connect to the network and reload the page."
   - no answer while online, or another status: "LOSAT Web could not load release `<id>` ..." and
     reload.
5. Add a `<link rel="stylesheet">` for each of `entry.styles` and wait for them (so that the
   application never shows unstyled), then a `<script type="module">` for `entry.script`; URLs are
   absolute, under `/r/<id>/`. A style sheet or the script that fails to load gives the
   "could not load" message.

The boot reads no research data and sends nothing but the two GETs of steps 2 and 4. Its messages
are text (`textContent`), in English; test IDs `boot-message` (with `data-kind`: `first-visit-offline`,
`no-index`, `not-served`, `offline-not-stored`, `failed`, `pin-not-saved`) and `boot-use-latest`.
`boot.css` styles only `#boot` and `#boot-noscript`, because the application shares the document.
The pure decisions are in `src/boot/plan.ts` (unit tests `tests/unit/boot.test.ts`).

**Layers** (ESLint, `web/app/eslint.config.js`; the rule's own test is `tests/unit/layers.test.ts`):
`src/boot/**` may import `src/domain/**` and `src/sw/protocol.ts` (the worker's URL and scope) only;
`src/sw/**` (the service worker) may import `src/domain/**` and its own files only; of the app's
layers, only `src/infra/**` may import from `src/sw/`, and only `src/sw/protocol.ts` (the names and
messages that the worker and the pages share). Nothing outside `src/boot/` imports it.

**Build.** The site build (`build/site.ts`) builds `src/boot/boot.ts` with a second Vite build
(`configFile: false`, one ES module, minified) into `boot.js`, and copies `src/boot/index.html`,
`src/boot/boot.css` and `src/boot/not-found.html` (as `404.html`). The release's entry comes from
the `index.html` that Vite builds for it; the build stops if that page has anything the boot page
would not load (an inline script or style, a link other than a style sheet or a module preload).

## Response headers (`web/app/public/_headers`)

`public/_headers` is the only source of the response headers (web/AGENTS.md rule 4). The build
copies it to the site root, where Cloudflare Pages reads it; the Vite servers use its `/*` rule,
and `build/site-server.ts`, the E2E harness server and the tests read the same file.

| Path | Headers |
|---|---|
| `/*` | `Cross-Origin-Opener-Policy: same-origin`, `Cross-Origin-Embedder-Policy: require-corp`, `Cross-Origin-Resource-Policy: same-origin`, the `Content-Security-Policy`, `X-Content-Type-Options: nosniff`, `Referrer-Policy: no-referrer` |
| `/r/*` | `Cache-Control: public, max-age=31536000, immutable`: a release's files never change (S16 item 7) |
| `/`, `/index.html`, `/boot.js`, `/boot.css`, `/releases.json`, `/404.html` | `Cache-Control: no-cache`: the site-level files are revalidated on every load |
| `/sw.js` | `Cache-Control: no-cache`: the service worker (site-level, not in the shell: the browser keeps the installed worker) |

The rules follow Cloudflare Pages (developers.cloudflare.com/pages/configuration/headers/, read
2026-10-11), as `build/headers.ts` (`headersFor(path)`) implements them:

- A rule is a URL pattern on a line of its own, followed by indented `Name: value` lines. At most
  100 rules; a line is at most 2,000 characters. (`build/headers.ts` also skips lines starting
  with `#`; the site's `_headers` has none.)
- A splat `*` (one per pattern) matches any characters, none included; `:name` (a letter, then
  letters, digits and `_`) matches the characters of one path segment (not `/`). Their values may be
  used in header values as `:splat` and `:name`. The request's port, protocol and query are not
  matched. (Pages also takes `https://host/...` patterns; this site does not use them, and
  `build/headers.ts` refuses them.)
- Every rule whose pattern matches applies. A header that more than one matching rule attaches is
  joined with `, `. `! Name` detaches a header that another rule attached.
- Headers are not applied to the responses of Pages Functions; a `_redirects` rule takes priority.

What Pages does besides (developers.cloudflare.com/pages/configuration/serving-pages/, read
2026-10-11), and what it means here:

- **404 and single-page mode.** "If your project does not include a top-level `404.html` file, Pages
  assumes that you are deploying a single-page application" and serves `/` for every path that is
  not a file, with status 200. The site therefore has a top-level `404.html`: a missing release file
  is then a 404, not the boot page. (The boot page also treats an answer that is not JSON as "not
  served", in case a host does answer every path with a page.) Pages looks for the nearest
  `404.html` up the directory tree.
- **`.html` redirects.** Pages redirects `/name.html` to `/name` and `/dir/index.html` to `/dir/`.
  Nothing links to `/index.html` or `/404.html`; their rules above are kept for a direct request.
- **ETag.** "Pages always sends `Etag` headers for `200 OK` responses" and answers a matching
  `If-None-Match` with `304 Not Modified`. Whether the `_headers` rules reach the 304 (WebKit judges
  a revalidated worker script by the headers of the 304; W3 decision 12) is checked on the preview
  (S16 item 7).
- **Default caching.** Without a rule, Pages sends `Cache-Control: public, max-age=0, must-revalidate`.
- The rules match the request's path, so a 404 under `/r/` also carries the immutable
  `Cache-Control` (if Pages applies `_headers` to its 404 answers; to be seen on the preview). A
  browser may then remember that a release is missing; release IDs are not reused, so this only
  matters for a release that is removed.

## The local server (`web/app/build/site-server.ts`)

`node build/site-server.ts <site directory> [--port 4173] [--host localhost]` serves a site as Pages
does; `npm run preview` serves `dist/` with it, and the E2E tests (`playwright.config.ts`) run against
it. `startSiteServer({ root, port?, host?, onRequest? })` starts it from a test; `onRequest` receives
every request (method, URL, path, headers, body) with the status it was answered with.

- Every response (200, 304, 308, 404, 405) has the headers of the `_headers` rules that match the
  request's path, read from the served directory's own `_headers`, and Pages' default
  `Cache-Control` where no rule sets one.
- Every 200 has a strong ETag (from the file's SHA-256); a matching `If-None-Match` gets 304.
- `/` and `/dir/` serve `index.html`; `/name` serves `name.html` if there is no file `name`;
  `/name.html` and `/dir/index.html` redirect (308) as on Pages.
- A missing file is 404 with the nearest `404.html` (a plain-text 404 if there is none: never a
  single-page fallback). `_headers`, `_redirects` and `_routes.json` are not served.
- Content types: `text/html`, `text/javascript`, `text/css` (with `charset=utf-8`), `application/json`,
  `application/wasm`; `application/octet-stream` for unknown extensions.
- GET and HEAD only (405 otherwise). No compression, no `_redirects`, no Functions.
- It reads the directory on every request, so a test may change the site while it runs (swap
  `releases.json`, add a release).

## Assembling a site with older releases (`web/app/build/assemble-site.ts`)

A build writes a site with one release. To deploy a site that still serves earlier releases (tabs
that run them keep loading their files; S16 item 3):

```bash
cd web/app
npm run build                                  # the new site in dist/
npm run assemble -- --site dist --previous <dir of an earlier site> [--previous ...] [--latest <id>]
```

- Every `r/<id>/` of each earlier site is checked against its manifest (every listed file present
  with its size and SHA-256, no file unlisted); a release that does not match stops the command
  before anything is copied. A release that the site already has must have the same manifest bytes.
- Each release is copied beside the releases, checked again, then moved into `r/`, so a release
  directory is complete or absent.
- `releases.json` is rewritten: every release (the new site's first, then each earlier site's in
  the order of its index), `latest` the new site's unless `--latest` names another, the shell of the
  new site, and the measurement endpoint of `LOSAT_WEB_ANALYTICS_ENDPOINT` if set, otherwise the new
  site's own setting.
- How many releases a deployment keeps is the Owner's decision (pending); the command keeps every
  release it is given.
