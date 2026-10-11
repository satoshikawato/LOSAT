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
removing the directory:

```
dist/
  _headers                       the response headers (copied from public/_headers, their only source)
  404.html                       the page of a missing file (src/boot/not-found.html)
  index.html                     the boot page (src/boot/index.html): no inline script or style
  boot.js                        the boot script (src/boot/boot.ts, one ES module)
  boot.css                       the style sheet of the boot page's messages and of 404.html
  releases.json                  the site index
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
- **Determinism.** No file holds a time; lists are sorted. Two builds of the same tree with the same
  engine modules write the same bytes.
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
| `app` | `{ "version", "build" }`: `package.json`'s version and the git short SHA of the build (`unknown` without git) |
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
The pure decisions are in `src/boot/plan.ts` (unit tests `tests/unit/boot.test.ts`); `src/boot/**`
may import `src/domain/**` only (ESLint).

**Build.** The site build (`build/site.ts`) builds `src/boot/boot.ts` with a second Vite build
(`configFile: false`, one ES module, minified) into `boot.js`, and copies `src/boot/index.html`,
`src/boot/boot.css` and `src/boot/not-found.html` (as `404.html`). The release's entry comes from
the `index.html` that Vite builds for it; the build stops if that page has anything the boot page
would not load (an inline script or style, a link other than a style sheet or a module preload).
