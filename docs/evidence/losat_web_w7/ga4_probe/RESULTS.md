# GA4 under COEP: variant (a) against variant (b) — probe results

Session S16 (W7), instruction item 4; plan §5.10; design §13.2. Probe run 2026-10-11 (UTC
03:10-03:17) on the S16 branch `s16-ga4` (from `f2f6af8c`). Nothing was configured in GA4, no real
measurement ID or API secret was used, and nothing was sent to Google. The only external requests
were the two header reads below (GET with the body discarded).

- **(a)** the design's plan: a measurement document on another origin (M), embedded as an iframe in
  the app (A), loads gtag.js and receives fixed events from A by `postMessage`.
- **(b)** the alternative: A posts one fixed event to its own `/api/visit`; a Cloudflare Pages
  Function forwards it to the GA4 Measurement Protocol. No document loads a Google script.

## How it was run

`./run.sh` (Playwright 1.61.1 from `web/app/node_modules`; `LOSAT_WEB_WEBKIT_EXECUTABLE` as in
`web/app/playwright.config.ts`), then `node summarize.mjs` for the tables. Files:
`server.mjs` (three local origins, records every request), `site/` (A's page, M's document, the
gtag stand-in), `probe.spec.ts`, `playwright.config.ts`, `origins.ts`, raw observations in
`results/<project>.json`, Google's headers in `google-headers.txt`.

| Stand-in | Origin | What it serves |
|---|---|---|
| A, the app | `http://localhost:4762` | A test page with the site's headers read from `web/app/public/_headers` (COOP same-origin, COEP require-corp, CORP same-origin, the CSP), plus `frame-src http://127.0.0.1:4862` where (a) needs it; `/api/visit` answers 204 (the Function's stand-in) |
| M, measurement document | `http://127.0.0.1:4862` | `m.html` with the variant's headers, a CSP allowing G and C, `frame-ancestors` A |
| G, googletagmanager.com | `http://localhost:4962/gtag/js` | A stand-in script **with the response headers recorded from googletagmanager.com**; it runs no Google code, marks that it ran and sends one hit to C by each gtag transport (sendBeacon, fetch keepalive no-cors, image) |
| C, the collect endpoint | `http://localhost:4962/g/collect` | Records what reaches it; answers 204 without CORP (the real endpoint's headers were not read) |

`localhost` and `127.0.0.1` are different sites, so M's frame is cross-site to A as on the real
hosts (all three are secure contexts). A CSP source cannot name an IPv6 literal, so G and C share
A's site; they are cross-site to M, where they are used, and A never loads them.

Browsers: Chromium 149.0.7827.55, Firefox 151.0, WebKit 26.5 (Playwright's builds). Playwright
changes two browsers' third-party storage rules, so two more projects restore the stock rules:
`chromium-stock-storage` (Playwright disables `ThirdPartyStoragePartitioning`; stock Chrome
partitions third-party storage) and `firefox-stock-storage` (Playwright sets
`network.cookie.cookieBehavior` 4; stock Firefox uses 5, Total Cookie Protection). Playwright's
WebKit is the WPE MiniBrowser on Linux, not Safari: Safari's tracking prevention is not in it.

### Google's response headers (2026-10-11T02:46Z, `google-headers.txt`)

| Header | `www.googletagmanager.com/gtag/js?id=G-XXXXXXXXXX` | `www.google-analytics.com/analytics.js` |
|---|---|---|
| `cross-origin-resource-policy` | `cross-origin` | `cross-origin` |
| `access-control-allow-origin` | `*` (with `access-control-allow-credentials: true`) | absent |
| `content-type` | `application/javascript; charset=UTF-8` | `text/javascript` (`x-content-type-options: nosniff`) |
| `cache-control` | `private, max-age=900` | `public, max-age=7200` |
| other | — | `cross-origin-opener-policy-report-only`, `content-security-policy-report-only` (report-only, no effect on a page that loads it) |

`CORP: cross-origin` is what lets gtag.js load in a document with `COEP: require-corp`. It is
Google's choice and may change; the probe cannot see the headers of the further scripts and
endpoints that the real gtag.js may use (reading them would send requests to Google).

## Results

Frame loading, script, hits and `postMessage` checks (A with `sandbox="allow-scripts allow-same-origin"` unless the row says otherwise):

| Variant | Chromium 149 | Firefox 151 | WebKit 26.5 |
|---|---|---|---|
| **b**: POST `/api/visit` (same origin), site CSP unchanged | 204; reached 1; body exact; no Referer; no Cookie; no CSP violation | same | same |
| a1: M without COEP/CORP | frame refused (COEP) | frame refused | frame refused |
| **a2**: M with COEP require-corp and CORP cross-origin | frame loads; gtag ran; C got beacon+fetch+image; checks pass | same | same |
| a3: M with COEP credentialless and CORP cross-origin | frame loads; gtag ran; C got beacon+fetch+image; checks pass | same | frame refused (no COEP credentialless) |
| a3n: M with COEP credentialless, no CORP | frame refused | frame refused | frame refused |
| a4: `<iframe credentialless>`, M without COEP/CORP | frame loads; gtag ran; C got all three | frame refused (attribute unknown) | frame refused (attribute unknown) |
| a2 under the site's CSP unchanged | frame refused by A's CSP (`frame-src`) | same | same |
| a2 without `sandbox` | as a2 | as a2 | as a2 |
| a2 with `sandbox="allow-scripts"` (opaque origin) | M's hello has origin `"null"`: A's origin check refuses it; a post to M's origin is dropped, only `'*'` reaches M; M's `localStorage` throws SecurityError | same | no hello reached A (not investigated) |
| A navigated again (commit) after a2 | commits | **does not commit in 15 s** (after b: commits) | commits |

In every loaded (a) case the hits reached C by all three transports, although COEP made the
page-side `fetch` reject and the image fail (COEP stops reading a response without CORP, not
sending the request). The hits carried M's origin as Referer; the real gtag.js also sends M's URL,
title, screen size, language and the browser's user agent, and Google sees the visitor's IP address.

The identifier kept by M's frame when A is opened again (a new tab, then a second tab; then M
opened top-level to see whether the frame's storage is partitioned from M's own):

| Variant | chromium (Playwright rules) | chromium-stock-storage | firefox (Playwright rules) | firefox-stock-storage | webkit |
|---|---|---|---|---|---|
| a2 | localStorage kept (unpartitioned); `SameSite=None; Secure` cookie kept; default (Lax) cookie not set | localStorage kept (partitioned); None cookie kept; Lax cookie not set | localStorage kept (unpartitioned); both cookies kept | localStorage kept (partitioned); both cookies kept | localStorage kept (unpartitioned); **no cookie set** |
| a3 | as a2 | as a2 | as a2 | as a2 | frame refused |
| a4 | **nothing kept** (a new ID on every load) | nothing kept | frame refused | frame refused | frame refused |

## Findings

1. **(b) works the same in all three browsers** under COOP/COEP and the unchanged CSP: the request
   reached the endpoint once with exactly the fixed body, without Referer or Cookie, with no CSP
   violation and no other request. It needs no header or CSP change.
2. **(a) needs M to send both COEP (`require-corp`, or `credentialless` outside WebKit) and
   `CORP: cross-origin`**, and A's CSP needs `frame-src <M's origin>`. Without either header the
   frame is refused in every browser. Only a2 loads in all three; WebKit has no COEP
   `credentialless`, and only Chromium has `<iframe credentialless>`, whose storage is new on each
   load (GA4 would count every visit as a new user).
3. **In a2, gtag.js loads only because Google sends `CORP: cross-origin` today.** M's own COEP then
   applies to everything gtag.js loads; the probe can check only the one script's headers.
4. **Counting users in (a) depends on third-party storage rules that differ by browser.** gtag keeps
   its client ID in the `_ga` cookie of M's document. In a cross-site frame, Chromium refuses a
   cookie without `SameSite=None` (gtag's default; `cookie_flags: 'SameSite=None;Secure'` would be
   needed), and Playwright's WebKit set no cookie at all, as Safari blocks all third-party cookies.
   With the stock rules Chromium and Firefox partition the frame's storage by top-level site, which
   still keeps an ID per site. Safari's own rules (third-party cookies blocked, third-party storage
   partitioned, script-writable storage capped at 7 days without interaction) were not testable here
   and are stated from their published descriptions, which this probe could not re-read (no external
   requests allowed). In (b) the client ID is first-party `localStorage` of the app itself, the same
   in every browser (Safari's 7-day cap applies to both variants alike).
5. **`postMessage` checks work only with `allow-same-origin`.** With
   `sandbox="allow-scripts allow-same-origin"` (or no sandbox) both sides see the right origin and
   source. Without `allow-same-origin` M's origin is opaque: A must refuse its messages, A cannot
   address it except by `'*'`, and M has no storage, so (a) would need `allow-same-origin`.
6. **Playwright's Firefox did not navigate A again after a cross-origin frame had been in it**
   (a1 to a3, with or without sandbox; a same-origin frame and the no-frame page navigated
   normally). Whether stock Firefox does the same was not established; it is a risk for (a) only.

### What a blocker does

Common filter lists (EasyPrivacy, uBlock Origin's defaults, Firefox's strict tracking protection,
Brave Shields) block `googletagmanager.com` and `google-analytics.com`. In (a) the frame on the
site's own measurement origin may load, but gtag.js and the hits are blocked, so a blocker stops
the measurement and leaves an empty frame. In (b) the request goes to the site's own path
`/api/visit`; we know of no common list rule for that path (generic rules target names such as
`analytics`, `collect`, `gtag`, `track`), but a list can add one at any time. In both variants a
blocked or failed request changes nothing else in the app (Phase 2's E2E tests abort the endpoint
and answer 500). Because (b) is harder for a blocker to see, the consent notice is its only gate:
it asks before anything is sent and keeps "Don't allow" as easy as "Allow".

### No research data in either path

Neither path carries research data. (b) sends exactly `{"v":1,"release":<release ID>,"client":<random UUID>}`,
the Function accepts only that text and forwards only the two IDs (no headers, no address); Google
sees Cloudflare's address, not the visitor's. In (a) A would post only fixed events to M, but gtag.js
in M also sends M's page URL and title, the screen, the language and the user agent from the
visitor's browser, and the parameters must be kept free of research data by configuration
(enhanced measurement off, no free-form dataLayer), which (b) does not need.

## Summary table (browser x variant; does the variant work end to end?)

| Variant | Chromium | Firefox | WebKit | CSP or header change in the app |
|---|---|---|---|---|
| (b) same-origin `/api/visit` | yes | yes | yes | none |
| (a) a1 M without COEP/CORP | no | no | no | `frame-src` |
| (a) a2 M COEP require-corp + CORP | yes (ID: cookie only with SameSite=None) | yes | loads; no third-party cookie | `frame-src`; depends on Google's CORP |
| (a) a3 M COEP credentialless + CORP | yes | yes | no | `frame-src` |
| (a) a4 `<iframe credentialless>` | loads; no lasting ID | no | no | `frame-src` |

## Recommendation: (b)

Adopted, Owner-delegated (standing instruction 2026-09-29, restated 2026-10-07):

- (b) works the same in Chromium, Firefox and WebKit under the site's COOP/COEP/CORP and the
  unchanged CSP; (a) works only with M sending COEP and CORP, a `frame-src` addition, and
  `allow-same-origin`, and in no variant identically in all three browsers.
- (b) loads no third-party script into any document, so the app's documents (which handle research
  data) stay free of Google code, as `web/AGENTS.md` rule 4 asks; (a) puts gtag.js in a frame of
  the app.
- (a) depends on Google's response headers (CORP on gtag.js and whatever it loads) and on
  third-party storage rules that differ by browser and change over time; (b) depends on neither: the
  client ID is the app's own first-party storage.
- (b) sends Google only the release ID and a random client ID, from Cloudflare's address; (a) sends
  Google the visitor's address, user agent and page details, and needs configuration to keep the
  rest out.
- Cost of (b): one Pages Function and two Cloudflare secrets (set by the maintainer), and GA4's
  geography reflects Cloudflare's data centre, not the visitor.
