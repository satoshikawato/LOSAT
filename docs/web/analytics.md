# LOSAT Web usage statistics (GA4)

Plan §5.10 (b), design §13.2, session S16 item 4. The choice between an iframe measurement document
(a) and a same-origin endpoint (b) was made from the probe in
[`docs/evidence/losat_web_w7/ga4_probe/RESULTS.md`](../evidence/losat_web_w7/ga4_probe/RESULTS.md):
(b), Owner-delegated (standing instruction 2026-09-29, restated 2026-10-07).

## What is sent, and when

Only after the visitor allows it, LOSAT Web posts one event to its own site:

```text
POST /api/visit
Content-Type: application/json

{"v":1,"release":"<release ID>","client":"<random client ID>"}
```

- `release` is the build's release ID (`build/release.ts`, for example `0123456789ab-1a2b3c4d`).
- `client` is a random UUID (`crypto.randomUUID()`) made at the moment of consent and kept in this
  browser, so that a returning visitor is counted once.
- One request per page load, and one at the moment of consent. No retry, no queue, no other event,
  no error reporting. The request omits credentials and the referrer (`credentials: 'omit'`,
  `referrerPolicy: 'no-referrer'`), uses `keepalive` and `cache: 'no-store'`. A blocked, failing or
  offline request is silent and changes nothing else.
- Never sent: sequences, headers, file names, results, notes, settings, sizes or counts of any of
  them, the page URL or title. No document of LOSAT Web loads a Google script.

The Pages Function forwards the event to the GA4 Measurement Protocol as

```json
{"client_id": "<client>", "non_personalized_ads": true,
 "events": [{"name": "app_visit", "params": {"release_id": "<release>"}}]}
```

and nothing of the visitor's request: no headers, cookies, referrer or address. Google receives the
request from Cloudflare, so GA4's geography shows Cloudflare's data centre, not the visitor.

## When measurement exists

Only when the site index `/releases.json` has `"analytics": {"endpoint": "/api/visit"}` (plan §5.9;
`src/domain/site-index.ts` accepts only a same-origin path under `/api/`). Without it, or when the
index cannot be read (offline), the page shows no consent text, makes no client ID and sends no
request. A build without a release ID measures nothing.

## Consent and withdrawal

- The notice ("Usage statistics", `ui/AnalyticsConsent.vue`) is in the page's flow above the tabs,
  non-modal, and takes no focus; it never covers the search form. It says that one event per visit
  (that LOSAT Web was opened, and its release ID) goes through this site to Google Analytics with a
  random ID kept in the browser, and that sequences, file names, results, notes and settings are
  never sent. Buttons: "Allow" and "Don't allow".
- The footer (`ui/AnalyticsSetting.vue`) shows "Usage statistics: on" or "off"; "Change" reopens
  the notice.
- The choice is kept in `localStorage['losat-web:analytics']`: `{"v":1,"consent":"granted","client":"<uuid>"}`
  or `{"v":1,"consent":"denied"}`. "Don't allow" (also to withdraw) stores the refusal, which removes
  the client ID; nothing more is sent. A new consent makes a new client ID.
- When storage fails (disabled, private window, full): nothing is sent without a stored consent, the
  footer says that the choice could not be kept, and the next visit asks again.

Code: `src/domain/analytics.ts` (the event and the stored choice), `src/ports/analytics.ts`,
`src/infra/browser/analytics.ts` (storage and the request), `src/application/analytics.ts`
(`Analytics`: reads the site index through a `siteIndex` function and the stored choice, sends,
allows, denies), wiring in `src/composition.ts`. Tests: `tests/unit/analytics.test.ts`,
`tests/unit/visit-function.test.ts`, `tests/e2e/analytics.spec.ts`.

## The Pages Function (`web/app/functions/api/visit.ts`)

- `POST` only (other methods: 405 with `Allow: POST`).
- 400 unless the `Origin` header equals the request's own origin, the content type is
  `application/json` (a `charset=utf-8` parameter is accepted), the body is at most 512 bytes of
  UTF-8, and it is exactly the text above: no other key, value form, order or white space
  (`parseVisitBody` in `src/domain/analytics.ts`, shared with the page). Nothing is forwarded for a
  refused request.
- 204 for an accepted event. It forwards (`context.waitUntil`, after the answer, failures dropped)
  only when both environment variables are set; otherwise it forwards nothing and still answers 204,
  so a preview behaves as production does without measuring.
- Pages does not apply `_headers` to a Function's responses, so the Function sets its own
  (`Cache-Control: no-store`, `Cross-Origin-Resource-Policy: same-origin`, `X-Content-Type-Options: nosniff`).
- No dependency: it declares the few types it uses and imports `src/domain/analytics.ts`, which
  Pages' build (wrangler) bundles with it. The file also exports `handleVisit` and
  `measurementRequest` for the unit tests; Pages makes routes only of `onRequest*` exports.
- Location: Pages reads `functions/` at the root of the Pages project. The site is deployed by
  direct upload (S16 decision 9, [`deploy_cloudflare.md`](deploy_cloudflare.md)):
  `wrangler pages deploy dist`, run from `web/app` so that `web/app/functions/` is uploaded with
  it. The static assets are served without running the Function (only `/api/visit` routes to it).
  `deploy_cloudflare.md` holds the Cloudflare details, read from Cloudflare's documentation on
  2026-10-11; where this document and that one differ, that one applies.
- The `Origin` check keeps other sites' pages from posting through a visitor's browser; it does not
  authenticate. Anyone can post a well-formed event with a made-up client ID, which can only add a
  visit count.

## What the maintainer sets

In GA4 (never in the repository):

1. A web data stream for LOSAT Web; its measurement ID (`G-…`) is `GA4_MEASUREMENT_ID`.
2. A Measurement Protocol API secret of that stream (the stream's "Measurement Protocol API
   secrets"); it is `GA4_API_SECRET`.
3. A custom dimension: scope Event, event parameter `release_id` (to see visits per release).
4. Data retention (Admin, Data settings): the shortest that serves visit counts (2 months is
   recommended; the maintainer decides).
5. Keep Google signals and advertising links off; enhanced measurement has no effect (no gtag), and
   can be turned off on the stream as well.

In Cloudflare Pages (S16 decision 9; the steps are in [`deploy_cloudflare.md`](deploy_cloudflare.md)):

1. Environment variables `GA4_MEASUREMENT_ID` and `GA4_API_SECRET`, the secret encrypted, in the
   production environment only. Previews have neither: their Function answers 204 and forwards
   nothing.
2. The site is built with `LOSAT_WEB_ANALYTICS_ENDPOINT=/api/visit`, so its index has
   `"analytics": {"endpoint": "/api/visit"}`, and the same `dist/` is deployed as a preview and
   then promoted to production. A preview therefore shows the consent notice and posts to its own
   `/api/visit`, which measures nothing there; the preview check exercises that path. (A site
   built without the variable measures nothing anywhere.)

GA4 counts `app_visit` events and users by client ID. Measurement Protocol events carry no session
or engagement data, so GA4's session and engagement reports stay empty or partial; read the event
count and the users per `release_id`.

## Why (b)

From the probe (RESULTS.md): (b) works the same in Chromium, Firefox and WebKit under the site's
COOP/COEP/CORP and the unchanged CSP. (a) needs the measurement document to send COEP and
`CORP: cross-origin`, a `frame-src` addition to the app's CSP and `allow-same-origin`; it relies on
Google keeping `CORP: cross-origin` on gtag.js, and its client ID depends on third-party cookie and
storage rules that differ by browser (Chromium refuses gtag's default cookie in a cross-site frame,
WebKit set no third-party cookie). (b) keeps every document of the app free of Google code
(`web/AGENTS.md` rule 4), and Google receives only the release ID and a random client ID, from
Cloudflare's address. Blockers that list Google's domains stop (a); a same-origin path is not on the
common lists we know of, so the consent notice is (b)'s only gate, and it asks before anything is
sent.
