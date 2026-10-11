// The site layout and the boot page (S16 item 1, docs/web/release_manifest.md): the boot page
// starts the pinned release from /r/<release-id>/ only; a pin to a release the site does not
// serve gives a message and a button, never a silent switch; a localStorage that throws does
// not stop the app; every response (200, 304, 404) has the site's isolation headers and CSP
// (public/_headers, served by build/site-server.ts as Cloudflare Pages serves it); release
// files are immutable in caches; a missing file is 404.
import { expect, test, type APIRequestContext, type Page, type Response } from '@playwright/test';
import { BUILD_HAS_ENGINE, runSearch } from './support/browser';
import { E2E_ORIGIN } from './support/origin';

const PIN_KEY = 'losat-web:release';
const SITE_FILES = new Set(['/', '/boot.js', '/boot.css', '/releases.json', '/sw.js']);
const ISOLATION = ['cross-origin-opener-policy', 'cross-origin-embedder-policy', 'cross-origin-resource-policy', 'content-security-policy'];

/** The site's latest release, from its index. */
async function latestRelease(request: APIRequestContext): Promise<string> {
  const index = (await (await request.get('/releases.json')).json()) as { latest: string };
  return index.latest;
}

/** Waits until the app has started (the search form is there). */
async function appStarted(page: Page): Promise<void> {
  await expect(page.getByTestId('add-to-queue')).toBeVisible();
}

const pathOf = (url: string) => new URL(url).pathname;
// Requests to the site itself (WebKit also reports the app's blob: URLs, which are not requests to the site).
const sameOrigin = (url: string) => url.startsWith(`${E2E_ORIGIN}/`);

test('the boot page pins the latest release on a first visit and loads only that release’s files', async ({ page, request, browserName }) => {
  const latest = await latestRelease(request);
  const paths: string[] = [];
  page.context().on('request', (r) => {
    if (sameOrigin(r.url())) paths.push(pathOf(r.url()));
  });
  await page.goto('/');
  await appStarted(page);
  expect(await page.evaluate((key) => localStorage.getItem(key), PIN_KEY)).toBe(latest);
  await runSearch(page, 1);

  const release = `/r/${latest}/`;
  expect(paths.filter((path) => !SITE_FILES.has(path) && !path.startsWith(release))).toEqual([]);
  expect(paths).toContain(`${release}release-manifest.json`);
  expect(paths.filter((path) => path.startsWith(`${release}assets/`) && path.endsWith('.css')).length).toBeGreaterThan(0);
  // The engine modules are files of the release (a build with the engine). Chromium reports the
  // requests of the Engine worker; the other browsers may not.
  if (BUILD_HAS_ENGINE && browserName === 'chromium') {
    expect(paths.some((path) => path.startsWith(`${release}engine/`) && path.endsWith('.wasm'))).toBe(true);
  }
});

test('a pinned release that the site no longer serves gives a message and a button, never a silent switch', async ({ page, request }) => {
  const latest = await latestRelease(request);
  await page.goto('/');
  await appStarted(page);
  await page.evaluate((key) => localStorage.setItem(key, 'gone-release'), PIN_KEY);
  await page.reload();

  await expect(page.getByTestId('boot-message')).toHaveText('This browser uses LOSAT Web release gone-release, which this site no longer serves.');
  await expect(page.getByTestId('boot-message')).toHaveAttribute('data-kind', 'not-served');
  await expect(page.getByTestId('add-to-queue')).toHaveCount(0);
  expect(await page.evaluate((key) => localStorage.getItem(key), PIN_KEY)).toBe('gone-release');

  await page.getByTestId('boot-use-latest').click();
  await appStarted(page);
  expect(await page.evaluate((key) => localStorage.getItem(key), PIN_KEY)).toBe(latest);
  await expect(page.getByTestId('boot-message')).toHaveCount(0);
});

test('a first visit without the site index says so, and starts nothing', async ({ page }) => {
  await page.route('**/releases.json', (route) => route.fulfill({ status: 503, body: 'unavailable' }));
  await page.goto('/');
  await expect(page.getByTestId('boot-message')).toHaveAttribute('data-kind', 'no-index');
  await expect(page.getByTestId('add-to-queue')).toHaveCount(0);
});

for (const [name, stub] of [
  ['blocked (reading localStorage throws)', 'access'],
  ['full (writing localStorage throws)', 'write'],
] as const) {
  test(`the app starts when localStorage is ${name}`, async ({ page }) => {
    await page.addInitScript((kind) => {
      if (kind === 'access') {
        Object.defineProperty(window, 'localStorage', {
          configurable: true,
          get() {
            throw new DOMException('The operation is insecure.', 'SecurityError');
          },
        });
      } else {
        Storage.prototype.setItem = () => {
          throw new DOMException('The quota has been exceeded.', 'QuotaExceededError');
        };
      }
    }, stub);
    await page.goto('/');
    await appStarted(page);
    await expect(page.getByTestId('boot-message')).toHaveCount(0);
  });
}

test('every response has the isolation headers and the CSP, on a reload too; release files are immutable', async ({ browser, request }) => {
  // Without the service worker, which would answer the reload itself (V-OFF, offline.spec.ts, covers it).
  const context = await browser.newContext({ baseURL: E2E_ORIGIN, serviceWorkers: 'block' });
  const page = await context.newPage();
  const responses: Response[] = [];
  page.context().on('response', (r) => {
    if (sameOrigin(r.url())) responses.push(r);
  });
  await page.goto('/');
  await appStarted(page);
  await page.reload();
  await appStarted(page);
  expect(await page.evaluate(() => globalThis.crossOriginIsolated)).toBe(true);

  const urls = new Set<string>();
  for (const response of responses) {
    const headers = await response.allHeaders();
    const path = pathOf(response.url());
    for (const name of ISOLATION) expect(headers[name], `${name} of ${path} (${response.status()})`).toBeTruthy();
    if (path.startsWith('/r/')) expect(headers['cache-control'], path).toBe('public, max-age=31536000, immutable');
    else expect(headers['cache-control'], path).toBe('no-cache');
    urls.add(response.url());
  }
  expect([...urls].some((url) => pathOf(url).startsWith('/r/'))).toBe(true);
  // The reload revalidates the site-level files (no-cache): the browsers report 304s for them.
  expect(responses.some((response) => response.status() === 304)).toBe(true);

  // The server's 304 for each file the page loaded: the same headers (WebKit judges a
  // revalidated worker script by the headers of the 304; W3 decision 12).
  for (const url of urls) {
    const first = await request.get(url);
    const etag = first.headers()['etag'];
    expect(etag, url).toBeTruthy();
    const revalidated = await request.get(url, { headers: { 'If-None-Match': etag! } });
    expect(revalidated.status(), url).toBe(304);
    for (const name of ISOLATION) expect(revalidated.headers()[name], `${name} of the 304 of ${url}`).toBe(first.headers()[name]);
  }
  await context.close();
});

test('a missing file is 404 with the site’s headers, and no page stands in for it', async ({ page, request }) => {
  const latest = await latestRelease(request);
  for (const path of [`/r/${latest}/assets/missing.js`, '/r/gone-release/release-manifest.json', '/search']) {
    const response = await request.get(path);
    expect(response.status(), path).toBe(404);
    for (const name of ISOLATION) expect(response.headers()[name], `${name} of ${path}`).toBeTruthy();
  }
  const navigation = await page.goto('/no-such-page');
  expect(navigation?.status()).toBe(404);
  await expect(page.getByRole('heading', { name: 'Not found' })).toBeVisible();
  await expect(page.getByTestId('add-to-queue')).toHaveCount(0);
});
