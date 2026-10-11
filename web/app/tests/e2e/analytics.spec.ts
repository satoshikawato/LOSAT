// Usage statistics (S16 item 4, plan §5.10 (b), docs/web/analytics.md) through the real
// application. The site index `/releases.json` turns measurement on (`analytics.endpoint`);
// the tests fulfil it and record or answer `/api/visit` with page.route, so they work with
// or without a served index. Before consent and after "Don't allow" nothing is sent; after
// "Allow" exactly one request per page load with exactly the fixed body; withdrawing stops
// it; storage that throws sends nothing; and a blocked, failing or offline endpoint changes
// nothing in a search, the results, an extraction, an export and a session save. No request
// goes to another origin.
import { expect, test, type BrowserContext, type Page, type Request } from '@playwright/test';
import { findReactors } from '../../build/reactors';
import { releaseId } from '../../build/release';
import { BUILD_HAS_ENGINE } from './support/browser';
import { E2E_ORIGIN } from './support/origin';
import { chooseRun, paste, program, showOutput, submit, waitStatus } from './support/search';

test.setTimeout(120_000);
// Without the site's service worker: once the first load's worker has claimed the page, WebKit
// under Playwright no longer routes the page's requests through context.route (the reload's
// /releases.json then reaches the real site, which does not measure), although the worker does
// not answer /releases.json or /api/visit itself (src/sw/serve.ts). These tests are not about
// offline use (V-OFF: offline.spec.ts). Found by the three-browser run of WP-GP (S16).
test.use({ serviceWorkers: 'block' });

/** The release ID of the build under test (vite.config.ts computes it from the same inputs). */
const RELEASE = releaseId(findReactors());
const STORAGE_KEY = 'losat-web:analytics';
const CLIENT_ID = /^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$/;

/** A site index naming this build, with or without the measurement endpoint. */
const siteIndex = (measured: boolean) => ({
  format: 1,
  latest: RELEASE,
  releases: [{ id: RELEASE, version: '0.0.0', build: 'e2e' }],
  ...(measured ? { analytics: { endpoint: '/api/visit' } } : {}),
});

type Endpoint = 'ok' | 'abort' | 'error';

interface Site {
  /** The requests that reached /api/visit, in order. */
  readonly visits: Request[];
  /** Every URL requested by the context's pages and workers. */
  readonly urls: string[];
  endpoint: Endpoint;
}

/** Serves the site index (or none) and answers /api/visit; records every request of the context. */
async function site(context: BrowserContext, index: 'measured' | 'unmeasured' | 'none'): Promise<Site> {
  const recorded: Site = { visits: [], urls: [], endpoint: 'ok' };
  context.on('request', (request) => recorded.urls.push(request.url()));
  if (index !== 'none') {
    await context.route('**/releases.json', (route) => route.fulfill({ json: siteIndex(index === 'measured') }));
  }
  await context.route('**/api/visit', async (route) => {
    recorded.visits.push(route.request());
    if (recorded.endpoint === 'abort') await route.abort('blockedbyclient');
    else await route.fulfill({ status: recorded.endpoint === 'error' ? 500 : 204 });
  });
  return recorded;
}

/** Opens the app and waits until it has read the site index. */
async function open(page: Page, reload = false): Promise<void> {
  const index = page.waitForResponse((response) => response.url().endsWith('/releases.json'));
  if (reload) await page.reload();
  else await page.goto('/');
  await index;
  await expect(page.getByTestId('add-to-queue')).toBeVisible();
  // Lets the application act on the index (it parses it and updates its state).
  await page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => setTimeout(resolve, 100))));
}

const stored = (page: Page) => page.evaluate((key) => localStorage.getItem(key), STORAGE_KEY);

/** Fails when a request of the context went to another origin (data: and blob: URLs are the page's own). */
function expectOnlyThisOrigin(recorded: Site): void {
  const other = recorded.urls.filter((url) => !url.startsWith(`${E2E_ORIGIN}/`) && !/^(data|blob):/.test(url));
  expect(other).toEqual([]);
}

/** Checks that a request is the fixed event: its body, method and headers. */
function expectVisit(request: Request): string {
  expect(request.method()).toBe('POST');
  expect(request.headers()['content-type']).toBe('application/json');
  const body = request.postData() ?? '';
  const parsed = JSON.parse(body) as { v: number; release: string; client: string };
  expect(body).toBe(JSON.stringify({ v: 1, release: RELEASE, client: parsed.client }));
  expect(parsed.client).toMatch(CLIENT_ID);
  return parsed.client;
}

test.describe('usage statistics', () => {
  for (const index of ['unmeasured', 'none'] as const) {
    test(`a site index ${index === 'none' ? 'that is not served' : 'without analytics'}: no consent text, no client ID, no request`, async ({ page, context }) => {
      const recorded = await site(context, index);
      await open(page);
      await expect(page.getByTestId('analytics-consent')).toHaveCount(0);
      await expect(page.getByTestId('analytics-setting')).toHaveCount(0);
      expect(await stored(page)).toBeNull();
      await open(page, true);
      expect(recorded.visits).toEqual([]);
      expect(recorded.urls.filter((url) => url.includes('/api/'))).toEqual([]);
      expectOnlyThisOrigin(recorded);
    });
  }

  test('before consent and after "Don\'t allow" nothing is sent, also after reloads', async ({ page, context }) => {
    const recorded = await site(context, 'measured');
    await open(page);
    const consent = page.getByTestId('analytics-consent');
    await expect(consent).toBeVisible();
    await expect(consent).toContainText('one event through this site to Google Analytics');
    await expect(consent).toContainText(`its release ID (${RELEASE})`);
    await expect(consent).toContainText('Your sequences, file names, results, notes and settings are never sent.');
    await expect(page.getByTestId('analytics-setting')).toHaveText('Usage statistics: off');
    await open(page, true);
    await expect(consent).toBeVisible();
    expect(recorded.visits).toEqual([]);

    await page.getByTestId('analytics-deny').click();
    await expect(consent).toHaveCount(0);
    await expect(page.getByTestId('analytics-setting')).toHaveText('Usage statistics: off');
    expect(await stored(page)).toBe('{"v":1,"consent":"denied"}');
    for (let n = 0; n < 2; n++) {
      await open(page, true);
      await expect(consent).toHaveCount(0);
      await expect(page.getByTestId('analytics-setting')).toHaveText('Usage statistics: off');
    }
    expect(recorded.visits).toEqual([]);
    expectOnlyThisOrigin(recorded);
  });

  test('after "Allow" exactly one request per page load, with exactly the fixed body and the client ID made at consent', async ({ page, context }) => {
    const recorded = await site(context, 'measured');
    await open(page);
    await page.getByTestId('analytics-allow').click();
    await expect(page.getByTestId('analytics-consent')).toHaveCount(0);
    await expect(page.getByTestId('analytics-setting')).toHaveText('Usage statistics: on');
    await expect.poll(() => recorded.visits.length).toBe(1);
    const client = expectVisit(recorded.visits[0]!);
    expect(JSON.parse((await stored(page))!)).toEqual({ v: 1, consent: 'granted', client });
    // Neither credentials nor a referrer go with it.
    expect(await recorded.visits[0]!.headerValue('cookie')).toBeNull();
    expect(await recorded.visits[0]!.headerValue('referer')).toBeNull();
    // The Function accepts only a request whose Origin is the site's own (functions/api/visit.ts).
    expect(await recorded.visits[0]!.headerValue('origin')).toBe(E2E_ORIGIN);

    for (const loads of [2, 3]) {
      await open(page, true);
      await expect(page.getByTestId('analytics-consent')).toHaveCount(0);
      await expect.poll(() => recorded.visits.length).toBe(loads);
      expect(expectVisit(recorded.visits[loads - 1]!)).toBe(client);
    }
    // No other request follows (no retry, no second event).
    await page.waitForTimeout(500);
    expect(recorded.visits).toHaveLength(3);
    expectOnlyThisOrigin(recorded);
  });

  test('withdrawing sends nothing more and removes the client ID', async ({ page, context }) => {
    const recorded = await site(context, 'measured');
    await open(page);
    await page.getByTestId('analytics-allow').click();
    await expect.poll(() => recorded.visits.length).toBe(1);
    await page.getByTestId('analytics-change').click();
    await expect(page.getByTestId('analytics-consent')).toBeVisible();
    await page.getByTestId('analytics-deny').click();
    await expect(page.getByTestId('analytics-setting')).toHaveText('Usage statistics: off');
    expect(await stored(page)).toBe('{"v":1,"consent":"denied"}');
    await open(page, true);
    await expect(page.getByTestId('analytics-consent')).toHaveCount(0);
    await page.waitForTimeout(500);
    expect(recorded.visits).toHaveLength(1);
    expectOnlyThisOrigin(recorded);
  });

  test('storage that throws: "Allow" sends nothing and says so, the next load asks again, and a search runs', async ({ page, context }) => {
    await context.addInitScript(() => {
      const refuse = () => {
        throw new DOMException('The operation is insecure.', 'SecurityError');
      };
      try {
        Object.defineProperty(window, 'localStorage', { configurable: true, get: refuse });
      } catch {
        // Where the property cannot be replaced, its methods throw.
      }
      for (const name of ['getItem', 'setItem', 'removeItem', 'clear', 'key'] as const) {
        Object.defineProperty(Storage.prototype, name, { configurable: true, value: refuse });
      }
    });
    const recorded = await site(context, 'measured');
    await open(page);
    await expect(page.getByTestId('analytics-consent')).toBeVisible();
    await page.getByTestId('analytics-allow').click();
    await expect(page.getByTestId('analytics-setting')).toHaveText('Usage statistics: off');
    await expect(page.getByTestId('analytics-unsaved')).toContainText('could not keep your choice');
    await program(page, 'blastn');
    await paste(page, 'query', '>q1\nACGTACGTACGTACGTACGT\n');
    await paste(page, 'subject', '>s1\nACGTACGTACGTACGTACGT\n');
    await submit(page);
    await waitStatus(page, 1, 'completed');
    await open(page, true);
    await expect(page.getByTestId('analytics-consent')).toBeVisible();
    expect(recorded.visits).toEqual([]);
    expectOnlyThisOrigin(recorded);
  });

  test.describe('on a phone', () => {
    test.use({ viewport: { width: 390, height: 844 }, isMobile: false });

    test('the notice is in the page\'s flow and covers no part of the search form', async ({ page, context }) => {
      await site(context, 'measured');
      await open(page);
      const consent = page.getByTestId('analytics-consent');
      await expect(consent).toBeVisible();
      expect(await consent.evaluate((element) => getComputedStyle(element).position)).toBe('static');
      const notice = (await consent.boundingBox())!;
      for (const testid of ['program-blastn', 'query-input', 'subject-input', 'add-to-queue']) {
        const box = (await page.getByTestId(testid).boundingBox())!;
        const overlaps = box.x < notice.x + notice.width && notice.x < box.x + box.width && box.y < notice.y + notice.height && notice.y < box.y + box.height;
        expect(overlaps, testid).toBe(false);
      }
      expect(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth)).toBe(true);
      if (process.env['LOSAT_WEB_ANALYTICS_SHOTS']) await page.screenshot({ path: `${process.env['LOSAT_WEB_ANALYTICS_SHOTS']}/phone-${test.info().project.name}.png`, fullPage: false });
    });
  });

  test.describe('a blocked, failing or offline endpoint changes nothing', () => {
    test.skip(BUILD_HAS_ENGINE, "the candidates follow the FakeEngine's HSPs");
    test.setTimeout(180_000);

    for (const mode of ['abort', 'error', 'offline'] as const) {
      test(`endpoint ${mode === 'abort' ? 'blocked' : mode === 'error' ? 'answering 500' : 'unreachable: the page is offline'}: a search, results, an extraction, an export and a session save work`, async ({
        page,
        context,
        browserName,
      }) => {
        const recorded = await site(context, 'measured');
        recorded.endpoint = mode === 'error' ? 'error' : 'abort';
        await open(page);
        if (mode !== 'offline') {
          // A stored consent, then a page load: its request fails as well.
          await page.getByTestId('analytics-allow').click();
          await expect.poll(() => recorded.visits.length).toBe(1);
          await open(page, true);
          await expect.poll(() => recorded.visits.length).toBe(2);
        }
        await program(page, 'blastn');
        await paste(page, 'query', '>q1 first query\nACGTACGTACGTACGTACGTACGTACGTACGTACGT\n');
        await paste(page, 'subject', '>s1 first subject\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n');
        if (mode === 'offline') {
          // In Playwright's WebKit, offline emulation (and a route that matches blob: URLs) also
          // fails the page's own reads of its input Blobs ("The I/O read operation failed"), so
          // there every HTTP request is aborted instead.
          if (browserName === 'webkit') await context.route((url) => url.protocol === 'http:' || url.protocol === 'https:', (route) => route.abort('internetdisconnected'));
          else await context.setOffline(true);
          await page.getByTestId('analytics-allow').click();
        }
        await expect(page.getByTestId('analytics-setting')).toHaveText('Usage statistics: on');
        await submit(page);
        await waitStatus(page, 1, 'completed');
        await chooseRun(page, 1);
        await page.getByTestId('results-view-hits').click();
        await page.getByTestId('descriptions-select-all').check();
        await page.getByTestId('descriptions-add-candidates').click();
        await expect(page.getByTestId('candidates-added')).toContainText('added to Candidates');
        await page.getByTestId('tab-candidates').click();
        let download = page.waitForEvent('download');
        await page.getByTestId('extract-download').click();
        expect((await download).suggestedFilename()).toBe('losat-candidates.fa');
        await showOutput(page, 1, 6);
        download = page.waitForEvent('download');
        await page.getByTestId('export-output').click();
        expect((await download).suggestedFilename()).toBe('losat-run1-blastn.outfmt6.txt');
        download = page.waitForEvent('download');
        await page.getByTestId('session-save').click();
        expect((await download).suggestedFilename()).toMatch(/^losat-session-\d{8}-\d{6}\.losat-session\.gz$/);
        await expect(page.getByTestId('session-message')).toContainText('Saved');
        await context.setOffline(false);
        expectOnlyThisOrigin(recorded);
      });
    }
  });
});
