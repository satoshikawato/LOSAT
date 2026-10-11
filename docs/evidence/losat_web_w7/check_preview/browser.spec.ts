// Checks 1, 3 and 4 of the preview check (README.md), in Chromium, Firefox and WebKit: the
// app starts through the boot page and is cross-origin isolated; a search cancelled while it
// runs is followed by a search that completes (WebKit: the Engine worker re-created after a
// cancel needs COEP/CORP on the 304 of its script); "Offline ready" appears and, in Chromium,
// an offline reload starts the app.
import { expect, test, type Page } from '@playwright/test';
import { fasta, openFiles, paste, program, submit, waitStatus } from '../../../../web/app/tests/e2e/support/search';
import { previewOrigin, readSite, record, request, type IndexView, type ManifestView } from './support';

const ORIGIN = previewOrigin();
const PIN_KEY = 'losat-web:release';
/** The proteome of the E2E tests' slow search: seconds in every browser, long enough to cancel (tests/e2e/search.spec.ts). */
const SLOW = 'NZ_CP006932.faa';

let site: { index: IndexView; manifest: ManifestView; swPresent: boolean };

test.beforeAll(async () => {
  const { index, manifest } = await readSite(ORIGIN);
  const sw = await request(ORIGIN, '/sw.js');
  site = { index, manifest, swPresent: sw.status === 200 && /javascript/i.test(sw.headers['content-type'] ?? '') };
});

async function appStarted(page: Page): Promise<void> {
  await expect(page.getByTestId('add-to-queue')).toBeVisible({ timeout: 60_000 });
}

test('the app starts through the boot page, isolated, on the latest release', async ({ page, browserName }, testInfo) => {
  const requested: string[] = [];
  page.on('request', (r) => {
    if (r.url().startsWith(`${ORIGIN}/`)) requested.push(new URL(r.url()).pathname);
  });
  await page.goto('/');
  await appStarted(page);
  const state = await page.evaluate((key) => ({
    isolated: crossOriginIsolated,
    sharedArrayBuffer: typeof SharedArrayBuffer === 'function',
    pin: localStorage.getItem(key),
    userAgent: navigator.userAgent,
  }), PIN_KEY);
  // The release shown by the app (its header, when this build shows it), else the pin the boot page wrote.
  const shownElement = page.getByTestId('release-id');
  const shown = (await shownElement.count()) > 0 ? ((await shownElement.first().textContent()) ?? '').trim() : state.pin;
  const latest = site.index.latest;
  await record(testInfo, {
    browser: browserName,
    latest,
    shown,
    shownFrom: (await shownElement.count()) > 0 ? 'release-id' : 'localStorage pin',
    ...state,
    entryRequested: requested.includes(`/r/${latest}/${site.manifest.entry.script}`),
    foreignPaths: requested.filter((path) => !['/', '/boot.js', '/boot.css', '/releases.json', '/sw.js', '/favicon.ico'].includes(path) && !path.startsWith(`/r/${latest}/`)),
  });
  expect(state.isolated, 'crossOriginIsolated').toBe(true);
  expect(state.sharedArrayBuffer, 'SharedArrayBuffer').toBe(true);
  expect(shown).toBe(latest);
  expect(requested, 'the entry script of the latest release').toContain(`/r/${latest}/${site.manifest.entry.script}`);
  await expect(page.getByTestId('boot-message')).toHaveCount(0);
});

test('a search cancelled while it runs, then a new search completes', async ({ page, browserName }, testInfo) => {
  test.skip(site.manifest.engine === null, 'a FakeEngine site: no search long enough to cancel');
  await page.goto('/');
  await appStarted(page);
  await program(page, 'blastp');
  await openFiles(page, 'query', [{ name: SLOW, text: fasta(SLOW) }]);
  await openFiles(page, 'subject', [{ name: SLOW, text: fasta(SLOW) }]);
  await submit(page);
  await expect(page.getByTestId('run-1-status')).toHaveText('running', { timeout: 60_000 });
  await page.getByTestId('run-1-cancel').click();
  await waitStatus(page, 1, 'cancelled');
  // The Engine worker is created again for the next run; the same subject, another query.
  await page.getByTestId('query-source-0-remove').click();
  await paste(page, 'query', '>small\nMKLVVLAAGGHHKLMKLVVLAAGGHHKL\n');
  await submit(page);
  await waitStatus(page, 2, 'completed');
  await record(testInfo, { browser: browserName, run1: 'cancelled', run2: 'completed' });
});

test('"Offline ready" appears; in Chromium an offline reload starts the app', async ({ page, context, browserName }, testInfo) => {
  test.skip(!site.swPresent, '/sw.js is not there: this site has no service worker');
  await page.goto('/');
  await appStarted(page);
  await expect(page.getByTestId('offline-status')).toHaveAttribute('data-state', 'ready', { timeout: 180_000 });
  const data: Record<string, unknown> = { browser: browserName, offlineReady: true };
  if (browserName === 'chromium') {
    // Playwright's offline switch is reliable only in Chromium (S16 decision 6).
    await context.setOffline(true);
    try {
      await page.reload();
      await appStarted(page);
      data['offlineReload'] = { started: true, isolated: await page.evaluate(() => crossOriginIsolated) };
      expect(data['offlineReload']).toEqual({ started: true, isolated: true });
    } finally {
      await context.setOffline(false);
      await record(testInfo, data);
    }
  } else {
    data['offlineReload'] = `not run in ${browserName}: only Chromium can go offline reliably`;
    await record(testInfo, data);
  }
});
