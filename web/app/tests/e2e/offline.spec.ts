// V-OFF (plan §6.2; S16 items 2, 3 and 7): offline start from the service worker's caches, the
// explicit update without mixing releases, an incomplete update that changes nothing, an older
// tab that keeps its data and its release, and a pin that the site no longer serves.
//
// Two releases, voff-a and voff-b, are built like `npm run build` (FakeEngine, or the engine of
// LOSAT_WEB_REACTORS in an engine build: then the Wasm modules are among the stored files) into
// temporary directories; B's deployment is assembled with build/assemble-site.ts so that it
// keeps A, and gets a changed worker script, as a new deployment does. Each site is served by
// build/site-server.ts behind one origin (support/front-server.ts) that can switch deployments,
// go offline (port closed: a real network failure for the page and its worker, in every
// browser; Chromium and Firefox also get Playwright's offline emulation, so that
// navigator.onLine is false - in WebKit that emulation fails navigations, probe 1a in
// logs/wp-u/progress.md), and cut a download in the middle.
import { expect, test, type BrowserContext, type Page } from '@playwright/test';
import { cpSync, mkdtempSync, rmSync } from 'node:fs';
import { readFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { startSiteServer, type SiteServer } from '../../build/site-server';
import { expectResult, heldSessionLocks, OUTFMT7_MARK, runSearch, sessionDirectories } from './support/browser';
import { startFrontServer, type FrontServer } from './support/front-server';
import { E2E_PORT } from './support/origin';
import { assemble, buildSite, changeServiceWorker, damage, expectedCaches, largestAsset, manifestOf } from './support/sites';

const A = 'voff-a';
const B = 'voff-b';
const PIN_KEY = 'losat-web:release';

let workDir = '';
let dirs: { a: string; ab: string; damaged: string };
let servers: SiteServer[] = [];
let origins: { a: string; ab: string; damaged: string };
let front: FrontServer;

test.describe.configure({ timeout: 180_000 });

test.beforeAll(async ({ browserName }, testInfo) => {
  testInfo.setTimeout(300_000);
  workDir = mkdtempSync(join(tmpdir(), `losat-voff-${browserName}-`));
  const a = join(workDir, 'a');
  const b = join(workDir, 'b');
  buildSite(a, A);
  buildSite(b, B);
  const ab = join(workDir, 'ab');
  cpSync(b, ab, { recursive: true });
  assemble(ab, [a]);
  changeServiceWorker(ab, 'the deployment of voff-b');
  const damaged = join(workDir, 'damaged');
  cpSync(ab, damaged, { recursive: true });
  damage(damaged, B, largestAsset(manifestOf(damaged, B)));
  dirs = { a, ab, damaged };
  servers = [await startSiteServer({ root: a }), await startSiteServer({ root: ab }), await startSiteServer({ root: damaged })];
  origins = { a: servers[0]!.origin, ab: servers[1]!.origin, damaged: servers[2]!.origin };
  // Another port for each worker that runs this file (common rule: the E2E port + 100, + 200).
  front = await startFrontServer(E2E_PORT + 100 * (testInfo.parallelIndex + 1), origins.a);
});

test.afterAll(async () => {
  await front?.close();
  for (const server of servers) await server.close();
  if (workDir !== '') rmSync(workDir, { recursive: true, force: true });
});

test.beforeEach(async () => {
  front.upstream = origins.a;
  front.cut.clear();
  await front.setOnline(true);
});

const open = (page: Page) => page.goto(`${front.origin}/`);

async function appStarted(page: Page, release: string): Promise<void> {
  await expect(page.getByTestId('release-id')).toHaveText(release);
  await expect(page.getByTestId('add-to-queue')).toBeVisible();
}

async function offlineReady(page: Page): Promise<void> {
  await expect(page.getByTestId('offline-status')).toHaveAttribute('data-state', 'ready', { timeout: 60_000 });
}

async function goOffline(context: BrowserContext, browserName: string): Promise<void> {
  await front.setOnline(false);
  if (browserName !== 'webkit') await context.setOffline(true);
}

const pinOf = (page: Page) => page.evaluate((key) => localStorage.getItem(key), PIN_KEY);
const isolated = (page: Page) => page.evaluate(() => globalThis.crossOriginIsolated);

/** Every cache of the origin: by name, the URL path of each entry and the SHA-256 of its body. */
function storedFiles(page: Page): Promise<Record<string, Record<string, string>>> {
  return page.evaluate(async () => {
    const hex = (bytes: ArrayBuffer) => Array.from(new Uint8Array(bytes), (byte) => byte.toString(16).padStart(2, '0')).join('');
    const stored: Record<string, Record<string, string>> = {};
    for (const name of await caches.keys()) {
      const cache = await caches.open(name);
      const files: Record<string, string> = {};
      for (const request of await cache.keys()) {
        const response = await cache.match(request);
        files[new URL(request.url).pathname] = hex(await crypto.subtle.digest('SHA-256', await response!.arrayBuffer()));
      }
      stored[name] = files;
    }
    return stored;
  });
}

/** The release files that this document loaded (from the network or from the worker). */
const loadedReleaseFiles = (page: Page) =>
  page.evaluate(() =>
    performance
      .getEntriesByType('resource')
      .map((entry) => new URL(entry.name).pathname)
      .filter((path) => path.startsWith('/r/')),
  );

/** The release files that reached the server since request `from`. */
const servedReleaseFiles = (from: number) => front.requests.slice(from).map((request) => request.path).filter((path) => path.startsWith('/r/'));

const workers = (page: Page) =>
  page.evaluate(async () => {
    const registration = await navigator.serviceWorker.getRegistration('/');
    return { active: Boolean(registration?.active), waiting: Boolean(registration?.waiting) };
  });

/** The file that clicking a button saves (the browser's download). */
async function saved(page: Page, testid: string): Promise<{ name: string; text: string }> {
  const download = page.waitForEvent('download');
  await page.getByTestId(testid).click();
  const file = await download;
  return { name: file.suggestedFilename(), text: await readFile((await file.path())!, 'latin1') };
}

/** "Review update", then "Download and check". */
async function reviewAndDownload(page: Page): Promise<void> {
  await expect(page.getByTestId('update-notice')).toHaveAttribute('data-state', 'available', { timeout: 30_000 });
  await page.getByTestId('update-review').click();
  await expect(page.getByTestId('update-impact')).toBeVisible();
  await page.getByTestId('update-download').click();
}

test('1. a first visit says "Offline ready" once every file of the release and the shell is stored and checked, and nothing else is stored', async ({ page }) => {
  await open(page);
  await appStarted(page, A);
  await offlineReady(page);
  expect(await storedFiles(page)).toEqual(expectedCaches(dirs.a, A));
});

test('1b. a first visit whose download breaks does not say "Offline ready" and stores nothing; "Try again" then stores the release', async ({ page }) => {
  // A file that the page does not load before a threaded search, so that the worker downloads it.
  const manifest = manifestOf(dirs.a, A);
  const file = (manifest.assets.find((asset) => asset.path.includes('thread-worker')) ?? manifest.assets[0]!).path;
  front.cut.add(`/r/${A}/${file}`);
  await open(page);
  await appStarted(page, A);
  await expect(page.getByTestId('offline-status')).toHaveAttribute('data-state', 'failed', { timeout: 60_000 });
  await expect(page.getByTestId('offline-failed')).toContainText(file);
  expect(await storedFiles(page)).toEqual({});
  front.cut.clear();
  await page.getByTestId('offline-retry').click();
  await offlineReady(page);
  expect(await storedFiles(page)).toEqual(expectedCaches(dirs.a, A));
});

test('2. offline, a reload and a new page start the stored release, cross-origin isolated; a search runs, its results open, an export downloads and a session saves', async ({
  page,
  context,
  browserName,
}) => {
  await open(page);
  await appStarted(page, A);
  await offlineReady(page);

  await goOffline(context, browserName);
  await page.reload();
  await appStarted(page, A);
  expect(await isolated(page)).toBe(true);
  const second = await context.newPage();
  await open(second);
  await appStarted(second, A);
  expect(await isolated(second)).toBe(true);
  // The release is stored whole: offline, the page still says so.
  await offlineReady(second);

  await runSearch(second, 1);
  await expectResult(second, 1);
  expect((await saved(second, 'export-output')).text).toContain(OUTFMT7_MARK);
  expect((await saved(second, 'session-save')).name).toMatch(/^losat-session-\d{8}-\d{6}\.losat-session\.gz$/);
});

test('3. an update changes nothing until the user chooses; next to the switch it says what the reload does now; the tab then runs B whole, and the new worker waits until then', async ({
  page,
}) => {
  await open(page);
  await appStarted(page, A);
  await offlineReady(page);

  front.upstream = origins.ab;
  const deployed = front.requests.length;
  await page.reload();
  await appStarted(page, A);
  await expect(page.getByTestId('update-notice')).toHaveAttribute('data-state', 'available', { timeout: 30_000 });
  await expect(page.getByTestId('update-notice')).toContainText(`This site now serves LOSAT Web release ${B}; this tab uses release ${A}.`);
  // The deployment's new worker script is installed and waits while this page uses the old one.
  await expect.poll(() => workers(page), { timeout: 30_000 }).toEqual({ active: true, waiting: true });
  expect(await pinOf(page)).toBe(A);
  const before = await loadedReleaseFiles(page);
  expect(before.length).toBeGreaterThan(0);
  expect(before.filter((path) => !path.startsWith(`/r/${A}/`))).toEqual([]);
  expect(servedReleaseFiles(deployed).filter((path) => !path.startsWith(`/r/${A}/`))).toEqual([]);

  await reviewAndDownload(page);
  await expect(page.getByTestId('update-downloaded')).toBeVisible({ timeout: 60_000 });
  // Downloaded and checked: still release A, pinned, and the new worker still waits.
  await expect(page.getByTestId('release-id')).toHaveText(A);
  expect(await pinOf(page)).toBe(A);
  expect(await workers(page)).toEqual({ active: true, waiting: true });
  const stored = await storedFiles(page);
  expect(stored[`losat-web-release-${B}`]).toEqual(expectedCaches(dirs.ab, B)[`losat-web-release-${B}`]);
  // Next to the switch, the screen says again what the reload does, counted now (plan §5.9):
  // a run that completes after the review is named there.
  const impact = page.getByTestId('update-notice').getByTestId('update-impact-saved');
  await expect(impact).toHaveText('No completed run or candidate is lost.');
  await runSearch(page, 1);
  await expect(page.getByTestId('update-downloaded')).toBeVisible();
  await expect(impact).toContainText('1 completed run, 0 candidates and 0 notes are not kept by the reload');

  const switching = front.requests.length;
  await page.getByTestId('update-switch').click();
  await appStarted(page, B);
  expect(await pinOf(page)).toBe(B);
  expect(await isolated(page)).toBe(true);
  expect(await workers(page)).toEqual({ active: true, waiting: false });
  const after = await loadedReleaseFiles(page);
  expect(after.length).toBeGreaterThan(0);
  expect(after.filter((path) => !path.startsWith(`/r/${B}/`))).toEqual([]);
  expect(servedReleaseFiles(switching).filter((path) => !path.startsWith(`/r/${B}/`))).toEqual([]);
});

test('3b. a release that one tab downloaded and checked survives another tab\'s cleanup until the switch, which then starts it offline', async ({
  context,
  browserName,
}) => {
  const first = await context.newPage();
  await open(first);
  await appStarted(first, A);
  await offlineReady(first);
  front.upstream = origins.ab;
  await first.reload();
  await appStarted(first, A);
  await reviewAndDownload(first);
  await expect(first.getByTestId('update-downloaded')).toBeVisible({ timeout: 60_000 });

  // A second tab starts A (the pin), is Offline ready, and has the worker clean up: B stays,
  // because the first tab holds B's lock from its download until it switches.
  const second = await context.newPage();
  await open(second);
  await appStarted(second, A);
  await offlineReady(second);
  await expect.poll(async () => Object.keys(await storedFiles(second)).sort()).toEqual([`losat-web-release-${A}`, `losat-web-release-${B}`, 'losat-web-shell']);

  await goOffline(context, browserName);
  await first.getByTestId('update-switch').click();
  await appStarted(first, B);
  expect(await pinOf(first)).toBe(B);
  expect(await isolated(first)).toBe(true);
});

for (const [name, breakB] of [
  ['a file of B with other bytes', () => (front.upstream = origins.damaged)],
  [
    'a connection cut in the middle of a file of B',
    () => {
      front.upstream = origins.ab;
      front.cut.add(`/r/${B}/${largestAsset(manifestOf(dirs.ab, B))}`);
    },
  ],
] as const) {
  test(`4. an incomplete update (${name}) fails with its message; the pin stays A, and A still starts offline`, async ({ page, context, browserName }) => {
    await open(page);
    await appStarted(page, A);
    await offlineReady(page);
    const before = await storedFiles(page);

    breakB();
    await page.reload();
    await appStarted(page, A);
    await reviewAndDownload(page);
    await expect(page.getByTestId('update-failed')).toHaveText(
      new RegExp(`The update could not be downloaded completely \\(.+\\)\\. LOSAT Web stays on release ${A}\\.`),
      { timeout: 60_000 },
    );
    expect(await pinOf(page)).toBe(A);
    await expect(page.getByTestId('release-id')).toHaveText(A);
    expect(await storedFiles(page)).toEqual(before);

    await goOffline(context, browserName);
    await page.reload();
    await appStarted(page, A);
    expect(await isolated(page)).toBe(true);
  });
}

test('5. an older tab keeps its results and its release while another tab switches to B; A goes once that tab is closed', async ({ context }) => {
  const first = await context.newPage();
  await open(first);
  await appStarted(first, A);
  await offlineReady(first);
  await runSearch(first, 1);
  const tokens = await heldSessionLocks(first);
  expect(tokens).toHaveLength(1);

  front.upstream = origins.ab;
  const second = await context.newPage();
  await open(second);
  await appStarted(second, A);
  await reviewAndDownload(second);
  await expect(second.getByTestId('update-downloaded')).toBeVisible({ timeout: 60_000 });
  await second.getByTestId('update-switch').click();
  await appStarted(second, B);
  // B's start-up has removed the temporary data that no tab owns, and stored B.
  await expect(second.getByTestId('storage-status')).toHaveAttribute('data-cleanup', 'done');
  await offlineReady(second);

  await expect(first.getByTestId('update-other-tab')).toHaveText(
    `LOSAT Web was switched to release ${B} in another tab. This tab keeps release ${A} until it is reloaded.`,
  );
  // Tab 1's tmp/<token>, under its lock losat-web:tmp:<token>, survived B's start-up (where the
  // results are in OPFS: a browser without it keeps them in the tab's memory, and the page cannot
  // list OPFS either).
  if ((await first.getByTestId('storage-status').getAttribute('data-backend')) === 'opfs') {
    expect(await sessionDirectories(second)).toContain(tokens[0]);
  }
  await expectResult(first, 1);
  expect((await saved(first, 'export-output')).text).toContain(OUTFMT7_MARK);
  // A's cache stays while tab 1 runs A ...
  expect(Object.keys(await storedFiles(second)).sort()).toEqual([`losat-web-release-${A}`, `losat-web-release-${B}`, 'losat-web-shell']);

  // ... and goes at B's next start-up once tab 1 is closed.
  await first.close();
  await second.reload();
  await appStarted(second, B);
  await offlineReady(second);
  await expect.poll(async () => Object.keys(await storedFiles(second)).sort()).toEqual([`losat-web-release-${B}`, 'losat-web-shell']);
});

test('6. a pin to a release that the site no longer serves and that is not stored gets a message and a button, never a silent switch', async ({
  page,
  context,
  browserName,
}) => {
  await open(page);
  await appStarted(page, A);
  await offlineReady(page);
  await page.evaluate((key) => localStorage.setItem(key, 'voff-gone'), PIN_KEY);
  await page.reload();
  const message = page.getByTestId('boot-message');
  await expect(message).toHaveAttribute('data-kind', 'not-served');
  await expect(message).toHaveText('This browser uses LOSAT Web release voff-gone, which this site no longer serves.');
  await expect(page.getByTestId('boot-use-latest')).toHaveText(`Use release ${A}`);
  await expect(page.getByTestId('add-to-queue')).toHaveCount(0);
  expect(await pinOf(page)).toBe('voff-gone');

  // Offline, a pinned release that is not stored (and that this browser never asked for: the
  // site's 404 under /r/ is cached as immutable) cannot be checked: the boot page says so, and
  // switches nothing.
  await page.evaluate((key) => localStorage.setItem(key, 'voff-never-seen'), PIN_KEY);
  await goOffline(context, browserName);
  await page.reload();
  // WebKit keeps navigator.onLine true (no emulation): the boot cannot tell that it is offline.
  await expect(message).toHaveAttribute('data-kind', browserName === 'webkit' ? 'failed' : 'offline-not-stored');
  await expect(page.getByTestId('add-to-queue')).toHaveCount(0);
  expect(await pinOf(page)).toBe('voff-never-seen');

  await context.setOffline(false);
  await front.setOnline(true);
  await page.reload();
  await page.getByTestId('boot-use-latest').click();
  await appStarted(page, A);
  expect(await pinOf(page)).toBe(A);
});

test('6c. a pinned release that the site did not serve once and serves again starts: the 404 that the HTTP cache kept is read again past it', async ({
  page,
}) => {
  await open(page);
  await appStarted(page, A);
  await offlineReady(page);
  // B is pinned (another tab's switch, say) while this deployment does not serve it: its
  // manifest is 404, which the rules for /r/* make immutable in the HTTP cache.
  await page.evaluate(([key, release]) => localStorage.setItem(key, release), [PIN_KEY, B] as const);
  await page.reload();
  await expect(page.getByTestId('boot-message')).toHaveAttribute('data-kind', 'not-served');
  // The deployment serves B again (a rollback undone): the boot does not keep the cached 404.
  front.upstream = origins.ab;
  await page.reload();
  await appStarted(page, B);
  expect(await pinOf(page)).toBe(B);
});

test('6b. without a pin and without the network, the stored boot page says that the first start needs the network', async ({ page, context, browserName }) => {
  await open(page);
  await appStarted(page, A);
  await offlineReady(page);
  await page.evaluate((key) => localStorage.removeItem(key), PIN_KEY);
  await goOffline(context, browserName);
  await page.reload();
  // WebKit keeps navigator.onLine true: the boot tries the site index, and says that it could not read it.
  await expect(page.getByTestId('boot-message')).toHaveAttribute('data-kind', browserName === 'webkit' ? 'no-index' : 'first-visit-offline');
  await expect(page.getByTestId('add-to-queue')).toHaveCount(0);
});
