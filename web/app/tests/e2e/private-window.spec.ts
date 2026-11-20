// A private window, a site whose data the browser blocks, or a storage that fails in a normal
// window: the site cannot store files (line WQA's sweep, D5, 2026-11-19). The application keeps
// results in the tab's memory and says why in LOSAT's words: the note quoted the browser's
// exception, "(Security error when calling GetDirectory)" in Firefox and "(The operation failed
// for an unknown transient reason (e.g. out of memory).)" in WebKit, which reads as if memory ran
// out. The browser's error, name and message, stays in the note's data-error attribute, for
// developers. The note says what differs from a tab that stores files: the results are held in
// memory, which keeps at most 512 MB of them. Losing them on a reload is not one of the
// differences, since LOSAT Web keeps nothing for the next visit in any browser (review of D5, M1).
//
// Private windows. Firefox: browser.privatebrowsing.autostart; it refuses the storage
// (SecurityError). WebKit: Playwright's context is an ephemeral session, as Safari's private
// windows are, and its getDirectory() rejects with the UnknownError that a failing storage gives
// too, so the page says only that the storage could not be opened and that a private window may
// have none (review storage-r2, M1). Chromium's incognito windows (Playwright's contexts are
// incognito) keep the site's storage: no note.
//
// A site whose data is blocked (review of D5, M2): the browser refuses the tab's Web Lock with
// SecurityError before any storage call, and getDirectory() too. Firefox with cookies blocked
// (network.cookie.cookieBehavior 2) runs here as it is; it also refuses `new BroadcastChannel()`
// there, and the page did not start until site-channel.ts. Chromium with cookies blocked answered
// "The request was denied." and "Storage directory access is denied." (a probe, 2026-11-19), but
// Playwright's Chromium (the headless shell) ignores a profile's cookie setting, so the Data
// worker's script is given those two refusals. Playwright's WebKit has no such setting (what
// Safari does when it blocks cookies is not known here); the same refusals stand in. No other tab
// was asked, so the page has no note about other tabs. The note names both causes of a refusal,
// and so what helps: in Chromium it is always blocked site data (review storage-r2, L1).
//
// A storage that fails in a normal window (review storage-r2, M1): a profile on disk whose folder
// for the sites' storage is a file, so the browser cannot make the site's. Each browser fails as
// when that folder cannot be written (a probe, 2026-11-19), and as root too (CI's container).
//
// A search where the site's data is blocked (line WQA's sweep, D8, 2026-11-20): in Firefox every search failed with
// "Run 1 failed: The operation is insecure. A failed run keeps no results." The engine's host opened BroadcastChannels
// on the page and in its workers to link its thread workers and the keeper of their memory to the page, and Firefox
// refuses them there; they are MessageChannels now (packages/losat/src/infra/engine-worker/link.ts).
import { rm, writeFile } from 'node:fs/promises';
import { join } from 'node:path';
import { expect, type Browser, type BrowserContext, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE, openApp } from './support/browser';
import { chooseRun, described, fasta, openFiles, paste, pickRun, program, submit, waitStatus } from './support/search';

/** The note when the browser refuses the site's storage (SecurityError). */
const REFUSED_NOTE =
  "Results are kept in this tab's memory (at most 512\u00a0MB), because the browser does not let this site store files (a private window, or site data blocked in the browser's settings).";
/** The note when the site's storage could not be opened (a lasting UnknownError, or another failure). */
const UNOPENED_NOTE =
  "Results are kept in this tab's memory (at most 512\u00a0MB), because the browser storage of this site could not be opened (a private window may have none).";
/** What each browser said, which the page must not repeat. */
const BROWSER_WORDS =
  /Security error|GetDirectory|transient|out of memory|access is denied|request was denied|not allowed in this context|create actor|state cached|interface object/i;
/** The note of a private window and the name of the browser's error in it. */
const PRIVATE = { firefox: [REFUSED_NOTE, 'SecurityError'], webkit: [UNOPENED_NOTE, 'UnknownError'] } as const;
/** The folder of a profile that holds the sites' storage, and the name of the error when it fails. */
const STORAGE_FOLDER = {
  chromium: ['Default/File System', 'InvalidStateError'],
  firefox: ['storage/default', 'UnknownError'],
  webkit: ['storage', 'UnknownError'],
} as const;

test("a private window keeps results in memory and says why in LOSAT words, not the browser's", async ({
  context,
  browserName,
  playwright,
}) => {
  let privateBrowser: Browser | undefined;
  let window: BrowserContext = context;
  if (browserName === 'firefox') {
    privateBrowser = await playwright.firefox.launch({ firefoxUserPrefs: { 'browser.privatebrowsing.autostart': true } });
    window = await privateBrowser.newContext({ baseURL: test.info().project.use.baseURL! });
  }
  try {
    const page = await openApp(window);
    if (browserName === 'chromium') {
      const status = page.getByTestId('storage-status');
      await expect(status).toHaveAttribute('data-backend', 'opfs');
      await expect(status).not.toHaveAttribute('data-cleanup', 'pending');
      await expect(status.locator('p')).toHaveCount(0);
    } else {
      const [note, errorName] = PRIVATE[browserName as keyof typeof PRIVATE];
      await expectMemoryNote(page, note, errorName);
    }
  } finally {
    await privateBrowser?.close();
  }
});

test('a site whose data the browser blocks keeps results in memory, with no note about other tabs', async ({
  context,
  browserName,
  playwright,
}) => {
  let blockingBrowser: Browser | undefined;
  let window: BrowserContext = context;
  if (browserName === 'firefox') {
    blockingBrowser = await playwright.firefox.launch({ firefoxUserPrefs: { 'network.cookie.cookieBehavior': 2 } });
    window = await blockingBrowser.newContext({ baseURL: test.info().project.use.baseURL! });
  } else {
    await context.route('**/assets/data-worker-*.js', async (route) => {
      const response = await route.fetch();
      const refuse =
        "LockManager.prototype.request = () => Promise.reject(new DOMException('The request was denied.', 'SecurityError'));\n" +
        'StorageManager.prototype.getDirectory = () => ' +
        "Promise.reject(new DOMException('Storage directory access is denied.', 'SecurityError'));\n";
      await route.fulfill({ response, body: refuse + (await response.text()) });
    });
  }
  try {
    await expectMemoryNote(await openApp(window), REFUSED_NOTE, 'SecurityError');
  } finally {
    await blockingBrowser?.close();
  }
});

// A number and its unit stay on one line (line WQA's sweep, D9, 2026-11-19): at 390 px the note broke
// "(at most 512" / "MB)". The size is joined to its unit by a no-break space, so that at no width of a
// phone's screen (the widths below, the note's text reflowed at each) do its two words fall on two lines.
test("the storage note's size and its unit never fall on two lines", async ({ context }) => {
  await context.route('**/assets/data-worker-*.js', async (route) => {
    const response = await route.fetch();
    const refuse =
      "LockManager.prototype.request = () => Promise.reject(new DOMException('The request was denied.', 'SecurityError'));\n" +
      'StorageManager.prototype.getDirectory = () => ' +
      "Promise.reject(new DOMException('Storage directory access is denied.', 'SecurityError'));\n";
    await route.fulfill({ response, body: refuse + (await response.text()) });
  });
  const page = await openApp(context);
  await expectMemoryNote(page, REFUSED_NOTE, 'SecurityError');
  const note = page.getByTestId('storage-memory-note');
  for (let width = 280; width <= 520; width += 4) {
    await page.setViewportSize({ width, height: 844 });
    const lines = await note.evaluate((element) => {
      const walker = document.createTreeWalker(element, NodeFilter.SHOW_TEXT);
      for (let node = walker.nextNode(); node !== null; node = walker.nextNode()) {
        const found = /512[\s\u00a0]MB/.exec(node.textContent ?? '');
        if (found === null) continue;
        const range = document.createRange();
        range.setStart(node, found.index);
        range.setEnd(node, found.index + found[0].length);
        return new Set([...range.getClientRects()].map((rect) => Math.round(rect.top))).size;
      }
      return 0;
    });
    expect(lines, `lines of "512 MB" at ${width} px`).toBe(1);
  }
});

test('a normal window whose storage fails keeps results in memory, and is not told it is a private window', async ({
  browserName,
  playwright,
}) => {
  const [folder, errorName] = STORAGE_FOLDER[browserName];
  const dir = test.info().outputPath('profile');
  // The project's launch options (WebKit's launcher) come with `playwright`, as in support/profile.ts.
  const options = { baseURL: test.info().project.use.baseURL! };
  await (await playwright[browserName].launchPersistentContext(dir, options)).close();
  await rm(join(dir, folder), { recursive: true, force: true });
  await writeFile(join(dir, folder), '');
  const window = await playwright[browserName].launchPersistentContext(dir, options);
  try {
    await expectMemoryNote(await openApp(window), UNOPENED_NOTE, errorName);
  } finally {
    await window.close();
  }
});

/**
 * What Firefox refuses where the site's data is blocked, on the page and in every worker, in its words (a probe,
 * 2026-11-20: losat-baselines/web-qa-20261119/fix-d8/probe/probe.log): BroadcastChannel, the Web Locks and the storage.
 * Chromium with cookies blocked refuses the last two and opens BroadcastChannels, but Playwright's Chromium (the headless
 * shell) ignores a profile's cookie setting, and Playwright's WebKit has no such setting: both get Firefox's refusals.
 */
const REFUSE_SITE_DATA = [
  "globalThis.BroadcastChannel = class { constructor() { throw new DOMException('The operation is insecure.', 'SecurityError'); } };",
  "if (typeof LockManager === 'function') LockManager.prototype.request = () => Promise.reject(new DOMException('LockManager.request: request() is not allowed in this context', 'SecurityError'));",
  "if (typeof StorageManager === 'function') StorageManager.prototype.getDirectory = () => Promise.reject(new DOMException('Security error when calling GetDirectory', 'SecurityError'));",
  '',
].join('\n');

/** Gives the page and the scripts of its workers Firefox's refusals of a site whose data is blocked. */
async function refuseSiteData(context: BrowserContext): Promise<void> {
  await context.addInitScript(REFUSE_SITE_DATA);
  await context.route(/\/assets\/[^/]*worker[^/]*\.js$/, async (route) => {
    const response = await route.fetch();
    await route.fulfill({ response, body: REFUSE_SITE_DATA + (await response.text()) });
  });
}

/** Counts the Engine workers that the page ends itself (`endedEngines.count`). */
async function countEndedEngines(context: BrowserContext): Promise<void> {
  await context.addInitScript(() => {
    const Native = window.Worker;
    const ended = { count: 0 };
    Object.assign(window, { endedEngines: ended });
    window.Worker = class extends Native {
      private readonly engine: boolean;
      constructor(url: string | URL, options?: WorkerOptions) {
        super(url, options);
        this.engine = options?.name === 'losat-engine';
      }
      override terminate(): void {
        if (this.engine) ended.count++;
        super.terminate();
      }
    };
  });
}

const PROTEIN = 'MKLVVLAAGGHHKLMKLVVLAAGGHHKLMSTNPKPQRKTKRNTNRRPQDVKFPGG';

test('a site whose data the browser blocks runs a threaded search to its results, and a cancelled one ends', async ({
  context,
  browserName,
  playwright,
}) => {
  test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine build has no Engine worker: the engine build runs this test');
  test.setTimeout(240_000);
  let blockingBrowser: Browser | undefined;
  let blocked: BrowserContext = context;
  if (browserName === 'firefox') {
    blockingBrowser = await playwright.firefox.launch({ firefoxUserPrefs: { 'network.cookie.cookieBehavior': 2 } });
    blocked = await blockingBrowser.newContext({ baseURL: test.info().project.use.baseURL! });
  } else {
    await refuseSiteData(context);
  }
  try {
    await countEndedEngines(blocked);
    const page = await openApp(blocked);
    await expectMemoryNote(page, REFUSED_NOTE, 'SecurityError');
    await described(page);
    // Two threads: the Engine worker links a thread worker and the keeper of their memory to the page.
    await program(page, 'blastp');
    await paste(page, 'query', `>q1\n${PROTEIN}\n`);
    await paste(page, 'subject', `>s1\n${PROTEIN}\n`);
    await page.getByTestId('threads').selectOption('2');
    await submit(page);
    await waitStatus(page, 1, 'completed');
    await chooseRun(page, 1);
    await expect(page.getByTestId('results-hits')).toContainText('s1');
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-record').locator('[data-detail="path"]')).toHaveText(/^\s*threaded, 2 threads\s*$/);

    // A threaded search cancelled while the Engine worker runs it: the page asks the keeper on its link to grow the
    // memory, and ends the worker on the keeper's answer.
    await page.getByTestId('tab-search').click();
    await program(page, 'tblastx');
    await openFiles(page, 'query', [{ name: 'AP027078.fasta', text: fasta('AP027078.fasta') }]);
    await openFiles(page, 'subject', [{ name: 'NZ_CP006932.fasta', text: fasta('NZ_CP006932.fasta') }]);
    await page.getByTestId('threads').selectOption('2');
    await submit(page);
    await expect(page.getByTestId('run-2-option')).toHaveAttribute('data-status', 'running', { timeout: 60_000 });
    const ended = () => page.evaluate(() => (window as unknown as { endedEngines: { count: number } }).endedEngines.count);
    const before = await ended();
    await pickRun(page, 2);
    await page.getByTestId('run-status-cancel').click();
    await waitStatus(page, 2, 'cancelled', 15_000);
    await expect.poll(ended, { timeout: 15_000 }).toBe(before + 1);
  } finally {
    await blockingBrowser?.close();
  }
});

/** The page keeps results in memory, says why in LOSAT's words, and, once the start-up cleanup is over, has no other note. */
async function expectMemoryNote(page: Page, text: string, errorName: string): Promise<void> {
  const status = page.getByTestId('storage-status');
  await expect(status).toHaveAttribute('data-backend', 'memory');
  await expect(status).not.toHaveAttribute('data-cleanup', 'pending');
  await expect(status.getByTestId('storage-cleanup-note')).toHaveCount(0);
  const note = status.getByTestId('storage-memory-note');
  await expect(note).toHaveText(text);
  await expect(note).toHaveAttribute('data-error', new RegExp(`^${errorName}: .`));
  expect(await status.innerText()).not.toMatch(BROWSER_WORDS);
}
