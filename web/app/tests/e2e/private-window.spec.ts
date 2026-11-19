// A private window, or a site whose data the browser blocks: the browser does not let the site
// store files (line WQA's sweep, D5, 2026-11-19). The application keeps results in the tab's
// memory and says why in LOSAT's words: the note quoted the browser's exception, "(Security error
// when calling GetDirectory)" in Firefox and "(The operation failed for an unknown transient
// reason (e.g. out of memory).)" in WebKit, which reads as if memory ran out. The browser's error,
// name and message, stays in the note's data-error attribute, for developers. The note says what
// differs from a tab that stores files: the results are held in memory, which keeps at most
// 512 MB of them. Losing them on a reload is not one of the differences, since LOSAT Web keeps
// nothing for the next visit in any browser (review of D5, M1).
//
// Private windows. Firefox: browser.privatebrowsing.autostart. WebKit: Playwright's context is an
// ephemeral session, as Safari's private windows are, and its getDirectory() rejects. Chromium's
// incognito windows (Playwright's contexts are incognito) keep the site's storage: no note.
//
// A site whose data is blocked (review of D5, M2): the browser refuses the tab's Web Lock with
// SecurityError before any storage call, and getDirectory() too. Firefox with cookies blocked
// (network.cookie.cookieBehavior 2) runs here as it is; it also refuses `new BroadcastChannel()`
// there, and the page did not start until site-channel.ts. Chromium with cookies blocked answered
// "The request was denied." and "Storage directory access is denied." (a probe, 2026-11-19), but
// Playwright's Chromium (the headless shell) ignores a profile's cookie setting, so the Data
// worker's script is given those two refusals. Playwright's WebKit has no such setting (what
// Safari does when it blocks cookies is not known here); the same refusals stand in. No other tab
// was asked, so the page has no note about other tabs.
import { expect, type Browser, type BrowserContext, type Page } from '@playwright/test';
import { test } from './support/test';
import { openApp } from './support/browser';

const NOTE =
  "Results are kept in this tab's memory (at most 512 MB), because the browser does not let this site store files here, as in a private window.";
/** What each browser said, which the page must not repeat. */
const BROWSER_WORDS = /Security error|GetDirectory|transient|out of memory|access is denied|request was denied|not allowed in this context/i;
const PRIVATE_ERROR = { firefox: 'SecurityError', webkit: 'UnknownError' } as const;

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
      await expectMemoryNote(page, PRIVATE_ERROR[browserName as keyof typeof PRIVATE_ERROR]);
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
    await expectMemoryNote(await openApp(window), 'SecurityError');
  } finally {
    await blockingBrowser?.close();
  }
});

/** The page keeps results in memory, says why in LOSAT's words, and, once the start-up cleanup is over, has no other note. */
async function expectMemoryNote(page: Page, errorName: string): Promise<void> {
  const status = page.getByTestId('storage-status');
  await expect(status).toHaveAttribute('data-backend', 'memory');
  await expect(status).not.toHaveAttribute('data-cleanup', 'pending');
  await expect(status.getByTestId('storage-cleanup-note')).toHaveCount(0);
  const note = status.getByTestId('storage-memory-note');
  await expect(note).toHaveText(NOTE);
  await expect(note).toHaveAttribute('data-error', new RegExp(`^${errorName}: .`));
  expect(await status.innerText()).not.toMatch(BROWSER_WORDS);
}
