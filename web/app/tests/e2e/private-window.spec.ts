// A private window, where the browser does not let the site store files (line WQA's sweep, D5,
// 2026-11-19). The application keeps results in the tab's memory and says why in LOSAT's words:
// the note quoted the browser's exception, "(Security error when calling GetDirectory)" in Firefox
// and "(The operation failed for an unknown transient reason (e.g. out of memory).)" in WebKit,
// which reads as if memory ran out. The name of the browser's error stays in the note's
// data-error attribute, for developers.
//
// Firefox: a private window (browser.privatebrowsing.autostart). WebKit: Playwright's context is
// an ephemeral session, as Safari's private windows are, and its getDirectory() rejects. Chromium's
// incognito windows keep OPFS, so there the Data worker's getDirectory() rejects as Chromium's does
// for a site that may not store data.
import { expect, type Browser, type BrowserContext } from '@playwright/test';
import { test } from './support/test';
import { openApp } from './support/browser';

const NOTE =
  "Results are kept in this tab's memory only, because the browser does not let this site store files here, as in a private window.";
/** What each browser said, which the page must not repeat. */
const BROWSER_WORDS = /Security error|GetDirectory|transient|out of memory|access is denied/i;
const ERROR_NAME = { chromium: 'SecurityError', firefox: 'SecurityError', webkit: 'UnknownError' } as const;

test('a private window keeps results in memory and says why in LOSAT words, not the browser\'s', async ({
  context,
  browserName,
  playwright,
}) => {
  let privateBrowser: Browser | undefined;
  let window: BrowserContext = context;
  if (browserName === 'firefox') {
    privateBrowser = await playwright.firefox.launch({ firefoxUserPrefs: { 'browser.privatebrowsing.autostart': true } });
    window = await privateBrowser.newContext({ baseURL: test.info().project.use.baseURL! });
  } else if (browserName === 'chromium') {
    await context.route('**/assets/data-worker-*.js', async (route) => {
      const response = await route.fetch();
      const refuse =
        'StorageManager.prototype.getDirectory = () => ' +
        "Promise.reject(new DOMException('Storage directory access is denied.', 'SecurityError'));\n";
      await route.fulfill({ response, body: refuse + (await response.text()) });
    });
  }
  try {
    const page = await openApp(window);
    const status = page.getByTestId('storage-status');
    await expect(status).toHaveAttribute('data-backend', 'memory');
    const note = status.getByTestId('storage-memory-note');
    await expect(note).toHaveText(NOTE);
    await expect(note).toHaveAttribute('data-error', ERROR_NAME[browserName as keyof typeof ERROR_NAME]);
    expect(await status.innerText()).not.toMatch(BROWSER_WORDS);
  } finally {
    await privateBrowser?.close();
  }
});
