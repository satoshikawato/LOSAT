// Saving a file in the E2E tests: a click on a button that saves, and the browser's download.
import type { Download, Page } from '@playwright/test';

/**
 * The least time between two downloads of one page. Chromium's renderer drops a download without a
 * word when its frame has started ten others and the first of them is less than a second old (a
 * burst limit in Blink; a click before each download does not lift it, and Firefox has none): in
 * fix round 5 a test that had saved ten files in 0.8 s never saw the eleventh, although the page
 * said it was saved. A person does not save eleven files in a second; a test that saves file after
 * file does. With this time between downloads, eleven of them take 1.5 s.
 */
export const DOWNLOAD_SPACING_MS = 150;

/** When each page's last download arrived, on a monotonic clock: the wall clock can jump (it did, by 1.6 s, while this was measured). */
const lastDownload = new WeakMap<Page, number>();

/**
 * Clicks the button `testid` and returns the download that it starts, at least DOWNLOAD_SPACING_MS
 * after the page's previous download arrived.
 */
export async function download(page: Page, testid: string): Promise<Download> {
  const wait = (lastDownload.get(page) ?? -Infinity) + DOWNLOAD_SPACING_MS - performance.now();
  if (wait > 0) await page.waitForTimeout(wait);
  const event = page.waitForEvent('download');
  await page.getByTestId(testid).click();
  const file = await event;
  lastDownload.set(page, performance.now());
  return file;
}
