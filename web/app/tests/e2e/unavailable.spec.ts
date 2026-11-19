// A part of LOSAT that could not be loaded (line WQA's sweep, D1, 2026-11-19;
// losat-baselines/web-qa-20261119/sweep/FINDINGS.md).
// D1: when the engine's serial module or the Data worker's script could not be loaded, every input card said "This input
// cannot be read: Failed to fetch" (the browser's words; Firefox and WebKit have others), Run LOSAT was refused for the
// pasted query, and nothing said for the page that LOSAT could not load a part of itself; a search whose threaded module
// or Engine worker could not be loaded failed with "Failed to fetch." or "The engine worker stopped: an error without a
// message". Now the page says once which part, why and the remedy, and the cards do not blame the inputs.
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { openSearch, paste, runMessage, submit, waitStatus } from './support/search';

const ENGINE_ONLY = 'the FakeEngine build has no engine modules nor Engine worker: the engine build runs this test';
const RELOAD = 'Reload the page.';
const DATA_NOT_STARTED = `LOSAT could not start the part that keeps the results. ${RELOAD}`;
const QUERY = '>q1\nACGTACGTACGTAGCTAGCTAGCATCGATCGACTAGCTAGCATCGACTAGC\n';
const SUBJECT = '>s1\nACGTACGTACGTAGCTAGCTAGCATCGATCGACTAGCTAGCATCGACTAGC\n';

/** Makes every request for a path that `pattern` matches fail: as a network error, or with an HTTP status. */
async function failRequests(page: Page, pattern: RegExp, how: 'network error' | number): Promise<void> {
  await page.context().route(
    (url) => pattern.test(url.pathname),
    (route) => (how === 'network error' ? route.abort('failed') : route.fulfill({ status: how, contentType: 'text/plain', body: 'error' })),
  );
}

/** Pastes both inputs and waits until each is read or not. */
async function pasteBoth(page: Page): Promise<void> {
  await paste(page, 'query', QUERY);
  await paste(page, 'subject', SUBJECT);
}

/** The cards do not blame the inputs: each says that it was not read. */
async function notRead(page: Page): Promise<void> {
  for (const role of ['query', 'subject'] as const) {
    await expect(page.getByTestId(`${role}-source-0-unread`)).toHaveText('Not read.');
    await expect(page.getByTestId(`${role}-source-0`)).not.toContainText('cannot be read');
  }
}

/** Run LOSAT is refused at once with `sentence`, not left "Starting…". */
async function refused(page: Page, sentence: string): Promise<void> {
  await page.getByTestId('run-losat').click();
  await expect(page.getByTestId('search-message')).toHaveText(sentence, { timeout: 15_000 });
  await expect(page.getByTestId('run-losat')).toHaveText('Run LOSAT');
}

for (const how of ['network error', 500] as const) {
  const why = how === 'network error' ? 'network error' : `HTTP ${how}`;
  test(`the page says that LOSAT could not download its engine file (${why}), and no input is blamed`, async ({ page }) => {
    test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
    await failRequests(page, /^\/engine\/losat-web-serial-[^/]*\.wasm$/, how);
    await openSearch(page);
    const sentence = new RegExp(
      `^LOSAT could not download its engine file /engine/losat-web-serial-[0-9a-f]+\\.wasm \\(${why}\\)\\. Reload the page\\.$`,
    );
    await expect(page.getByTestId('losat-unavailable')).toHaveText(sentence);
    await expect(page.getByTestId('losat-unavailable')).toHaveAttribute('role', 'alert');
    await pasteBoth(page);
    await notRead(page);
    const line = (await page.getByTestId('losat-unavailable').textContent())!.trim();
    await refused(page, line);
  });
}

test('the page says that LOSAT could not start the part that keeps the results when its script cannot be loaded', async ({ page }) => {
  await failRequests(page, /^\/assets\/data-worker-[^/]*\.js$/, 'network error');
  await openSearch(page);
  await expect(page.getByTestId('losat-unavailable')).toHaveText(DATA_NOT_STARTED);
  await pasteBoth(page);
  await notRead(page);
  await refused(page, DATA_NOT_STARTED);
});

test('a search whose threaded module cannot be downloaded fails with the file and the remedy', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  await failRequests(page, /^\/engine\/losat-web-threads-[^/]*\.wasm$/, 'network error');
  await openSearch(page);
  await pasteBoth(page);
  await submit(page);
  await waitStatus(page, 1, 'failed', 30_000);
  expect(await runMessage(page, 1)).toMatch(
    /^Run 1 failed: LOSAT could not download its engine file \/engine\/losat-web-threads-[0-9a-f]+\.wasm \(network error\)\. Reload the page\. A failed run keeps no results\.$/,
  );
});

test('a search whose Engine worker cannot be loaded fails with what could not start and the remedy', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  await failRequests(page, /^\/assets\/engine-worker-[^/]*\.js$/, 'network error');
  await openSearch(page);
  await pasteBoth(page);
  await submit(page);
  await waitStatus(page, 1, 'failed', 30_000);
  expect(await runMessage(page, 1)).toBe(
    `Run 1 failed: LOSAT could not start the part that runs searches. ${RELOAD} A failed run keeps no results.`,
  );
});
