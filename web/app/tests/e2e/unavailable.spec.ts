// A part of LOSAT that could not be loaded, or that ended (line WQA's sweep, D1 and D2, 2026-11-19;
// losat-baselines/web-qa-20261119/sweep/FINDINGS.md).
// D1: when the engine's serial module or the Data worker's script could not be loaded, every input card said "This input
// cannot be read: Failed to fetch" (the browser's words; Firefox and WebKit have others), Run LOSAT was refused for the
// pasted query, and nothing said for the page that LOSAT could not load a part of itself; a search whose threaded module
// or Engine worker could not be loaded failed with "Failed to fetch." or "The engine worker stopped: an error without a
// message". Now the page says once which part, why and the remedy, and the cards do not blame the inputs.
// D2: a worker that ends without an error event, as one that the browser ends does, was never noticed: a run stayed
// "Organizing results" (the Data worker) or "Searching" (the Engine worker) for good, Run LOSAT "Starting…", and a run
// stuck while organizing its results had no Cancel. The tests end a worker with Worker.terminate(), which gives the
// page no event, as the browser's end of a worker does.
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { fasta, openFiles, openSearch, paste, program, runMessage, submit, waitStatus } from './support/search';

const ENGINE_ONLY = 'the FakeEngine build has no engine modules nor Engine worker: the engine build runs this test';
const RELOAD = 'Reload the page.';
const DATA_NOT_STARTED = `LOSAT could not start the part that keeps the results. ${RELOAD}`;
const DATA_ENDED = `The part of LOSAT that keeps the results stopped. ${RELOAD}`;
const QUERY = '>q1\nACGTACGTACGTAGCTAGCTAGCATCGATCGACTAGCTAGCATCGACTAGC\n';
const SUBJECT = '>s1\nACGTACGTACGTAGCTAGCTAGCATCGATCGACTAGCTAGCATCGACTAGC\n';

/** Keeps the workers that the page starts, by their names ('losat-data', 'losat-engine'), so that a test can end one. */
async function keepWorkers(page: Page): Promise<void> {
  await page.addInitScript(() => {
    const Native = window.Worker;
    const kept: Record<string, Worker[]> = {};
    (window as unknown as { keptWorkers: Record<string, Worker[]> }).keptWorkers = kept;
    window.Worker = class extends Native {
      constructor(url: string | URL, options?: WorkerOptions) {
        super(url, options);
        (kept[options?.name ?? ''] ??= []).push(this);
      }
    };
  });
}

/** Ends the latest worker of `name` from outside, without an error event, as a browser that ends it does. */
async function endWorker(page: Page, name: 'losat-data' | 'losat-engine'): Promise<void> {
  const ended = await page.evaluate((which) => {
    const worker = (window as unknown as { keptWorkers: Record<string, Worker[]> }).keptWorkers[which]?.at(-1);
    worker?.terminate();
    return worker !== undefined;
  }, name);
  expect(ended, `a worker named ${name}`).toBe(true);
}

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

test('the page says when the part that keeps the results ended, and nothing waits for it for good', async ({ page }) => {
  await keepWorkers(page);
  await openSearch(page);
  await pasteBoth(page);
  await endWorker(page, 'losat-data');
  await expect(page.getByTestId('losat-unavailable')).toHaveText(DATA_ENDED, { timeout: 15_000 });
  // Before, Run LOSAT stayed "Starting…" and a new text "Reading the records…" for good.
  await refused(page, DATA_ENDED);
  await page.getByTestId('query-input').fill(`${QUERY}>q2\nACGT\n`);
  await expect(page.getByTestId('query-source-0-unread')).toHaveText('Not read.', { timeout: 15_000 });
});

/** A search of several seconds (TBLASTX of two genomes, 2 threads): 5 s in Chromium on the sweep's machine. */
async function startLongSearch(page: Page): Promise<void> {
  await openSearch(page);
  await program(page, 'tblastx');
  await openFiles(page, 'query', [{ name: 'AP027078.fasta', text: fasta('AP027078.fasta') }]);
  await openFiles(page, 'subject', [{ name: 'NZ_CP006932.fasta', text: fasta('NZ_CP006932.fasta') }]);
  await page.getByTestId('threads').selectOption('2');
  await submit(page);
  await expect(page.getByTestId('run-1-option')).toHaveAttribute('data-status', 'running', { timeout: 60_000 });
}

test('a search whose Data worker ended fails with what stopped, instead of organizing its results for good', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  test.setTimeout(180_000);
  await keepWorkers(page);
  await startLongSearch(page);
  await endWorker(page, 'losat-data');
  await waitStatus(page, 1, 'failed', 30_000);
  expect(await runMessage(page, 1)).toBe(`Run 1 failed: ${DATA_ENDED} A failed run keeps no results.`);
  await expect(page.getByTestId('losat-unavailable')).toHaveText(DATA_ENDED);
});

test('a search whose Engine worker ended fails with what stopped, and the next search runs', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  test.setTimeout(180_000);
  await keepWorkers(page);
  await startLongSearch(page);
  await endWorker(page, 'losat-engine');
  await waitStatus(page, 1, 'failed', 30_000);
  expect(await runMessage(page, 1)).toBe('Run 1 failed: The part of LOSAT that runs searches stopped. A failed run keeps no results.');
  await expect(page.getByTestId('losat-unavailable')).toHaveCount(0);
  await page.getByTestId('tab-search').click();
  await submit(page);
  await waitStatus(page, 2, 'completed', 120_000);
});

test('the next search runs after the Engine worker ended between searches', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  test.setTimeout(120_000);
  await keepWorkers(page);
  await openSearch(page);
  await pasteBoth(page);
  await submit(page);
  await waitStatus(page, 1, 'completed', 60_000);
  await endWorker(page, 'losat-engine');
  // Before, run 2 stayed "Preparing" for good.
  await page.getByTestId('tab-search').click();
  await submit(page);
  await waitStatus(page, 2, 'completed', 30_000);
});

test('Cancel ends a search whose results are not organized', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  test.setTimeout(180_000);
  // Without Web Locks on the page, the page cannot learn that the Data worker ended (its workers keep theirs): the run
  // waits for the commit of its results, as it does for a Data worker that does not answer.
  await page.addInitScript(() => Object.defineProperty(Navigator.prototype, 'locks', { configurable: true, get: () => undefined }));
  await keepWorkers(page);
  await startLongSearch(page);
  await endWorker(page, 'losat-data');
  await expect(page.getByTestId('run-1-option')).toHaveAttribute('data-status', 'finalizing', { timeout: 60_000 });
  // Before, a run that organizes its results had no Cancel.
  await expect(page.getByTestId('run-status-cancel')).toBeVisible({ timeout: 5_000 });
  await page.getByTestId('run-status-cancel').click();
  await waitStatus(page, 1, 'cancelled', 15_000);
});
