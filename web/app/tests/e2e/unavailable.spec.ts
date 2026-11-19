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
// The review of D1 and D2 (2026-11-19): the first search of a new Engine worker was not watched (H1), and an error
// during it read as a worker that could not start (M1); an Example with LOSAT unavailable downloaded its whole set and
// then blamed the file (M2); "Task Selection" and "Algorithm parameters" showed controls with no choices or values
// (L1); a search went on after the Data worker ended (L2). Its second review: one failed download of the engine file
// was said for the whole session (M1); a threaded runtime whose worker ended was never finished in Firefox (L1); "The part
// of LOSAT that runs searches stopped." did not say what to do (L2).
import { gzipSync } from 'node:zlib';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { exampleById } from '../../src/domain/examples';
import { BUILD_HAS_ENGINE } from './support/browser';
import { E2E_ORIGIN } from './support/origin';
import { fasta, openFiles, openParameters, openSearch, paste, pickRun, program, runMessage, submit, waitStatus } from './support/search';

const ENGINE_ONLY = 'the FakeEngine build has no engine modules nor Engine worker: the engine build runs this test';
const RELOAD = 'Reload the page.';
const DATA_NOT_STARTED = `LOSAT could not start the part that keeps the results. ${RELOAD}`;
const DATA_ENDED = `The part of LOSAT that keeps the results stopped. ${RELOAD}`;
// The next search starts another Engine worker (the second review of D1 and D2, L2).
const ENGINE_ENDED = 'The part of LOSAT that runs searches stopped. Run the search again.';
const QUERY = '>q1\nACGTACGTACGTAGCTAGCTAGCATCGATCGACTAGCTAGCATCGACTAGC\n';
const SUBJECT = '>s1\nACGTACGTACGTAGCTAGCTAGCATCGATCGACTAGCTAGCATCGACTAGC\n';

/**
 * Keeps the workers that the page starts, by their names ('losat-data', 'losat-engine'), so that a test can end one,
 * and counts those that the page ended itself (`endedWorkers`).
 */
async function keepWorkers(page: Page): Promise<void> {
  await page.addInitScript(() => {
    const Native = window.Worker;
    const kept: Record<string, Worker[]> = {};
    const ended: Record<string, number> = {};
    Object.assign(window, { keptWorkers: kept, endedWorkers: ended });
    window.Worker = class extends Native {
      private readonly kind: string;
      constructor(url: string | URL, options?: WorkerOptions) {
        super(url, options);
        this.kind = options?.name ?? '';
        (kept[this.kind] ??= []).push(this);
      }
      override terminate(): void {
        ended[this.kind] = (ended[this.kind] ?? 0) + 1;
        super.terminate();
      }
    };
  });
}

/** Counts the page's open BroadcastChannels by name (`openChannels`). */
async function countChannels(page: Page): Promise<void> {
  await page.addInitScript(() => {
    const Native = window.BroadcastChannel;
    const open: Record<string, number> = {};
    Object.assign(window, { openChannels: open });
    window.BroadcastChannel = class extends Native {
      constructor(name: string) {
        super(name);
        open[name] = (open[name] ?? 0) + 1;
      }
      override close(): void {
        open[this.name] = (open[this.name] ?? 0) - 1;
        super.close();
      }
    };
  });
}

/** How many of the page's channels to the keepers of the engine's shared memory are open. */
async function openKeeperChannels(page: Page): Promise<number> {
  return page.evaluate(() =>
    Object.entries((window as unknown as { openChannels: Record<string, number> }).openChannels)
      .filter(([name]) => name.startsWith('losat-web:engine-keeper:'))
      .reduce((sum, [, count]) => sum + count, 0),
  );
}

/** How many workers of `name` the page ended itself. */
async function endedByPage(page: Page, name: 'losat-data' | 'losat-engine'): Promise<number> {
  return page.evaluate((which) => (window as unknown as { endedWorkers: Record<string, number> }).endedWorkers[which] ?? 0, name);
}

/**
 * Ends the latest worker of `name` from outside, without an error event, as a browser that ends it does; or, with
 * `error`, gives the page the error event of a worker that died of an uncaught error with that message.
 */
async function endWorker(page: Page, name: 'losat-data' | 'losat-engine', error?: string): Promise<void> {
  const ended = await page.evaluate(
    ([which, message]) => {
      const worker = (window as unknown as { keptWorkers: Record<string, Worker[]> }).keptWorkers[which]?.at(-1);
      if (worker === undefined) return false;
      // Not counted as ended by the page: Native's terminate, as the browser ends it.
      if (message === undefined) Object.getPrototypeOf(Object.getPrototypeOf(worker)).terminate.call(worker);
      else worker.dispatchEvent(new ErrorEvent('error', { message, cancelable: true }));
      return true;
    },
    [name, error] as const,
  );
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

const SIXTEEN_S = exampleById('refseq-16s-bacteria');

/** Answers the download of the 16S Example set with a small set; returns the requests that leave the site. */
async function stubExampleHost(page: Page): Promise<string[]> {
  const away: string[] = [];
  page.on('request', (request) => {
    if (/^https?:/.test(request.url()) && !request.url().startsWith(E2E_ORIGIN)) away.push(request.url());
  });
  // A build without the engine has no decoder, and gets the set decoded on the way (examples.spec.ts).
  const body = BUILD_HAS_ENGINE ? gzipSync(SUBJECT) : Buffer.from(SUBJECT);
  await page.route(
    (url) => url.href === SIXTEEN_S.url,
    (route) => route.fulfill({ status: 200, contentType: 'application/x-gzip', headers: { 'Access-Control-Allow-Origin': '*' }, body }),
  );
  return away;
}

/**
 * An Example is refused at once with `sentence`, as Run LOSAT is, and downloads nothing. Before, the set was downloaded
 * whole, and then the file was blamed: "The file from ftp.ncbi.nlm.nih.gov could not be decompressed (LOSAT could not
 * download its engine file …)" (the review of D1, M2).
 */
async function exampleRefused(page: Page, sentence: string, away: readonly string[]): Promise<void> {
  await page.getByTestId('example-refseq-16s-bacteria').click();
  await expect(page.getByTestId('example-error')).toHaveText(sentence, { timeout: 15_000 });
  expect(away).toEqual([]);
}

/**
 * The options were not read: "Task Selection" and "Algorithm parameters" say so, as the cards do, with no control.
 * Before, "Optimize for -task" had no choices, and every number field was blank (the review of D1, L1).
 */
async function optionsNotRead(page: Page): Promise<void> {
  const tasks = page.getByTestId('program-selection');
  await expect(tasks.getByTestId('param-task-unread')).toHaveText('Not read.');
  await expect(tasks.locator('input, select')).toHaveCount(0);
  await openParameters(page);
  await expect(page.getByTestId('parameters-unread')).toHaveText('Not read.');
  await expect(page.getByTestId('parameter-form').locator('input, select')).toHaveCount(0);
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
    const away = await stubExampleHost(page);
    await pasteBoth(page);
    await notRead(page);
    const line = (await page.getByTestId('losat-unavailable').textContent())!.trim();
    await refused(page, line);
    await exampleRefused(page, line, away);
  });
}

test('once the engine file downloads after a failed download, the page no longer says that LOSAT is unavailable', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  test.setTimeout(120_000);
  // The first request for the serial module fails, as a network glitch at load does; the next ones are answered.
  let requests = 0;
  await page.context().route(
    (url) => /^\/engine\/losat-web-serial-[^/]*\.wasm$/.test(url.pathname),
    (route) => (++requests === 1 ? route.abort('failed') : route.continue()),
  );
  const away = await stubExampleHost(page);
  await openSearch(page);
  await pasteBoth(page);
  // Before, the page said for the whole session that LOSAT could not download its engine file, and to reload it,
  // losing the results; the options said "Not read." and the Examples were refused.
  await expect(page.getByTestId('losat-unavailable')).toHaveCount(0, { timeout: 15_000 });
  expect(requests).toBeGreaterThan(1);
  await expect(page.getByTestId('program-selection').locator('input')).not.toHaveCount(0);
  await submit(page);
  await waitStatus(page, 1, 'completed', 60_000);
  await page.getByTestId('tab-search').click();
  await page.getByTestId('example-refseq-16s-bacteria').click();
  await expect(page.getByTestId('subject-source-0')).toContainText('bacteria.16SrRNA.fna', { timeout: 30_000 });
  await expect(page.getByTestId('example-error')).toHaveCount(0);
  expect(away).toEqual([SIXTEEN_S.url]);
  await expect(page.getByTestId('losat-unavailable')).toHaveCount(0);
});

test('the options say that they were not read when the engine file cannot be downloaded', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  await failRequests(page, /^\/engine\/losat-web-serial-[^/]*\.wasm$/, 'network error');
  await openSearch(page);
  await expect(page.getByTestId('losat-unavailable')).toBeVisible();
  await optionsNotRead(page);
});

test('the page says that LOSAT could not start the part that keeps the results when its script cannot be loaded', async ({ page }) => {
  await failRequests(page, /^\/assets\/data-worker-[^/]*\.js$/, 'network error');
  const away = await stubExampleHost(page);
  await openSearch(page);
  await expect(page.getByTestId('losat-unavailable')).toHaveText(DATA_NOT_STARTED);
  await pasteBoth(page);
  await notRead(page);
  await refused(page, DATA_NOT_STARTED);
  await exampleRefused(page, DATA_NOT_STARTED, away);
  // The FakeEngine describes the options itself; the engine build asks the Data worker.
  if (BUILD_HAS_ENGINE) await optionsNotRead(page);
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

/** The search form with a search of several seconds (TBLASTX of two genomes): 5 s with 2 threads in Chromium on the sweep's machine. */
async function longSearchForm(page: Page, threads: '1' | '2'): Promise<void> {
  await openSearch(page);
  await program(page, 'tblastx');
  await openFiles(page, 'query', [{ name: 'AP027078.fasta', text: fasta('AP027078.fasta') }]);
  await openFiles(page, 'subject', [{ name: 'NZ_CP006932.fasta', text: fasta('NZ_CP006932.fasta') }]);
  await page.getByTestId('threads').selectOption(threads);
}

/** Starts that search with 2 threads. */
async function startLongSearch(page: Page): Promise<void> {
  await longSearchForm(page, '2');
  await submit(page);
  await expect(page.getByTestId('run-1-option')).toHaveAttribute('data-status', 'running', { timeout: 60_000 });
}

test('a search whose Data worker ended fails at once with what stopped, and its search ends', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  test.setTimeout(180_000);
  await keepWorkers(page);
  await startLongSearch(page);
  const engines = await endedByPage(page, 'losat-engine');
  await endWorker(page, 'losat-data');
  // The page ends the search's Engine worker. Before, the search went on, and the run failed only when it ended, at the
  // commit of its results, while the page already said that the part that keeps the results stopped (the review of D2, L2).
  await expect.poll(() => endedByPage(page, 'losat-engine'), { timeout: 15_000 }).toBe(engines + 1);
  await waitStatus(page, 1, 'failed', 30_000);
  expect(await runMessage(page, 1)).toBe(`Run 1 failed: ${DATA_ENDED} A failed run keeps no results.`);
  await expect(page.getByTestId('losat-unavailable')).toHaveText(DATA_ENDED);
});

test('a search whose Engine worker ended fails with what stopped, and the next search runs', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  test.setTimeout(180_000);
  await keepWorkers(page);
  await countChannels(page);
  await startLongSearch(page);
  await endWorker(page, 'losat-engine');
  await waitStatus(page, 1, 'failed', 30_000);
  expect(await runMessage(page, 1)).toBe(`Run 1 failed: ${ENGINE_ENDED} A failed run keeps no results.`);
  // The worker's threads and its keeper ended with it: the page ends its runtime at once. Before, the page waited for the
  // keeper's answer for good, with the keeper's channel open (here in Chromium, Firefox and WebKit).
  await expect.poll(() => openKeeperChannels(page), { timeout: 15_000 }).toBe(0);
  await expect(page.getByTestId('losat-unavailable')).toHaveCount(0);
  await page.getByTestId('tab-search').click();
  await submit(page);
  await waitStatus(page, 2, 'completed', 120_000);
});

// The first search of a new Engine worker is posted together with its `init`: after a cancel (here), before the warm-up,
// after a renewal. Before, the worker said that it started only once it held its lock, which came after that search had
// started, so the page did not watch it: a worker ended during that search left the run "Searching" for good (the review
// of D2, H1), and an error then read "LOSAT could not start the part that runs searches. Reload the page." (M1).
for (const [end, sentence] of [
  ['ends', ENGINE_ENDED],
  ['fails', 'The engine worker stopped: out of memory.'],
] as const) {
  test(`the first search of a new Engine worker fails with what stopped when the worker ${end} during it`, async ({ page }) => {
    test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
    test.setTimeout(180_000);
    await keepWorkers(page);
    // 1 thread: the search runs in the Engine worker itself, which takes no message until it ends.
    await longSearchForm(page, '1');
    await submit(page);
    await expect(page.getByTestId('run-1-option')).toBeAttached({ timeout: 30_000 });
    await page.getByTestId('tab-search').click();
    await submit(page);
    await expect(page.getByTestId('run-2-option')).toHaveAttribute('data-status', 'queued', { timeout: 30_000 });
    await expect(page.getByTestId('run-1-option')).toHaveAttribute('data-status', 'running', { timeout: 60_000 });
    await pickRun(page, 1);
    await page.getByTestId('run-status-cancel').click();
    await waitStatus(page, 1, 'cancelled', 15_000);
    await expect(page.getByTestId('run-2-option')).toHaveAttribute('data-status', 'running', { timeout: 60_000 });
    await endWorker(page, 'losat-engine', end === 'fails' ? 'out of memory' : undefined);
    await waitStatus(page, 2, 'failed', 30_000);
    expect(await runMessage(page, 2)).toBe(`Run 2 failed: ${sentence} A failed run keeps no results.`);
  });
}

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
