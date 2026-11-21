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
// of LOSAT that runs searches stopped." did not say what to do (L2). Its third review: Run LOSAT and the Examples asked
// the engine nothing again, and were refused once its file downloaded (M1); a runtime already waiting for its keeper
// was never finished when its worker ended (L2). Its fourth review: an Example that waited while LOSAT was asked again
// showed nothing and could not be ended (M1); Run LOSAT did not ask again what was in flight at the press (L1). Its fifth
// review: an Example pressed again while LOSAT was asked again was refused at once (M1); the waiting row offered the
// download from NCBI (L1); a row refused while LOSAT was unavailable said "Reload the page" once LOSAT worked (L2). Its
// sixth review: the same for a set whose inflation waited for the engine file, which said meanwhile that it downloaded.
import { gzipSync } from 'node:zlib';
import { expect, type Page, type Route } from '@playwright/test';
import { test } from './support/test';
import { exampleById } from '../../src/domain/examples';
import { BUILD_HAS_ENGINE } from './support/browser';
import { navigate } from './support/navigate';
import { E2E_ORIGIN } from './support/origin';
import { fasta, openFiles, openParameters, openSearch, paste, pickRun, program, runMessage, submit, waitStatus } from './support/search';

const ENGINE_ONLY = 'the FakeEngine build has no engine modules nor Engine worker: the engine build runs this test';
const RELOAD = 'Reload the page.';
/** The words of an Example's row while LOSAT is asked again, before any request. */
const WAITING = 'Waiting for LOSAT…';
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

/**
 * Keeps the page's ends of the links to the keepers of the engine's shared memory (link.ts), in their order, with the
 * Engine worker that posted each and whether it is open (`keeperLinks`). With `holdKeeper`, the page does not hear a
 * keeper answer that it grew the memory, as while its grow takes long (up to 2 s on a busy machine, gateway.ts).
 */
async function watchKeepers(page: Page, { holdKeeper = false } = {}): Promise<void> {
  await page.addInitScript((hold) => {
    const keepers: Array<{ readonly worker: Worker; open: boolean }> = [];
    Object.assign(window, { keeperLinks: keepers });
    // Worker.prototype's own onmessage, also under a subclass of Worker (keepWorkers).
    let prototype: object = Worker.prototype;
    while (Object.getOwnPropertyDescriptor(prototype, 'onmessage') === undefined) prototype = Object.getPrototypeOf(prototype);
    const workerHandler = Object.getOwnPropertyDescriptor(prototype, 'onmessage')!;
    const portHandler = Object.getOwnPropertyDescriptor(MessagePort.prototype, 'onmessage')!;
    const keep = (worker: Worker, port: MessagePort) => {
      const entry = { worker, open: true };
      keepers.push(entry);
      const close = port.close.bind(port);
      port.close = () => {
        entry.open = false;
        close();
      };
      if (!hold) return;
      Object.defineProperty(port, 'onmessage', {
        get: () => portHandler.get!.call(port),
        set: (listener: ((event: MessageEvent) => void) | null) =>
          portHandler.set!.call(port, listener === null ? null : (event: MessageEvent) => event.data?.type !== 'grown' && listener.call(port, event)),
      });
    };
    Object.defineProperty(prototype, 'onmessage', {
      configurable: true,
      get(this: Worker) {
        return workerHandler.get!.call(this);
      },
      set(this: Worker, listener: ((event: MessageEvent) => void) | null) {
        workerHandler.set!.call(
          this,
          listener === null
            ? null
            : (event: MessageEvent) => {
                if (event.data?.type === 'link' && event.data.role === 'keeper') keep(this, event.data.port);
                listener.call(this, event);
              },
        );
      },
    });
  }, holdKeeper);
}

/** Whether the page's link to the keeper of the `index`th Engine worker it started is open. */
async function keeperOpen(page: Page, index: number): Promise<boolean> {
  return page.evaluate((at) => {
    const { keeperLinks, keptWorkers } = window as unknown as {
      keeperLinks: Array<{ worker: Worker; open: boolean }>;
      keptWorkers: Record<string, Worker[]>;
    };
    const worker = keptWorkers['losat-engine']?.[at];
    return keeperLinks.some((link) => link.worker === worker && link.open);
  }, index);
}

/** How many workers of `name` the page started. */
async function startedWorkers(page: Page, name: 'losat-data' | 'losat-engine'): Promise<number> {
  return page.evaluate((which) => (window as unknown as { keptWorkers: Record<string, Worker[]> }).keptWorkers[which]?.length ?? 0, name);
}

/** How many of the page's links to the keepers of the engine's shared memory are open. */
async function openKeeperLinks(page: Page): Promise<number> {
  return page.evaluate(() => (window as unknown as { keeperLinks: Array<{ open: boolean }> }).keeperLinks.filter((link) => link.open).length);
}

/** How many workers of `name` the page ended itself. */
async function endedByPage(page: Page, name: 'losat-data' | 'losat-engine'): Promise<number> {
  return page.evaluate((which) => (window as unknown as { endedWorkers: Record<string, number> }).endedWorkers[which] ?? 0, name);
}

/**
 * Ends the latest worker of `name` (or the `at`th) from outside, without an error event, as a browser that ends it does;
 * or, with `error`, gives the page the error event of a worker that died of an uncaught error with that message.
 */
async function endWorker(page: Page, name: 'losat-data' | 'losat-engine', error?: string, at = -1): Promise<void> {
  const ended = await page.evaluate(
    ([which, message, index]) => {
      const worker = (window as unknown as { keptWorkers: Record<string, Worker[]> }).keptWorkers[which]?.at(index);
      if (worker === undefined) return false;
      // Not counted as ended by the page: Native's terminate, as the browser ends it.
      if (message === undefined) Object.getPrototypeOf(Object.getPrototypeOf(worker)).terminate.call(worker);
      else worker.dispatchEvent(new ErrorEvent('error', { message, cancelable: true }));
      return true;
    },
    [name, error, at] as const,
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
const ITS = exampleById('refseq-its-fungi');

/** Answers the download of the 16S Example set with a small set; returns the requests that leave the site. */
async function stubExampleHost(page: Page): Promise<string[]> {
  const away: string[] = [];
  page.on('request', (request) => {
    if (/^https?:/.test(request.url()) && !request.url().startsWith(E2E_ORIGIN)) away.push(request.url());
  });
  // A build without the engine has no decoder, and gets the set decoded on the way (examples.spec.ts).
  const body = BUILD_HAS_ENGINE ? gzipSync(SUBJECT) : Buffer.from(SUBJECT);
  await page.route(
    (url) => url.href === SIXTEEN_S.url || url.href === ITS.url,
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

// The third review of D1 and D2, M1 (2026-11-19): with the inputs pasted while the engine file could not be downloaded,
// Run LOSAT and the Examples were refused at once with "Reload the page", and asked nothing, also once the file could
// be downloaded: a reload would have lost the pasted inputs, though one more request read them.
/** The serial engine file: its requests fail, wait (`held`, until the test answers them) or download, by `mode`. */
async function engineFile(page: Page) {
  const file = { mode: 'fail' as 'fail' | 'hold' | 'pass', refused: 0, passed: 0, held: [] as Route[] };
  await page.context().route(
    (url) => /^\/engine\/losat-web-serial-[^/]*\.wasm$/.test(url.pathname),
    async (route) => {
      if (file.mode === 'hold') file.held.push(route);
      else if (file.mode === 'pass') {
        file.passed++;
        await route.continue();
      } else {
        file.refused++;
        await route.abort('failed');
      }
    },
  );
  return file;
}

/**
 * Both inputs pasted while the engine file cannot be downloaded, until everything that the form asked has failed and
 * nothing more is asked: a question still pending when the connection came back would read the inputs again by itself.
 */
async function pastedWhileUnavailable(page: Page, file: Awaited<ReturnType<typeof engineFile>>): Promise<void> {
  await openSearch(page);
  await pasteBoth(page);
  await notRead(page);
  await expect(page.getByTestId('losat-unavailable')).toBeVisible();
  await expect(page.getByTestId('search-panel')).toHaveAttribute('data-describe', 'failed');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'unavailable');
  let last = -1;
  await expect
    .poll(() => {
      const settled = file.refused === last;
      last = file.refused;
      return settled;
    }, { intervals: [1000], timeout: 15_000 })
    .toBe(true);
}

for (const press of ['Run LOSAT', 'an Example'] as const) {
  test(`${press} asks the engine again once its file can be downloaded, and is not refused`, async ({ page }) => {
    test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
    test.setTimeout(120_000);
    const file = await engineFile(page);
    const away = await stubExampleHost(page);
    await pastedWhileUnavailable(page, file);
    // The connection is back; nothing is edited.
    file.mode = 'pass';
    if (press === 'Run LOSAT') {
      await page.getByTestId('run-losat').click();
      await waitStatus(page, 1, 'completed', 60_000);
      await page.getByTestId('tab-search').click();
    } else {
      await page.getByTestId('example-refseq-16s-bacteria').click();
      await expect(page.getByTestId('subject-source-0')).toContainText('bacteria.16SrRNA.fna', { timeout: 30_000 });
      await expect(page.getByTestId('example-error')).toHaveCount(0);
      expect(away).toEqual([SIXTEEN_S.url]);
    }
    expect(file.passed).toBeGreaterThan(0);
    await expect(page.getByTestId('losat-unavailable')).toHaveCount(0);
  });
}

// The fourth review of D1 and D2, M1 (2026-11-19): while LOSAT was asked again, the Example's load had no state: the row
// showed no Cancel, and a mode change did not end it, so 16S (8 MB) downloaded and was saved for a mode the user left.
for (const end of ['a mode change', 'Cancel'] as const) {
  test(`an Example that waits while LOSAT is asked again has Cancel, and ${end} ends it before any request`, async ({ page }) => {
    test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
    test.setTimeout(120_000);
    const file = await engineFile(page);
    const away = await stubExampleHost(page);
    await pastedWhileUnavailable(page, file);
    // The connection is back, and the engine file is on its way.
    file.mode = 'hold';
    const row = 'example-refseq-16s-bacteria';
    const offer = 'Downloads 8.0 MB from ftp.ncbi.nlm.nih.gov';
    await page.getByTestId(row).click();
    // Busy from the press, with Cancel; the row says that it waits for LOSAT, not that it downloads (the fifth review, L1).
    await expect(page.getByTestId(`${row}-cancel`)).toBeVisible();
    await expect(page.getByTestId(`${row}-note`)).toHaveText(WAITING);
    await expect.poll(() => file.held.length).toBeGreaterThan(0);
    if (end === 'Cancel') await page.getByTestId(`${row}-cancel`).click();
    else await page.getByTestId('program-blastp').check();
    file.mode = 'pass';
    for (const route of file.held.splice(0)) await route.continue();
    if (end === 'a mode change') await program(page, 'blastn');
    await expect(page.getByTestId('losat-unavailable')).toHaveCount(0, { timeout: 30_000 });
    await page.waitForTimeout(1000);
    // Nothing was downloaded nor kept: the row offers the download as before.
    await expect(page.getByTestId(`${row}-note`)).toHaveText(offer);
    await expect(page.getByTestId(`${row}-cancel`)).toHaveCount(0);
    expect(away).toEqual([]);
  });
}

// W-3 of line WQA (2026-11-21): a mode pressed as LOSAT was asked again was not chosen (Firefox in CI, run 44767268395).
// A mode change asks LOSAT again, and the alert above the tabs goes a moment later; when that came between the next
// press and its release, the form moved up under the pointer by as much, and the release chose nothing. The page now
// says what it said until the press ends (ui/press-hold.ts).
test('a mode pressed as LOSAT is asked again is chosen, and the alert goes after the press', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  test.setTimeout(120_000);
  const file = await engineFile(page);
  await stubExampleHost(page);
  await pastedWhileUnavailable(page, file);
  // The connection is back, and the engine file is on its way: an Example asks LOSAT again.
  file.mode = 'hold';
  await page.getByTestId('example-refseq-16s-bacteria').click();
  await expect.poll(() => file.held.length).toBeGreaterThan(0);
  await expect(page.locator('main')).toHaveAttribute('data-losat', 'unavailable');
  await expect(page.getByTestId('losat-unavailable')).toBeVisible();
  // A press on BLASTN at the top of the page (as in CI: no scroll can keep it under the pointer).
  await page.evaluate(() => window.scrollTo(0, 0));
  const mode = page.getByTestId('program-blastn');
  const pressed = (await mode.boundingBox())!;
  await page.mouse.move(pressed.x + pressed.width / 2, pressed.y + pressed.height / 2);
  await page.mouse.down();
  // The change to BLASTP that came before the press in CI comes during it here, and so does what it ends: the page
  // no longer reads LOSAT as unavailable before the release.
  await page.getByTestId('program-blastp').evaluate((input: HTMLInputElement) => input.click());
  await expect(page.locator('main')).toHaveAttribute('data-losat', 'available', { timeout: 30_000 });
  expect.soft((await mode.boundingBox())!.y, 'where BLASTN is at the release').toBe(pressed.y);
  await page.mouse.up();
  await expect(mode).toBeChecked();
  await expect(page.getByTestId('losat-unavailable')).toHaveCount(0);
  file.mode = 'pass';
  for (const route of file.held.splice(0)) await route.continue();
});

// The fifth review of D1 and D2, M1 (2026-11-19): what the form asked again read as asked while it was on its way, so a
// second press (after Cancel, or of the other set) asked nothing: it was refused at once with "Reload the page", which
// stayed once the engine file downloaded, and nothing was downloaded.
for (const again of ['Cancel and the same set', 'the other set'] as const) {
  test(`an Example pressed again while LOSAT is asked again (${again}) waits for the answers and loads`, async ({ page }) => {
    test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
    test.setTimeout(120_000);
    const file = await engineFile(page);
    const away = await stubExampleHost(page);
    await pastedWhileUnavailable(page, file);
    // The connection is back, and the engine file is on its way.
    file.mode = 'hold';
    const first = 'example-refseq-16s-bacteria';
    await page.getByTestId(first).click();
    await expect.poll(() => file.held.length).toBeGreaterThan(0);
    // The validation was asked again too: it no longer reads as unavailable.
    await expect(page.getByTestId('argv-validation')).not.toHaveAttribute('data-state', 'unavailable');
    const second = again === 'the other set' ? 'example-refseq-its-fungi' : first;
    if (again === 'Cancel and the same set') await page.getByTestId(`${first}-cancel`).click();
    await page.getByTestId(second).click();
    await page.waitForTimeout(1000);
    await expect(page.getByTestId('example-error')).toHaveCount(0);
    await expect(page.getByTestId(`${second}-note`)).toHaveText(WAITING);
    await expect(page.getByTestId(`${second}-cancel`)).toBeVisible();
    expect(away).toEqual([]);
    file.mode = 'pass';
    for (const route of file.held.splice(0)) await route.continue();
    const set = again === 'the other set' ? ITS : SIXTEEN_S;
    await expect(page.getByTestId('subject-source-0')).toContainText(set.fileName, { timeout: 30_000 });
    await expect(page.getByTestId('example-error')).toHaveCount(0);
    await expect(page.getByTestId('losat-unavailable')).toHaveCount(0);
    expect(away).toEqual([set.url]);
  });
}

// The fifth review of D1 and D2, L2 (2026-11-19): a row refused while the engine file could not be downloaded kept
// "Reload the page" once Run LOSAT had brought LOSAT back.
test('an Example refused while LOSAT was unavailable no longer says to reload once LOSAT works', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  test.setTimeout(120_000);
  const file = await engineFile(page);
  const away = await stubExampleHost(page);
  await pastedWhileUnavailable(page, file);
  await page.getByTestId('example-refseq-16s-bacteria').click();
  await expect(page.getByTestId('example-error')).toHaveText(/^LOSAT could not download its engine file .*\. Reload the page\.$/, { timeout: 15_000 });
  file.mode = 'pass';
  await page.getByTestId('run-losat').click();
  await waitStatus(page, 1, 'completed', 60_000);
  await page.getByTestId('tab-search').click();
  await expect(page.getByTestId('losat-unavailable')).toHaveCount(0);
  await expect(page.getByTestId('example-error')).toHaveCount(0);
  await expect(page.getByTestId('example-refseq-16s-bacteria-note')).toHaveText('Downloads 8.0 MB from ftp.ncbi.nlm.nih.gov');
  expect(away).toEqual([]);
});

// The sixth review of D1 and D2 (2026-11-20): a set pressed while the engine file was still on its way at load downloaded
// at once (nothing read as unavailable yet), and its inflation waited for that file: the row said "Downloading 43 bytes
// of 43 bytes from ftp.ncbi.nlm.nih.gov…" meanwhile, and once the file failed, "Reload the page" also after LOSAT worked.
test('an Example whose inflation waits for the engine file says so, and no longer says to reload once LOSAT works', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  test.setTimeout(120_000);
  const file = await engineFile(page);
  file.mode = 'hold';
  const away = await stubExampleHost(page);
  await navigate(page, '/');
  const row = 'example-refseq-16s-bacteria';
  await expect(page.getByTestId(row)).toBeVisible();
  await expect.poll(() => file.held.length).toBeGreaterThan(0);
  await page.getByTestId(row).click();
  await expect.poll(() => away.length, { timeout: 15_000 }).toBe(1);
  // Soft: the rest says whether the refusal goes too.
  await expect.soft(page.getByTestId(`${row}-note`)).toHaveText(WAITING);
  await expect(page.getByTestId(`${row}-cancel`)).toBeVisible();
  file.mode = 'fail';
  for (const route of file.held.splice(0)) await route.abort('failed');
  await expect(page.getByTestId('example-error')).toHaveText(/^LOSAT could not download its engine file .*\. Reload the page\.$/, { timeout: 15_000 });
  // The connection is back; a pasted query reads the engine file again.
  file.mode = 'pass';
  await paste(page, 'query', QUERY);
  await expect(page.getByTestId('losat-unavailable')).toHaveCount(0, { timeout: 30_000 });
  await expect(page.getByTestId('example-error')).toHaveCount(0);
  expect(away).toEqual([SIXTEEN_S.url]);
});

// The fourth review of D1 and D2, L1 (2026-11-19): a question in flight at the press that failed after it was not asked
// again, so Run LOSAT could say "Reload the page" though the next press succeeded. In the browser a later question
// downloads the engine file again and its answer reads the inputs again, so this run started before the fix too: the
// repeat itself is tested where nothing else asks (availability.test.ts). This guards the press as the user makes it.
test('Run LOSAT starts the search when a read in flight at the press fails after it', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  test.setTimeout(120_000);
  const file = await engineFile(page);
  await pastedWhileUnavailable(page, file);
  // An edit of the query is read while the engine file is on its way.
  file.mode = 'hold';
  await page.getByTestId('query-input').fill(`${QUERY}ACGT\n`);
  await expect.poll(() => file.held.length, { timeout: 15_000 }).toBeGreaterThan(0);
  await page.getByTestId('run-losat').click();
  await page.waitForTimeout(1000);
  // The connection is back, but the requests made before fail.
  file.mode = 'pass';
  for (const route of file.held.splice(0)) await route.abort('failed');
  await waitStatus(page, 1, 'completed', 60_000);
  expect(file.passed).toBeGreaterThan(0);
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
    /^Run 1 failed: LOSAT could not download its engine file \/engine\/losat-web-threads-[0-9a-f]+\.wasm \(network error\)\. Reload the page\. It keeps no outputs\.$/,
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
    `Run 1 failed: LOSAT could not start the part that runs searches. ${RELOAD} It keeps no outputs.`,
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
  expect(await runMessage(page, 1)).toBe(`Run 1 failed: ${DATA_ENDED} It keeps no outputs.`);
  await expect(page.getByTestId('losat-unavailable')).toHaveText(DATA_ENDED);
});

test('a search whose Engine worker ended fails with what stopped, and the next search runs', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  test.setTimeout(180_000);
  await keepWorkers(page);
  await watchKeepers(page);
  await startLongSearch(page);
  expect(await openKeeperLinks(page)).toBe(1);
  await endWorker(page, 'losat-engine');
  await waitStatus(page, 1, 'failed', 30_000);
  expect(await runMessage(page, 1)).toBe(`Run 1 failed: ${ENGINE_ENDED} It keeps no outputs.`);
  // The worker's threads and its keeper ended with it: the page ends its runtime at once. Before, the page waited for the
  // keeper's answer for good, with the keeper's channel open (here in Chromium, Firefox and WebKit).
  await expect.poll(() => openKeeperLinks(page), { timeout: 15_000 }).toBe(0);
  await expect(page.getByTestId('losat-unavailable')).toHaveCount(0);
  await page.getByTestId('tab-search').click();
  await submit(page);
  await waitStatus(page, 2, 'completed', 120_000);
});

// The third review of D1 and D2, L2 (2026-11-19): a runtime that waits for its keeper after a Cancel (or a renewal, a
// thread fault, an error) no longer watched its worker: a worker that ended before the keeper answered left the wait
// for good, with the keeper's channel open. The page watches the worker until the runtime has ended.
test('a cancelled threaded search whose Engine worker ends before its keeper answers ends at once', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, ENGINE_ONLY);
  test.setTimeout(180_000);
  await keepWorkers(page);
  await watchKeepers(page, { holdKeeper: true });
  await startLongSearch(page);
  const engine = (await startedWorkers(page, 'losat-engine')) - 1;
  expect(await keeperOpen(page, engine)).toBe(true);
  await pickRun(page, 1);
  await page.getByTestId('run-status-cancel').click();
  await waitStatus(page, 1, 'cancelled', 15_000);
  // The threads were in use: the page waits for the keeper, which has not answered.
  await page.waitForTimeout(1000);
  expect(await keeperOpen(page, engine)).toBe(true);
  await endWorker(page, 'losat-engine', undefined, engine);
  await expect.poll(() => keeperOpen(page, engine), { timeout: 15_000 }).toBe(false);
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
    expect(await runMessage(page, 2)).toBe(`Run 2 failed: ${sentence} It keeps no outputs.`);
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
