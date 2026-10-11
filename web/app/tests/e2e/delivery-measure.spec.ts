// Measurements of delivery (S16, W7: the release layout, offline use and onboarding; plan §6.2,
// V-PERF style: one warm-up, then three runs, with the median and the range). Not a test: it
// runs only with LOSAT_WEB_MEASURE containing `delivery` (or `all`), in the engine build, and
// writes its records to LOSAT_WEB_EVIDENCE (`delivery-visits-<browser>.json` and
// `delivery-onboarding-<browser>.json`). Run it with --workers=1, one browser at a time.
//
// 1. `visits`: in a new browser context each time (nothing stored, an empty HTTP cache), a first
//    visit, then a second visit in the same page (the service worker answers from what it stored),
//    then an offline visit (the server stopped, as V-OFF goes offline: offline.spec.ts; Chromium
//    and Firefox also get Playwright's offline emulation). The site is the engine build that
//    Playwright serves, behind the front server of V-OFF (support/front-server.ts), which can stop
//    and counts the bytes of every body it sends. For each visit, from the start of the navigation
//    (the page's `performance.now()`): `app`, the search form is on the screen (the boot page has
//    started the release, and the user can paste, open files and press "Run LOSAT"); `checked`, the
//    engine has checked the form's default options (its module is loaded and answers);
//    `preparing` and `ready`, when the header's offline status first says "Preparing offline use"
//    and "Offline ready" (`prepareMs` between the two). An observer in the page (an init script)
//    takes each time when the change reaches the document, before it is painted. `bytes`: what the
//    server sent from the navigation until "Offline ready" (`network`), the part of it sent while
//    the preparation ran (`preparation`: the requests that reached the server after the page said
//    "Preparing offline use", the files the worker downloaded; files that the page had already
//    loaded come from the browser's HTTP cache), and what the worker stores (`stored`:
//    every file of the release manifest, the manifest and the shell, as the site index and the
//    manifest list them).
// 2. `onboarding`: in one page of a first visit, after "Offline ready": the self-test ("Test this
//    browser", in "Getting started"), its time from the click until its summary is on the screen
//    (`ms`) and the time of each case as the panel reports it (the app's own clock around each
//    search); then, for each program that LOSAT Web searches, "Try an example" (from the click until
//    both inputs are read and checked: `load`), "Run LOSAT" (from the click until the queue says
//    `completed`: `search`), and "Open results" (from the click until the first HSP of the run is
//    read in the results: `open`); `firstResults` = `search` + `open`. The first warm-up of the
//    page also loads the threaded engine (the first search of a page): it is kept in the record.
//
// In-page clocks (`timed`) start at the click in the page and stop at the first animation frame in
// which the screen shows the result; their resolution is a frame (about 17 ms).
import { mkdirSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, test, type BrowserContext, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NO_ENGINE_REASON } from './support/engine';
import { startFrontServer, type FrontRequest, type FrontServer } from './support/front-server';
import { E2E_ORIGIN, E2E_PORT } from './support/origin';
import { program } from './support/search';

const MEASURE = new Set((process.env['LOSAT_WEB_MEASURE'] ?? '').split(',').filter((name) => name !== ''));
const EVIDENCE = process.env['LOSAT_WEB_EVIDENCE'] || undefined;
const REPETITIONS = 3;
const RUNNABLE = ['blastn', 'blastp', 'tblastn', 'tblastx'] as const;

test.skip(!MEASURE.has('delivery') && !MEASURE.has('all'), 'measurements run only with LOSAT_WEB_MEASURE');
test.skip(!BUILD_HAS_ENGINE, NO_ENGINE_REASON);
test.setTimeout(1_800_000);
test.use({ actionTimeout: 120_000 });

function save(part: string, browser: string, record: unknown): void {
  if (EVIDENCE === undefined) return;
  mkdirSync(EVIDENCE, { recursive: true });
  writeFileSync(join(EVIDENCE, `delivery-${part}-${browser}.json`), `${JSON.stringify(record, null, 2)}\n`);
}

const message = (error: unknown) => (error instanceof Error ? error.message : String(error));
const median = (values: readonly number[]) => [...values].sort((a, b) => a - b)[Math.floor(values.length / 2)]!;

/** Median, minimum and maximum of every number of the samples (nested records included). */
function summarize(samples: readonly Record<string, unknown>[]): Record<string, unknown> {
  const out: Record<string, unknown> = {};
  for (const key of Object.keys(samples[0] ?? {})) {
    const values = samples.map((sample) => sample[key]);
    if (values.every((value) => typeof value === 'number')) {
      const numbers = values as number[];
      out[key] = { median: median(numbers), min: Math.min(...numbers), max: Math.max(...numbers) };
    } else if (values.every((value) => typeof value === 'object' && value !== null && !Array.isArray(value))) {
      out[key] = summarize(values as Record<string, unknown>[]);
    }
  }
  return out;
}

// --- the clock in the page -----------------------------------------------------------------------

/** What the screen shows when the step is done: the element exists, and its attribute or text matches. */
interface Until {
  readonly testid: string;
  readonly attr?: string;
  /** A regular expression source for the attribute (with `attr`) or the text. */
  readonly matches?: string;
}

/** Clicks `testid` in the page, and returns the milliseconds until a frame in which every `until` holds. */
function timed(page: Page, testid: string, until: readonly Until[], timeoutMs = 600_000): Promise<number> {
  return page.evaluate(
    ({ testid, until, timeoutMs }) =>
      new Promise<number>((resolve, reject) => {
        const holds = () =>
          until.every((u) => {
            const element = document.querySelector(`[data-testid="${u.testid}"]`);
            if (element === null) return false;
            if (u.matches === undefined) return true;
            const value = u.attr === undefined ? element.textContent : element.getAttribute(u.attr);
            return new RegExp(u.matches).test(value ?? '');
          });
        const target = document.querySelector<HTMLElement>(`[data-testid="${testid}"]`);
        if (target === null) {
          reject(new Error(`nothing to click: ${testid}`));
          return;
        }
        const start = performance.now();
        target.click();
        const tick = () => {
          const now = performance.now();
          if (holds()) resolve(Math.round(now - start));
          else if (now - start > timeoutMs) reject(new Error(`${testid}: the screen did not show ${JSON.stringify(until)}`));
          else requestAnimationFrame(tick);
        };
        requestAnimationFrame(tick);
      }),
    { testid, until, timeoutMs },
  );
}

// --- measurement 1: visits -------------------------------------------------------------------------

/** The marks that the init script takes, from the navigation's start. */
interface Marks {
  app?: number;
  checked?: number;
  'offline:preparing'?: number;
  'offline:ready'?: number;
  'offline:failed'?: number;
  'offline:unsupported'?: number;
  boot?: number;
}

/** The test's side of the marks: the server's request count when the page first said "Preparing offline use". */
interface ServerMarks {
  preparing?: number;
}

/**
 * Runs in every document of the context: records when the app's states first reach the document,
 * and tells the test (the binding `__losatServerMark`) when the preparation begins, so that the
 * server's requests from then on are counted as the preparation's (by their order at the server,
 * not by the clocks of two processes).
 */
function installMarks(): void {
  const marks: Record<string, number> = {};
  (globalThis as unknown as { __losatMarks: Record<string, number> }).__losatMarks = marks;
  const final = ['offline:ready', 'offline:failed', 'offline:unsupported'];
  const observer = new MutationObserver(() => {
    const now = performance.now();
    const element = (testid: string) => document.querySelector(`[data-testid="${testid}"]`);
    if (marks['app'] === undefined && element('add-to-queue') !== null) marks['app'] = now;
    if (marks['boot'] === undefined && element('boot-message') !== null) marks['boot'] = now;
    const argv = element('argv-validation')?.getAttribute('data-state');
    if (marks['checked'] === undefined && (argv === 'ok' || argv === 'invalid')) marks['checked'] = now;
    const offline = element('offline-status')?.getAttribute('data-state');
    if (offline !== undefined && offline !== null && marks[`offline:${offline}`] === undefined) {
      marks[`offline:${offline}`] = now;
      if (offline === 'preparing') void (globalThis as unknown as { __losatServerMark?: (name: string) => Promise<void> }).__losatServerMark?.('preparing');
    }
    if (marks['app'] !== undefined && marks['checked'] !== undefined && final.some((name) => marks[name] !== undefined)) observer.disconnect();
  });
  observer.observe(document, { subtree: true, childList: true, attributes: true, attributeFilter: ['data-state'] });
}

interface Visit {
  app: number;
  checked: number;
  preparing: number;
  ready: number;
  prepareMs: number;
  /** The navigation's own times (PerformanceNavigationTiming): the boot page's response and load. */
  navigation: { responseEnd: number; domContentLoaded: number; load: number };
  bytes: { network: number; preparation: number };
  requests: number;
}

/** Opens the site in `page`, and waits until the app has started, checked its options and said "Offline ready". */
async function visit(page: Page, front: FrontServer, at: ServerMarks): Promise<Visit> {
  const from = front.requests.length;
  delete at.preparing;
  await page.goto(`${front.origin}/`);
  await page.waitForFunction(
    () => {
      const m = (globalThis as unknown as { __losatMarks?: Record<string, number> }).__losatMarks;
      if (m === undefined) return false;
      if (m['boot'] !== undefined || m['offline:failed'] !== undefined || m['offline:unsupported'] !== undefined) return true;
      return m['app'] !== undefined && m['checked'] !== undefined && m['offline:ready'] !== undefined;
    },
    undefined,
    { timeout: 300_000, polling: 100 },
  );
  const seen = await page.evaluate(() => {
    const entry = performance.getEntriesByType('navigation')[0] as PerformanceNavigationTiming | undefined;
    return {
      marks: (globalThis as unknown as { __losatMarks: Record<string, number> }).__losatMarks,
      navigation: {
        responseEnd: Math.round(entry?.responseEnd ?? -1),
        domContentLoaded: Math.round(entry?.domContentLoadedEventEnd ?? -1),
        load: Math.round(entry?.loadEventEnd ?? -1),
      },
    };
  });
  const marks = seen.marks as Marks;
  if (marks.boot !== undefined) throw new Error(`the boot page showed a message: ${await page.getByTestId('boot-message').textContent()}`);
  if (marks['offline:ready'] === undefined) throw new Error(`no "Offline ready": ${await page.getByTestId('offline-status').textContent()}`);
  const preparing = marks['offline:preparing'] ?? marks['offline:ready'];
  const ready = marks['offline:ready'];
  // Every request since the navigation: nothing else asks the server after "Offline ready" (the
  // worker's cleanup reads only its caches).
  const requests: FrontRequest[] = front.requests.slice(from);
  const sum = (list: readonly FrontRequest[]) => list.reduce((total, request) => total + request.bytes, 0);
  return {
    app: Math.round(marks.app!),
    checked: Math.round(marks.checked!),
    preparing: Math.round(preparing),
    ready: Math.round(ready),
    prepareMs: Math.round(ready - preparing),
    navigation: seen.navigation,
    bytes: {
      network: sum(requests),
      preparation: sum(front.requests.slice(at.preparing ?? front.requests.length)),
    },
    requests: requests.length,
  };
}

async function goOffline(front: FrontServer, context: BrowserContext, browserName: string): Promise<void> {
  await front.setOnline(false);
  // WebKit fails navigations under Playwright's offline emulation (offline.spec.ts): the stopped server alone.
  if (browserName !== 'webkit') await context.setOffline(true);
}

async function goOnline(front: FrontServer, context: BrowserContext, browserName: string): Promise<void> {
  if (browserName !== 'webkit') await context.setOffline(false);
  await front.setOnline(true);
}

test('visits: first, second (from the service worker) and offline; the preparation until "Offline ready"', async ({ browser, browserName, request }, testInfo) => {
  const record: Record<string, unknown> = { origin: 'front server before the E2E site (engine build)' };
  const index = (await (await request.get('/releases.json')).json()) as { latest: string; shell: { size: number }[] };
  const manifestResponse = await request.get(`/r/${index.latest}/release-manifest.json`);
  const manifest = (await manifestResponse.json()) as { assets: { size: number }[] };
  record['release'] = index.latest;
  record['stored'] = {
    files: index.shell.length + manifest.assets.length + 1,
    bytes:
      index.shell.reduce((total, file) => total + file.size, 0) +
      manifest.assets.reduce((total, file) => total + file.size, 0) +
      (await manifestResponse.body()).length,
  };
  const front = await startFrontServer(E2E_PORT + 100 * (testInfo.parallelIndex + 1), E2E_ORIGIN);
  const samples: Record<string, unknown>[] = [];
  const at: ServerMarks = {};
  try {
    for (let i = 0; i <= REPETITIONS; i++) {
      const context = await browser.newContext();
      try {
        await context.exposeBinding('__losatServerMark', (_source, name: string) => {
          if (name === 'preparing') at.preparing ??= front.requests.length;
        });
        await context.addInitScript(installMarks);
        const page = await context.newPage();
        const first = await visit(page, front, at);
        const second = await visit(page, front, at);
        await goOffline(front, context, browserName);
        let offline: Visit;
        try {
          offline = await visit(page, front, at);
        } finally {
          await goOnline(front, context, browserName);
        }
        const sample = { first, second, offline };
        if (i === 0) record['warmup'] = sample;
        else samples.push(sample);
        console.log(`${browserName} visit ${i === 0 ? 'warm-up' : i}: ${JSON.stringify(sample)}`);
      } finally {
        await context.close();
      }
    }
    record['samples'] = samples;
    record['summary'] = summarize(samples);
  } catch (error) {
    record['error'] = message(error);
    throw error;
  } finally {
    save('visits', browserName, record);
    await front.close();
  }
});

// --- measurement 2: the self-test and the examples ---------------------------------------------------

/** Removes every input source (the search form keeps the sources of the last search). */
async function clearInputs(page: Page): Promise<void> {
  for (const role of ['query', 'subject'] as const) {
    const sources = page.locator(`[data-testid^="${role}-source-"][data-status]`);
    while ((await sources.count()) > 0) await page.getByTestId(`${role}-source-0-remove`).click();
  }
}

interface SelfTestRun {
  ms: number;
  cases: Record<string, number>;
}

async function selfTest(page: Page): Promise<SelfTestRun> {
  const ms = await timed(page, 'self-test-run', [
    { testid: 'self-test-result', attr: 'data-status', matches: '^finished$' },
    { testid: 'self-test-summary' },
  ]);
  const summary = (await page.getByTestId('self-test-summary').textContent()) ?? '';
  if (!/^All \d+ cases passed/.test(summary)) throw new Error(`the self-test did not pass: ${summary}`);
  const cases: Record<string, number> = {};
  for (const item of await page.locator('[data-testid^="self-test-case-"]').all()) {
    const id = (await item.getAttribute('data-testid'))!.replace(/^self-test-case-/, '');
    const time = /, ([\d,]+) ms/.exec((await item.textContent()) ?? '');
    if (time === null) throw new Error(`no time for the self-test case ${id}`);
    cases[id] = Number(time[1]!.replace(/,/g, ''));
  }
  return { ms, cases };
}

interface ExampleRun {
  load: number;
  search: number;
  open: number;
  firstResults: number;
}

async function example(page: Page, id: string, number: number): Promise<ExampleRun> {
  await page.getByTestId('tab-search').click();
  await clearInputs(page);
  await program(page, id);
  await expect(page.getByTestId('try-example')).toBeEnabled();
  const load = await timed(page, 'try-example', [
    { testid: 'example-loaded' },
    { testid: 'query-source-0', attr: 'data-status', matches: '^ready$' },
    { testid: 'subject-source-0', attr: 'data-status', matches: '^ready$' },
    { testid: 'query-source-0-check', attr: 'data-check', matches: '^(?!pending$)' },
    { testid: 'subject-source-0-check', attr: 'data-check', matches: '^(?!pending$)' },
    { testid: 'argv-validation', attr: 'data-state', matches: '^ok$' },
  ]);
  const search = await timed(page, 'add-to-queue', [{ testid: `run-${number}-status`, matches: '^(completed|failed|cancelled)$' }]);
  const status = await page.getByTestId(`run-${number}-status`).textContent();
  if (status !== 'completed') throw new Error(`run ${number} (${id}) ended ${status}`);
  const open = await timed(page, `run-${number}-open`, [
    { testid: 'results-hits', attr: 'data-run', matches: `^${number}$` },
    { testid: 'hsp-detail', attr: 'data-state', matches: '^ready$' },
  ]);
  return { load, search, open, firstResults: search + open };
}

test('onboarding: the self-test, and "Try an example" then "Run LOSAT" to the first results', async ({ page, browserName }) => {
  const record: Record<string, unknown> = {};
  try {
    // A first visit: "Getting started", where the self-test is, is open.
    await page.goto('/');
    await expect(page.getByTestId('offline-status')).toHaveAttribute('data-state', 'ready', { timeout: 300_000 });
    await expect(page.getByTestId('getting-started')).toBeVisible();

    await page.getByTestId('self-test').locator('summary').click();
    const tests: SelfTestRun[] = [];
    for (let i = 0; i <= REPETITIONS; i++) {
      const run = await selfTest(page);
      if (i === 0) record['selfTestWarmup'] = run;
      else tests.push(run);
      console.log(`${browserName} self-test ${i === 0 ? 'warm-up' : i}: ${JSON.stringify(run)}`);
    }
    record['selfTest'] = { samples: tests, summary: summarize(tests as unknown as Record<string, unknown>[]) };
    save('onboarding', browserName, record);

    let number = 0;
    const examples: Record<string, unknown> = {};
    record['examples'] = examples;
    for (const id of RUNNABLE) {
      const samples: ExampleRun[] = [];
      const entry: Record<string, unknown> = {};
      examples[id] = entry;
      for (let i = 0; i <= REPETITIONS; i++) {
        const run = await example(page, id, ++number);
        if (i === 0) entry['warmup'] = run;
        else samples.push(run);
        console.log(`${browserName} ${id} ${i === 0 ? 'warm-up' : i}: ${JSON.stringify(run)}`);
      }
      entry['samples'] = samples;
      entry['summary'] = summarize(samples as unknown as Record<string, unknown>[]);
      save('onboarding', browserName, record);
    }
  } catch (error) {
    record['error'] = message(error);
    throw error;
  } finally {
    save('onboarding', browserName, record);
  }
});
