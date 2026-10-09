// Measurements of the results screen with many queries and many HSPs (S13; plan §2.3 "図と一覧の
// 集約表示": the maintainer sets the drawing-count threshold from these numbers; W3's handoff:
// "the results screen's query picker and lists must handle a run of 100,000 queries,
// virtualised"). Not a test: it runs only with LOSAT_WEB_MEASURE containing `results` (or
// `all`), with the engine build, and writes its records to LOSAT_WEB_EVIDENCE
// (`results-<browser>.json`; docs/evidence/losat_web_w4/ reports the results). The query counts
// come from LOSAT_WEB_MEASURE_QUERIES (default 10000,100000), the copies of the repeated unit of
// the dot-plot case from LOSAT_WEB_MEASURE_COPIES (default 1500,3000,4000: in the try runs of
// S13, 3000 was the largest that completed in Chromium and Firefox, WebKit's in-memory store
// refused 3000, and the engine ran out of memory at 4000; a run that fails is recorded with its
// error, and the measurement goes on).
//
// Two measurements, each through the application with a real search:
//
// 1. `queries`: for each count N, a BLASTN search of N queries of 100 letters (most with hits:
//    70% one subject, 20% three, every 1000th "wide" query with 200 subjects, 10% random without
//    hits; tests/e2e/support/synthetic.ts) against 203 subject records. It records the input
//    (indexing and checking the query file, after the subject's), the search (from "Add to queue"
//    until it completes; the run's elapsed time as the queue shows it, and the runtime path), the
//    number of HSPs (the lines of the run's outfmt 6 in the Outputs view, once), and, in one
//    warm-up and three measured repetitions, in the same page: opening the run from the queue
//    until the hits view is ready (the HSP records read and transferred, the index built), and
//    again until the first query's first HSP is read as well, the rows the query picker draws,
//    scrolling it to its end, finding the last query but one with `query-filter` (which selects
//    it), clearing the filter, selecting the query before it by a click, "With hits only",
//    selecting the second query with 200 subjects by a click and sorting its subject list by
//    each column. In Chromium also the JS heap (after a garbage
//    collection) before opening, after the first opening and after the repetitions, and the
//    memory of the page and its workers where the browser exposes it.
// 2. `pair`: a query that repeats a unit of 20 letters (with 1% differences between the copies)
//    searched against itself: every pair of copies is a diagonal, so one query-subject pair has
//    about two HSPs per copy. It records the search, the segments on the dot plot, and, in one
//    warm-up and three repetitions, the time to open the run, to show the dot plot, to zoom in,
//    out and to the selected HSP, to select an HSP by `n`, by a click on its segment and in the
//    list, to sort and to scroll the HSP list.
//
// Times are milliseconds, taken in the page with `performance.now()` from the action to the first
// animation frame in which the screen shows the result (`ms`) and to the frame after it
// (`paintMs`); the resolution is a frame (about 17 ms). An action that the test performs through
// Playwright's mouse or keyboard (selecting an HSP by `n` or by a click on the dot plot) starts its
// clock at the event in the page, and the clock stops at the first frame that shows the result
// after the test starts looking, so these times include Playwright's round trip (upper bounds).
import { mkdirSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, test, type CDPSession, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NO_ENGINE_REASON } from './support/engine';
import { program, settled, submit } from './support/search';
import { manyQueries, queryId, repeats } from './support/synthetic';

const MEASURE = new Set((process.env['LOSAT_WEB_MEASURE'] ?? '').split(',').filter((name) => name !== ''));
const EVIDENCE = process.env['LOSAT_WEB_EVIDENCE'] || undefined;
const COUNTS = (process.env['LOSAT_WEB_MEASURE_QUERIES'] ?? '10000,100000').split(',').map(Number);
const COPIES = (process.env['LOSAT_WEB_MEASURE_COPIES'] ?? '1500,3000,4000').split(',').map(Number);
const REPETITIONS = 3;
const VERBOSE = process.env['LOSAT_WEB_MEASURE_VERBOSE'] !== undefined;
const note = (text: string) => VERBOSE && console.log(`${new Date().toISOString().slice(11, 19)} ${text}`);
/** A search that takes longer than this is a failure of the measurement (the record says so). */
const SEARCH_MINUTES = Number(process.env['LOSAT_WEB_MEASURE_SEARCH_MINUTES'] ?? 30);

test.skip(!MEASURE.has('results') && !MEASURE.has('all'), 'measurements run only with LOSAT_WEB_MEASURE');
test.skip(!BUILD_HAS_ENGINE, NO_ENGINE_REASON);
test.setTimeout(7_200_000);
// A click that cannot happen fails after two minutes, not at the end of the test.
test.use({ actionTimeout: 120_000 });

// --- the clock in the page ---------------------------------------------------------------------

/** What the page does at time zero. */
type Act =
  | { readonly kind: 'click'; readonly testid: string }
  | { readonly kind: 'fill'; readonly testid: string; readonly value: string }
  | {
      readonly kind: 'scroll';
      readonly testid: string;
      readonly to: 'start' | 'end';
    }
  /** The clock starts at the event of an `arm`ed element (Playwright's own mouse or keyboard acts). */
  | { readonly kind: 'armed' };

/** What the screen shows when the step is done (all must hold). */
interface Cond {
  /** The element's test ID, or (with `css`) ignored. */
  readonly testid: string;
  /** A CSS selector for the first matching element, in place of `testid`. */
  readonly css?: string;
  /** The parent of the element, not the element. */
  readonly parent?: boolean;
  readonly absent?: boolean;
  readonly attr?: string;
  readonly equals?: string;
  /** The attribute is something else than this (for an act that happened before the step began). */
  readonly differs?: string;
  /** A regular expression source. */
  readonly matches?: string;
  /** The attribute differs from what it was at time zero. */
  readonly changed?: boolean;
  /** A regular expression source for the element's text. */
  readonly text?: string;
}

interface Timing {
  readonly ms: number;
  readonly paintMs: number;
}

async function arm(page: Page, testid: string, event: string): Promise<void> {
  await page.evaluate(
    ({ testid: id, event: name }) => {
      const element = document.querySelector(`[data-testid="${id}"]`)!;
      element.addEventListener(name, () => ((window as unknown as { __t0: number }).__t0 = performance.now()), {
        once: true,
        capture: true,
      });
    },
    { testid, event },
  );
}

/** Performs `act` in the page and waits (in the page) until the screen shows `until`. */
async function step(page: Page, act: Act, until: readonly Cond[], timeoutMs = 600_000): Promise<Timing> {
  return page.evaluate(
    async ({ act: action, until: conditions, timeoutMs: limit }) => {
      const find = (id: string) => document.querySelector(`[data-testid="${id}"]`) as HTMLElement | null;
      const target = (c: { testid: string; css?: string; parent?: boolean }) => {
        const element = c.css === undefined ? find(c.testid) : (document.querySelector(c.css) as HTMLElement | null);
        return c.parent && element ? element.parentElement : element;
      };
      const initial = conditions.map((c) => (c.attr === undefined ? null : (target(c)?.getAttribute(c.attr) ?? null)));
      let t0 = performance.now();
      if (action.kind === 'armed') t0 = (window as unknown as { __t0: number }).__t0;
      else if (action.kind === 'click') find(action.testid)!.click();
      else if (action.kind === 'scroll') {
        const element = find(action.testid)!;
        element.scrollTop = action.to === 'end' ? element.scrollHeight : 0;
      } else {
        const input = find(action.testid) as HTMLInputElement;
        input.value = action.value;
        input.dispatchEvent(new Event('input', { bubbles: true }));
      }
      const holds = () =>
        conditions.every((c, i) => {
          const element = target(c);
          if (c.absent) return element === null;
          if (element === null) return false;
          if (c.attr !== undefined) {
            const value = element.getAttribute(c.attr);
            if (c.changed) return value !== initial[i];
            if (c.equals !== undefined) return value === c.equals;
            if (c.differs !== undefined) return value !== c.differs;
            if (c.matches !== undefined) return value !== null && new RegExp(c.matches).test(value);
            return value !== null;
          }
          if (c.text !== undefined) return new RegExp(c.text).test(element.textContent ?? '');
          return true;
        });
      return new Promise<{ ms: number; paintMs: number }>((resolve, reject) => {
        const poll = () => {
          if (holds()) {
            const ms = performance.now() - t0;
            requestAnimationFrame(() => resolve({ ms, paintMs: performance.now() - t0 }));
          } else if (performance.now() - t0 > limit) {
            reject(new Error(`timed out after ${limit} ms waiting for ${JSON.stringify(conditions)}`));
          } else requestAnimationFrame(poll);
        };
        // The first look is after the next frame, when the page has rendered what the act changed.
        requestAnimationFrame(poll);
      });
    },
    { act, until, timeoutMs },
  );
}

const rounded = (timing: Timing): Timing => ({
  ms: Math.round(timing.ms),
  paintMs: Math.round(timing.paintMs),
});
const median = (values: readonly number[]) => [...values].sort((a, b) => a - b)[Math.floor(values.length / 2)]!;

/** Median, minimum and maximum of every number of the samples (nested records included). */
function summarize(samples: readonly Record<string, unknown>[]): Record<string, unknown> {
  const out: Record<string, unknown> = {};
  for (const key of Object.keys(samples[0] ?? {})) {
    const values = samples.map((sample) => sample[key]);
    if (values.every((value) => typeof value === 'number')) {
      const numbers = values as number[];
      out[key] = {
        median: median(numbers),
        min: Math.min(...numbers),
        max: Math.max(...numbers),
      };
    } else if (values.every((value) => typeof value === 'object' && value !== null)) {
      out[key] = summarize(values as Record<string, unknown>[]);
    }
  }
  return out;
}

// --- the records -------------------------------------------------------------------------------

interface Records {
  queries: unknown[];
  pair: unknown[];
}
const records: Records = { queries: [], pair: [] };
let browser = 'unknown';

function save(): void {
  if (EVIDENCE === undefined) return;
  mkdirSync(EVIDENCE, { recursive: true });
  writeFileSync(join(EVIDENCE, `results-${browser}.json`), `${JSON.stringify(records, null, 2)}\n`);
}

const message = (error: unknown) => (error instanceof Error ? error.message : String(error));

// --- Chromium's memory ---------------------------------------------------------------------------

interface Heap {
  readonly usedMB: number;
  readonly totalMB: number;
  /** `performance.measureUserAgentSpecificMemory()`: the page and its workers, where the browser has it. */
  readonly agentMB?: number;
}

async function heap(page: Page, session: CDPSession | undefined): Promise<Heap | undefined> {
  if (session === undefined) return undefined;
  await session.send('HeapProfiler.collectGarbage');
  const usage = (await session.send('Runtime.getHeapUsage')) as {
    usedSize: number;
    totalSize: number;
  };
  const agent = await page
    .evaluate(async () => {
      const measure = (
        performance as unknown as {
          measureUserAgentSpecificMemory?: () => Promise<{ bytes: number }>;
        }
      ).measureUserAgentSpecificMemory;
      if (measure === undefined || !crossOriginIsolated) return undefined;
      return (await measure.call(performance)).bytes;
    })
    .catch(() => undefined);
  return {
    usedMB: Math.round(usage.usedSize / 1e5) / 10,
    totalMB: Math.round(usage.totalSize / 1e5) / 10,
    ...(agent === undefined ? {} : { agentMB: Math.round(agent / 1e5) / 10 }),
  };
}

// --- the search -----------------------------------------------------------------------------------

async function clearInputs(page: Page): Promise<void> {
  for (const role of ['query', 'subject'] as const) {
    const sources = page.locator(`[data-testid^="${role}-source-"][data-status]`);
    while ((await sources.count()) > 0) await page.getByTestId(`${role}-source-0-remove`).click();
  }
}

interface Search {
  /** Choosing the query file until its record table is shown, and until the engine's check answers. */
  readonly indexMs: number;
  readonly checkMs: number;
  /** "Add to queue" until the run is in the queue (the run input, its SHA-256, the engine's validation). */
  readonly enqueueMs: number;
  /** The run in the queue until it completes, seen from the page. */
  readonly searchMs: number;
  /** What the queue shows for the run. */
  readonly elapsed: string;
  readonly path: string;
}

/** Queues a search of the files with the form's defaults and waits until run `number` completes. */
async function search(
  page: Page,
  number: number,
  query: { name: string; buffer: Buffer },
  subject: { name: string; buffer: Buffer },
): Promise<Search> {
  await page.getByTestId('tab-search').click();
  await clearInputs(page);
  await program(page, 'blastn');
  await page.getByTestId('subject-files').setInputFiles({
    name: subject.name,
    mimeType: 'text/plain',
    buffer: subject.buffer,
  });
  await settled(page, 'subject');
  // The clock of the input starts after the subject is read and checked: it times the query file.
  const t0 = Date.now();
  await page.getByTestId('query-files').setInputFiles({
    name: query.name,
    mimeType: 'text/plain',
    buffer: query.buffer,
  });
  const source = page.getByTestId('query-source-0');
  await expect(source).toHaveAttribute('data-status', 'ready', {
    timeout: 900_000,
  });
  const t1 = Date.now();
  note(`run ${number}: query file read`);
  await expect(page.getByTestId('query-source-0-check')).toHaveAttribute('data-check', 'ok', { timeout: 900_000 });
  const t2 = Date.now();
  note(`run ${number}: query checked`);
  await expect(page.getByTestId('argv-validation')).not.toHaveAttribute('data-state', 'checking');
  const tSubmit = Date.now();
  await submit(page);
  await expect(page.getByTestId(`run-${number}`)).toBeVisible({
    timeout: 600_000,
  });
  const t3 = Date.now();
  note(`run ${number}: in the queue`);
  const status = page.getByTestId(`run-${number}-status`);
  // A long search is reported once a minute, with the phase the queue shows.
  for (let minutes = 1; ; minutes++) {
    try {
      await expect(status).toHaveText(/^(completed|cancelled|failed)$/, { timeout: 60_000 });
      break;
    } catch (error) {
      if (minutes >= SEARCH_MINUTES) throw error;
      const phase =
        (await page
          .getByTestId(`run-${number}-phase`)
          .textContent()
          .catch(() => null)) ?? '?';
      const elapsed =
        (await page
          .getByTestId(`run-${number}-elapsed`)
          .textContent()
          .catch(() => null)) ?? '';
      console.log(`run ${number}: ${minutes} min, ${phase.trim()}, ${elapsed.trim()}`);
    }
  }
  const t4 = Date.now();
  if ((await status.textContent()) !== 'completed') {
    throw new Error(
      `run ${number} ended ${await status.textContent()}: ${
        (await page
          .getByTestId(`run-${number}-error`)
          .textContent()
          .catch(() => '')) ?? ''
      }`,
    );
  }
  const elapsed =
    (await page
      .getByTestId(`run-${number}-elapsed`)
      .textContent()
      .catch(() => '')) ?? '';
  await page.getByTestId(`run-${number}`).locator('summary', { hasText: 'Details' }).click();
  const details = page.getByTestId(`run-${number}-details`);
  const path =
    (await details
      .locator('[data-detail="path"]')
      .textContent()
      .catch(() => '')) ?? '';
  return {
    indexMs: t1 - t0,
    checkMs: t2 - t1,
    enqueueMs: t3 - tSubmit,
    searchMs: t4 - t3,
    elapsed: elapsed.trim(),
    path: path.trim().replace(/\s+/g, ' '),
  };
}

// --- measurement 1: many queries ------------------------------------------------------------------

const countText = (n: number) => n.toLocaleString('en-US');
const rowsDrawn = (page: Page, prefix: string) => page.locator(`[data-testid^="${prefix}"]`).count();

interface Repetition {
  open: Timing;
  firstHsp: Timing;
  queryRowsDrawn: number;
  scrollToEnd: Timing;
  rowsDrawnAtEnd: number;
  scrollToStart: Timing;
  findFar: Timing;
  selectFar: Timing;
  clearFilter: Timing;
  hitsOnly: Timing;
  hitsOnlyQueries: number;
  hitsOnlyOff: Timing;
  selectWide: Timing;
  wideSubjects: number;
  subjectRowsDrawn: number;
  sortSubjects: Record<string, Timing>;
}

/**
 * Opens run `number` from the queue while run `other` is shown, twice: `open` until its hits view
 * is ready, and `firstHsp` until the first HSP of its first query with hits (selected when the run
 * opens) is also read from outfmt 0.
 */
async function openRun(page: Page, number: number, other: number): Promise<{ open: Timing; firstHsp: Timing }> {
  const open = await step(page, { kind: 'click', testid: `run-${number}-open` }, [
    { testid: 'results-hits', attr: 'data-run', equals: String(number) },
    { testid: 'query-list' },
  ]);
  await page.getByTestId(`run-${other}-open`).click();
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', String(other), { timeout: 120_000 });
  const firstHsp = await step(page, { kind: 'click', testid: `run-${number}-open` }, [
    { testid: 'results-hits', attr: 'data-run', equals: String(number) },
    { testid: 'hsp-detail', attr: 'data-state', equals: 'ready' },
  ]);
  return { open: rounded(open), firstHsp: rounded(firstHsp) };
}

async function repetition(page: Page, count: number): Promise<Repetition> {
  const far = count - 2;
  const total = countText(count);
  // The small run (run 2) first, so that opening run 1 is a change of run.
  await page.getByTestId('run-2-open').click();
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '2', { timeout: 120_000 });
  const { open, firstHsp } = await openRun(page, 1, 2);
  await expect(page.getByTestId('query-count')).toHaveText(`${total} of ${total} queries`);
  const queryRowsDrawn = await rowsDrawn(page, 'query-row-');

  const scrollToEnd = await step(page, { kind: 'scroll', testid: 'query-list', to: 'end' }, [{ testid: `query-row-${count - 1}` }]);
  const rowsDrawnAtEnd = await rowsDrawn(page, 'query-row-');
  const scrollToStart = await step(page, { kind: 'scroll', testid: 'query-list', to: 'start' }, [{ testid: 'query-row-0' }]);

  // A query far down, by its ID: the filter lists it alone and moves the selection to it (the
  // query selected before is hidden; decision 10), whose first HSP is read.
  const findFar = await step(page, { kind: 'fill', testid: 'query-filter', value: queryId(far) }, [
    { testid: 'query-count', text: `^1 of ${total} queries$` },
    { testid: `query-row-${far}` },
    { testid: 'hsp-detail', attr: 'data-hsp', matches: `^${far}:` },
    { testid: 'hsp-detail', attr: 'data-state', equals: 'ready' },
  ]);
  const clearFilter = await step(page, { kind: 'fill', testid: 'query-filter', value: '' }, [
    { testid: 'query-count', text: `^${total} of ${total} queries$` },
  ]);
  // A click on the query before it (a query with hits; the list keeps the selected query in view).
  const near = far - 1;
  await expect(page.getByTestId(`query-row-${near}`)).toHaveCount(1);
  const selectFar = await step(page, { kind: 'click', testid: `query-row-${near}` }, [
    { testid: 'hsp-detail', attr: 'data-hsp', matches: `^${near}:` },
    { testid: 'hsp-detail', attr: 'data-state', equals: 'ready' },
  ]);

  // "With hits only" (the selected query has hits).
  const hitsOnly = await step(page, { kind: 'click', testid: 'filter-hits-only' }, [{ testid: 'query-count', text: `^(?!${total} of)` }]);
  const [, withHits] = /^([\d,]+) of/.exec((await page.getByTestId('query-count').textContent()) ?? '')!;
  const hitsOnlyOff = await step(page, { kind: 'click', testid: 'filter-hits-only' }, [
    { testid: 'query-count', text: `^${total} of ${total} queries$` },
  ]);

  // The queries with 200 subjects: the filter selects the first of them (decision 10); a click on
  // the second, whose subject list is then sorted by each column.
  await page.getByTestId('query-filter').fill('wide');
  await expect(page.getByTestId('query-count')).toHaveText(new RegExp(`^[\\d,]+ of ${total} queries$`));
  const wideRows = page.locator('[data-testid^="query-row-"]');
  const rowIndex = async (n: number) => Number((await wideRows.nth(n).getAttribute('data-testid'))!.replace('query-row-', ''));
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-hsp', new RegExp(`^${await rowIndex(0)}:`));
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-state', 'ready');
  const wideIdx = await rowIndex(1);
  const selectWide = await step(page, { kind: 'click', testid: `query-row-${wideIdx}` }, [
    { testid: 'hsp-detail', attr: 'data-hsp', matches: `^${wideIdx}:` },
    { testid: 'hsp-detail', attr: 'data-state', equals: 'ready' },
    { testid: 'subject-list', attr: 'data-count', matches: '^[1-9]\\d{2,}$' },
  ]);
  const wideSubjects = Number(await page.getByTestId('subject-list').getAttribute('data-count'));
  const subjectRowsDrawn = await rowsDrawn(page, 'subject-row-');
  const sortSubjects: Record<string, Timing> = {};
  for (const key of ['length', 'bitScore', 'eValue', 'hsps', 'order']) {
    sortSubjects[key] = rounded(
      await step(page, { kind: 'click', testid: `subject-sort-${key}` }, [
        {
          testid: `subject-sort-${key}`,
          parent: true,
          attr: 'aria-sort',
          changed: true,
        },
      ]),
    );
  }
  // The view filters of the next repetition start empty, whatever opening a run does with them.
  await page.getByTestId('query-filter').fill('');
  await expect(page.getByTestId('query-count')).toHaveText(`${total} of ${total} queries`);
  await page.getByTestId('filter-clear').click();

  return {
    open,
    firstHsp,
    queryRowsDrawn,
    scrollToEnd: rounded(scrollToEnd),
    rowsDrawnAtEnd,
    scrollToStart: rounded(scrollToStart),
    findFar: rounded(findFar),
    selectFar: rounded(selectFar),
    clearFilter: rounded(clearFilter),
    hitsOnly: rounded(hitsOnly),
    hitsOnlyQueries: Number(withHits!.replace(/,/g, '')),
    hitsOnlyOff: rounded(hitsOnlyOff),
    selectWide: rounded(selectWide),
    wideSubjects,
    subjectRowsDrawn,
    sortSubjects,
  };
}

// One test, with its own page and storage, for each count.
for (const count of COUNTS) {
  test(`${countText(count)} queries: opening the run, the query picker, finding and sorting`, async ({ page, browserName }) => {
    browser = browserName;
    const session = browserName === 'chromium' ? await page.context().newCDPSession(page) : undefined;
    const record: Record<string, unknown> = { count };
    records.queries.push(record);
    try {
      await page.goto('/');
      const inputs = manyQueries(count);
      Object.assign(record, {
        queryBytes: inputs.queries.length,
        subjectBytes: inputs.subjects.length,
        kinds: inputs.kinds,
      });
      // Run 1: the many queries; run 2: a small one, to change runs with.
      record['search'] = await search(
        page,
        1,
        { name: `queries-${count}.fna`, buffer: inputs.queries },
        { name: 'subjects.fna', buffer: inputs.subjects },
      );
      console.log(`${browserName} ${count} queries: search ${JSON.stringify(record['search'])}`);
      const small = manyQueries(20);
      record['smallSearch'] = await search(
        page,
        2,
        { name: 'queries-20.fna', buffer: small.queries },
        { name: 'subjects.fna', buffer: inputs.subjects },
      );
      save();

      const before = await heap(page, session);
      await page.getByTestId('tab-results').click();
      const warmup = await repetition(page, count);
      record['warmup'] = warmup;
      const afterFirst = await heap(page, session);
      const samples: Repetition[] = [];
      for (let i = 0; i < REPETITIONS; i++) samples.push(await repetition(page, count));
      record['samples'] = samples;
      record['summary'] = summarize(samples as unknown as Record<string, unknown>[]);
      const afterAll = await heap(page, session);
      if (before !== undefined)
        record['heap'] = {
          beforeOpening: before,
          afterFirstOpening: afterFirst,
          afterRepetitions: afterAll,
        };
      save();

      // The number of HSPs: the lines of outfmt 6 as the Outputs view shows them (once, after the repetitions).
      await page.getByTestId('run-1-open').click();
      await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1', { timeout: 120_000 });
      await page.getByTestId('results-view-outputs').click();
      const outputs = await step(
        page,
        { kind: 'click', testid: 'format-6' },
        [{ testid: 'result-output', attr: 'data-shown', equals: '1:6' }],
        900_000,
      ).catch(async (error) => {
        record['outfmt6Error'] = message(error);
        return undefined;
      });
      if (outputs !== undefined) {
        const lines = await page.getByTestId('result-output').evaluate((element) => {
          const text = element.textContent ?? '';
          return {
            bytes: text.length,
            hsps: text.split('\n').filter((line) => line !== '' && !line.startsWith('#')).length,
          };
        });
        record['outfmt6'] = {
          ...lines,
          showMs: Math.round(outputs.ms),
          showPaintMs: Math.round(outputs.paintMs),
        };
      }
      console.log(
        `${browserName} ${count} queries: ${JSON.stringify({ summary: record['summary'], heap: record['heap'], outfmt6: record['outfmt6'] })}`,
      );
    } catch (error) {
      record['error'] = message(error);
      console.log(`${browserName} ${count} queries: FAILED ${message(error)}`);
    }
    save();
  });
}

// --- measurement 2: many HSPs in one query-subject pair ---------------------------------------------

interface PairRepetition {
  open: Timing;
  hsps: number;
  hspRowsDrawn: number;
  showDotPlot: Timing;
  segments: number;
  coloured: number;
  zoomIn: Timing;
  zoomOut: Timing;
  zoomToHsp: Timing;
  reset: Timing;
  selectByKey: Timing;
  selectByClick: Timing;
  clickHitTarget: number;
  selectInList: Timing;
  sortHsps: Timing;
  scrollHspList: Timing;
}

async function pairRepetition(page: Page, number: number, other: number): Promise<PairRepetition> {
  await page.getByTestId(`run-${other}-open`).click();
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', String(other), { timeout: 120_000 });
  const open = await step(page, { kind: 'click', testid: `run-${number}-open` }, [
    { testid: 'results-hits', attr: 'data-run', equals: String(number) },
    { testid: 'hsp-detail', attr: 'data-state', equals: 'ready' },
  ]);
  await page.getByTestId('pane-alignment').click();
  const hsps = Number(await page.getByTestId('hsp-list').getAttribute('data-count'));
  const hspRowsDrawn = await rowsDrawn(page, 'hsp-row-');

  const showDotPlot = await step(page, { kind: 'click', testid: 'pane-dotplot' }, [
    { testid: 'dotplot-canvas', attr: 'data-segments', equals: String(hsps) },
  ]);
  const canvas = page.getByTestId('dotplot-canvas');
  const segments = Number(await canvas.getAttribute('data-segments'));
  // Pixels in the colours of the HSPs (forward blue, reverse orange): the plot is drawn.
  const coloured = await canvas.evaluate((element) => {
    const c = element as HTMLCanvasElement;
    const { data } = c.getContext('2d')!.getImageData(0, 0, c.width, c.height);
    let n = 0;
    for (let i = 0; i < data.length; i += 4) {
      const [r, g, b] = [data[i]!, data[i + 1]!, data[i + 2]!];
      if (Math.abs(r - 31) + Math.abs(g - 95) + Math.abs(b - 191) < 60 || Math.abs(r - 194) + Math.abs(g - 65) + Math.abs(b - 12) < 60) n++;
    }
    return n;
  });
  const whole = (await canvas.getAttribute('data-view'))!;
  const zoomIn = await step(page, { kind: 'click', testid: 'dotplot-zoom-in' }, [
    { testid: 'dotplot-canvas', attr: 'data-view', changed: true },
  ]);
  const zoomOut = await step(page, { kind: 'click', testid: 'dotplot-zoom-out' }, [
    { testid: 'dotplot-canvas', attr: 'data-view', equals: whole },
  ]);

  // Selection: the key n, a click on a segment (the first 200 are offered to tests), a click in the list.
  const selected = (await canvas.getAttribute('data-selected'))!;
  await canvas.focus();
  await arm(page, 'dotplot-canvas', 'keydown');
  await page.keyboard.press('n');
  const selectByKey = await step(
    page,
    { kind: 'armed' },
    [{ testid: 'dotplot-canvas', attr: 'data-selected', differs: selected }],
    30_000,
  ).catch(() => {
    throw new Error(`n did not change the selected HSP (${selected})`);
  });
  const targets = JSON.parse((await canvas.getAttribute('data-targets'))!) as {
    hsp: string;
    x: number;
    y: number;
  }[];
  const now = (await canvas.getAttribute('data-selected'))!;
  const target = targets.find((t) => t.hsp !== now && t.hsp !== selected)!;
  await arm(page, 'dotplot-canvas', 'pointerdown');
  await canvas.click({ position: { x: target.x, y: target.y } });
  const selectByClick = await step(page, { kind: 'armed' }, [{ testid: 'dotplot-canvas', attr: 'data-selected', differs: now }], 30_000);
  // Where segments lie closer together than the click's reach, the click picks the nearest of them, not always the one aimed at.
  const clickHitTarget = (await canvas.getAttribute('data-selected')) === target.hsp ? 1 : 0;

  // "Zoom to HSP" frames the selected HSP. An HSP that spans nearly the whole of both sequences hardly changes the view, so the clock stops at the next frame.
  const zoomToHsp = await step(page, { kind: 'click', testid: 'dotplot-zoom-hsp' }, [{ testid: 'dotplot-zoom-hsp' }]);
  const reset = await step(page, { kind: 'click', testid: 'dotplot-reset' }, [
    { testid: 'dotplot-canvas', attr: 'data-view', equals: whole },
  ]);

  await page.getByTestId('pane-alignment').click();
  const row = page.getByTestId('hsp-list').locator('[data-testid^="hsp-row-"][aria-pressed="false"]').nth(2);
  const rowId = (await row.getAttribute('data-testid'))!.replace('hsp-row-', '').replace('-', ':');
  const selectInList = await step(page, { kind: 'click', testid: `hsp-row-${rowId.replace(':', '-')}` }, [
    { testid: 'hsp-detail', attr: 'data-hsp', equals: rowId },
    { testid: 'hsp-detail', attr: 'data-state', equals: 'ready' },
  ]);
  const sortHsps = await step(page, { kind: 'click', testid: 'hsp-sort-bitScore' }, [
    {
      testid: 'hsp-sort-bitScore',
      parent: true,
      attr: 'aria-sort',
      changed: true,
    },
  ]);
  const firstRow = '[data-testid="hsp-list"] [data-testid^="hsp-row-"]';
  const lastRow = await step(
    page,
    { kind: 'scroll', testid: 'hsp-list', to: 'end' },
    [{ testid: '', css: firstRow, attr: 'data-testid', changed: true }],
    60_000,
  );
  await step(
    page,
    { kind: 'scroll', testid: 'hsp-list', to: 'start' },
    [{ testid: '', css: firstRow, attr: 'data-testid', changed: true }],
    60_000,
  );
  return {
    open: rounded(open),
    hsps,
    hspRowsDrawn,
    showDotPlot: rounded(showDotPlot),
    segments,
    coloured,
    zoomIn: rounded(zoomIn),
    zoomOut: rounded(zoomOut),
    zoomToHsp: rounded(zoomToHsp),
    reset: rounded(reset),
    selectByKey: rounded(selectByKey),
    selectByClick: rounded(selectByClick),
    clickHitTarget,
    selectInList: rounded(selectInList),
    sortHsps: rounded(sortHsps),
    scrollHspList: rounded(lastRow),
  };
}

test('many HSPs in one query-subject pair: the dot plot and the HSP list', async ({ page, browserName }) => {
  browser = browserName;
  await page.goto('/');
  // Run 1 is a small search to change runs with; the repeats are run 2, 3, ...
  const small = manyQueries(20);
  const first = await search(
    page,
    1,
    { name: 'queries-20.fna', buffer: small.queries },
    { name: 'subjects.fna', buffer: manyQueries(20).subjects },
  );
  void first;
  let number = 1;
  for (const copies of COPIES) {
    number++;
    const record: Record<string, unknown> = { unit: 20, copies };
    records.pair.push(record);
    try {
      const file = {
        name: `repeat-20x${copies}.fna`,
        buffer: repeats(20, copies),
      };
      Object.assign(record, { bytes: file.buffer.length });
      record['search'] = await search(page, number, file, file);
      console.log(`${browserName} ${copies} copies: search ${JSON.stringify(record['search'])}`);
      save();
      await page.getByTestId('tab-results').click();
      const warmup = await pairRepetition(page, number, 1);
      record['warmup'] = warmup;
      const samples: PairRepetition[] = [];
      for (let i = 0; i < REPETITIONS; i++) samples.push(await pairRepetition(page, number, 1));
      record['samples'] = samples;
      record['summary'] = summarize(samples as unknown as Record<string, unknown>[]);
      console.log(`${browserName} ${copies} copies: ${warmup.hsps} HSPs ${JSON.stringify(record['summary'])}`);
    } catch (error) {
      record['error'] = message(error);
      console.log(`${browserName} ${copies} copies: FAILED ${message(error)}`);
    }
    save();
  }
});
