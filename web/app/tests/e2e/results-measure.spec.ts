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
//    selecting the second query with 200 subjects by a click, showing its Graphic Summary (W4b)
//    and sorting its subject list (the Descriptions) by each column. Since W4b the results screen
//    has tabs: the selection and the HSP's detail are timed in the Alignments tab (the detail is
//    there), the subject list in the Descriptions. In Chromium also the JS heap (after a garbage
//    collection) before opening, after the first opening and after the repetitions, and the
//    memory of the page and its workers where the browser exposes it.
// 2. `pair`: a query that repeats a unit of 20 letters (with 1% differences between the copies)
//    searched against itself: every pair of copies is a diagonal, so one query-subject pair has
//    about two HSPs per copy. It records the search, the segments on the dot plot, and, in one
//    warm-up and three repetitions, the time to open the run, to show the dot plot, to zoom in,
//    out and to the selected HSP, to select an HSP by `n`, by a click on its segment and in the
//    list, to open the selected HSP's popup with Enter (W4b), to sort and to scroll the HSP list.
//
// Since W5 (S14: "抽出と候補トレイは、この規模で同じ操作の速さを保つ"), each measurement ends with the
// candidate tray filled from it, after its own records, which keep their names and their order so
// that they stay comparable with W4 and W4b. In one warm-up and three repetitions, from an empty
// tray: `queries` marks every subject of the first query with 200 subjects in the Descriptions
// ("select all") and adds their HSPs ("Add to candidates"); `pair`, for each count of copies whose
// measurement completed, adds every HSP of the pair in the Alignments ("Add all matches to
// candidates"). Then, the same in both (`tray` in the record): opening the Candidates tab, scrolling
// its list to the end, sorting it by Subject, extracting the hit regions of its first 100
// candidates (the clock stops at the frame that shows the extraction's summary; `downloadMs`, taken
// by the test from before the act to Playwright's download event, is an upper bound), and removing
// all of them.
//
// Since W6 (S15: "書き出しとセッションは、この規模で画面を止めない"), each measurement then ends with the
// files of its runs, after the records above (which keep their names and order): from the Outputs
// tab, one at a time, CSV, JSON (aligned rows on) and the report of the whole run (the scope that
// holds every HSP; no view filter is set), and the export of the stored outfmt 6; for `pair` also
// the dot plot as SVG; and last the session file: the tray filled again (`queries`: the HSPs of the
// first query with 200 subjects; `pair`: every HSP of the largest pair measured), "Save session",
// then the measured page closed and the file opened in a new page of the same browser. Each file
// records its size, the time from the act until the screen shows its summary or message (`ms`,
// `paintMs`; the outfmt 6 export and the SVG show none), the time until Playwright's download event
// (`downloadMs`, taken by the test: an upper bound), and the frames that the page drew meanwhile with
// the longest time between two of them (`maxFrameGapMs`): the page stays responsive while a large
// file is written or read, whatever the file's total time.
//
// Times are milliseconds, taken in the page with `performance.now()` from the action to the first
// animation frame in which the screen shows the result (`ms`) and to the frame after it
// (`paintMs`); the resolution is a frame (about 17 ms). The canvases (the dot plot, the Graphic
// Summary) count the frames they have drawn (`data-drawn`) and the dot plot the HSP drawn as
// selected (`data-drawn-selected`): their clocks wait for the drawing of their own change (W4b F1;
// W4 decision 14), and the page is first looked at after the act's own updates have asked for
// their frames. An action that the test performs through
// Playwright's mouse or keyboard (selecting an HSP by `n` or by a click on the dot plot) starts its
// clock at the event in the page, and the clock stops at the first frame that shows the result
// after the test starts looking, so these times include Playwright's round trip (upper bounds).
import { mkdirSync, statSync, writeFileSync } from 'node:fs';
import { readFile } from 'node:fs/promises';
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
  /** Chooses an option of a select element (its change event). */
  | { readonly kind: 'select'; readonly testid: string; readonly value: string }
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
  /** The element is laid out: neither it nor a parent has `display: none` (the tabs that stay mounted). */
  readonly visible?: boolean;
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

/** The distance in pixels from a point to the line segment `[ax, ay, bx, by]`. */
function segmentDistance(x: number, y: number, [ax, ay, bx, by]: readonly number[]): number {
  const lx = bx! - ax!;
  const ly = by! - ay!;
  const length = lx * lx + ly * ly;
  const t = length === 0 ? 0 : Math.max(0, Math.min(1, ((x - ax!) * lx + (y - ay!) * ly) / length));
  return Math.hypot(x - (ax! + t * lx), y - (ay! + t * ly));
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
      } else if (action.kind === 'select') {
        const select = find(action.testid) as HTMLSelectElement;
        select.value = action.value;
        select.dispatchEvent(new Event('change', { bubbles: true }));
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
          if (c.visible && element.getClientRects().length === 0) return false;
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
      // The act's updates (Vue's, in microtasks) mount components and ask for the frames that draw
      // them; the first look is registered after them, so that it comes after their drawing in the
      // next frame, and does not stop the clock in a frame that has not drawn yet (W4b F1).
      for (let i = 0; i < 5; i++) await Promise.resolve();
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

/**
 * Waits until the page has drawn 5 frames in a row less than 50 ms apart: what the act before
 * left to do (the Alignments' sections that arrive and are laid out after a selection) is done,
 * so that the next clock times its own act (W4 decision 14; W4b F1).
 */
async function quiet(page: Page): Promise<void> {
  await page.evaluate(
    () =>
      new Promise<void>((resolve) => {
        let last = performance.now();
        let calm = 0;
        const tick = (now: number) => {
          calm = now - last < 50 ? calm + 1 : 0;
          last = now;
          if (calm >= 5) resolve();
          else requestAnimationFrame(tick);
        };
        requestAnimationFrame(tick);
      }),
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
  showGraphic: Timing;
  graphicRows: number;
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
  // The small run (run 2) first, so that opening run 1 is a change of run. The Alignments tab,
  // which the next run opens on, shows the selected HSP's detail.
  await page.getByTestId('run-2-open').click();
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '2', { timeout: 120_000 });
  await page.getByTestId('results-view-hits').click();
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
    { testid: 'alignments-subject' },
  ]);
  // Its Graphic Summary: the first 100 of its subjects drawn on the canvas.
  const showGraphic = await step(page, { kind: 'click', testid: 'results-view-hits' }, [
    { testid: 'graphic-canvas', attr: 'data-rows', matches: '^[1-9]\\d{2,}$' },
    { testid: 'graphic-canvas', attr: 'data-drawn' },
  ]);
  const graphicRows = Number(await page.getByTestId('graphic-canvas').getAttribute('data-rows'));
  await page.getByTestId('results-view-hits').click();
  await listWhole(page);
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
    showGraphic: rounded(showGraphic),
    graphicRows,
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
      save();

      // The candidate tray (W5), after the records above.
      const tray: Record<string, unknown> = {};
      record['tray'] = tray;
      try {
        const subjects = await wideQuery(page, count);
        await trayRepetitions(tray, () => queriesTrayRepetition(page, subjects));
        console.log(`${browserName} ${count} queries: tray ${JSON.stringify(tray['summary'])}`);
      } catch (error) {
        tray['error'] = message(error);
        console.log(`${browserName} ${count} queries: tray FAILED ${message(error)}`);
      }
      save();

      // The run's files and the session file (W6), after the records above.
      const files: Record<string, unknown> = {};
      record['files'] = files;
      try {
        await page.getByTestId('run-1-open').click();
        await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1', { timeout: 120_000 });
        Object.assign(files, await runFilesOf(page, 1));
        console.log(`${browserName} ${count} queries: files ${JSON.stringify(files)}`);
      } catch (error) {
        files['error'] = message(error);
        console.log(`${browserName} ${count} queries: files FAILED ${message(error)}`);
      }
      save();
      const sessionFile: Record<string, unknown> = {};
      record['session'] = sessionFile;
      try {
        // The tray: the HSPs of the first query with 200 subjects.
        await page.getByTestId('tab-results').click();
        await wideQuery(page, count);
        await page.getByTestId('descriptions-select-all').check();
        await page.getByTestId('descriptions-add-candidates').click();
        await expect(page.getByTestId('tab-candidates-count')).toHaveText(/^[1-9][\d,]*$/);
        await sessionRoundTrip(page, sessionFile);
        console.log(`${browserName} ${count} queries: session ${JSON.stringify(sessionFile)}`);
      } catch (error) {
        sessionFile['error'] = message(error);
        console.log(`${browserName} ${count} queries: session FAILED ${message(error)}`);
      }
    } catch (error) {
      record['error'] = message(error);
      console.log(`${browserName} ${count} queries: FAILED ${message(error)}`);
    }
    save();
  });
}

// --- the candidate tray (W5), filled from either measurement ------------------------------------------

/** The candidates whose hit regions are extracted. */
const EXTRACTED = 100;
const candidateRows = (page: Page) => page.getByTestId('candidate-list').locator('[data-key]').count();

interface TrayRepetition {
  candidates: number;
  openTray: Timing;
  scrollToEnd: Timing;
  rowsDrawnAtEnd: number;
  sortBySubject: Timing;
  extract: Timing;
  /** From before the act until Playwright's download event, taken by the test: an upper bound. */
  downloadMs: number;
  sequences: number;
  bytes: number;
  removeAll: Timing;
}

/**
 * The tray's operations on the `count` candidates just added (all selected, in the order added),
 * while the results screen is shown; the tray is empty again at the end.
 */
async function trayOperations(page: Page, count: number): Promise<TrayRepetition> {
  const total = countText(count);
  // The sort's clock looks for the first row, which is drawn at the end of the list only for a short list.
  expect(count).toBeGreaterThan(EXTRACTED);
  // The tray stays mounted, hidden, while the results are shown (its rows are drawn as candidates
  // are added): the tab shows it in place of the results screen.
  await quiet(page);
  const openTray = await step(page, { kind: 'click', testid: 'tab-candidates' }, [
    { testid: 'candidates', visible: true },
    { testid: 'candidate-list', attr: 'data-count', equals: String(count) },
    { testid: 'candidate-1' },
  ]);
  await quiet(page);
  const scrollToEnd = await step(page, { kind: 'scroll', testid: 'candidate-list', to: 'end' }, [{ testid: `candidate-${count}` }]);
  const rowsDrawnAtEnd = await candidateRows(page);
  // A sort shows the list's first rows again in the same update: row 1, not drawn at the end, is drawn in the new order.
  await quiet(page);
  const sortBySubject = await step(page, { kind: 'select', testid: 'candidates-sort', value: 'subject' }, [{ testid: 'candidate-1' }]);
  await expect(page.getByTestId('candidates-sort')).toHaveValue('subject');

  // The first 100 in that order: none selected, then a click on each one's mark (each scrolled into view).
  await page.getByTestId('candidates-select-all').uncheck();
  await expect(page.getByTestId('candidates-selected')).toHaveText(`0 of ${total} candidates selected`);
  for (let n = 1; n <= EXTRACTED; n++) await page.getByTestId(`candidate-mark-${n}`).check();
  await expect(page.getByTestId('candidates-selected')).toHaveText(`${EXTRACTED} of ${total} candidates selected`);
  await expect(page.getByTestId('extract-region-hit')).toBeChecked();
  await expect(page.getByTestId('extract-join-separate')).toBeChecked();
  // The summary of the extraction before stays until the next one replaces it: the clock waits for a new element.
  await page.evaluate(() => document.querySelector('[data-testid="extract-summary"]')?.setAttribute('data-measured', 'before'));
  await quiet(page);
  const arrival = page.waitForEvent('download', { timeout: 600_000 }).then((download) => ({ download, at: Date.now() }));
  const t0 = Date.now();
  const extract = await step(page, { kind: 'click', testid: 'extract-download' }, [
    { testid: 'extract-summary', attr: 'data-measured', differs: 'before' },
    { testid: 'extract-summary', text: `Saved losat-candidates\\.fa: ${EXTRACTED} sequences\\s` },
  ]);
  const { download, at } = await arrival;
  const text = await readFile((await download.path())!, 'utf8');
  const sequences = text.split('\n').filter((line) => line.startsWith('>')).length;
  expect(sequences).toBe(EXTRACTED);

  await page.getByTestId('candidates-select-all').check();
  await expect(page.getByTestId('candidates-selected')).toHaveText(`${total} of ${total} candidates selected`);
  await quiet(page);
  const removeAll = await step(page, { kind: 'click', testid: 'candidates-remove-selected' }, [
    { testid: 'candidates-empty' },
    { testid: 'tab-candidates-count', text: '^0$' },
  ]);
  return {
    candidates: count,
    openTray: rounded(openTray),
    scrollToEnd: rounded(scrollToEnd),
    rowsDrawnAtEnd,
    sortBySubject: rounded(sortBySubject),
    extract: rounded(extract),
    downloadMs: at - t0,
    sequences,
    bytes: Buffer.byteLength(text),
    removeAll: rounded(removeAll),
  };
}

/**
 * The Descriptions of the selected query listed whole, as W5 measured them: since 2026-10-10 the
 * one page lists the first 100 subjects until "Show all", which stays for the query once pressed.
 */
async function listWhole(page: Page): Promise<void> {
  const [, shown] = /^([\d,]+) shown$/.exec(((await page.getByTestId('subject-count').textContent()) ?? '').trim())!;
  const subjects = shown!.replace(/,/g, '');
  if ((await page.getByTestId('subject-list').getAttribute('data-count')) !== subjects) await page.getByTestId('descriptions-show-all').click();
  await expect(page.getByTestId('subject-list')).toHaveAttribute('data-count', subjects);
}

/** One warm-up and the repetitions, written into `into` as they come (a failure keeps what was measured). */
async function trayRepetitions<T extends object>(into: Record<string, unknown>, repeat: () => Promise<T>): Promise<void> {
  into['warmup'] = await repeat();
  const samples: T[] = [];
  into['samples'] = samples;
  for (let i = 0; i < REPETITIONS; i++) samples.push(await repeat());
  into['summary'] = summarize(samples as unknown as Record<string, unknown>[]);
}

/**
 * Selects the first query with 200 subjects of run 1 (`query-filter` selects it; decision 10) and
 * shows its Descriptions, with the filter cleared again; returns the number of its subjects.
 */
async function wideQuery(page: Page, count: number): Promise<number> {
  const total = countText(count);
  await page.getByTestId('results-view-hits').click();
  await page.getByTestId('query-filter').fill('wide');
  await expect(page.getByTestId('query-count')).toHaveText(new RegExp(`^[\\d,]+ of ${total} queries$`));
  await expect(page.locator('[data-testid^="query-row-"]').first()).toHaveAttribute('aria-pressed', 'true');
  await page.getByTestId('query-filter').fill('');
  await expect(page.getByTestId('query-count')).toHaveText(`${total} of ${total} queries`);
  await listWhole(page);
  return Number(await page.getByTestId('subject-list').getAttribute('data-count'));
}

interface QueriesTrayRepetition extends TrayRepetition {
  markAll: Timing;
  addMarked: Timing;
}

/** Marks every subject of the selected query ("select all"), adds them, and the tray's operations. */
async function queriesTrayRepetition(page: Page, subjects: number): Promise<QueriesTrayRepetition> {
  await page.getByTestId('tab-results').click();
  await expect(page.getByTestId('subject-list')).toHaveAttribute('data-count', String(subjects));
  await expect(page.getByTestId('descriptions-selected')).toHaveText('0 sequences selected');
  await expect(page.getByTestId('tab-candidates-count')).toHaveText('0');
  await quiet(page);
  const markAll = await step(page, { kind: 'click', testid: 'descriptions-select-all' }, [
    { testid: 'descriptions-selected', text: `^${countText(subjects)} sequences selected$` },
  ]);
  await quiet(page);
  const addMarked = await step(page, { kind: 'click', testid: 'descriptions-add-candidates' }, [
    { testid: 'tab-candidates-count', text: '^[1-9][\\d,]*$' },
  ]);
  const added = Number((await page.getByTestId('tab-candidates-count').textContent())!.replace(/,/g, ''));
  // The marks stay after the addition; the next repetition marks the subjects again.
  await page.getByTestId('descriptions-select-all').uncheck();
  await expect(page.getByTestId('descriptions-selected')).toHaveText('0 sequences selected');
  return { markAll: rounded(markAll), addMarked: rounded(addMarked), ...(await trayOperations(page, added)) };
}

/** Shows run `number` (a pair) in the Alignments, opened after run 1 so that the clocks see a change of run. */
async function openPair(page: Page, number: number): Promise<void> {
  for (const run of [1, number]) {
    await page.getByTestId(`run-${run}-open`).click();
    await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', String(run), { timeout: 120_000 });
  }
  await page.getByTestId('results-view-hits').click();
}

interface PairTrayRepetition extends TrayRepetition {
  addAll: Timing;
}

/** Adds every HSP of the pair shown ("Add all matches to candidates"), and the tray's operations. */
async function pairTrayRepetition(page: Page): Promise<PairTrayRepetition> {
  await page.getByTestId('tab-results').click();
  await expect(page.getByTestId('alignments-add-subject')).toHaveText('Add all matches to candidates');
  await expect(page.getByTestId('tab-candidates-count')).toHaveText('0');
  const hsps = Number(await page.getByTestId('hsp-list').getAttribute('data-count'));
  await quiet(page);
  const addAll = await step(page, { kind: 'click', testid: 'alignments-add-subject' }, [
    { testid: 'tab-candidates-count', text: `^${countText(hsps)}$` },
    { testid: 'alignments-add-subject', text: '^\\s*All matches in candidates\\s*$' },
  ]);
  return { addAll: rounded(addAll), ...(await trayOperations(page, hsps)) };
}

// --- the run's files and the session file (W6) ------------------------------------------------------

/** A large file is written or read within this time, or the measurement of it fails. */
const FILE_TIMEOUT_MS = 1_800_000;

interface Frames {
  /** Frames drawn from the act until the file arrived, and the longest time between two of them. */
  readonly frames: number;
  readonly maxFrameGapMs: number;
  /** When the longest gap ended, from the act (where in the work the page stopped drawing). */
  readonly maxFrameGapEndMs: number;
  /**
   * The longest delay of a 10 ms timer on the page's main thread: a long task of the page shows
   * here and as a frame gap; a frame gap without it is the browser not drawing (not the page's script).
   */
  readonly maxTaskLagMs: number;
}

/** Starts counting the frames that the page draws and the longest gap between two of them (a long task shows as a long gap). */
async function watchFrames(page: Page): Promise<void> {
  await page.evaluate(() => {
    const start = performance.now();
    const frames = { on: true, count: 0, maxGap: 0, maxGapEnd: 0, start, last: start, maxLag: 0 };
    (window as unknown as { __frames: typeof frames }).__frames = frames;
    const TIMER_MS = 10;
    const timer = (due: number) => {
      if (!frames.on) return;
      const now = performance.now();
      frames.maxLag = Math.max(frames.maxLag, now - due);
      // When the timer is due: taken now, when it is set (fix round 2: taken when it fired, the
      // delay was always about 0, and only the first timer counted).
      const next = now + TIMER_MS;
      setTimeout(() => timer(next), TIMER_MS);
    };
    setTimeout(() => timer(start + TIMER_MS), TIMER_MS);
    const tick = (now: number) => {
      if (!frames.on) return;
      if (now - frames.last > frames.maxGap) {
        frames.maxGap = now - frames.last;
        frames.maxGapEnd = now - frames.start;
      }
      frames.last = now;
      frames.count++;
      requestAnimationFrame(tick);
    };
    requestAnimationFrame(tick);
  });
}

async function framesSeen(page: Page): Promise<Frames> {
  return page.evaluate(() => {
    const frames = (
      window as unknown as { __frames: { on: boolean; count: number; maxGap: number; maxGapEnd: number; start: number; last: number; maxLag: number } }
    ).__frames;
    frames.on = false;
    const now = performance.now();
    if (now - frames.last > frames.maxGap) {
      frames.maxGap = now - frames.last;
      frames.maxGapEnd = now - frames.start;
    }
    return {
      frames: frames.count,
      maxFrameGapMs: Math.round(frames.maxGap),
      maxFrameGapEndMs: Math.round(frames.maxGapEnd),
      maxTaskLagMs: Math.round(frames.maxLag),
    };
  });
}

interface FileRecord extends Frames {
  readonly name: string;
  readonly bytes: number;
  /** From the act until the summary or message shows (none for the outfmt 6 export and the SVG). */
  readonly ms?: number;
  readonly paintMs?: number;
  /** From before the act until Playwright's download event, taken by the test: an upper bound. */
  readonly downloadMs: number;
}

/**
 * Clicks a button that saves a file and waits until the screen shows `until` (none: the download
 * only) and the file has arrived. A summary element `marker` that an earlier file left is marked,
 * so that the clock waits for the new one. The file is deleted, or kept at `keep`.
 */
async function timedFile(page: Page, testid: string, until: readonly Cond[], marker?: string, keep?: string): Promise<FileRecord> {
  if (marker !== undefined) await page.evaluate((id) => document.querySelector(`[data-testid="${id}"]`)?.setAttribute('data-measured', 'before'), marker);
  await quiet(page);
  const arrival = page.waitForEvent('download', { timeout: FILE_TIMEOUT_MS }).then((download) => ({ download, at: Date.now() }));
  await watchFrames(page);
  const t0 = Date.now();
  let timing: Timing | undefined;
  if (until.length > 0) timing = rounded(await step(page, { kind: 'click', testid }, until, FILE_TIMEOUT_MS));
  else await page.evaluate((id) => (document.querySelector(`[data-testid="${id}"]`) as HTMLElement).click(), testid);
  const { download, at } = await arrival;
  const path = await download.path();
  const frames = await framesSeen(page);
  const name = download.suggestedFilename();
  const bytes = statSync(path).size;
  if (keep !== undefined) await download.saveAs(keep);
  await download.delete();
  return { name, bytes, ...(timing ?? {}), downloadMs: at - t0, ...frames };
}

/** The run's own files of the whole run and its outfmt 6 export, from the Outputs tab of run `number` (shown), one at a time. */
async function runFilesOf(page: Page, number: number): Promise<Record<string, unknown>> {
  await page.getByTestId('tab-results').click();
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', String(number));
  const shown = Date.now();
  await page.getByTestId('results-view-outputs').click();
  await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', `${number}:6`, { timeout: FILE_TIMEOUT_MS });
  const outputsShownMs = Date.now() - shown;
  const hsps = Number(await page.getByTestId('export-scope-all').getAttribute('data-count'));
  await expect(page.getByTestId('export-scope-all')).toBeChecked();
  await expect(page.getByTestId('export-json-aligned')).toBeChecked();
  const files: Record<string, unknown> = {
    scope: 'all',
    hsps,
    filteredHsps: Number(await page.getByTestId('export-scope-filtered').getAttribute('data-count')),
    outputsShownMs,
  };
  for (const format of ['csv', 'json', 'report'] as const) {
    files[format] = await timedFile(
      page,
      `export-${format}`,
      [
        { testid: 'export-summary', attr: 'data-measured', differs: 'before' },
        { testid: 'export-summary', attr: 'data-format', equals: format },
        { testid: 'export-summary', attr: 'data-hsps', equals: String(hsps) },
      ],
      'export-summary',
    );
    note(`run ${number}: ${format} ${JSON.stringify(files[format])}`);
  }
  files['outfmt6'] = await timedFile(page, 'export-output', []);
  return files;
}

/** The dot plot of the pair shown, saved as SVG (the whole plot, every HSP). */
async function dotPlotSvg(page: Page): Promise<FileRecord & { segments: number }> {
  await page.getByTestId('pane-dotplot').click();
  const hsps = Number(await page.getByTestId('hsp-list').getAttribute('data-count'));
  await expect(page.getByTestId('dotplot-canvas')).toHaveAttribute('data-segments', String(hsps), { timeout: 120_000 });
  const file = await timedFile(page, 'dotplot-svg', []);
  await expect(page.getByTestId('dotplot-svg-error')).toHaveCount(0);
  return { ...file, segments: hsps };
}

/**
 * Saves the session file of the page's completed runs with the tray's candidates, closes the
 * page, and opens the file in a new page of the same browser (a new working session).
 */
async function sessionRoundTrip(page: Page, into: Record<string, unknown>): Promise<void> {
  into['candidates'] = Number((await page.getByTestId('tab-candidates-count').textContent())!.replace(/,/g, ''));
  into['runs'] = await page.getByTestId('queue').locator('li[data-status="completed"]').count();
  await expect(page.getByTestId('session-include-candidates')).toBeChecked();
  const kept = join(test.info().outputPath('session'), 'measured.losat-session.gz');
  mkdirSync(join(kept, '..'), { recursive: true });
  into['save'] = await timedFile(page, 'session-save', [{ testid: 'session-message', text: '^Saved ' }], undefined, kept);
  into['saveMessage'] = await page.getByTestId('session-message').textContent();
  save();
  const context = page.context();
  await page.close();
  const fresh = await context.newPage();
  try {
    await fresh.goto('/');
    await expect(fresh.getByTestId('storage-status')).toBeVisible();
    await quiet(fresh);
    await arm(fresh, 'session-file', 'change');
    await watchFrames(fresh);
    await fresh.getByTestId('session-file').setInputFiles(kept);
    const open = await step(fresh, { kind: 'armed' }, [{ testid: 'session-message', text: '^(Opened |.* was not opened)' }], FILE_TIMEOUT_MS);
    const frames = await framesSeen(fresh);
    const text = (await fresh.getByTestId('session-message').textContent())!.trim();
    into['open'] = { ...rounded(open), ...frames };
    into['openMessage'] = text;
    if (!text.startsWith('Opened ')) throw new Error(text);
    into['openedRuns'] = await fresh.getByTestId('queue').locator('li[data-status="completed"]').count();
    into['openedCandidates'] = Number((await fresh.getByTestId('tab-candidates-count').textContent())!.replace(/,/g, ''));
  } finally {
    await fresh.close();
  }
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
  popupByKey: Timing;
  selectInList: Timing;
  sortHsps: Timing;
  scrollHspList: Timing;
}

async function pairRepetition(page: Page, number: number, other: number): Promise<PairRepetition> {
  await page.getByTestId(`run-${other}-open`).click();
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', String(other), { timeout: 120_000 });
  // The run opens on the Alignments tab, which shows the selected HSP's detail (W4b).
  await page.getByTestId('results-view-hits').click();
  const open = await step(page, { kind: 'click', testid: `run-${number}-open` }, [
    { testid: 'results-hits', attr: 'data-run', equals: String(number) },
    { testid: 'hsp-detail', attr: 'data-state', equals: 'ready' },
  ]);
  await page.getByTestId('results-view-hits').click();
  const hsps = Number(await page.getByTestId('hsp-list').getAttribute('data-count'));
  const hspRowsDrawn = await rowsDrawn(page, 'hsp-row-');

  // Each clock of the dot plot stops at the frame that drew its change (`data-drawn`).
  const drew: Cond = { testid: 'dotplot-canvas', attr: 'data-drawn', changed: true };
  const showDotPlot = await step(page, { kind: 'click', testid: 'pane-dotplot' }, [
    { testid: 'dotplot-canvas', attr: 'data-segments', equals: String(hsps) },
    drew,
  ]);
  const canvas = page.getByTestId('dotplot-canvas');
  const segments = Number(await canvas.getAttribute('data-segments'));
  // Pixels in the colours of the HSPs (blast2dotplot.py's #1f77b4 and #ff7f0e at the opacities of
  // the identity classes, over the white plot): the plot is drawn.
  const coloured = await canvas.evaluate((element) => {
    const c = element as HTMLCanvasElement;
    const { data } = c.getContext('2d')!.getImageData(0, 0, c.width, c.height);
    const tints = [
      [31, 119, 180],
      [255, 127, 14],
    ].flatMap((rgb) => [1, 0.8, 0.6, 0.4].map((alpha) => rgb.map((v) => Math.round(255 - alpha * (255 - v)))));
    let n = 0;
    for (let i = 0; i < data.length; i += 4) {
      const [r, g, b] = [data[i]!, data[i + 1]!, data[i + 2]!];
      if (tints.some(([tr, tg, tb]) => Math.abs(r - tr!) + Math.abs(g - tg!) + Math.abs(b - tb!) < 30)) n++;
    }
    return n;
  });
  const whole = (await canvas.getAttribute('data-view'))!;
  const zoomIn = await step(page, { kind: 'click', testid: 'dotplot-zoom-in' }, [
    { testid: 'dotplot-canvas', attr: 'data-view', changed: true },
    drew,
  ]);
  const zoomOut = await step(page, { kind: 'click', testid: 'dotplot-zoom-out' }, [
    { testid: 'dotplot-canvas', attr: 'data-view', equals: whole },
    drew,
  ]);

  // Selection: the key n, a click on a segment (the first 200 are offered to tests), a click in the list.
  const selected = (await canvas.getAttribute('data-selected'))!;
  await canvas.focus();
  await arm(page, 'dotplot-canvas', 'keydown');
  await page.keyboard.press('n');
  // The act happened before the clock looks: it waits until the canvas has drawn another HSP as selected.
  const selectByKey = await step(
    page,
    { kind: 'armed' },
    [{ testid: 'dotplot-canvas', attr: 'data-drawn-selected', differs: selected }],
    30_000,
  ).catch(() => {
    throw new Error(`n did not change the selected HSP (${selected})`);
  });
  const targets = JSON.parse((await canvas.getAttribute('data-targets'))!) as {
    hsp: string;
    x: number;
    y: number;
    ends: [number, number, number, number];
  }[];
  const now = (await canvas.getAttribute('data-selected'))!;
  // The click selects the nearest line, so it aims at the HSP whose midpoint lies farthest from the
  // selected HSP's line: where lines lie a pixel apart (the repeats), a target next to the selected
  // line would leave the selection unchanged.
  const current = targets.find((t) => t.hsp === now);
  const away = (t: { x: number; y: number }) => (current === undefined ? 0 : segmentDistance(t.x, t.y, current.ends));
  const target = targets
    .filter((t) => t.hsp !== now && t.hsp !== selected)
    .reduce((best, t) => (away(t) > away(best) ? t : best));
  await arm(page, 'dotplot-canvas', 'pointerdown');
  await canvas.click({ position: { x: target.x, y: target.y } });
  const selectByClick = await step(page, { kind: 'armed' }, [{ testid: 'dotplot-canvas', attr: 'data-drawn-selected', differs: now }], 30_000);
  // Where segments lie closer together than the click's reach, the click picks the nearest of them, not always the one aimed at.
  const clickHitTarget = (await canvas.getAttribute('data-selected')) === target.hsp ? 1 : 0;

  // The popup of the selected HSP (W4b): the click opened it; closed, Enter on the plot opens it again.
  const popup = page.getByTestId('dotplot-popup');
  await expect(popup).toHaveCount(1);
  await canvas.focus();
  await page.keyboard.press('Escape');
  await expect(popup).toHaveCount(0);
  const chosen = (await canvas.getAttribute('data-selected'))!;
  const drawnBefore = (await canvas.getAttribute('data-drawn'))!;
  await arm(page, 'dotplot-canvas', 'keydown');
  await page.keyboard.press('Enter');
  // The popup is placed by the frame that the key asked for.
  const popupByKey = await step(
    page,
    { kind: 'armed' },
    [
      { testid: 'dotplot-popup', attr: 'data-hsp', equals: chosen },
      { testid: 'dotplot-canvas', attr: 'data-drawn', differs: drawnBefore },
    ],
    30_000,
  );
  await page.keyboard.press('Escape');
  await expect(popup).toHaveCount(0);

  // "Zoom to HSP" frames the selected HSP. An HSP that spans nearly the whole of both sequences hardly changes the view, so the clock stops at the frame drawn after the click.
  const zoomToHsp = await step(page, { kind: 'click', testid: 'dotplot-zoom-hsp' }, [drew]);
  const reset = await step(page, { kind: 'click', testid: 'dotplot-reset' }, [
    { testid: 'dotplot-canvas', attr: 'data-view', equals: whole },
    drew,
  ]);

  await page.getByTestId('results-view-hits').click();
  const row = page.getByTestId('hsp-list').locator('[data-testid^="hsp-row-"][aria-pressed="false"]').nth(2);
  await quiet(page);
  const rowId = (await row.getAttribute('data-testid'))!.replace('hsp-row-', '').replace('-', ':');
  const selectInList = await step(page, { kind: 'click', testid: `hsp-row-${rowId.replace(':', '-')}` }, [
    { testid: 'hsp-detail', attr: 'data-hsp', equals: rowId },
    { testid: 'hsp-detail', attr: 'data-state', equals: 'ready' },
  ]);
  await quiet(page);
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
    popupByKey: rounded(popupByKey),
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
  /** The pairs whose measurement completed, for the candidate tray. */
  const measured: { number: number; copies: number; record: Record<string, unknown> }[] = [];
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
      measured.push({ number, copies, record });
    } catch (error) {
      record['error'] = message(error);
      console.log(`${browserName} ${copies} copies: FAILED ${message(error)}`);
    }
    save();
  }

  // The candidate tray (W5) of each pair, after the records above (which keep the order of W4b).
  for (const { number: run, copies, record } of measured) {
    const tray: Record<string, unknown> = {};
    record['tray'] = tray;
    try {
      await openPair(page, run);
      await trayRepetitions(tray, () => pairTrayRepetition(page));
      console.log(`${browserName} ${copies} copies: tray ${JSON.stringify(tray['summary'])}`);
    } catch (error) {
      tray['error'] = message(error);
      console.log(`${browserName} ${copies} copies: tray FAILED ${message(error)}`);
    }
    save();
  }

  // The files of each pair and its dot plot as SVG (W6), after the records above.
  for (const { number: run, copies, record } of measured) {
    const files: Record<string, unknown> = {};
    record['files'] = files;
    try {
      await openPair(page, run);
      Object.assign(files, await runFilesOf(page, run));
      files['svg'] = await dotPlotSvg(page);
      console.log(`${browserName} ${copies} copies: files ${JSON.stringify(files)}`);
    } catch (error) {
      files['error'] = message(error);
      console.log(`${browserName} ${copies} copies: files FAILED ${message(error)}`);
    }
    save();
  }

  // The session file of the runs, with every HSP of the largest pair measured in the tray (W6).
  const largest = measured.at(-1);
  if (largest !== undefined) {
    const session: Record<string, unknown> = {};
    largest.record['session'] = session;
    try {
      await openPair(page, largest.number);
      await expect(page.getByTestId('alignments-add-subject')).toHaveText('Add all matches to candidates');
      await page.getByTestId('alignments-add-subject').click();
      await expect(page.getByTestId('alignments-add-subject')).toHaveText(/All matches in candidates/, { timeout: 120_000 });
      await sessionRoundTrip(page, session);
      console.log(`${browserName} ${largest.copies} copies: session ${JSON.stringify(session)}`);
    } catch (error) {
      session['error'] = message(error);
      console.log(`${browserName} ${largest.copies} copies: session FAILED ${message(error)}`);
    }
    save();
  }
});
