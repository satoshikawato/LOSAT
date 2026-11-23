// Runs walks of walk/plan.ts with the steps of walk/steps.ts, and shrinks a failing walk to a
// small reproduction (walk/shrink.ts) by replaying parts of it, each in a new browser context (no
// storage, worker or lock left by the walk or an earlier replay), within a time limit.
// A walk's steps and their draws are attached to the test as walk-steps.json, and a shrunk
// reproduction as walk-shrunk.json; LOSAT_WEB_WALK_REPLAY=<one of them> replays it
// (random-walk.spec.ts).
import { readFileSync, writeFileSync } from 'node:fs';
import type { Browser, BrowserContextOptions, Page, TestInfo } from '@playwright/test';
import { BUILD_HAS_ENGINE } from '../support/browser';
import { NATIVE } from '../support/native';
import { standInWakeLock } from '../support/test';
import { makeFixtures, type FixtureLibrary } from './fixtures';
import { createRandom, type PlannedStep, type WalkKind } from './plan';
import { shrinkWalk } from './shrink';
import { check, Draws, open, STEPS, StepFailure, watch, type RecordedStep, type Walker } from './steps';

/** LOSAT_WEB_WALK_TRACE=1 prints each step as it ends. */
const TRACE = process.env['LOSAT_WEB_WALK_TRACE'] === '1';

export interface WalkFailure {
  readonly index: number;
  readonly kind: WalkKind;
  readonly signature: string;
  readonly message: string;
}

export interface WalkResult {
  readonly steps: readonly RecordedStep[];
  readonly notes: readonly string[];
  readonly counts?: Walker['counts'];
  readonly failure?: WalkFailure;
}

/** A walk as a file: what walk-steps.json and walk-shrunk.json hold. */
export interface WalkFile {
  readonly seed: number;
  readonly walk: number;
  readonly browser: string;
  readonly steps: readonly RecordedStep[];
  readonly failure?: WalkFailure;
  /** A recorded walk's least number of files with records left out that the page must offer and the CLI run on (line WQA's D22). */
  readonly minOfferedLeftOut?: number;
}

/**
 * Runs `steps` in `page` from a fresh start. Steps without draws draw them from the walk's
 * generator (`seed`); recorded draws are replayed. Stops at the first failure.
 */
export async function runSteps(
  page: Page,
  steps: ReadonlyArray<PlannedStep | RecordedStep>,
  { seed, fixtures, outDir }: { seed: number; fixtures: FixtureLibrary; outDir: string },
): Promise<WalkResult> {
  const random = createRandom(seed);
  const walker = await watch(page, fixtures, outDir);
  const recorded: RecordedStep[] = [];
  const fail = (index: number, kind: WalkKind, error: unknown): WalkResult => {
    const failure: WalkFailure =
      error instanceof StepFailure
        ? { index, kind, signature: error.signature, message: error.message }
        : { index, kind, signature: errorSignature(error), message: error instanceof Error ? (error.stack ?? error.message) : String(error) };
    return { steps: recorded, notes: walker.notes, counts: walker.counts, failure };
  };
  try {
    await open(walker);
    await check(walker);
  } catch (error) {
    return fail(-1, 'reload', error);
  }
  for (const [index, step] of steps.entries()) {
    const draws = new Draws('args' in step ? step.args : undefined, random);
    let note = '';
    const started = Date.now();
    try {
      note = await STEPS[step.kind](walker, draws);
      if (TRACE) console.log(`step ${index + 1} ${step.kind}: ${note} (${Date.now() - started} ms)`);
      recorded.push({ kind: step.kind, args: [...draws.taken], note });
      await check(walker);
    } catch (error) {
      if (recorded.length === index) recorded.push({ kind: step.kind, args: [...draws.taken], note: 'failed' });
      return fail(index, step.kind, error);
    }
  }
  return { steps: recorded, notes: walker.notes, counts: walker.counts };
}

/**
 * An error that is not a step's own failure, named by its first line without numbers (a timeout's
 * milliseconds, a run's number), so that a replay reproduces the same error and not just any.
 */
function errorSignature(error: unknown): string {
  const line = (error instanceof Error ? error.message : String(error)).split('\n', 1)[0]!;
  return `error:${line.replace(/\d+/g, 'N').slice(0, 100)}`;
}

/** The options of the test's own context, for the new contexts of a shrink's replays. */
const CONTEXT_OPTIONS = ['baseURL', 'viewport', 'userAgent', 'deviceScaleFactor', 'isMobile', 'hasTouch', 'locale', 'timezoneId'] as const;

/** Replays a candidate in a new context: whether it fails at its last step as `failure` did. */
async function reproduces(
  browser: Browser,
  contextOptions: BrowserContextOptions,
  steps: readonly RecordedStep[],
  failure: WalkFailure,
  options: { seed: number; fixtures: FixtureLibrary; outDir: string },
): Promise<boolean> {
  const context = await browser.newContext(contextOptions);
  try {
    if (browser.browserType().name() === 'webkit') await context.addInitScript(standInWakeLock);
    const page = await context.newPage();
    const result = await runSteps(page, steps, options);
    return result.failure !== undefined && result.failure.index === steps.length - 1 && result.failure.signature === failure.signature;
  } finally {
    await context.close().catch(() => undefined);
  }
}

export interface ShrinkLimits {
  /** Replays per failing walk. */
  readonly replays: number;
  /** The time that a failing walk's shrinking may take, in ms. */
  readonly ms: number;
}

/** The shrink limits of a run: LOSAT_WEB_WALK_SHRINK_REPLAYS, else the default. */
export function shrinkLimits(defaults: ShrinkLimits): ShrinkLimits {
  const replays = Number(process.env['LOSAT_WEB_WALK_SHRINK_REPLAYS'] ?? defaults.replays);
  return { replays: Number.isInteger(replays) && replays >= 0 ? replays : defaults.replays, ms: defaults.ms };
}

/**
 * Runs one walk of a plan as a test body: the walk, then, if it failed, its shrinking; attaches
 * walk-steps.json (and walk-shrunk.json), and fails with the failure and the replay's command.
 */
export async function walkTest(
  { page, browser, browserName }: { page: Page; browser: Browser; browserName: string },
  testInfo: TestInfo,
  { seed, walk, steps, shrink, minOfferedLeftOut }: { seed: number; walk: number; steps: ReadonlyArray<PlannedStep | RecordedStep>; shrink: ShrinkLimits; minOfferedLeftOut?: number },
): Promise<void> {
  const fixtures = makeFixtures(seed);
  const outDir = testInfo.outputPath('native');
  const options = { seed: seed + walk, fixtures, outDir };
  const result = await runSteps(page, steps, options);
  // A walk that never searched, or never compared with the CLI it was given, checked too little:
  // each walk opens with a run of the inputs it starts with (plan.ts OPENING).
  if (result.failure === undefined && result.counts !== undefined && BUILD_HAS_ENGINE) {
    if (result.counts.completed === 0) throw new Error(`walk ${walk} (seed ${seed}) completed no run in ${result.steps.length} steps`);
    if (NATIVE !== undefined && result.counts.compared === 0) throw new Error(`walk ${walk} (seed ${seed}) compared no output with the native CLI`);
    // A recorded walk's steps are draws: a fixture added to the walk's files moves them, and the walk may then reach no file with records left out.
    if (NATIVE !== undefined && minOfferedLeftOut !== undefined && result.counts.offeredLeftOut < minOfferedLeftOut) {
      throw new Error(`walk ${walk} (seed ${seed}) ran the CLI on ${result.counts.offeredLeftOut} files with records left out, fewer than the ${minOfferedLeftOut} that it was recorded for`);
    }
  }
  const file: WalkFile = { seed, walk, browser: browserName, steps: result.steps, ...(result.failure ? { failure: result.failure } : {}) };
  // In the test's output directory as well, where a passing walk's steps can be read.
  writeFileSync(testInfo.outputPath('walk-steps.json'), JSON.stringify(file, null, 1));
  await testInfo.attach('walk-steps.json', { body: JSON.stringify(file, null, 1), contentType: 'application/json' });
  for (const note of result.notes) testInfo.annotations.push({ type: 'walk note', description: note });
  if (result.counts !== undefined) {
    const { compared, completed, failed, cancelled, refused, unmatched, offered, offeredLeftOut } = result.counts;
    const line = `runs: ${completed} completed, ${failed} failed, ${cancelled} cancelled, ${refused} presses refused; outputs equal to the native CLI's: ${compared}; completed runs not compared (inputs of one name): ${unmatched}; files of LOSAT's that the page offered and the CLI ran on (records left out, several files, pasted text): ${offered}, with records left out: ${offeredLeftOut}`;
    testInfo.annotations.push({ type: 'walk runs', description: line });
    if (TRACE) console.log(line);
  }
  const failure = result.failure;
  if (failure === undefined) return;
  let shrunk: readonly RecordedStep[] = result.steps;
  let shrinkNote = 'not shrunk';
  if (failure.index > 0 && shrink.replays > 0) {
    await page.close().catch(() => undefined);
    const prefix = [...Array(failure.index).keys()];
    const use = testInfo.project.use as Record<string, unknown>;
    const contextOptions = Object.fromEntries(CONTEXT_OPTIONS.flatMap((key) => (use[key] === undefined ? [] : [[key, use[key]]]))) as BrowserContextOptions;
    const outcome = await shrinkWalk({
      prefix,
      failing: failure.index,
      budget: shrink.replays,
      deadline: Date.now() + shrink.ms,
      reproduces: (indices) =>
        reproduces(
          browser,
          contextOptions,
          indices.map((i) => result.steps[i]!),
          failure,
          options,
        ),
    });
    shrunk = outcome.indices.map((i) => result.steps[i]!);
    shrinkNote = `${outcome.confirmed ? 'shrunk' : 'not reproduced in'} ${outcome.replays} replays to ${shrunk.length} of ${failure.index + 1} steps${outcome.finished ? '' : ' (replays or time used up)'}`;
    const shrunkFile: WalkFile = { ...file, steps: shrunk };
    await testInfo.attach('walk-shrunk.json', { body: JSON.stringify(shrunkFile, null, 1), contentType: 'application/json' });
  }
  const lines = shrunk.map((step, i) => `  ${i + 1}. ${step.kind}${step.note ? `: ${step.note}` : ''}`);
  throw new Error(
    `walk ${walk} (seed ${seed}) failed at step ${failure.index + 1} (${failure.kind}), ${failure.signature}:\n${failure.message}\n` +
      `${shrinkNote}:\n${lines.join('\n')}\n` +
      `Replay: LOSAT_WEB_WALK_REPLAY=<walk-shrunk.json or walk-steps.json of this test> npx playwright test random-walk --project=${browserName}`,
  );
}

/** The walk of LOSAT_WEB_WALK_REPLAY, if set. */
export function replayFile(): WalkFile | undefined {
  const path = process.env['LOSAT_WEB_WALK_REPLAY'];
  if (path === undefined || path === '') return undefined;
  return JSON.parse(readFileSync(path, 'utf8')) as WalkFile;
}
