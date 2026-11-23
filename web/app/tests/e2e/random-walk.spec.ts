// The random walk of the search form and the results (line WQA, after gbdraw's random walk): the
// seeded plan of walk/plan.ts, whose walks hold every step kind and every ordered pair of kinds,
// each walk in a fresh page, each step checked as walk/steps.ts says (no page or console error, no
// request to another origin, no script of the inputs run, no page wider than its window, the form
// as the walk left it; each completed run's stored outputs are the bytes that its command line
// writes in the native CLI of the same commit, with LOSAT_WEB_NATIVE). A failing walk is shrunk to
// a small reproduction (walk/shrink.ts), attached as walk-shrunk.json with the whole walk
// (walk-steps.json).
//
// LOSAT_WEB_WALK_LEVEL=full runs the long plan (every ordered pair also as adjacent steps),
// LOSAT_WEB_WALK_SEED another seed, LOSAT_WEB_WALK_REPLAY=<walk file> one recorded walk, and
// LOSAT_WEB_WALK_SHRINK_REPLAYS the replays a failing walk may use to shrink.
import { readFileSync, readdirSync } from 'node:fs';
import { join } from 'node:path';
import { expect } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NATIVE } from './support/native';
import { test } from './support/test';
import { coverageOf, planWalk, PR_SEED, WALK_KINDS, type PlanLevel } from './walk/plan';
import { replayFile, shrinkLimits, walkTest, type WalkFile } from './walk/runner';

const level = (process.env['LOSAT_WEB_WALK_LEVEL'] ?? 'pr') as PlanLevel;
const seed = Number(process.env['LOSAT_WEB_WALK_SEED'] ?? PR_SEED);
const plan = planWalk(seed, { level });
/**
 * A failing walk is shrunk for at most 3 min at the pull-request level (the CI jobs of web.yml end
 * at 25 and 35 min, and their artifacts are uploaded only if the job ends), 15 min at the full one.
 */
const SHRINK = shrinkLimits({ replays: level === 'pr' ? 8 : 16, ms: (level === 'pr' ? 3 : 15) * 60_000 });
/** A walk's time: most steps take well under a second, a run's step some seconds; 4 s a step, and 2 min. */
const STEP_MS = 4_000;
const WALK_MS = (steps: number) => 120_000 + STEP_MS * steps;

test.describe.configure({ retries: 0 });
// An action on a control that never comes fails in seconds, as a step of the walk, not at the
// walk's own time limit.
test.use({ actionTimeout: 15_000, navigationTimeout: 30_000 });

test('the plan holds every kind and every ordered pair of kinds', () => {
  expect(coverageOf(plan, level).missing).toEqual([]);
  const drawn = new Set(plan.walks.flatMap(({ steps }) => steps.map(({ kind }) => kind)));
  expect(WALK_KINDS.filter((kind) => !drawn.has(kind))).toEqual([]);
});

const replay = replayFile();
if (replay !== undefined) {
  test(`random walk replay: walk ${replay.walk} of seed ${replay.seed}, ${replay.steps.length} steps`, async ({ page, browser, browserName }, testInfo) => {
    test.setTimeout(WALK_MS(replay.steps.length));
    await walkTest({ page, browser, browserName }, testInfo, { seed: replay.seed, walk: replay.walk, steps: replay.steps, shrink: { replays: 0, ms: 0 } });
  });
} else {
  for (const [index, walk] of plan.walks.entries()) {
    test(`random walk ${index + 1} of ${plan.walks.length} (seed ${seed}, ${walk.steps.length} steps)`, async ({ page, browser, browserName }, testInfo) => {
      test.setTimeout(WALK_MS(walk.steps.length) + SHRINK.ms + 60_000);
      await walkTest({ page, browser, browserName }, testInfo, { seed, walk: index, steps: walk.steps, shrink: SHRINK });
    });
  }
}

// Recorded walks that the drawn ones rarely reach (walk/recorded/*.json, replays of walk files): an input that is not the user's
// file as given (a file with records left out, line WQA's D22) must be named and offered as a file of its own, and the CLI must
// give the stored bytes on the shown command with the offered file. Run in every browser, whatever the seed, but only with the engine
// and the native CLI: without them (the app jobs) no such case is checked.
if (replay === undefined && BUILD_HAS_ENGINE && NATIVE !== undefined) {
  const directory = join(process.cwd(), 'tests/e2e/walk/recorded');
  for (const name of readdirSync(directory).filter((file) => file.endsWith('.json')).sort()) {
    const recorded = JSON.parse(readFileSync(join(directory, name), 'utf8')) as WalkFile;
    test(`recorded walk ${name.replace(/\.json$/, '')} (${recorded.steps.length} steps)`, async ({ page, browser, browserName }, testInfo) => {
      // a recorded walk has up to three runs of the engine, some tens of seconds each on a loaded machine
      test.setTimeout(WALK_MS(recorded.steps.length) + 180_000);
      await walkTest({ page, browser, browserName }, testInfo, { seed: recorded.seed, walk: recorded.walk, steps: recorded.steps, shrink: { replays: 0, ms: 0 }, ...(recorded.minOfferedLeftOut === undefined ? {} : { minOfferedLeftOut: recorded.minOfferedLeftOut }) });
    });
  }
}
