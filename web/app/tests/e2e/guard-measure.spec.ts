// Traps of the threaded engine with and without the shared-memory guard (release v0.2.0
// DECISIONS 52). Not a test: it runs only with LOSAT_WEB_MEASURE containing `guard` (or
// `all`) and writes its records to LOSAT_WEB_EVIDENCE.
//
// Two harness builds in one browser: the threaded module as the build serves it (guarded
// by build/reactors.ts) and the reactor artifact of LOSAT_WEB_REACTORS as it is. Each round
// runs the same BLASTP searches in both, in a new page each, the order alternating by round
// (ABBA). Every search starts a new runtime (renewal after one search), so the shared memory
// grows from the module's initial size in each. Several pairs run at the same time, as
// gbdraw runs its pair workers (gbdraw OV-2369: 4 searches of 4 threads).
import { createHash } from 'node:crypto';
import { readFileSync, writeFileSync } from 'node:fs';
import { isAbsolute, join } from 'node:path';
import { chromium, expect, test, type BrowserContext } from '@playwright/test';
import type { Reactors } from '../../build/reactors';
import type { SearchCase, TimedRun } from './harness/engine';
import { buildHarness } from './support/browser';
import { ENGINE, NO_ENGINE_REASON } from './support/engine';
import { openHarness, REPOSITORY, startHarnessServer, type HarnessServer } from './support/harness-server';

const MEASURE = new Set((process.env['LOSAT_WEB_MEASURE'] ?? '').split(',').filter((name) => name !== ''));
const EVIDENCE = process.env['LOSAT_WEB_EVIDENCE'] || undefined;
const env = (name: string, fallback: string) => process.env[name] || fallback;
const ROUNDS = Number(env('LOSAT_WEB_GUARD_ROUNDS', '10'));
const SEARCHES = Number(env('LOSAT_WEB_GUARD_SEARCHES', '1'));
const THREADS = Number(env('LOSAT_WEB_GUARD_THREADS', '8'));
// The searches of a round, run at the same time in one page, each in its own runtime:
// `query=subject` pairs separated by `;`, each side a `,` list of FASTA files to join
// (absolute, or relative to the repository).
const PAIRS = env(
  'LOSAT_WEB_GUARD_PAIRS',
  'LOSAT/tests/fasta/AP027078.faa,LOSAT/tests/fasta/AP027131.faa=LOSAT/tests/fasta/AP027132.faa,LOSAT/tests/fasta/AP027133.faa',
)
  .split(';')
  .map((pair) => pair.split('=').map((side) => side.split(',')) as [string[], string[]]);
const ARMS = env('LOSAT_WEB_GUARD_ARMS', 'guarded,unguarded').split(',') as Arm[];
// Another Chromium build than Playwright's own (for example the one that gbdraw measured with).
const CHROMIUM = process.env['LOSAT_WEB_GUARD_CHROMIUM'] || undefined;

test.skip(!MEASURE.has('guard') && !MEASURE.has('all'), 'measurements run only with LOSAT_WEB_MEASURE');
test.skip(ENGINE === undefined, NO_ENGINE_REASON);
test.setTimeout(24 * 3_600_000);

type Arm = 'guarded' | 'unguarded';

/** The reactors with the threaded artifact served as it is, without the guard. */
function unguarded(reactors: Reactors): Reactors {
  const bytes = readFileSync(join(reactors.dir, 'losat-web-threads.wasm'));
  const sha256 = createHash('sha256').update(bytes).digest('hex');
  return { ...reactors, threads: { name: `engine/losat-web-threads-${sha256.slice(0, 16)}.wasm`, bytes, sha256, artifactSha256: sha256 } };
}

const concat = (paths: readonly string[]) => Buffer.concat(paths.map((path) => readFileSync(isAbsolute(path) ? path : join(REPOSITORY, path))));
const extra = new Map<string, Uint8Array>(
  PAIRS.flatMap(([query, subject], i) => [
    [`guard_query_${i}.faa`, concat(query)],
    [`guard_subject_${i}.faa`, concat(subject)],
  ]),
);

const search = (id: string, i: number): SearchCase => ({
  id,
  program: 'blastp',
  options: [],
  query: { url: `/__extra/guard_query_${i}.faa`, name: `guard_query_${i}.faa` },
  subject: { url: `/__extra/guard_subject_${i}.faa`, name: `guard_subject_${i}.faa` },
  threads: THREADS,
});

interface Sample {
  readonly round: number;
  readonly arm: Arm;
  readonly runs: readonly TimedRun[];
}

async function round(context: BrowserContext, server: HarnessServer, label: string): Promise<TimedRun[]> {
  const page = await openHarness(context, server);
  try {
    // One timed() per pair, at the same time: each has its own Data worker and runtime.
    const lists = PAIRS.map((_, i) => Array.from({ length: SEARCHES }, (_, j) => search(`${label} p${i} ${j}`, i)));
    return await page.evaluate(
      async ([lists]) =>
        (await Promise.all(lists.map((list) => window.losatHarness!.engine.timed(list, { renewal: { maxRuns: 1 }, reindex: true })))).flat(),
      [lists] as const,
    );
  } finally {
    await page.close();
  }
}

test('guard: traps of the threaded engine with and without the shared-memory guard', async ({ context: own, browserName, browser: ownBrowser }) => {
  const browser = CHROMIUM === undefined ? ownBrowser : await chromium.launch({ executablePath: CHROMIUM });
  const context = CHROMIUM === undefined ? own : await browser.newContext();
  const servers = new Map<Arm, { server: HarnessServer; reactors: Reactors }>();
  for (const arm of ARMS) {
    const reactors = arm === 'guarded' ? ENGINE! : unguarded(ENGINE!);
    servers.set(arm, { server: await startHarnessServer(await buildHarness(reactors), { extra, noStore: PAIRS.length > 1 }), reactors });
  }
  const samples: Sample[] = [];
  try {
    for (let r = 0; r < ROUNDS; r++) {
      const order = r % 2 === 0 ? ARMS : [...ARMS].reverse();
      for (const arm of order) {
        const { server } = servers.get(arm)!;
        const runs = await round(context, server, `${arm} r${r}`);
        samples.push({ round: r, arm, runs });
        const traps = runs.filter((run) => run.error !== undefined).length;
        console.log(`round ${r} ${arm}: ${traps}/${runs.length} stopped; ${runs.map((run) => run.error ?? `${run.runtime?.path} n=${run.runtime?.threads} ${Math.round(run.runMs ?? 0)} ms`).join('; ')}`);
        if (EVIDENCE !== undefined) {
          writeFileSync(join(EVIDENCE, `guard-${browserName}.json`), `${JSON.stringify({ browser: browser.version(), threads: THREADS, pairs: PAIRS, samples }, null, 2)}\n`);
        }
      }
    }
  } finally {
    for (const { server } of servers.values()) await server.close();
    if (CHROMIUM !== undefined) await browser.close();
  }
  for (const sample of samples) {
    for (const run of sample.runs) if (run.error === undefined) expect(run.runtime?.path).toBe('threaded');
  }
});
