// Measurements of the browser runtime (S09; plan §2.3, §4.6, §5.5, DW-8). Not a test: it
// runs only with LOSAT_WEB_MEASURE (a comma list of dw8, memory, auto, cancel, or all) and
// writes its records to LOSAT_WEB_EVIDENCE. docs/evidence/losat_web_w1/README.md reports
// the results and the choices made from them.
//
// - dw8: the share of a warm search that preprocessing determined by the subject alone can
//   take. Every search of the same subject reuses it (R1), so the search phase of a minimal
//   query (a few residues of the real one) against the same subject contains all of that
//   work, plus the scan of the subject and the fixed cost of a run: its ratio to the search
//   phase of the real query is an upper bound of the share. One warm-up and three measured
//   searches each (AGENTS.md). The engine's own timing (LOSAT_TIMING) is not available in
//   the reactors (BLASTP compiles it out on wasm32; TBLASTN has none).
// - memory: the linear memory of the engine instance over repeated searches (G12,
//   docs/wasm_reactor_memory_followup_20260914.md), inputs indexed again for every search
//   as the application does.
// - auto: the search phase at 1, 2 and 4 threads by input size (the Auto rule).
// - cancel: the time to cancel, and to prepare the next search in a new runtime (worker,
//   instance, registration of the subject again) against a warm one (the cooperative
//   cancel threshold, plan §2.3).
import { readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, test } from '@playwright/test';
import type { SearchCase, TimedOptions, TimedRun } from './harness/engine';
import { buildHarness, type HarnessFiles } from './support/browser';
import { ENGINE, NO_ENGINE_REASON } from './support/engine';
import { openHarness, REPOSITORY, startHarnessServer, type HarnessServer } from './support/harness-server';

const MEASURE = new Set((process.env['LOSAT_WEB_MEASURE'] ?? '').split(',').filter((name) => name !== ''));
const wants = (name: string) => MEASURE.has(name) || MEASURE.has('all');
const EVIDENCE = process.env['LOSAT_WEB_EVIDENCE'] || undefined;

test.skip(MEASURE.size === 0, 'measurements run only with LOSAT_WEB_MEASURE');
test.skip(ENGINE === undefined, NO_ENGINE_REASON);
test.setTimeout(3_600_000);

// --- inputs ---------------------------------------------------------------------------------

interface Fasta {
  readonly header: string;
  readonly sequence: string;
}

function readFasta(path: string): Fasta[] {
  const records: Fasta[] = [];
  for (const block of readFileSync(join(REPOSITORY, path), 'utf8').split(/^>/m).slice(1)) {
    const [header = '', ...lines] = block.split('\n');
    records.push({ header, sequence: lines.join('').replace(/\s+/g, '') });
  }
  return records;
}

function writeFasta(records: readonly Fasta[]): Uint8Array {
  const text = records.map(({ header, sequence }) => `>${header}\n${sequence.replace(/(.{70})/g, '$1\n').replace(/\n$/, '')}\n`).join('');
  return new TextEncoder().encode(text);
}

/** The first `residues` residues of the first record. */
const firstResidues = (path: string, residues: number) => {
  const [first] = readFasta(path);
  return writeFasta([{ header: first!.header, sequence: first!.sequence.slice(0, residues) }]);
};

/** Whole records from the start, until their residues reach `residues`. */
const firstRecords = (path: string, residues: number) => {
  const kept: Fasta[] = [];
  let total = 0;
  for (const record of readFasta(path)) {
    if (total >= residues) break;
    kept.push(record);
    total += record.sequence.length;
  }
  return writeFasta(kept);
};

const extra = new Map<string, Uint8Array>();
const file = (path: string) => ({ url: `/__files/${path}`, name: path.split('/').pop()! });
const generated = (name: string, bytes: Uint8Array) => {
  extra.set(name, bytes);
  return { url: `/__extra/${name}`, name };
};

interface MeasuredCase {
  readonly id: string;
  readonly program: SearchCase['program'];
  readonly options: readonly string[];
  readonly query: SearchCase['query'];
  readonly subject: SearchCase['subject'];
  /** A minimal query: the first residues of the first query record. */
  readonly qmin: SearchCase['query'];
}

/** One representative fixture per program: the V-PERF fixtures (docs/evidence/losat_web_e1a/measure_perf.py). */
function perfCases(): MeasuredCase[] {
  return [
    {
      id: 'blastp SicyWSV/PajaWSV',
      program: 'blastp',
      options: ['-max_hsps', '1'],
      query: file('LOSAT/tests/fasta/SicyWSV.faa'),
      subject: file('LOSAT/tests/fasta/PajaWSV.faa'),
      qmin: generated('qmin-blastp.faa', firstResidues('LOSAT/tests/fasta/SicyWSV.faa', 30)),
    },
    {
      id: 'tblastn AvCLPV protein/AvCLPV',
      program: 'tblastn',
      options: ['-task', 'tblastn'],
      query: file('docs/evidence/tlosan_stage_g/benchmark/first_AvCLPV_protein.faa'),
      subject: file('LOSAT/tests/fasta/AvCLPV.fasta'),
      qmin: generated('qmin-tblastn.faa', firstResidues('docs/evidence/tlosan_stage_g/benchmark/first_AvCLPV_protein.faa', 30)),
    },
    {
      id: 'blastn megablast AP027152/LC738884',
      program: 'blastn',
      options: ['-task', 'megablast'],
      query: file('LOSAT/tests/fasta/AP027152.fasta'),
      subject: file('LOSAT/tests/fasta/LC738884.fasta'),
      qmin: generated('qmin-blastn.fasta', firstResidues('LOSAT/tests/fasta/AP027152.fasta', 60)),
    },
    {
      id: 'tblastx LC738884/LC741431',
      program: 'tblastx',
      options: [],
      query: file('LOSAT/tests/fasta/LC738884.fasta'),
      subject: file('LOSAT/tests/fasta/LC741431.fasta'),
      qmin: generated('qmin-tblastx.fasta', firstResidues('LOSAT/tests/fasta/LC738884.fasta', 90)),
    },
  ];
}

/** A large BLASTN pair: two E. coli genomes (native: 1.3 s, 244 MB). */
const LARGE: Omit<MeasuredCase, 'qmin'> = {
  id: 'blastn megablast Sakai/MG1655',
  program: 'blastn',
  options: ['-task', 'megablast'],
  query: file('LOSAT/tests/fasta/Sakai.fna'),
  subject: file('LOSAT/tests/fasta/MG1655.fna'),
};

const search = (c: Omit<MeasuredCase, 'qmin'>, threads: number, id = c.id, query = c.query): SearchCase => ({
  id,
  program: c.program,
  options: c.options,
  query,
  subject: c.subject,
  threads,
});

/** Inputs of the Auto measurement: query and subject of about `total / 2` residues each. */
function autoCases(total: number): Array<Omit<MeasuredCase, 'qmin'>> {
  const half = Math.round(total / 2);
  const pair = (program: SearchCase['program'], options: string[], q: string, s: string, cut: typeof firstResidues) => ({
    id: `${program} ${total}`,
    program,
    options,
    query: generated(`auto-${program}-${total}-q.fa`, cut(q, half)),
    subject: generated(`auto-${program}-${total}-s.fa`, cut(s, half)),
  });
  return [
    pair('blastn', ['-task', 'megablast'], 'LOSAT/tests/fasta/Sakai.fna', 'LOSAT/tests/fasta/MG1655.fna', firstResidues),
    pair('tblastx', [], 'LOSAT/tests/fasta/LC738884.fasta', 'LOSAT/tests/fasta/LC741431.fasta', firstResidues),
    pair('blastp', [], 'LOSAT/tests/fasta/NZ_CP006932.faa', 'LOSAT/tests/fasta/NZ_CP006932.faa', firstRecords),
    {
      id: `tblastn ${total}`,
      program: 'tblastn',
      options: [],
      query: generated(`auto-tblastn-${total}-q.faa`, firstRecords('LOSAT/tests/fasta/NZ_CP006932.faa', Math.round(total / 8))),
      subject: generated(`auto-tblastn-${total}-s.fa`, firstResidues('LOSAT/tests/fasta/NZ_CP006932.fasta', Math.round((total * 7) / 8))),
    },
  ];
}

// 500,000 residues already takes the whole NZ_CP006932 proteome for BLASTP.
const AUTO_TOTALS = [20_000, 60_000, 200_000, 500_000];

// --- the runs --------------------------------------------------------------------------------

let harness: HarnessFiles;
let server: HarnessServer;

test.beforeAll(async () => {
  perfCases();
  for (const total of AUTO_TOTALS) autoCases(total);
  harness = await buildHarness();
  server = await startHarnessServer(harness, { extra });
});

test.afterAll(async () => {
  await server?.close();
});

function record(name: string, browser: string, value: unknown): void {
  if (EVIDENCE !== undefined) writeFileSync(join(EVIDENCE, `measure-${name}-${browser}.json`), `${JSON.stringify(value, null, 2)}\n`);
}

const median = (values: readonly number[]) => {
  const sorted = [...values].sort((a, b) => a - b);
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2 === 1 ? sorted[middle]! : (sorted[middle - 1]! + sorted[middle]!) / 2;
};

async function timed(page: import('@playwright/test').Page, list: SearchCase[], options: TimedOptions = {}): Promise<TimedRun[]> {
  const results = await page.evaluate(([searches, opts]) => window.losatHarness!.engine.timed(searches, opts), [list, options] as const);
  for (const result of results) expect.soft(result.error, `${result.id}`).toBeUndefined();
  return results;
}

test('dw8: the share of a warm search that subject-only preprocessing can take', async ({ context, browserName }) => {
  test.skip(!wants('dw8'));
  const rows = [];
  for (const threads of [1, 4]) {
    for (const c of perfCases()) {
      const page = await openHarness(context, server);
      // One warm-up and three measured searches of each query; the subject is retained.
      const list = [
        ...[0, 1, 2, 3].map((i) => search(c, threads, `${c.id} query ${i}`)),
        ...[0, 1, 2, 3].map((i) => search(c, threads, `${c.id} qmin ${i}`, c.qmin)),
      ];
      const results = await timed(page, list);
      const full = results.slice(1, 4).map((r) => r.runMs!);
      const minimal = results.slice(5, 8).map((r) => r.runMs!);
      rows.push({
        id: c.id,
        threads,
        runMs: full.map(Math.round),
        qminRunMs: minimal.map(Math.round),
        medianRunMs: Math.round(median(full)),
        medianQminRunMs: Math.round(median(minimal)),
        upperBoundShare: Number((median(minimal) / median(full)).toFixed(4)),
        retained: results.map((r) => r.runtime?.subjectRetained),
      });
      console.log(JSON.stringify(rows[rows.length - 1]));
      await page.close();
    }
  }
  record('dw8', browserName, rows);
});

test('memory: the linear memory of one instance over repeated searches', async ({ context, browserName }) => {
  test.skip(!wants('memory'));
  const rows = [];
  for (const threads of [1, 4]) {
    for (const c of [...perfCases(), LARGE]) {
      const page = await openHarness(context, server);
      const results = await timed(page, Array.from({ length: 20 }, (_, i) => search(c, threads, `${c.id} ${i}`)), { reindex: true });
      rows.push({
        id: c.id,
        threads,
        linearBytesAfter: results.map((r) => r.runtime?.memory?.linearBytesAfter),
        runMs: results.map((r) => Math.round(r.runMs ?? NaN)),
        outputBytes: results[0]?.outputBytes,
        generations: [...new Set(results.map((r) => r.runtime?.runtimeGeneration))],
      });
      console.log(JSON.stringify({ id: c.id, threads, last: rows[rows.length - 1]!.linearBytesAfter.slice(-3) }));
      await page.close();
    }
    // Every program in turn in one instance.
    const page = await openHarness(context, server);
    const mixed = Array.from({ length: 6 }, (_, round) => [...perfCases(), LARGE].map((c) => search(c, threads, `${c.id} round ${round}`))).flat();
    const results = await timed(page, mixed, { reindex: true });
    rows.push({
      id: 'mixed',
      threads,
      order: mixed.map((s) => s.id),
      linearBytesAfter: results.map((r) => r.runtime?.memory?.linearBytesAfter),
      runMs: results.map((r) => Math.round(r.runMs ?? NaN)),
    });
    await page.close();
  }
  record('memory', browserName, rows);
});

test('auto: the search phase at 1, 2 and 4 threads by input size', async ({ context, browserName }) => {
  test.skip(!wants('auto'));
  const rows = [];
  for (const total of AUTO_TOTALS) {
    for (const c of autoCases(total)) {
      const page = await openHarness(context, server);
      const list = [1, 2, 4].flatMap((threads) => [0, 1, 2, 3].map((i) => search(c, threads, `${c.id} n=${threads} ${i}`)));
      const results = await timed(page, list);
      const byThreads = Object.fromEntries(
        [1, 2, 4].map((threads, k) => [threads, results.slice(k * 4 + 1, k * 4 + 4).map((r) => r.runMs!)]),
      ) as Record<number, number[]>;
      rows.push({
        id: c.id,
        program: c.program,
        totalResidues: total,
        bytes: (extra.get(c.query.name)?.length ?? 0) + (extra.get(c.subject.name)?.length ?? 0),
        runMs: byThreads,
        medianRunMs: Object.fromEntries(Object.entries(byThreads).map(([n, values]) => [n, Math.round(median(values))])),
      });
      console.log(JSON.stringify(rows[rows.length - 1]!.medianRunMs), c.id);
      await page.close();
    }
  }
  record('auto', browserName, rows);
});

test('cancel: the time to cancel and to prepare the next search in a new runtime', async ({ context, browserName }) => {
  test.skip(!wants('cancel'));
  const rows = [];
  const long = perfCases()[3]!;
  for (const threads of [1, 4]) {
    for (const next of [perfCases()[2]!, LARGE]) {
      const page = await openHarness(context, server);
      // A warm runtime holds `next`'s subject; the long search (another subject) is cancelled;
      // the next search runs in a new runtime and registers its subject again.
      const list = [
        search(next, threads, `${next.id} first`),
        search(next, threads, `${next.id} warm`),
        search(long, threads, `${long.id} cancelled`),
        search(next, threads, `${next.id} after cancel`),
        search(next, threads, `${next.id} warm again`),
      ];
      const results = await page.evaluate(
        ([searches, opts]) => window.losatHarness!.engine.timed(searches, opts),
        [list, { cancel: { index: 2, afterMs: 500 } }] as const,
      );
      rows.push({
        next: next.id,
        threads,
        runs: results.map((r) => ({
          id: r.id,
          error: r.error,
          prepareMs: Math.round(r.prepareMs ?? NaN),
          runMs: Math.round(r.runMs ?? NaN),
          cancelMs: r.cancelMs === undefined ? undefined : Math.round(r.cancelMs),
          generation: r.runtime?.runtimeGeneration,
          retained: r.runtime?.subjectRetained,
        })),
      });
      console.log(JSON.stringify(rows[rows.length - 1]));
      await page.close();
    }
  }
  record('cancel', browserName, rows);
});
