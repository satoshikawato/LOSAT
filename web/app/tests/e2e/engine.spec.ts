// The real engine in the browsers (S09, plan §6.2 V-BR and the S10 contracts): the S10
// contract suites against the adapter's serial reactor and the Engine worker; V-BR,
// searches through the application compared with the native CLI of the same commit and
// with NCBI's frozen bytes (TD-5) at 1, 2 and 4 threads; the retained subject (R1); cancel
// and recovery; renewal of the runtime; and the serial fallback. Every project (Chromium,
// Firefox, WebKit) runs them, on a harness server (support/harness-server.ts). The
// storage-full cases need the Chrome DevTools Protocol and run in contracts.spec.ts.
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { expect, test } from '@playwright/test';
import type { OutputFormat } from '../../src/domain/output-format';
import type { RetentionResult, SearchCase, SearchResult } from './harness/engine';
import { buildHarness, type HarnessFiles } from './support/browser';
import { ENGINE, NO_ENGINE_REASON } from './support/engine';
import { openHarness, REPOSITORY, startHarnessServer, type HarnessServer } from './support/harness-server';
import { manifestCase, NATIVE, nativeExpectation, NO_NATIVE_REASON, searchOf, vbrCases, type VbrCase } from './support/native';

test.skip(ENGINE === undefined, NO_ENGINE_REASON);
test.setTimeout(900_000);

/** Where the tests write their records for the gate record (docs/evidence/losat_web_w1). */
const EVIDENCE = process.env['LOSAT_WEB_EVIDENCE'] || undefined;
const FORMATS: readonly OutputFormat[] = [0, 6, 7];

let harness: HarnessFiles;
let isolated: HarnessServer;
let notIsolated: HarnessServer;

/**
 * A TBLASTN search of several query batches (62,500 residues of NZ_CP006932's proteins
 * against 437,500 nt of its genome), as files in a directory of their own: the query is
 * searched in batches of 20,000 residues, each with its own thread pool. The directory is
 * made by the worker that runs these tests and removed after them.
 */
let batchesDir: string | undefined;
function writeBatchInputs(): Map<string, Uint8Array> {
  batchesDir = mkdtempSync(join(tmpdir(), 'losat-web-batches-'));
  const read = (path: string) => readFileSync(join(REPOSITORY, path), 'utf8');
  let query = '';
  let residues = 0;
  for (const record of read('LOSAT/tests/fasta/NZ_CP006932.faa').split(/^>/m).slice(1)) {
    if (residues >= 62_500) break;
    query += `>${record}`;
    residues += record.split('\n').slice(1).join('').replace(/\s/g, '').length;
  }
  const [header, ...lines] = read('LOSAT/tests/fasta/NZ_CP006932.fasta').split('\n');
  const subject = `${header}\n${lines.join('').slice(0, 437_500).replace(/(.{70})/g, '$1\n')}\n`;
  const files = new Map([
    ['batches_query.faa', new TextEncoder().encode(query)],
    ['batches_subject.fasta', new TextEncoder().encode(subject)],
  ]);
  for (const [name, bytes] of files) writeFileSync(join(batchesDir, name), bytes);
  return files;
}

test.beforeAll(async () => {
  harness = await buildHarness();
  isolated = await startHarnessServer(harness, { extra: writeBatchInputs() });
  notIsolated = await startHarnessServer(harness, { isolated: false });
});

test.afterAll(async () => {
  await isolated?.close();
  await notIsolated?.close();
  if (batchesDir !== undefined) rmSync(batchesDir, { recursive: true, force: true });
});

const failures = (results: ReadonlyArray<{ ok: boolean }>) => results.filter((result) => !result.ok);

function record(name: string, browser: string, value: unknown): void {
  if (EVIDENCE !== undefined) writeFileSync(join(EVIDENCE, `${name}-${browser}.json`), `${JSON.stringify(value, null, 2)}\n`);
}

interface Outputs {
  readonly id: string;
  readonly argv?: readonly string[];
  readonly sha256: Partial<Record<OutputFormat, string>>;
  readonly diagnosticsSha256?: string;
}

/** Compares one search's stored outputs with the native CLI and NCBI's frozen bytes. */
function expectOutputs(vbr: VbrCase, result: Outputs) {
  const native = nativeExpectation(result.argv!, vbr.cwd);
  const frozen: Partial<Record<OutputFormat, boolean>> = {};
  for (const format of FORMATS) {
    expect.soft(result.sha256[format], `${result.id} outfmt ${format} = native CLI`).toBe(native.sha256[format]);
    const fixed = vbr.frozen[format];
    if (fixed === undefined) continue;
    frozen[format] = result.sha256[format] === fixed.sha256;
    expect.soft(result.sha256[format], `${result.id} outfmt ${format} = NCBI frozen`).toBe(fixed.sha256);
  }
  expect.soft(result.diagnosticsSha256, `${result.id} diagnostics = native CLI standard error`).toBe(native.stderrSha256);
  return { native: native.sha256, frozen };
}

test('record scanner contract: the serial reactor', async ({ context }) => {
  const page = await openHarness(context, isolated);
  const results = await page.evaluate(() => window.losatHarness!.engine.recordScanner());
  expect(failures(results)).toEqual([]);
  expect(results.length).toBeGreaterThan(20);
});

for (const threads of [1, 4]) {
  test(`engine input contract: the real EngineGateway at ${threads} thread(s)`, async ({ context }) => {
    const page = await openHarness(context, isolated);
    const results = await page.evaluate(
      (n) =>
        window.losatHarness!.engine.engineInput({
          argv: ['blastn', '-query', 'query.fa', '-subject', 'subject.fa', '-task', 'blastn'],
          query: '>q1 first\nACGTTGCAAGGCTTAACCGGTTAAACCCGGGTTTAAACGTACGTAGC\n>q2\nTAGCTAGGATCCGATCGATTAGCATGCATGCAAATTTGGGCCC\n',
          subject:
            '>s1\nACGTTGCAAGGCTTAACCGGTTAAACCCGGGTTTAAACGTACGTAGCTAGCTAGGATCCGATCGATTAGCATGCATGCAAATTTGGGCCC\n' +
            '>s2 second\nTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA\n',
          queryRecords: [
            { id: 'q1', length: 47 },
            { id: 'q2', length: 43 },
          ],
          subjectRecords: [
            { id: 's1', length: 90 },
            { id: 's2', length: 32 },
          ],
          threads: n,
        }),
      threads,
    );
    expect(failures(results)).toEqual([]);
    expect(results).toHaveLength(4);
  });
}

test('run output contract: the writer in the Engine worker, the real Data worker (without storage full)', async ({ context }) => {
  const page = await openHarness(context, isolated);
  await page.exposeFunction('setStorageQuota', async () => undefined);
  const { results } = await page.evaluate(() => window.losatHarness!.engine.runOutputEngine());
  const withoutQuota = results.filter((result) => !result.name.startsWith('storage full'));
  expect(failures(withoutQuota)).toEqual([]);
  expect(withoutQuota).toHaveLength(10);
});

test('V-BR: searches through the application give the expected bytes at 1, 2 and 4 threads', async ({ context, browserName }) => {
  test.skip(NATIVE === undefined, NO_NATIVE_REASON);
  const page = await openHarness(context, isolated);
  const cases = vbrCases();
  const list: SearchCase[] = cases.flatMap((vbr) => [1, 2, 4].map((threads) => searchOf(vbr, threads)));
  const results: SearchResult[] = await page.evaluate((searches) => window.losatHarness!.engine.searches(searches), list);
  const rows = [];
  for (const [i, result] of results.entries()) {
    const vbr = cases[Math.floor(i / 3)]!;
    const threads = list[i]!.threads as number;
    expect.soft(result.status, `${result.id}: ${result.record.error ?? ''}`).toBe('completed');
    if (result.status !== 'completed') continue;
    expect.soft(result.record.runtimePath, `${result.id} ${result.record.fallbackReason ?? ''}`).toBe(threads === 1 ? 'serial' : 'threaded');
    expect.soft(result.record.threads, result.id).toBe(threads);
    expect.soft(result.hitsMatchRows, `${result.id} HSP records`).toBe(true);
    const compared = expectOutputs(vbr, result);
    rows.push({
      id: result.id,
      argv: result.argv,
      runtimePath: result.record.runtimePath,
      threads: result.record.threads,
      engineBuild: result.record.engineBuild,
      sha256: result.sha256,
      native: compared.native,
      frozenMatch: compared.frozen,
      hits: result.hits,
      elapsedMs: Math.round(result.elapsedMs),
      memory: result.record.memory,
    });
  }
  record('v-br', browserName, rows);
});

test('R1: searches of one subject with other queries and options reuse it and give the expected bytes', async ({ context, browserName }) => {
  test.skip(NATIVE === undefined, NO_NATIVE_REASON);
  const page = await openHarness(context, isolated);
  // The same subject file and dataset revision; the query or the options change.
  const steps: Array<readonly [VbrCase, number, boolean]> = [
    [manifestCase('multi.blastn', ['multi.blastn']), 1, false],
    [manifestCase('multi.megablast', ['multi.megablast']), 1, true],
    [manifestCase('iupac.beyond_table.blastn', ['scoring.beyond_table_iupac.blastn']), 1, true],
    [manifestCase('rna.megablast', ['input.rna.megablast']), 1, true],
    [manifestCase('multi.mts3.blastn', ['multi.mts3.blastn']), 1, true],
    // The threaded module is another instance: it registers the subject once.
    [manifestCase('multi.ws7_e1000.blastn', ['multi.ws7_e1000.blastn']), 4, false],
    [manifestCase('multi.besthit.blastn', ['multi.besthit.blastn']), 4, true],
    // Another program reads the subject on its own.
    [manifestCase('tblastx.multi', ['tblastx.multi.0', 'tblastx.multi.7']), 1, false],
    [manifestCase('tblastx.batch', ['tblastx.batch.0', 'tblastx.batch.7']), 1, true],
    [manifestCase('tblastx.mts1', ['tblastx.mts1.0', 'tblastx.mts1.7']), 4, false],
    [manifestCase('tblastx.thrwin', ['tblastx.thrwin.0']), 4, true],
  ];
  const list = steps.map(([vbr, threads]) => searchOf(vbr, threads));
  const results: RetentionResult[] = await page.evaluate((searches) => window.losatHarness!.engine.retention(searches), list);
  const rows = [];
  for (const [i, result] of results.entries()) {
    const [vbr, threads, retained] = steps[i]!;
    expect.soft(result.error, result.id).toBeUndefined();
    if (result.error !== undefined) continue;
    expect.soft(result.runtime?.subjectRetained, `${result.id} retained`).toBe(retained);
    expect.soft(result.runtime?.threads, result.id).toBe(threads);
    const compared = expectOutputs(vbr, result);
    rows.push({
      id: result.id,
      argv: result.argv,
      retained: result.runtime?.subjectRetained,
      sha256: result.sha256,
      native: compared.native,
      frozenMatch: compared.frozen,
      memory: result.runtime?.memory,
    });
  }
  record('r1', browserName, rows);
});

test('a search that builds one thread pool after another runs threaded from a new runtime', async ({ context }) => {
  // The first threaded search of a runtime: each query batch starts its pool while the
  // threads of the previous pool are still returning (the ThreadHost's spare thread workers).
  const page = await openHarness(context, isolated);
  const batches: VbrCase = {
    id: 'tblastn.query_batches',
    program: 'tblastn',
    options: [],
    cwd: batchesDir!,
    query: 'batches_query.faa',
    subject: 'batches_subject.fasta',
    frozen: {},
  };
  const file = (name: string) => ({ url: `/__extra/${name}`, name });
  for (const threads of [2, 4]) {
    const search: SearchCase = { ...searchOf(batches, threads), query: file(batches.query), subject: file(batches.subject) };
    // A new application (a session with its own options), so that this is the first search of its runtime.
    const [result] = await page.evaluate((s) => window.losatHarness!.engine.searches([s], { renewal: {} }), search);
    expect(result!.status, result!.record.error).toBe('completed');
    expect(result!.record.runtimePath).toBe('threaded');
    expect(result!.record.threads).toBe(threads);
    if (NATIVE !== undefined) expectOutputs(batches, result!);
  }
});

test('a cancelled search ends its runtime; the next search runs in a new one with the expected bytes', async ({ context, browserName }) => {
  const page = await openHarness(context, isolated);
  const next = vbrCases()[0]!;
  const long: SearchCase = {
    id: 'long tblastx',
    program: 'tblastx',
    options: [],
    query: { url: '/__files/LOSAT/tests/fasta/MjeNMV.fasta', name: 'MjeNMV.fasta' },
    subject: { url: '/__files/LOSAT/tests/fasta/MelaMJNV.fasta', name: 'MelaMJNV.fasta' },
    threads: 4,
  };
  const rows = [];
  for (const threads of [1, 4]) {
    const result = await page.evaluate(
      ([a, b]) => window.losatHarness!.engine.cancelThenSearch(a, b),
      [{ ...long, threads }, searchOf(next, threads)] as const,
    );
    expect(result.cancelled.status).toBe('cancelled');
    expect(result.next.status, result.next.record.error).toBe('completed');
    expect(result.next.sha256[0]).toBe(next.frozen[0]!.sha256);
    expect(result.next.record.runtimeGeneration).toBe(2);
    expect(result.next.record.runtimePath).toBe(threads === 1 ? 'serial' : 'threaded');
    if (NATIVE !== undefined) expectOutputs(next, result.next);
    rows.push({ threads, cancelMs: Math.round(result.cancelMs), nextPrepareMs: Math.round(result.nextPrepareMs) });
    console.log(`cancel n=${threads}: cancelled in ${Math.round(result.cancelMs)} ms; next search prepared in ${Math.round(result.nextPrepareMs)} ms`);
  }
  record('cancel', browserName, rows);
});

test('a runtime renewed after every search gives the same bytes and registers the subject again', async ({ context }) => {
  const page = await openHarness(context, isolated);
  const vbr = vbrCases()[1]!;
  // The same subject revision, so that a warm runtime would reuse it (R1); 1, 1, 4, 4
  // threads, so that without renewal the second search of each module reuses it.
  const list: SearchCase[] = [1, 1, 4, 4].map((threads, i) => ({ ...searchOf(vbr, threads), id: `renewed ${i} n=${threads}` }));
  const warm: RetentionResult[] = await page.evaluate((cases) => window.losatHarness!.engine.retention(cases), list);
  expect(warm.map((result) => result.runtime?.subjectRetained)).toEqual([false, true, false, true]);
  expect(warm.map((result) => result.runtime?.runtimeGeneration)).toEqual([1, 1, 1, 1]);
  const renewed: RetentionResult[] = await page.evaluate(
    (cases) => window.losatHarness!.engine.retention(cases, { renewal: { maxRuns: 1 } }),
    list,
  );
  expect(renewed.map((result) => result.error)).toEqual([undefined, undefined, undefined, undefined]);
  expect(renewed.map((result) => result.runtime?.runtimeGeneration)).toEqual([1, 2, 3, 4]);
  expect(renewed.map((result) => result.runtime?.memory?.instanceRuns)).toEqual([1, 1, 1, 1]);
  expect(renewed.map((result) => result.runtime?.subjectRetained)).toEqual([false, false, false, false]);
  for (const result of [...warm, ...renewed]) {
    expect(result.sha256[0], result.id).toBe(vbr.frozen[0]!.sha256);
    if (NATIVE !== undefined) expectOutputs(vbr, result);
  }
});

test('without cross-origin isolation, a threaded search runs on the serial module with the reason', async ({ context }) => {
  const page = await openHarness(context, notIsolated);
  expect(await page.evaluate(() => globalThis.crossOriginIsolated)).toBe(false);
  const vbr = vbrCases()[0]!;
  const [result] = await page.evaluate((search) => window.losatHarness!.engine.searches([search]), searchOf(vbr, 4));
  expect(result!.status, result!.record.error).toBe('completed');
  expect(result!.sha256[0]).toBe(vbr.frozen[0]!.sha256);
  expect(result!.record.runtimePath).toBe('serial');
  expect(result!.record.threads).toBe(1);
  expect(result!.record.fallbackReason).toMatch(/not cross-origin isolated/);
});
