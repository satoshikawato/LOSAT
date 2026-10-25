// The shared memory of the threaded engine in the browsers (plan TD-7): each of Chromium,
// Firefox and WebKit provides the maximum that the module declares (65,535 pages), and one
// search per program gives the same bytes on a memory that starts above 2 GiB as on the
// ordinary one, at 4 threads (docs/web/abi_v2.md §6: pointers are unsigned), and the bytes
// of the native CLI of the same commit (and NCBI's frozen bytes where the manifest has
// them). Reactors built before the 1 GiB limit was lifted (plan DW-28) cannot hold a heap
// above 2 GiB; the test skips them.
import { writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, test } from '@playwright/test';
import { OUTPUT_FORMATS } from '../../src/domain/output-format';
import { inspectReactor } from '../../src/infra/reactor/artifact';
import { buildHarness, type HarnessFiles } from './support/browser';
import { ENGINE, NO_ENGINE_REASON } from './support/engine';
import { openHarness, startHarnessServer, type HarnessServer } from './support/harness-server';
import { manifestCase, NATIVE, nativeExpectation, NO_NATIVE_REASON, searchOf, type VbrCase } from './support/native';

const TWO_GIB = 2 ** 31;
const DECLARED = ENGINE === undefined ? undefined : inspectReactor(ENGINE.threads.bytes).memory.maximum;
const EVIDENCE = process.env['LOSAT_WEB_EVIDENCE'] || undefined;

/** One search per program of LOSAT Web: V-BR cells of LOSAT/tests/outfmt0_manifest.tsv. */
const CASES: readonly VbrCase[] = [
  manifestCase('compact.blastn', ['compact.blastn']),
  manifestCase('width.blastp', ['width.blastp']),
  manifestCase('width.tblastn', ['width.tblastn']),
  manifestCase('tblastx.multi', ['tblastx.multi.0', 'tblastx.multi.7']),
];
const argvOf = (c: VbrCase) => [c.program, '-query', c.query, '-subject', c.subject, ...c.options];

test.skip(ENGINE === undefined, NO_ENGINE_REASON);
test.skip(NATIVE === undefined, NO_NATIVE_REASON);
test.skip(DECLARED !== undefined && DECLARED * 65_536 <= TWO_GIB, 'these reactors declare a shared memory of at most 2 GiB');
test.setTimeout(300_000);

let harness: HarnessFiles;
let server: HarnessServer;

test.beforeAll(async () => {
  harness = await buildHarness();
  server = await startHarnessServer(harness);
});

test.afterAll(async () => {
  await server?.close();
});

test('the browser provides the declared maximum, and searches above 2 GiB give the same bytes', async ({ context, browserName }) => {
  const page = await openHarness(context, server);
  const result = await page.evaluate(
    (list) => window.losatHarness!.engine.highMemory(list, 4),
    CASES.map((c) => {
      const search = searchOf(c, 4);
      return { program: c.program, argv: argvOf(c), query: search.query, subject: search.subject };
    }),
  );
  if (EVIDENCE !== undefined) writeFileSync(join(EVIDENCE, `engine-memory-${browserName}.json`), `${JSON.stringify(result, null, 2)}\n`);

  // The whole declared maximum, in both modes.
  expect(result.declaredPages).toBe(DECLARED);
  expect([result.normal.maximumPages, result.high.maximumPages]).toEqual([DECLARED, DECLARED]);
  // The high mode really ran above 2 GiB.
  expect(result.normal.initialBytes).toBeLessThan(TWO_GIB);
  expect(result.high.initialBytes).toBeGreaterThan(TWO_GIB);
  expect(result.high.firstHighPointer).toBeGreaterThanOrEqual(TWO_GIB);
  expect(result.high.linearBytes).toBeGreaterThanOrEqual(result.high.initialBytes);

  expect(result.high.outputs.map((o) => o.program)).toEqual(CASES.map((c) => c.program));
  for (const [index, testCase] of CASES.entries()) {
    const normal = result.normal.outputs[index]!;
    const high = result.high.outputs[index]!;
    // Every stream (the formats, the HSP records and the diagnostics) is the same.
    expect(high, `${testCase.id}: high memory`).toEqual(normal);
    const native = nativeExpectation(argvOf(testCase), testCase.cwd);
    for (const format of OUTPUT_FORMATS) {
      expect({ format, sha256: high.sha256[format] }, `${testCase.id} outfmt ${format}: native`).toEqual({ format, sha256: native.sha256[format] });
      const frozen = testCase.frozen[format];
      if (frozen !== undefined) expect(high.sha256[format], `${testCase.id} outfmt ${format}: NCBI`).toBe(frozen.sha256);
    }
  }
});
