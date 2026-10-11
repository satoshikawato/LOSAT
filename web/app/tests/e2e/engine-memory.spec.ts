// The shared memory of the threaded engine in the browsers (plan TD-7): each of Chromium,
// Firefox and WebKit provides the maximum that the module declares (65,535 pages), and the
// self-test searches give the same bytes on a memory that starts above 2 GiB as on the
// ordinary one, at 4 threads (docs/web/abi_v2.md §6: pointers are unsigned). Reactors
// built before the 1 GiB limit was lifted (plan DW-28) cannot hold a heap above 2 GiB; the
// test skips them.
import { readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, test } from '@playwright/test';
import { parseSelfTestTable } from '../../src/domain/examples';
import { OUTPUT_FORMATS } from '../../src/domain/output-format';
import { programById } from '../../src/domain/programs';
import { inspectReactor } from '../../src/infra/reactor/artifact';
import { buildHarness, type HarnessFiles } from './support/browser';
import { ENGINE, NO_ENGINE_REASON } from './support/engine';
import { openHarness, startHarnessServer, type HarnessServer } from './support/harness-server';

const TWO_GIB = 2 ** 31;
const DECLARED = ENGINE === undefined ? undefined : inspectReactor(ENGINE.threads.bytes).memory.maximum;
const EVIDENCE = process.env['LOSAT_WEB_EVIDENCE'] || undefined;
const TABLE: unknown = JSON.parse(readFileSync(new URL('../../src/examples/self-test.json', import.meta.url), 'utf8'));

test.skip(ENGINE === undefined, NO_ENGINE_REASON);
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
  const table = parseSelfTestTable(TABLE);
  const cases = table.cases.filter((c) => programById(c.program).unavailable === undefined);
  const page = await openHarness(context, server);
  const result = await page.evaluate(
    (list) => window.losatHarness!.engine.highMemory(list, 4),
    cases.map((c) => ({
      program: c.program,
      argv: c.argv,
      query: { url: `/__files/web/app/src/examples/${c.inputs.query.file}`, name: c.inputs.query.file },
      subject: { url: `/__files/web/app/src/examples/${c.inputs.subject.file}`, name: c.inputs.subject.file },
    })),
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

  expect(result.high.outputs.map((o) => o.program)).toEqual(cases.map((c) => c.program));
  for (const [index, testCase] of cases.entries()) {
    const normal = result.normal.outputs[index]!;
    const high = result.high.outputs[index]!;
    // Every stream (the formats, the HSP records and the diagnostics) is the same.
    expect(high, `${testCase.program}: high memory`).toEqual(normal);
    for (const format of OUTPUT_FORMATS) {
      expect({ format, bytes: high.bytes[format], sha256: high.sha256[format] }, `${testCase.program} outfmt ${format}`).toEqual({
        format,
        bytes: testCase.outputs[format].bytes,
        sha256: testCase.outputs[format].sha256,
      });
    }
  }
});
