// The self-test's expected outputs (src/examples/self-test.json) are the engine's: the serial
// reactor runs each example of a program that LOSAT Web searches in Node, as the Engine worker
// does (register both inputs, the form's argv with -num_threads 1), and the SHA-256 of streams
// 0, 6 and 7 are the table's. It needs LOSAT_WEB_REACTORS (CI's browser-engine job and the gate
// run it): after an engine change that alters these outputs it fails until the table is made
// again with `web/app/tools/examples.py expected` from the native CLI of the new engine commit.
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { findReactors } from '../../build/reactors';
import { parseSelfTestTable } from '../../src/domain/examples';
import { OUTPUT_FORMATS } from '../../src/domain/output-format';
import { programById } from '../../src/domain/programs';
import { ROLE_QUERY, ROLE_SUBJECT } from '../../src/infra/reactor/abi';
import { instantiateSerial } from '../../src/infra/reactor/instance';
import TABLE from '../../src/examples/self-test.json';

const EXAMPLES = resolve(import.meta.dirname, '../../src/examples');
const reactors = findReactors();
const table = parseSelfTestTable(TABLE);
const sha256 = (bytes: Uint8Array) => createHash('sha256').update(bytes).digest('hex');

describe.skipIf(reactors === undefined)("the self-test's expected outputs are the serial reactor's", () => {
  const module = reactors === undefined ? undefined : new WebAssembly.Module(reactors.serial.bytes as BufferSource);
  for (const testCase of table.cases.filter((c) => programById(c.program).unavailable === undefined)) {
    it(`${testCase.program}: outfmt 0, 6 and 7`, async () => {
      const { abi } = await instantiateSerial(module!);
      const read = (role: 'query' | 'subject') => readFileSync(resolve(EXAMPLES, testCase.inputs[role].file));
      const query = abi.register(testCase.program, ROLE_QUERY, read('query'));
      const subject = abi.register(testCase.program, ROLE_SUBJECT, read('subject'));
      const streams = new Map<number, Uint8Array[]>();
      try {
        abi.run([...testCase.argv, '-num_threads', '1'], query.handle, subject.handle, (stream, bytes) =>
          streams.set(stream, [...(streams.get(stream) ?? []), bytes.slice()]),
        );
      } finally {
        abi.release(query.handle);
        abi.release(subject.handle);
      }
      for (const format of OUTPUT_FORMATS) {
        const bytes = Buffer.concat(streams.get(format) ?? []);
        const expected = testCase.outputs[format];
        expect({ format, bytes: bytes.length, sha256: sha256(bytes) }, `${testCase.program} outfmt ${format}`).toEqual({
          format,
          bytes: expected.bytes,
          sha256: expected.sha256,
        });
      }
      // Each example finds hits: the HSP records of stream 1.
      expect(Buffer.concat(streams.get(1) ?? []).toString('utf8').trim().split('\n').length).toBeGreaterThan(1);
    });
  }
});
