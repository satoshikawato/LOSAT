// Runs the contract suites (tests/contract) against the implementations that work in Node.
// The OPFS store and the cross-worker channel run the same suites in a browser
// (tests/e2e/contracts.spec.ts).
import { describe, expect, it } from 'vitest';
import { RunOutputWriter, sha256Hex } from 'losat';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import { BLOCK_STORE_CASES } from '../contract/block-store.contract';
import { ENGINE_INPUT_CASES } from '../contract/engine-input.contract';
import { RECORD_SCANNER_CASES } from '../contract/record-scanner.contract';
import { RUN_OUTPUT_CASES, type RemoteWriter } from '../contract/run-output.contract';
import { NodeInflater } from './support/node-inflater';

describe('BlockStore contract: memory', () => {
  for (const contractCase of BLOCK_STORE_CASES) {
    it(contractCase.name, async () => {
      const store = new MemoryBlockStore();
      await contractCase.run({
        store,
        exhaust: async () => store.setCapacity(store.usage()),
        restore: async () => store.setCapacity(Number.POSITIVE_INFINITY),
      });
    });
  }
});

describe('RecordScanner contract: FakeScanner', () => {
  for (const contractCase of RECORD_SCANNER_CASES) {
    it(contractCase.name, () => contractCase.run({ scanner: new FakeScanner() }));
  }
});

/** The engine side in the same thread as the data layer. */
function localWriter(port: MessagePort): RemoteWriter {
  const writer = new RunOutputWriter(port);
  return {
    write: async (stream, bytes) => writer.write(stream, bytes),
    end: async () => writer.end(),
    post: async (message) => port.postMessage(message),
  };
}

describe('Run output contract: DataService with the memory store, in one thread', () => {
  for (const contractCase of RUN_OUTPUT_CASES) {
    it(contractCase.name, async () => {
      const store = new MemoryBlockStore();
      let token = 0;
      const data = new DataService({
        store,
        scanner: new FakeScanner(),
        inflater: new NodeInflater(),
        checker: new FakeInputChecker(),
        digest: sha256Hex,
        newToken: () => `token-${++token}`,
        cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
      });
      await contractCase.run({
        data,
        writer: async (port) => localWriter(port),
        usage: async () => store.usage(),
        exhaust: async () => store.setCapacity(store.usage()),
        restore: async () => store.setCapacity(Number.POSITIVE_INFINITY),
      });
    });
  }
});

describe('Engine input contract: FakeEngine', () => {
  const encoder = new TextEncoder();
  for (const contractCase of ENGINE_INPUT_CASES) {
    it(contractCase.name, () =>
      contractCase.run({
        engine: new FakeEngine(),
        argv: ['blastn', '-query', 'query.fa', '-subject', 'subject.fa'],
        query: encoder.encode('>q1 first\nACGTACGT\n>q2\nGGCC\n'),
        subject: encoder.encode('>s1\nACGTACGTAA\n>s2 second\nTTTT\n'),
        queryRecords: [
          { id: 'q1', length: 8 },
          { id: 'q2', length: 4 },
        ],
        subjectRecords: [
          { id: 's1', length: 10 },
          { id: 's2', length: 4 },
        ],
      }),
    );
  }
});

// The fake's argument error names an unknown option as NCBI's reader does: without its first dash
// and without a value joined by `=` (line RG item K4i-2b round 2, T-L2).
describe('FakeEngine argument errors', () => {
  it('names -foo=3 as "foo", as NCBI splits the word at its first =', async () => {
    const engine = new FakeEngine();
    for (const [word, named] of [['-foo=3', 'foo'], ['-foo', 'foo'], ['--foo=3', '-foo']] as const) {
      const result = await engine.validate(['blastp', '-evalue', '1', word]);
      expect(result.ok, word).toBe(false);
      expect(result.ok ? '' : result.message, word).toContain(`Error: Unknown argument: "${named}"\nError:  Unknown argument: "${named}"\n`);
    }
  });
});
