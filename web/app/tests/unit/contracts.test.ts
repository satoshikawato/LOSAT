// Runs the contract suites (tests/contract) against the implementations that work in Node.
// The OPFS store and the cross-worker channel run the same suites in a browser
// (tests/e2e/contracts.spec.ts).
import { describe, it } from 'vitest';
import { sha256Hex } from '../../src/infra/browser/platform';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeScanner } from '../../src/infra/fake/fake-fasta';
import { RunOutputWriter } from '../../src/infra/run-output/writer';
import { BLOCK_STORE_CASES } from '../contract/block-store.contract';
import { ENGINE_INPUT_CASES } from '../contract/engine-input.contract';
import { RECORD_SCANNER_CASES } from '../contract/record-scanner.contract';
import { RUN_OUTPUT_CASES, type RemoteWriter } from '../contract/run-output.contract';

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
