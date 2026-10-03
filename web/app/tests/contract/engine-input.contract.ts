// Engine input contract (plan §5.4; ports/engine.ts EngineInput): after `register`, the
// engine compares the records that its parser read with the record table of the request
// and stops the run, before searching and before any output, if they differ. The
// FakeEngine runs it in Vitest; S09 runs it against the real EngineGateway and Engine worker
// with real FASTA inputs.
import type { RecordKey } from '../../src/domain/dataset';
import type { EngineGateway, EngineRunRequest } from '../../src/ports/engine';
import type { RunOutputMessage } from '../../src/ports/run-output';
import { check, rejects, same, type ContractCase } from './contract';

export interface EngineInputEnv {
  readonly engine: EngineGateway;
  /** A valid argv for `query` and `subject`, without `-num_threads`. */
  readonly argv: readonly string[];
  /** FASTA inputs of two or more records each. */
  readonly query: Uint8Array;
  readonly subject: Uint8Array;
  /** What the engine's own parser reads from `query` and `subject`. */
  readonly queryRecords: readonly RecordKey[];
  readonly subjectRecords: readonly RecordKey[];
  /** Threads of the runs; 1 by default. */
  readonly threads?: number;
}

interface Observed {
  readonly port: MessagePort;
  readonly messages: RunOutputMessage[];
  close(): void;
}

/** A port whose messages the case can see, in place of the data layer. */
function observedPort(): Observed {
  const channel = new MessageChannel();
  const messages: RunOutputMessage[] = [];
  channel.port2.onmessage = (event: MessageEvent<RunOutputMessage>) => messages.push(event.data);
  return { port: channel.port1, messages, close: () => channel.port2.close() };
}

async function waitFor(condition: () => boolean, ms: number): Promise<boolean> {
  const deadline = Date.now() + ms;
  while (!condition()) {
    if (Date.now() > deadline) return false;
    await new Promise((resolve) => setTimeout(resolve, 10));
  }
  return true;
}

async function sha256(bytes: Uint8Array): Promise<string> {
  const digest = await crypto.subtle.digest('SHA-256', bytes as BufferSource);
  return Array.from(new Uint8Array(digest), (byte) => byte.toString(16).padStart(2, '0')).join('');
}

async function request(
  env: EngineInputEnv,
  runId: string,
  change: Partial<Record<'query' | 'subject', readonly RecordKey[]>>,
): Promise<EngineRunRequest> {
  return {
    runId,
    argv: env.argv,
    query: { bytes: env.query, sha256: await sha256(env.query), revisionIds: ['query'], records: change.query ?? env.queryRecords },
    subject: {
      bytes: env.subject,
      sha256: await sha256(env.subject),
      revisionIds: ['subject'],
      records: change.subject ?? env.subjectRecords,
    },
    requestedThreads: env.threads ?? 1,
  };
}

async function expectMismatch(
  env: EngineInputEnv,
  runId: string,
  change: Partial<Record<'query' | 'subject', readonly RecordKey[]>>,
  role: 'query' | 'subject',
): Promise<void> {
  const observed = observedPort();
  try {
    const error = await rejects(
      env.engine.run(await request(env, runId, change), observed.port, () => undefined),
      new RegExp(`^InputMismatchError: The ${role} records`),
      'run',
    );
    check(error.message.includes('record table'), 'the message names the record table');
    await new Promise((resolve) => setTimeout(resolve, 300));
    same(observed.messages.length, 0, 'messages sent by a run stopped at register');
  } finally {
    observed.close();
  }
}

export const ENGINE_INPUT_CASES: readonly ContractCase<EngineInputEnv>[] = [
  {
    name: 'register: the run proceeds when the record tables agree',
    async run(env) {
      const observed = observedPort();
      try {
        await env.engine.run(await request(env, 'input-ok', {}), observed.port, () => undefined);
        check(await waitFor(() => observed.messages.some((m) => m.type === 'end'), 2000), 'the output ends');
      } finally {
        observed.close();
      }
    },
  },
  {
    name: 'register: a query record with another ID stops the run before any output',
    async run(env) {
      const [first, ...rest] = env.queryRecords;
      check(first !== undefined, 'the query has a record');
      await expectMismatch(env, 'input-id', { query: [{ ...first, id: `${first.id}-other` }, ...rest] }, 'query');
    },
  },
  {
    name: 'register: a subject record with another length stops the run before any output',
    async run(env) {
      const records = [...env.subjectRecords];
      const last = records.pop();
      check(last !== undefined, 'the subject has a record');
      await expectMismatch(env, 'input-length', { subject: [...records, { ...last, length: last.length + 1 }] }, 'subject');
    },
  },
  {
    name: 'register: a missing record stops the run before any output',
    async run(env) {
      await expectMismatch(env, 'input-count', { query: env.queryRecords.slice(1) }, 'query');
    },
  },
];
