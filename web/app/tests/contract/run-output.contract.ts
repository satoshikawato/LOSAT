// Run output channel contract (ports/run-output.ts, ports/data.ts RunStore): an engine-side
// writer sends the outputs of a run over the MessagePort of `openRun`, and the data layer
// stages, commits, discards and reads them. Vitest runs it in one thread; Playwright runs
// it with the writer in a separate worker and the real Data worker (OPFS). S09 runs it with
// the writer inside the real Engine worker.
import type { RunStore } from '../../src/ports/data';
import { DIAGNOSTICS_STREAM, HITS_STREAM, type OutputStream } from '../../src/ports/run-output';
import {
  check,
  concatBytes,
  MiB,
  pattern,
  rejects,
  same,
  sameBytes,
  settlesWithin,
  type ContractCase,
} from './contract';

/** The engine side of one run; it may live in another thread. */
export interface RemoteWriter {
  write(stream: OutputStream, bytes: Uint8Array): Promise<void>;
  end(): Promise<void>;
  /** Posts a raw message on the port, as a broken engine might. */
  post(message: unknown): Promise<void>;
}

export interface RunOutputEnv {
  /** The data layer; the cases use distinct run ids, so one data layer may serve them all. */
  readonly data: RunStore;
  /** An engine-side writer on `port`; the environment may transfer the port to a worker. */
  writer(port: MessagePort): Promise<RemoteWriter>;
  /** Bytes that the data layer holds in temporary storage. */
  usage(): Promise<number>;
  /** Makes the storage refuse more data; `restore` undoes it. */
  exhaust(): Promise<void>;
  restore(): Promise<void>;
}

const encoder = new TextEncoder();
const bytes = (value: string) => encoder.encode(value);

function hitLine(index: number): string {
  return `${JSON.stringify({ index, q_idx: 0, s_idx: 0, rank: index, out6: null, out0: null, out0_subject: null })}\n`;
}

async function openWriter(env: RunOutputEnv, runId: string): Promise<RemoteWriter> {
  return env.writer(await env.data.openRun(runId));
}

export const RUN_OUTPUT_CASES: readonly ContractCase<RunOutputEnv>[] = [
  {
    name: 'every stream arrives in order and the run commits',
    async run(env) {
      const writer = await openWriter(env, 'run-a');
      const written: Record<number, Uint8Array[]> = { 0: [], 6: [], 7: [], 1: [], 3: [] };
      const send = async (stream: OutputStream, value: Uint8Array) => {
        written[stream]!.push(value);
        await writer.write(stream, value);
      };
      await send(0, bytes('Query= q1\n'));
      await send(6, bytes('q1\ts1\t100.00\n'));
      await send(7, bytes('# BLASTN\n'));
      await send(6, bytes('q1\ts2\t99.00\n'));
      await send(HITS_STREAM, bytes(hitLine(0)));
      await send(HITS_STREAM, bytes(hitLine(1).slice(0, 10)));
      await send(HITS_STREAM, bytes(hitLine(1).slice(10)));
      await send(DIAGNOSTICS_STREAM, bytes('Warning: one\n'));
      await send(0, bytes('Lambda\n'));
      await send(7, bytes('# 2 hits found\n'));
      await writer.end();
      const result = await env.data.commitRun('run-a');
      same(result.runId, 'run-a', 'runId');
      same(result.hitCount, 2, 'hitCount');
      for (const format of [0, 6, 7] as const) {
        const expected = concatBytes(written[format]!);
        same(result.byteLengths[format], expected.length, `byteLengths[${format}]`);
        sameBytes(await env.data.readOutput('run-a', format), expected, `outfmt ${format}`);
      }
      const hits = await env.data.readHits('run-a');
      same(hits.length, 2, 'HSP records');
      same(hits[1]?.index, 1, 'second HSP record');
      same(await env.data.readDiagnostics('run-a'), 'Warning: one\n', 'diagnostics');
    },
  },
  {
    name: 'many small chunks and one large chunk arrive unchanged',
    async run(env) {
      const writer = await openWriter(env, 'run-b');
      const parts: Uint8Array[] = [];
      for (let i = 0; i < 500; i++) parts.push(pattern(1, i));
      parts.push(pattern(8 * MiB, 3));
      for (const part of parts) await writer.write(6, part);
      await writer.end();
      const result = await env.data.commitRun('run-b');
      same(result.byteLengths[6], 500 + 8 * MiB, 'length');
      sameBytes(await env.data.readOutput('run-b', 6), concatBytes(parts), 'outfmt 6');
    },
  },
  {
    name: 'a run without output commits with empty outputs',
    async run(env) {
      const writer = await openWriter(env, 'run-c');
      await writer.end();
      const result = await env.data.commitRun('run-c');
      same(result.hitCount, 0, 'hitCount');
      for (const format of [0, 6, 7] as const) {
        same((await env.data.readOutput('run-c', format)).length, 0, `outfmt ${format}`);
      }
      same((await env.data.readHits('run-c')).length, 0, 'HSP records');
      same(await env.data.readDiagnostics('run-c'), '', 'diagnostics');
    },
  },
  {
    name: 'staging: a run cannot be read before it is committed, and the commit waits for end',
    async run(env) {
      const writer = await openWriter(env, 'run-d');
      await writer.write(6, bytes('row\n'));
      await rejects(env.data.readOutput('run-d', 6), /no committed result/, 'read before commit');
      const commit = env.data.commitRun('run-d');
      check(!(await settlesWithin(commit, 300)), 'the commit must wait for end');
      await writer.end();
      same((await commit).byteLengths[6], 4, 'length after end');
    },
  },
  {
    name: 'an end whose totals disagree fails the commit and discards the run',
    async run(env) {
      const writer = await openWriter(env, 'run-e');
      await writer.write(6, bytes('row\n'));
      await writer.post({ type: 'end', chunks: 2, bytes: 99 });
      await rejects(env.data.commitRun('run-e'), /incomplete/, 'commit');
      await rejects(env.data.readOutput('run-e', 6), /no committed result/, 'read after the failed commit');
    },
  },
  {
    name: 'a message outside the protocol fails the commit',
    async run(env) {
      const writer = await openWriter(env, 'run-f');
      await writer.post({ type: 'chunk', stream: 2, bytes: bytes('{}') });
      await writer.end();
      await rejects(env.data.commitRun('run-f'), /incomplete/, 'commit');
    },
  },
  {
    name: 'discard drops a staged run and ignores its late chunks',
    async run(env) {
      const baseline = await env.usage();
      const writer = await openWriter(env, 'run-g');
      await writer.write(6, bytes('first\n'));
      await env.data.discardRun('run-g');
      await writer.write(6, bytes('late\n'));
      await writer.end();
      await new Promise((resolve) => setTimeout(resolve, 300));
      same(await env.usage(), baseline, 'bytes held after the discard and the late chunks');
      await rejects(env.data.readOutput('run-g', 6), /no committed result/, 'read after discard');
      await rejects(env.data.commitRun('run-g'), /not staged/, 'commit after discard');
      const again = await openWriter(env, 'run-g');
      await again.write(6, bytes('again\n'));
      await again.end();
      same((await env.data.commitRun('run-g')).byteLengths[6], 6, 'a new run with the same id');
      sameBytes(await env.data.readOutput('run-g', 6), bytes('again\n'), 'its output');
    },
  },
  {
    name: 'a committed run is kept by discard and removed by delete',
    async run(env) {
      const writer = await openWriter(env, 'run-h');
      await writer.write(0, bytes('kept\n'));
      await writer.end();
      await env.data.commitRun('run-h');
      await env.data.discardRun('run-h');
      sameBytes(await env.data.readOutput('run-h', 0), bytes('kept\n'), 'after discard');
      await env.data.deleteRun('run-h');
      await rejects(env.data.readOutput('run-h', 0), /no committed result/, 'after delete');
    },
  },
  {
    name: 'two staged runs keep their outputs apart',
    async run(env) {
      const one = await openWriter(env, 'run-i');
      const two = await openWriter(env, 'run-j');
      await one.write(6, bytes('one-1\n'));
      await two.write(6, bytes('two-1\n'));
      await one.write(6, bytes('one-2\n'));
      await two.end();
      await one.end();
      await env.data.commitRun('run-j');
      await env.data.commitRun('run-i');
      sameBytes(await env.data.readOutput('run-i', 6), bytes('one-1\none-2\n'), 'first run');
      sameBytes(await env.data.readOutput('run-j', 6), bytes('two-1\n'), 'second run');
    },
  },
  {
    name: 'the writer refuses to write after end',
    async run(env) {
      const writer = await openWriter(env, 'run-k');
      await writer.end();
      await rejects(writer.write(6, bytes('late\n')), /already ended/, 'write after end');
      await rejects(writer.end(), /already ended/, 'a second end');
    },
  },
  {
    name: 'storage full: the run fails with the reason, and committed runs stay readable',
    async run(env) {
      const earlier = await openWriter(env, 'run-l');
      await earlier.write(6, bytes('earlier result\n'));
      await earlier.end();
      await env.data.commitRun('run-l');
      const writer = await openWriter(env, 'run-m');
      await env.exhaust();
      try {
        await writer.write(0, pattern(4 * MiB, 5));
        await writer.write(6, bytes('row\n'));
        await writer.end();
        await rejects(env.data.commitRun('run-m'), /StorageFullError: Not enough temporary storage/, 'commit');
        await rejects(env.data.readOutput('run-m', 6), /no committed result/, 'the failed run');
        sameBytes(await env.data.readOutput('run-l', 6), bytes('earlier result\n'), 'the earlier run');
      } finally {
        await env.restore();
      }
      const after = await openWriter(env, 'run-n');
      await after.write(6, bytes('after\n'));
      await after.end();
      same((await env.data.commitRun('run-n')).byteLengths[6], 6, 'a run after the storage is restored');
    },
  },
];
