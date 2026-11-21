// Run output channel contract (losat/ports run-output.ts, ports/data.ts RunStore): an engine-side
// writer sends the outputs of a run over the MessagePort of `openRun`, and the data layer
// stages, commits, discards and reads them. Vitest runs it in one thread; Playwright runs
// it with the writer in a separate worker and the real Data worker (OPFS). S09 runs it with
// the writer inside the real Engine worker.
import { CUSTOM_OUTPUT, type RunOutput } from 'losat/domain';
import { DIAGNOSTICS_STREAM, HITS_STREAM, type OutputStream } from 'losat/ports';
import type { RunStore } from '../../src/ports/data';
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
/** A ZIP archive of stored entries, as the engine writes the archive of outfmt 13 and 14 (docs/web/abi_v2.md §5). */
function zip(entries: ReadonlyArray<readonly [string, string]>): Uint8Array {
  const parts: Uint8Array[] = [];
  const directory: Uint8Array[] = [];
  let offset = 0;
  const u16 = (value: number) => Uint8Array.of(value & 255, value >> 8);
  const u32 = (value: number) => Uint8Array.of(value & 255, (value >> 8) & 255, (value >> 16) & 255, (value >>> 24) & 255);
  for (const [name, text] of entries) {
    const data = new TextEncoder().encode(text);
    const label = new TextEncoder().encode(name);
    const local = concatBytes([u32(0x04034b50), u16(10), u16(0), u16(0), u16(0), u16(0x21), u32(0), u32(data.length), u32(data.length), u16(label.length), u16(0), label, data]);
    directory.push(concatBytes([u32(0x02014b50), u16(10), u16(10), u16(0), u16(0), u16(0), u16(0x21), u32(0), u32(data.length), u32(data.length), u16(label.length), u16(0), u16(0), u16(0), u16(0), u32(0), u32(offset), label]));
    parts.push(local);
    offset += local.length;
  }
  const central = concatBytes(directory);
  return concatBytes([...parts, central, u32(0x06054b50), u16(0), u16(0), u16(entries.length), u16(entries.length), u32(central.length), u32(offset), u16(0)]);
}

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

async function openWriter(env: RunOutputEnv, runId: string, outputs?: readonly RunOutput[]): Promise<RemoteWriter> {
  return env.writer(await env.data.openRun(runId, outputs));
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
    name: "a discard drops what arrived of one output, which stays empty; the others and the end's totals are kept",
    async run(env) {
      const baseline = await env.usage();
      const writer = await openWriter(env, 'run-o');
      await writer.write(6, pattern(2 * MiB, 1));
      await writer.write(0, bytes('Query= q1\n'));
      await writer.write(6, pattern(MiB, 2));
      // The engine's signal of stream 9 (ABI v2 §5): outfmt 6 failed alone after these chunks, its error on stream 3.
      await writer.post({ type: 'discard', stream: 6, error: 'diagnostics' });
      await writer.write(7, bytes('# BLASTN\n'));
      await writer.write(DIAGNOSTICS_STREAM, bytes('Error: refused\n'));
      await writer.write(0, bytes('Lambda\n'));
      await writer.end();
      const result = await env.data.commitRun('run-o');
      same(result.byteLengths[6], 0, 'outfmt 6');
      same((await env.data.readOutput('run-o', 6)).length, 0, 'outfmt 6 bytes');
      sameBytes(await env.data.readOutput('run-o', 0), bytes('Query= q1\nLambda\n'), 'outfmt 0');
      sameBytes(await env.data.readOutput('run-o', 7), bytes('# BLASTN\n'), 'outfmt 7');
      same(await env.data.readDiagnostics('run-o'), 'Error: refused\n', 'diagnostics');
      same(result.hitCount, 0, 'hitCount');
      same(JSON.stringify(result.failed), JSON.stringify([{ output: 6, error: 'diagnostics' }]), 'failed');
      // The 3 MiB of outfmt 6 are not held: only the 41 bytes of the other streams.
      same(await env.usage(), baseline + 41, 'bytes held after the commit');
    },
  },
  {
    name: 'a discard of a stream that is not an output fails the commit',
    async run(env) {
      const writer = await openWriter(env, 'run-p');
      await writer.write(HITS_STREAM, bytes(hitLine(0)));
      await writer.post({ type: 'discard', stream: HITS_STREAM, error: 'diagnostics' });
      await writer.end();
      await rejects(env.data.commitRun('run-p'), /incomplete/, 'commit');
    },
  },
  {
    // Line WQA's sweep, D3 review (2026-11-19) and D15 (2026-11-20): a run that the engine failed after it began to
    // write ends with `fail`; the data layer waits for it, then commits the run with the outputs that the engine did
    // not discard (the bytes that the CLI leaves before its error) and reads its stream 3.
    name: "a failed run is kept with the outputs that the engine did not discard, and its warnings are read after its fail",
    async run(env) {
      const writer = await openWriter(env, 'run-q');
      const warnings = [bytes('Warning: [blastn] Subject_1 empty record: Subject sequence contains no data\n'), bytes('Warning: caf\u00e9\n')];
      await writer.write(6, pattern(MiB, 4));
      await writer.write(6, bytes('last row\n'));
      // What the engine wrote of an output that it then drops stays in the totals of `fail`, and is dropped.
      await writer.write(0, bytes('dropped\n'));
      await writer.post({ type: 'discard', stream: 0, error: 'run' });
      await writer.post({ type: 'discard', stream: 7, error: 'run' });
      await writer.write(DIAGNOSTICS_STREAM, warnings[0]!);
      // Its UTF-8 split across two chunks.
      await writer.write(DIAGNOSTICS_STREAM, warnings[1]!.subarray(0, 13));
      await writer.write(DIAGNOSTICS_STREAM, warnings[1]!.subarray(13));
      const kept = env.data.keepFailedRun('run-q');
      check(!(await settlesWithin(kept, 300)), 'the keep must wait for fail');
      await writer.post({ type: 'fail', chunks: 6, bytes: MiB + 9 + 8 + warnings[0]!.length + warnings[1]!.length });
      const run = await kept;
      check(run !== undefined, 'the failed run is kept');
      same(run!.diagnostics, 'Warning: [blastn] Subject_1 empty record: Subject sequence contains no data\nWarning: caf\u00e9\n', 'warnings');
      same(JSON.stringify(run!.result!.byteLengths), JSON.stringify({ 0: 0, 6: MiB + 9, 7: 0 }), 'byte lengths');
      same(run!.result!.hitCount, 0, 'no HSP records');
      same(JSON.stringify(run!.result!.failed), JSON.stringify([{ output: 0, error: 'run' }, { output: 7, error: 'run' }]), 'dropped outputs');
      sameBytes(await env.data.readOutput('run-q', 6), concatBytes([pattern(MiB, 4), bytes('last row\n')]), 'its output');
      same((await env.data.readOutput('run-q', 0)).length, 0, 'a dropped output is empty');
      same(await env.data.readDiagnostics('run-q'), run!.diagnostics, 'diagnostics as a completed run reads them');
      await rejects(env.data.commitRun('run-q'), /not staged/, 'commit after the keep');
      same(await env.data.keepFailedRun('run-q'), undefined, 'a second keep');
      await env.data.deleteRun('run-q');
      await rejects(env.data.readOutput('run-q', 6), /no committed result/, 'read after delete');

      // A fail whose totals disagree keeps nothing; a commit of a failed run fails.
      const baseline = await env.usage();
      const broken = await openWriter(env, 'run-r');
      await broken.write(DIAGNOSTICS_STREAM, bytes('Warning: one\n'));
      await broken.post({ type: 'fail', chunks: 2, bytes: 99 });
      same(await env.data.keepFailedRun('run-r'), undefined, 'incomplete');
      const failed = await openWriter(env, 'run-s');
      await failed.write(6, bytes('row\n'));
      await failed.post({ type: 'fail', chunks: 1, bytes: 4 });
      await rejects(env.data.commitRun('run-s'), /failed/, 'commit of a failed run');
      same(await env.usage(), baseline, 'bytes held at the end');
    },
  },
  {
    // D15 round 2 (2026-11-21): a failed run that wrote nothing to the outputs it keeps (an empty file, an archive of
    // empty entries) keeps only its warnings, and holds no blocks; an archive with data in an entry is kept.
    name: 'a failed run that kept no bytes is dropped, and keeps its warnings',
    async run(env) {
      const baseline = await env.usage();
      const empty = await openWriter(env, 'run-u');
      await empty.write(DIAGNOSTICS_STREAM, bytes('Warning: one\n'));
      await empty.post({ type: 'fail', chunks: 1, bytes: 13 });
      const kept = await env.data.keepFailedRun('run-u');
      same(kept?.diagnostics, 'Warning: one\n', 'warnings of a run without output');
      same(kept?.result, undefined, 'no result');
      await rejects(env.data.readOutput('run-u', 6), /no committed result/, 'the run is dropped');
      same(await env.usage(), baseline, 'bytes held at the end');

      const emptyArchive = zip([['results.json', '']]);
      const hollow = await openWriter(env, 'run-v', [CUSTOM_OUTPUT]);
      await hollow.write(CUSTOM_OUTPUT, emptyArchive);
      await hollow.write(DIAGNOSTICS_STREAM, bytes('Warning: two\n'));
      await hollow.post({ type: 'fail', chunks: 2, bytes: emptyArchive.length + 13 });
      const dropped = await env.data.keepFailedRun('run-v');
      same(dropped?.result, undefined, 'an archive of empty entries is no output');
      same(dropped?.diagnostics, 'Warning: two\n', 'its warnings');

      const full = zip([['results_1.json', '{"a":1}'], ['results.json', '']]);
      const filled = await openWriter(env, 'run-w', [CUSTOM_OUTPUT]);
      await filled.write(CUSTOM_OUTPUT, full);
      await filled.post({ type: 'fail', chunks: 1, bytes: full.length });
      const keptArchive = await env.data.keepFailedRun('run-w');
      check(keptArchive?.result !== undefined, 'an archive with data in an entry is kept');
      sameBytes(await env.data.readOutput('run-w', CUSTOM_OUTPUT), full, 'its archive');
      await env.data.deleteRun('run-w');
      same(await env.usage(), baseline, 'bytes held at the end of the archives');
    },
  },
  {
    name: 'a discard without a place of its error fails the commit',
    async run(env) {
      const writer = await openWriter(env, 'run-t');
      await writer.write(6, bytes('row\n'));
      await writer.post({ type: 'discard', stream: 6 });
      await writer.end();
      await rejects(env.data.commitRun('run-t'), /incomplete/, 'commit');
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
