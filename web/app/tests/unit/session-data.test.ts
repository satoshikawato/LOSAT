// What session files need of the infrastructure: the Data worker's description of a run input,
// its reading of a run's streams and its count of a staged run's bytes (infra/data/data-service.ts),
// and the gzip sink and reader of the browser (infra/browser/compression.ts).
import { createHash } from 'node:crypto';
import { gunzipSync, gzipSync } from 'node:zlib';
import { describe, expect, it } from 'vitest';
import { sha256Hex } from '../../src/infra/browser/platform';
import { browserCompression } from '../../src/infra/browser/compression';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import { RunOutputWriter } from '../../src/infra/run-output/writer';
import { DIAGNOSTICS_STREAM, HITS_STREAM } from '../../src/ports/run-output';

const encoder = new TextEncoder();
const sha256 = (text: string) => createHash('sha256').update(text).digest('hex');

function service(readChunkBytes?: number) {
  let token = 0;
  return new DataService({
    ...(readChunkBytes === undefined ? {} : { readChunkBytes }),
    store: new MemoryBlockStore(),
    scanner: new FakeScanner(),
    checker: new FakeInputChecker(),
    digest: sha256Hex,
    newToken: () => `token-${++token}`,
    cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
  });
}

const tick = () => new Promise((resolve) => setTimeout(resolve, 0));

describe('DataService for session files', () => {
  it('describes a run input: its records in order, with their SHA-256, and its sources with the exclusions', async () => {
    const data = service();
    const a = '>a1 x\nACGT\n>a2\nGG\n>a3\nTTT\n';
    const b = '>b1\nCCCC\n';
    const ra = await data.indexSource((await data.addSource(new File([a], 'a.fa'))).sourceId, 1);
    const revised = await data.reviseDataset(ra.revisionId, [1]);
    const rb = await data.indexSource((await data.addSource(new File([b], 'b.fa'))).sourceId, 1);
    expect(await data.describeRunInput([revised.revisionId, rb.revisionId])).toEqual({
      reader: 1,
      records: { id: ['a1', 'a3', 'b1'], length: [4, 3, 4], sha256: [sha256('>a1 x\nACGT\n'), sha256('>a3\nTTT\n'), sha256(b)] },
      sources: [
        { name: 'a.fa', size: a.length, records: 3, excluded: [1] },
        { name: 'b.fa', size: b.length, records: 1, excluded: [] },
      ],
    });
    const protein = await data.indexSource((await data.addSource(new File(['>p\nMK\n'], 'p.fa'))).sourceId, 2);
    await expect(data.describeRunInput([rb.revisionId, protein.revisionId])).rejects.toThrow(/different reader kinds/);
  });

  it('reads the streams of a committed run in ranges, and gives their lengths', async () => {
    const data = service();
    const port = await data.openRun('r');
    const writer = new RunOutputWriter(port);
    writer.write(0, encoder.encode('outfmt 0 text'));
    writer.write(HITS_STREAM, encoder.encode('{"index":0}\n'));
    writer.write(DIAGNOSTICS_STREAM, encoder.encode('Warning\n'));
    writer.end();
    await data.commitRun('r');
    expect(await data.runBlockLengths('r')).toEqual({ 0: 13, 6: 0, 7: 0, [HITS_STREAM]: 12, [DIAGNOSTICS_STREAM]: 8 });
    expect(new TextDecoder().decode(await data.readRunBlock('r', HITS_STREAM, 2, 7))).toBe('index');
    expect(new TextDecoder().decode(await data.readRunBlock('r', DIAGNOSTICS_STREAM, 0, 8))).toBe('Warning\n');
    await expect(data.readRunBlock('r', 6, 0, 1)).rejects.toThrow(RangeError);
    await expect(data.readRunBlock('r', 2 as never, 0, 0)).rejects.toThrow(/not a stream/);
    await expect(data.runBlockLengths('other')).rejects.toThrow(/no committed result/);
  });

  it('checks the HSP records of a committed run as the JSON of their lines, read in bounded ranges', async () => {
    const data = service(64);
    const record = (index: number, change: Record<string, unknown> = {}) =>
      JSON.stringify({
        index,
        q_idx: 0,
        s_idx: 0,
        rank: index,
        raw_score: 1,
        bit_score: 2,
        e_value: 0.5,
        q_start: 1,
        q_end: 4,
        s_start: 1,
        s_end: 4,
        query_frame: null,
        subject_frame: null,
        subject_length: 4,
        query_aligned: 'ACGT',
        subject_aligned: 'ACGT',
        out6: null,
        out0: null,
        out0_subject: null,
        ...change,
      });
    const committed = async (runId: string, lines: readonly string[]) => {
      const writer = new RunOutputWriter(await data.openRun(runId));
      writer.write(HITS_STREAM, encoder.encode(`${lines.join('\n')}\n`));
      writer.end();
      await data.commitRun(runId);
    };
    const bounds = { count: 3, queries: 1, subjects: 1, out0: 0, out6: 0 };
    await committed('good', [record(0), record(1), record(2)]);
    expect(await data.checkHspRecords('good', bounds)).toBeUndefined();
    // The records that come after a refused one are not read.
    await committed('bad', [record(0), record(1, { s_idx: null }), 'not JSON']);
    expect(await data.checkHspRecords('bad', bounds)).toBe('HSP record 2 has s_idx null, not a whole number of 0 or more');
    await committed('broken', [record(0), 'not JSON', record(2)]);
    await expect(data.checkHspRecords('broken', bounds)).rejects.toThrow(/^HSP record 2 is not JSON/);
    expect(await data.checkHspRecords('good', { ...bounds, count: 4 })).toBe('there are 3 HSP records, but the manifest gives 4');
    // The lines found for the check serve readHspRecords too.
    expect((await data.readHspRecords('good', [2]))[0]!.rank).toBe(2);
  });

  it('tells a writer how many bytes of a staged run are stored, once they reach a count', async () => {
    const data = service();
    const port = await data.openRun('r');
    const writer = new RunOutputWriter(port);
    let stored: number | undefined;
    void data.stagedBytes('r', 10).then((n) => (stored = n));
    writer.write(0, encoder.encode('12345'));
    await tick();
    await tick();
    expect(stored).toBeUndefined();
    writer.write(6, encoder.encode('678901'));
    await expect(data.stagedBytes('r', 10)).resolves.toBe(11);
    expect(stored).toBe(11);
    expect(await data.stagedBytes('r', 3)).toBe(11);
    // A dropped run, and one that is not staged, keep no writer waiting.
    const waiting = data.stagedBytes('r', 1000);
    await data.deleteRun('r');
    await expect(waiting).resolves.toBe(11);
    expect(await data.stagedBytes('unknown', 5)).toBe(0);
  });
});

describe('browserCompression', () => {
  it('gzips blocks as they are written and hands the compressed blocks on in order', async () => {
    const out: Uint8Array[] = [];
    const sink = browserCompression.gzip(async (bytes) => void out.push(bytes.slice()));
    const text = 'LOSAT-WEB-SESSION 1\n'.repeat(5000);
    const block = encoder.encode(text);
    for (let at = 0; at < block.length; at += 1000) await sink.write(block.subarray(at, at + 1000));
    await sink.close();
    expect(new TextDecoder().decode(gunzipSync(Buffer.concat(out)))).toBe(text);
    await expect(sink.write(block)).rejects.toThrow(/no longer open/);
    // Data that does not compress leaves in blocks of about 1 MiB, not as one file.
    const blocks: Uint8Array[] = [];
    const big = browserCompression.gzip(async (bytes) => void blocks.push(bytes.slice()));
    let state = 1;
    const noise = new Uint8Array(3 << 20).map(() => (state = (Math.imul(state, 1103515245) + 12345) >>> 0) >>> 24);
    for (let at = 0; at < noise.length; at += 1 << 16) await big.write(noise.subarray(at, at + (1 << 16)));
    await big.close();
    expect(blocks.length).toBeGreaterThanOrEqual(3);
    expect(Math.max(...blocks.map((b) => b.length))).toBeLessThan(2 << 20);
    expect(Buffer.from(gunzipSync(Buffer.concat(blocks))).equals(Buffer.from(noise))).toBe(true);
  });

  it('fails the writes when the compressed blocks cannot be handed on, and aborts without output', async () => {
    const failing = browserCompression.gzip(async () => {
      throw new Error('disk full');
    });
    const big = new Uint8Array(4 << 20).map((_, i) => (i * 7919) % 251);
    await expect(
      (async () => {
        for (let i = 0; i < 8; i++) await failing.write(big);
        await failing.close();
      })(),
    ).rejects.toThrow('disk full');
    const out: Uint8Array[] = [];
    const aborted = browserCompression.gzip(async (bytes) => void out.push(bytes));
    await aborted.write(encoder.encode('x'));
    aborted.abort();
    aborted.abort();
    await expect(aborted.close()).rejects.toThrow(/no longer open/);
  });

  it('reads gzip data in chunks, and throws on damaged data', async () => {
    const text = 'ACGT'.repeat(100_000);
    const chunks: Uint8Array[] = [];
    for await (const chunk of browserCompression.gunzip(new Blob([gzipSync(text)]))) chunks.push(chunk);
    expect(new TextDecoder().decode(Buffer.concat(chunks))).toBe(text);
    const damaged = gzipSync(text);
    damaged[damaged.length - 6]! ^= 1;
    await expect(
      (async () => {
        for await (const chunk of browserCompression.gunzip(new Blob([damaged]))) void chunk;
      })(),
    ).rejects.toThrow();
    // Stopping early stops the reading.
    for await (const chunk of browserCompression.gunzip(new Blob([gzipSync(text)]))) {
      expect(chunk.length).toBeGreaterThan(0);
      break;
    }
  });
});
