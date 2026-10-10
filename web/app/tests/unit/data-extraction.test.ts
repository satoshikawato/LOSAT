// The Data worker's reads for extraction (DataService.readResidues and readHspRecords): records
// that the test writes with a record table it computes itself (support/fasta-writer.ts), given
// to the service by a stub RecordScanner, and HSP records of a committed run.
import { createHash } from 'node:crypto';
import { describe, expect, it } from 'vitest';
import type { FastaParserKind, IndexedRecord } from '../../src/domain/dataset';
import { sha256Hex } from '../../src/infra/browser/platform';
import { DataService, type DataServiceDeps } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeInputChecker } from '../../src/infra/fake/fake-fasta';
import { RunOutputWriter } from '../../src/infra/run-output/writer';
import type { HspRecord } from '../../src/ports/engine';
import type { RecordScanner } from '../../src/ports/scan';
import { FastaWriter, indexedRecord, residues, seeded, type ReaderKind } from './support/fasta-writer';

const encoder = new TextEncoder();
const latin1 = new TextDecoder('latin1');
const sha256 = (bytes: Uint8Array | string) => createHash('sha256').update(bytes).digest('hex');

/** A scanner that returns the record table that the test computed for each text it was given. */
class TableScanner implements RecordScanner {
  private readonly tables = new Map<string, readonly IndexedRecord[]>();

  add(text: string, records: readonly IndexedRecord[]): void {
    this.tables.set(text, records);
  }

  async scan(_parser: FastaParserKind, chunks: AsyncIterable<Uint8Array>): Promise<{ records: readonly IndexedRecord[] }> {
    const parts: Uint8Array[] = [];
    for await (const chunk of chunks) parts.push(chunk);
    const text = latin1.decode(Buffer.concat(parts));
    const records = this.tables.get(text);
    if (records === undefined) throw new Error('the test gave no record table for this text');
    return { records };
  }
}

/** A File that records the size of every slice read from it. */
class CountingFile extends File {
  readonly slices: number[] = [];
  override slice(start?: number, end?: number, contentType?: string): Blob {
    const blob = super.slice(start, end, contentType);
    this.slices.push(blob.size);
    return blob;
  }
}

function service(scanner: RecordScanner, overrides: Partial<DataServiceDeps> = {}) {
  let token = 0;
  const store = new MemoryBlockStore();
  const data = new DataService({
    store,
    scanner,
    checker: new FakeInputChecker(),
    digest: sha256Hex,
    newToken: () => `token-${++token}`,
    cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
    ...overrides,
  });
  return { data, store };
}

/** Writes a source, gives its record table to the scanner, and indexes it. */
async function source(data: DataService, scanner: TableScanner, writer: FastaWriter, kind: ReaderKind, name: string, file?: (text: string) => File) {
  const text = writer.toString();
  scanner.add(text, writer.records.map((record, i) => indexedRecord(record, i, kind, 5)));
  const ref = await data.addSource(file?.(text) ?? new File([text], name));
  return data.indexSource(ref.sourceId, kind as FastaParserKind);
}

const read = (bytes: Uint8Array) => latin1.decode(bytes);

describe('DataService.readResidues', () => {
  it('reads the record at a position of the run input, in buildRunInput order, with exclusions and duplicate IDs', async () => {
    const scanner = new TableScanner();
    const { data } = service(scanner, { readChunkBytes: 16 });
    const random = seeded(11);
    const first = new FastaWriter(random, 1).raw('; leading comment\n');
    const a = first.record('dup first', residues(random, 37, 'ACGTU', 0.5), { kind: 'ragged', minWidth: 2, maxWidth: 9, eol: '\n', noise: true });
    const b = first.record('dup second', residues(random, 23, 'acgtn'), { kind: 'uniform', width: 6, eol: '\r\n' });
    const c = first.record('c', residues(random, 12, 'ACGT'), { kind: 'uniform', width: 5, eol: '  \n' });
    const second = new FastaWriter(random, 1);
    const d = second.record('dup', residues(random, 15, 'ACGT', 0.2), { kind: 'ragged', minWidth: 1, maxWidth: 4, eol: '\r\n', noise: true });
    const r1 = await source(data, scanner, first, 1, 'first.fa');
    const r2 = await source(data, scanner, second, 1, 'second.fa');

    const all = [r1.revisionId, r2.revisionId];
    const input = await data.buildRunInput(all);
    expect(input.records.map((record) => record.id)).toEqual(['dup', 'dup', 'c', 'dup']);
    for (const [position, written] of [a, b, c, d].entries()) {
      const got = await data.readResidues(all, position, [{ from: 1, to: written.letters.length }, { from: 2, to: 3 }]);
      expect(got.residues.map(read)).toEqual([written.letters, written.letters.slice(1, 3)]);
    }
    const origin = (await data.readResidues(all, 3, [])).origin;
    expect(origin).toEqual({
      sourceId: r2.sourceId,
      sourceName: 'second.fa',
      revisionId: r2.revisionId,
      recordIndex: 0,
      id: 'dup',
      length: 15,
      sha256: sha256(second.toString().slice(d.headerOffset, d.endOffset)),
    });

    // Leaving out the first record shifts the positions, as it shifts q_idx and s_idx.
    const shifted = await data.reviseDataset(r1.revisionId, [0]);
    const revisions = [shifted.revisionId, r2.revisionId];
    expect((await data.buildRunInput(revisions)).records.map((record) => record.id)).toEqual(['dup', 'c', 'dup']);
    const at0 = await data.readResidues(revisions, 0, [{ from: 4, to: 20 }]);
    expect(at0.origin.recordIndex).toBe(1);
    expect(read(at0.residues[0]!)).toBe(b.letters.slice(3, 20));
    expect(read((await data.readResidues(revisions, 2, [{ from: 15, to: 15 }])).residues[0]!)).toBe(d.letters[14]);
  });

  it('reads a first record without a defline (offsets 0 and 0): from the formula, and after comment lines from its checkpoints', async () => {
    const scanner = new TableScanner();
    const { data } = service(scanner, { readChunkBytes: 16 });
    const random = seeded(29);
    // Residues from the first byte, in uniform lines: the formula counts from offset 0.
    const plain = new FastaWriter(random, 1);
    const a = plain.record(null, residues(random, 47, 'ACGTU', 0.3), { kind: 'uniform', width: 10, eol: '\n' });
    const b = plain.record('b second', residues(random, 20, 'ACGT'), { kind: 'uniform', width: 7, eol: '\r\n' });
    // After a comment and a blank line, the record still starts the input (its bytes include them).
    const commented = new FastaWriter(random, 2).raw('; a comment\n\n');
    const c = commented.record(null, residues(random, 30, 'ACDEFGHIKLMNPQRSTVWY*', 0.2), { kind: 'ragged', minWidth: 3, maxWidth: 8, eol: '\n', noise: true });
    const r1 = await source(data, scanner, plain, 1, 'plain.fa');
    const r2 = await source(data, scanner, commented, 2, 'commented.fa');
    expect(r1.records.map((r) => [r.id, r.header_offset, r.sequence_offset, r.line_layout.kind])).toEqual([
      ['', 0, 0, 'uniform'],
      ['b', b.headerOffset, b.sequenceOffset, 'uniform'],
    ]);
    expect(r2.records.map((r) => [r.id, r.header_offset, r.sequence_offset, r.line_layout.kind])).toEqual([['', 0, 0, 'checkpoints']]);

    for (const [revision, written, text] of [
      [r1, a, plain.toString()],
      [r2, c, commented.toString()],
    ] as const) {
      const length = written.letters.length;
      const got = await data.readResidues([revision.revisionId], 0, [{ from: 1, to: length }, { from: 1, to: 1 }, { from: 12, to: length - 3 }]);
      expect(got.residues.map(read)).toEqual([written.letters, written.letters[0], written.letters.slice(11, length - 3)]);
      expect(got.origin).toMatchObject({ recordIndex: 0, id: '', length, sha256: sha256(text.slice(0, written.endOffset)) });
    }
    expect(read((await data.readResidues([r1.revisionId], 1, [{ from: 1, to: 20 }])).residues[0]!)).toBe(b.letters);
    // Left out, the record without a defline shifts the positions as any other.
    const withoutFirst = await data.reviseDataset(r1.revisionId, [0]);
    const shifted = await data.readResidues([withoutFirst.revisionId], 0, [{ from: 3, to: 9 }]);
    expect([shifted.origin.id, shifted.origin.recordIndex, read(shifted.residues[0]!)]).toEqual(['b', 1, b.letters.slice(2, 9)]);
  });

  it('reads a large record in bounded slices of its File, from its checkpoints', async () => {
    const scanner = new TableScanner();
    const { data } = service(scanner, { readChunkBytes: 4096 });
    const random = seeded(5);
    const writer = new FastaWriter(random, 2);
    const big = writer.record('big protein', residues(random, 100_000, 'ACDEFGHIKLMNPQRSTVWY*', 0.1), { kind: 'ragged', minWidth: 40, maxWidth: 80, eol: '\n', noise: true });
    let file: CountingFile | undefined;
    const revision = await source(data, scanner, writer, 2, 'big.fa', (text) => (file = new CountingFile([text], 'big.fa')));
    file!.slices.length = 0;
    const got = await data.readResidues([revision.revisionId], 0, [{ from: 1, to: 100_000 }, { from: 99_990, to: 100_000 }, { from: 7, to: 7 }]);
    expect(got.residues.map(read)).toEqual([big.letters, big.letters.slice(99_989), big.letters[6]]);
    expect(Math.max(...file!.slices)).toBeLessThanOrEqual(4096);
    expect(file!.slices.length).toBeGreaterThan(100_000 / 4096);
    // Each array owns its whole buffer, so the Data worker transfers it (rpc.ts).
    for (const bytes of got.residues) expect(bytes.byteLength).toBe(bytes.buffer.byteLength);
  });

  it('refuses intervals outside the record and positions outside the run input', async () => {
    const scanner = new TableScanner();
    const { data } = service(scanner);
    const writer = new FastaWriter(seeded(2), 1);
    writer.record('x', 'ACGTACGTAC', { kind: 'uniform', width: 4, eol: '\n' });
    const revision = await source(data, scanner, writer, 1, 'x.fa');
    const ids = [revision.revisionId];
    await expect(data.readResidues(ids, 0, [{ from: 0, to: 3 }])).rejects.toThrow(RangeError);
    await expect(data.readResidues(ids, 0, [{ from: 5, to: 11 }])).rejects.toThrow(/not an interval of record 1 \("x", 10 letters\)/);
    await expect(data.readResidues(ids, 0, [{ from: 5, to: 4 }])).rejects.toThrow(RangeError);
    await expect(data.readResidues(ids, 1, [{ from: 1, to: 1 }])).rejects.toThrow(/has 1 records, so it has no record 2/);
    await expect(data.readResidues(ids, -1, [])).rejects.toThrow(RangeError);
    await expect(data.readResidues([], 0, [])).rejects.toThrow(/at least one dataset revision/);
  });

  it('refuses a source that no longer matches its record table', async () => {
    const scanner = new TableScanner();
    const { data } = service(scanner);
    const text = '>x\nACGTAC\nGT\n';
    const record: IndexedRecord = {
      index: 0,
      id: 'x',
      header_offset: 0,
      sequence_offset: 3,
      end_offset: text.length,
      length: 8,
      line_layout: { kind: 'checkpoints', every: 65_536, offsets: [3] },
      residue_counts: { A: 2, C: 2, G: 2, T: 2 },
    };
    // The table counts other residues than the file has.
    scanner.add(text, [{ ...record, residue_counts: { A: 3, C: 1, G: 2, T: 2 } }]);
    const changed = await data.indexSource((await data.addSource(new File([text], 'x.fa'))).sourceId, 1);
    await expect(data.readResidues([changed.revisionId], 0, [{ from: 1, to: 8 }])).rejects.toThrow(
      'The source "x.fa" no longer matches its record table: record 1 ("x"), read for 1-8: its residues are not those counted in the record table. Add the file again.',
    );
    // A part of the record is not checked against the counts.
    expect(read((await data.readResidues([changed.revisionId], 0, [{ from: 2, to: 7 }])).residues[0]!)).toBe('CGTACG');
    // The table has more residues than the file.
    const longer = `${text}\n`;
    scanner.add(longer, [{ ...record, end_offset: longer.length, length: 9 }]);
    const short = await data.indexSource((await data.addSource(new File([longer], 'y.fa'))).sourceId, 1);
    await expect(data.readResidues([short.revisionId], 0, [{ from: 5, to: 9 }])).rejects.toThrow(/only 4 of its 5 residues were found/);
  });

  it('refuses records indexed with another reader kind than 1 or 2', async () => {
    const scanner = new TableScanner();
    const { data } = service(scanner);
    const writer = new FastaWriter(seeded(3), 1);
    writer.record('x', 'ACGT', { kind: 'uniform', width: 4, eol: '\n' });
    const revision = await source(data, scanner, writer, 0 as ReaderKind, 'x.fa');
    await expect(data.readResidues([revision.revisionId], 0, [{ from: 1, to: 4 }])).rejects.toThrow(/kind 0 cannot be extracted/);
  });
});

describe('DataService.readHspRecords', () => {
  const hit = (index: number, aligned: string): HspRecord => ({
    index,
    q_idx: 0,
    s_idx: index % 3,
    rank: index,
    raw_score: 10,
    bit_score: 20.5,
    e_value: 1e-5,
    q_start: 1,
    q_end: aligned.length,
    s_start: aligned.length,
    s_end: 1,
    query_frame: null,
    subject_frame: null,
    subject_length: 100,
    query_aligned: aligned,
    subject_aligned: aligned.toLowerCase(),
    out6: [index * 10, index * 10 + 10],
    out0: null,
    out0_subject: null,
  });

  async function committed(lines: readonly string[], chunkBytes: number) {
    const store = new MemoryBlockStore();
    const reads: number[] = [];
    const read = store.read.bind(store);
    store.read = async (path, offset, length) => {
      reads.push(length);
      return read(path, offset, length);
    };
    const { data } = service(new TableScanner(), { store, readChunkBytes: chunkBytes });
    const writer = new RunOutputWriter(await data.openRun('run-1'));
    // Split anywhere, with blank lines, as the engine's chunks may come.
    const bytes = encoder.encode(lines.join('\n'));
    for (let at = 0; at < bytes.length; at += 7) writer.write(1, bytes.subarray(at, at + 7));
    writer.end();
    await data.commitRun('run-1');
    return { data, reads };
  }

  it('returns the records of the indices asked, in that order, with their aligned rows, reading in bounded ranges', async () => {
    const records = Array.from({ length: 50 }, (_, i) => hit(i, `ACGT-${'A'.repeat(i)}`));
    const lines = records.flatMap((record, i) => (i % 10 === 3 ? ['', '  ', JSON.stringify(record)] : [JSON.stringify(record)]));
    const { data, reads } = await committed(lines, 64);
    const got = await data.readHspRecords('run-1', [42, 0, 7, 42]);
    expect(got).toEqual([records[42], records[0], records[7], records[42]]);
    expect(Math.max(...reads)).toBeLessThanOrEqual(Math.max(64, ...lines.map((line) => line.length)));
    reads.length = 0;
    expect(await data.readHspRecords('run-1', [49])).toEqual([records[49]]);
    expect(reads).toEqual([JSON.stringify(records[49]).length]);
    await expect(data.readHspRecords('run-1', [50])).rejects.toThrow('run run-1 has no HSP 50 (it has 50 HSP records)');
    await expect(data.readHspRecords('run-1', [-1])).rejects.toThrow(RangeError);
    await expect(data.readHspRecords('run-2', [0])).rejects.toThrow('run run-2 has no committed result');
    expect(await data.readHspRecords('run-1', [])).toEqual([]);
  });

  it('finds the records by their index where the lines are not in index order', async () => {
    const records = Array.from({ length: 12 }, (_, i) => hit(i, 'MKV'.repeat(i + 1)));
    const shuffled = [...records].reverse();
    const { data } = await committed(shuffled.map((record) => JSON.stringify(record)), 50);
    const got = await data.readHspRecords('run-1', [3, 11, 0]);
    expect(got.map((record) => record.index)).toEqual([3, 11, 0]);
    expect(got).toEqual([records[3], records[11], records[0]]);
  });
});
