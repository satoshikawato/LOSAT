import { createHash } from 'node:crypto';
import { describe, expect, it } from 'vitest';
import type { FastaParserKind, IndexedRecord } from '../../src/domain/dataset';
import { sha256Hex } from '../../src/infra/browser/platform';
import { MEMORY_FULL_MESSAGE } from '../../src/infra/data/block-store';
import { DataService, type DataServiceDeps } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { lineInMessage, ncbiLineStart } from '../../src/infra/data/message-line';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import { RunOutputWriter } from '../../src/infra/run-output/writer';
import type { InputChecker } from '../../src/ports/input-check';
import type { RecordScanner } from '../../src/ports/scan';

const encoder = new TextEncoder();
const decoder = new TextDecoder();
const sha256 = (bytes: Uint8Array | string) => createHash('sha256').update(bytes).digest('hex');

function service(overrides: Partial<DataServiceDeps> & { store?: MemoryBlockStore } = {}) {
  let token = 0;
  const store = overrides.store ?? new MemoryBlockStore();
  const data = new DataService({
    store,
    scanner: new FakeScanner(),
    checker: new FakeInputChecker(),
    digest: sha256Hex,
    newToken: () => `token-${++token}`,
    cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
    ...overrides,
  });
  return { data, store };
}

const THREE = '>a one\nACGT\nAC\n>b\nGG\n>c three\nTTT';

describe('DataService sources and record tables', () => {
  it('keeps a File reference and builds its record table with per-record SHA-256', async () => {
    const { data } = service({ readChunkBytes: 5 });
    const file = new File([THREE], 'three.fa');
    const source = await data.addSource(file);
    expect(source).toEqual({ sourceId: 'token-1', name: 'three.fa', size: file.size });
    const revision = await data.indexSource(source.sourceId, 1);
    expect(Object.isFrozen(revision)).toBe(true);
    expect(revision.excluded).toEqual([]);
    expect(revision.records.map((r) => [r.id, r.length])).toEqual([
      ['a', 6],
      ['b', 2],
      ['c', 3],
    ]);
    for (const record of revision.records) {
      expect(record.sha256).toBe(sha256(THREE.slice(record.header_offset, record.end_offset)));
    }
  });

  it('keeps records with the same ID apart and keeps the original headers and case (REQ-05)', async () => {
    const { data } = service();
    const text = '>dup First Header\nacgtACGT\n>dup second\nGGCC\n';
    const source = await data.addSource(new File([text], 'dup.fa'));
    const revision = await data.indexSource(source.sourceId, 1);
    expect(revision.records.map((r) => [r.index, r.id, r.length])).toEqual([
      [0, 'dup', 8],
      [1, 'dup', 4],
    ]);
    expect(revision.records[0]!.sha256).not.toBe(revision.records[1]!.sha256);
    // The residues are counted upper-cased; the run input keeps the original bytes.
    expect(revision.records[0]!.residue_counts).toEqual({ A: 2, C: 2, G: 2, T: 2 });
    const second = await data.reviseDataset(revision.revisionId, [0]);
    const input = await data.buildRunInput([second.revisionId]);
    expect(decoder.decode(input.bytes)).toBe('>dup second\nGGCC\n');
  });

  it('builds the run input and checks an input of more records than a call can take as arguments', async () => {
    const { data } = service();
    const count = 150_000;
    const text = Array.from({ length: count }, (_, i) => `>r${i}\nA\n`).join('');
    const source = await data.addSource(new File([text], 'many.fa'));
    const revision = await data.indexSource(source.sourceId, 1);
    const excluded = await data.reviseDataset(revision.revisionId, [0]);
    const input = await data.buildRunInput([excluded.revisionId]);
    expect(input.records).toHaveLength(count - 1);
    expect(input.records[0]).toEqual({ id: 'r1', length: 1 });
    expect(await data.checkInput('blastn', 'query', [excluded.revisionId])).toMatchObject({ ok: true });
  }, 60_000);

  it("rejects a source with the index scan's error", async () => {
    const message = "CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.";
    const failing: RecordScanner = {
      async scan() {
        throw new Error(message);
      },
    };
    const { data } = service({ scanner: failing });
    const source = await data.addSource(new File(['>a\n1234567890\n'], 'digits.fa'));
    await expect(data.indexSource(source.sourceId, 1)).rejects.toThrow(message);
    // The app has the kinds of the engine's reader only (kind 0 stays the adapter's).
    const fake = service().data;
    const other = await fake.addSource(new File(['>a\nAC\n'], 'a.fa'));
    await expect(fake.indexSource(other.sourceId, 0 as unknown as FastaParserKind)).rejects.toThrow('parser kind 0');
  });

  it('rejects a record table whose offsets are not ordered inside the source', async () => {
    const broken: RecordScanner = {
      async scan() {
        const record: IndexedRecord = {
          index: 0,
          id: 'x',
          header_offset: 0,
          sequence_offset: 3,
          end_offset: 999,
          length: 1,
          line_layout: { kind: 'uniform', width: 1, eol: 1 },
          residue_counts: { A: 1 },
        };
        return { records: [record] };
      },
    };
    const { data } = service({ scanner: broken });
    const source = await data.addSource(new File(['>x\nA\n'], 'x.fa'));
    await expect(data.indexSource(source.sourceId, 1)).rejects.toThrow('invalid record table');
  });

  it('makes a new revision for a new selection and refuses unknown records', async () => {
    const { data } = service();
    const source = await data.addSource(new File([THREE], 'three.fa'));
    const base = await data.indexSource(source.sourceId, 1);
    const revised = await data.reviseDataset(base.revisionId, [2, 0, 2]);
    expect(revised.revisionId).not.toBe(base.revisionId);
    expect(revised.excluded).toEqual([0, 2]);
    expect(base.excluded).toEqual([]);
    await expect(data.reviseDataset(base.revisionId, [3])).rejects.toThrow('not in the record table');
  });

  it('gives the source unchanged when every record is included', async () => {
    const { data } = service();
    const text = `${THREE}\n\n`;
    const source = await data.addSource(new File([text], 'three.fa'));
    const revision = await data.indexSource(source.sourceId, 1);
    const input = await data.buildRunInput([revision.revisionId]);
    expect(decoder.decode(input.bytes)).toBe(text);
    expect(input.sha256).toBe(sha256(text));
    expect(input.records).toEqual([
      { id: 'a', length: 6 },
      { id: 'b', length: 2 },
      { id: 'c', length: 3 },
    ]);
  });

  it('gives the included original records in order when records are left out', async () => {
    const { data } = service();
    const source = await data.addSource(new File([THREE], 'three.fa'));
    const revision = await data.reviseDataset((await data.indexSource(source.sourceId, 1)).revisionId, [1]);
    const input = await data.buildRunInput([revision.revisionId]);
    expect(decoder.decode(input.bytes)).toBe('>a one\nACGT\nAC\n>c three\nTTT');
    expect(input.records).toEqual([
      { id: 'a', length: 6 },
      { id: 'c', length: 3 },
    ]);
  });

  it('adds one newline after a source without a final newline when another follows', async () => {
    const { data } = service();
    const first = await data.indexSource((await data.addSource(new File(['>x\nAC'], 'x.fa'))).sourceId, 1);
    const second = await data.indexSource((await data.addSource(new File(['>y\nGG\n'], 'y.fa'))).sourceId, 1);
    const input = await data.buildRunInput([first.revisionId, second.revisionId, first.revisionId]);
    expect(decoder.decode(input.bytes)).toBe('>x\nAC\n>y\nGG\n>x\nAC');
    expect(input.records.map((r) => r.id)).toEqual(['x', 'y', 'x']);
    // A source whose records are all left out adds nothing, not even the newline.
    const none = await data.reviseDataset(second.revisionId, [0]);
    expect(decoder.decode((await data.buildRunInput([first.revisionId, none.revisionId])).bytes)).toBe('>x\nAC');
  });

  it('keeps a first record without a defline: its bytes start the source, and it can be left out', async () => {
    const { data } = service();
    const text = '#c\nACGT\nAC\n>b\nGG\n';
    const revision = await data.indexSource((await data.addSource(new File([text], 'plain.fa'))).sourceId, 1);
    expect(revision.records.map(({ id, header_offset, sequence_offset, end_offset, length }) => [id, header_offset, sequence_offset, end_offset, length])).toEqual([
      ['', 0, 0, 11, 6],
      ['b', 11, 14, 17, 2],
    ]);
    expect(revision.records[0]!.sha256).toBe(sha256('#c\nACGT\nAC\n'));
    expect(decoder.decode((await data.buildRunInput([revision.revisionId])).bytes)).toBe(text);
    const withoutB = await data.reviseDataset(revision.revisionId, [1]);
    expect(decoder.decode((await data.buildRunInput([withoutB.revisionId])).bytes)).toBe('#c\nACGT\nAC\n');
    const withoutFirst = await data.reviseDataset(revision.revisionId, [0]);
    const input = await data.buildRunInput([withoutFirst.revisionId]);
    expect(decoder.decode(input.bytes)).toBe('>b\nGG\n');
    expect(input.records).toEqual([{ id: 'b', length: 2 }]);
  });

  it('keeps a first record without a defline at the start of a combined input only', async () => {
    const { data } = service();
    const index = async (name: string, text: string) =>
      (await data.indexSource((await data.addSource(new File([text], name))).sourceId, 1)).revisionId;
    const plain = await index('plain.fa', 'ACGT\n>b\nGG\n');
    const fasta = await index('x.fa', '>x\nAC\n');
    const comments = await index('notes.fa', '# no records\n\n');
    expect(decoder.decode((await data.buildRunInput([plain, fasta])).bytes)).toBe('ACGT\n>b\nGG\n>x\nAC\n');
    // After a source without records (blank and comment lines), it is still the first record.
    expect((await data.buildRunInput([comments, plain])).records).toEqual([
      { id: '', length: 4 },
      { id: 'b', length: 2 },
    ]);
    // After another record, the engine would read its residues as part of that record.
    await expect(data.buildRunInput([fasta, plain])).rejects.toThrow(
      'the first record of plain.fa has no defline (a line that begins with ">"), so it cannot follow other records in one search input',
    );
    await expect(data.checkInput('blastn', 'query', [fasta, plain])).rejects.toThrow('the first record of plain.fa has no defline');
  });

  it('accepts a source without records: white space, blank and comment lines', async () => {
    const { data } = service();
    const text = ' \n\n# a comment\n;another\n';
    const revision = await data.indexSource((await data.addSource(new File([text], 'empty.fa'))).sourceId, 1);
    expect(revision.records).toEqual([]);
    const input = await data.buildRunInput([revision.revisionId]);
    expect(decoder.decode(input.bytes)).toBe(text);
    expect(input.records).toEqual([]);
    expect(await data.checkInput('blastn', 'subject', [revision.revisionId])).toEqual({ ok: true, records: [] });
  });

  it('accepts a record without a defline as the first record of a table only', async () => {
    const record = (index: number, header: number, sequence: number, end: number): IndexedRecord => ({
      index,
      id: '',
      header_offset: header,
      sequence_offset: sequence,
      end_offset: end,
      length: 1,
      line_layout: { kind: 'uniform', width: 1, eol: 1 },
      residue_counts: { A: 1 },
    });
    const tables: ReadonlyArray<readonly [IndexedRecord[], boolean]> = [
      [[record(0, 0, 0, 2), record(1, 2, 5, 7)], true],
      [[record(0, 2, 5, 7), record(1, 7, 7, 12)], false],
    ];
    for (const [records, valid] of tables) {
      const { data } = service({ scanner: { scan: async () => ({ records }) } });
      const source = await data.addSource(new File(['A\n>b\nA\n>c\nA\n'], 'x.fa'));
      const indexing = data.indexSource(source.sourceId, 1);
      if (valid) await expect(indexing).resolves.toMatchObject({ records: [{ index: 0 }, { index: 1 }] });
      else await expect(indexing).rejects.toThrow('invalid record table (record 2)');
    }
  });
});

describe('DataService checkInput: the record that a refusal names by its line', () => {
  const near = (line: number) =>
    `BLAST query error: CFastaReader: Near line ${line}, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.`;
  function setup() {
    let message = '';
    const refusing: InputChecker = { check: async () => ({ ok: false, message }) };
    const { data } = service({ checker: refusing });
    const index = async (name: string, text: string) =>
      (await data.indexSource((await data.addSource(new File([text], name))).sourceId, 1)).revisionId;
    const positionOf = async (refusal: string, revisionIds: readonly string[]) => {
      message = refusal;
      const check = await data.checkInput('blastn', 'query', revisionIds);
      expect(check).toMatchObject({ ok: false, message: refusal });
      return check.ok ? undefined : check.recordPosition;
    };
    return { data, index, positionOf };
  }

  it('finds the record of a line with LF and with CR LF line ends; blank lines count', async () => {
    const { index, positionOf } = setup();
    const lf = await index('lf.fa', '>a\nACGT\n>b\nACGT\n\n>c\nACGT\n');
    expect(await positionOf(near(2), [lf])).toBe(0);
    expect(await positionOf(near(4), [lf])).toBe(1);
    expect(await positionOf(near(5), [lf])).toBe(1);
    expect(await positionOf(near(6), [lf])).toBe(2);
    expect(await positionOf(near(8), [lf])).toBeUndefined();
    const crlf = await index('crlf.fa', '>a\r\nACGT\r\n>b\r\nACGT\r\n\r\n>c\r\nACGT\r\n');
    expect(await positionOf(near(5), [crlf])).toBe(1);
    expect(await positionOf(near(6), [crlf])).toBe(2);
  });

  it("counts the lines as NCBI's line reader does where it joins two lines of the file", async () => {
    const { index, positionOf } = setup();
    // A lone CR in a file of LF line ends: "C" and ">b" are one line (3), so line 5 is ">c".
    const lone = await index('lone-cr.fa', '>a\nA\rC\n>b\nGG\n>c\nTT\n');
    expect(await positionOf(near(5), [lone])).toBe(2);
    // An LF in a file of CR line ends: "GG" and ">b" are one line (3), so line 5 is ">c".
    const cr = await index('cr.fa', '>a\rAC\nGG\r>b\rTT\r>c\rCC\r');
    expect(await positionOf(near(5), [cr])).toBe(2);
  });

  it('finds a first record without a defline, and no record for a line before the first record', async () => {
    const { index, positionOf } = setup();
    const plain = await index('plain.fa', '#c\nACGT\n>b\nAC\n');
    expect(await positionOf(near(2), [plain])).toBe(0);
    expect(await positionOf(near(4), [plain])).toBe(1);
    const seqId =
      'the first line ("AB123456") is not a defline and may be a sequence identifier that NCBI BLAST+ fetches through a data loader (from GenBank or a BLAST database), which is not supported by LOSAT Web (start the input with a \'>\' defline)';
    expect(await positionOf(seqId, [plain])).toBe(0);
    const commented = await index('commented.fa', '#c\n\n>a\nAC\n');
    expect(await positionOf(near(2), [commented])).toBeUndefined();
    expect(await positionOf(near(4), [commented])).toBe(0);
  });

  it('finds the record in a run input of several revisions, with an added newline and exclusions', async () => {
    const { data, index, positionOf } = setup();
    const x = await index('x.fa', '>x\nAC');
    const yzw = await index('yzw.fa', '>y\nGG\n>z\nTT\n>w\nCC\n');
    // The run input is ">x\nAC\n>y\nGG\n>z\nTT\n>w\nCC\n": the added newline ends line 2.
    expect(await positionOf(near(2), [x, yzw])).toBe(0);
    expect(await positionOf(near(4), [x, yzw])).toBe(1);
    const gap =
      "line 5 is a gap line ('>?'), whose residues have no bytes in the input for LOSAT Web's index to locate; gap lines are not supported by LOSAT Web";
    expect(await positionOf(gap, [x, yzw])).toBe(2);
    // Without z, line 5 is ">w", the third record of the input.
    const withoutZ = (await data.reviseDataset(yzw, [1])).revisionId;
    expect(await positionOf(near(5), [x, withoutZ])).toBe(2);
    expect(await positionOf(near(7), [x, withoutZ])).toBeUndefined();
  });

  it("numbers the lines as NCBI's line reader: every line counts, CR LF and a lone CR end one", () => {
    const starts = (text: string) => {
      const bytes = encoder.encode(text);
      const out: number[] = [];
      for (let line = 1; ; line++) {
        const start = ncbiLineStart(bytes, line);
        if (start === undefined) return out;
        out.push(start);
      }
    };
    expect(starts('')).toEqual([]);
    expect(starts('>a\nAC')).toEqual([0, 3]);
    expect(starts('\n\n>a\n')).toEqual([0, 1, 2]);
    expect(starts('>a\r\nAC\r\n\r\nGG')).toEqual([0, 4, 8, 10]);
    expect(starts('>a\rAC\r\rGG\r')).toEqual([0, 3, 6, 7]);
    // CR LF after LF line ends is still one line end; a lone CR there joins the next line.
    expect(starts('>a\nAC\r\nGG\n')).toEqual([0, 3, 7]);
    expect(starts('>a\nA\rC\nGG\n>b\n')).toEqual([0, 3, 5, 10]);
    expect(ncbiLineStart(encoder.encode('>a\n'), 0)).toBeUndefined();
    expect(lineInMessage('BLAST query error: CFastaReader: Near line 12, there is a line')).toBe(12);
    expect(lineInMessage("line 3 is a gap line ('>?')")).toBe(3);
    expect(lineInMessage('the first line ("AB123456") is not a defline')).toBe(1);
    expect(lineInMessage('record 2 ("b"): NCBI BLAST+ reads two of its lines as one')).toBeUndefined();
  });

  it('gives no position for a message that names no line', async () => {
    const { index, positionOf } = setup();
    const source = await index('a.fa', '>a\nAC\n>b\nGG\n');
    expect(await positionOf('BLAST engine error: Empty CBlastQueryVector', [source])).toBeUndefined();
    expect(await positionOf('query record 2 (b) has no residues', [source])).toBeUndefined();
  });
});

describe('DataService runs and storage status', () => {
  it('reports the storage it uses, the session bytes and the cleanup state', async () => {
    let settle!: (value: { state: 'done'; removedSessions: number }) => void;
    const { data } = service({
      fallbackReason: 'this browser has no Origin Private File System',
      cleanup: new Promise((resolve) => {
        settle = resolve;
      }),
      estimate: async () => ({ usage: 10, quota: 1000 }),
    });
    expect(await data.storageInfo()).toEqual({
      backend: 'memory',
      fallbackReason: 'this browser has no Origin Private File System',
      sessionBytes: 0,
      estimate: { usage: 10, quota: 1000 },
      cleanup: { state: 'pending' },
    });
    settle({ state: 'done', removedSessions: 2 });
    const port = await data.openRun('run-1');
    const writer = new RunOutputWriter(port);
    writer.write(6, encoder.encode('12345'));
    writer.end();
    await data.commitRun('run-1');
    const info = await data.storageInfo();
    expect(info.sessionBytes).toBe(5);
    expect(info.cleanup).toEqual({ state: 'done', removedSessions: 2 });
  });

  it('refuses a second run with the same id and a commit of an unknown run', async () => {
    const { data } = service();
    await data.openRun('run-1');
    await expect(data.openRun('run-1')).rejects.toThrow('already exists');
    await expect(data.commitRun('run-2')).rejects.toThrow('not staged');
  });

  it('reports the storage failure even if removing the failed run fails too', async () => {
    const store = new MemoryBlockStore();
    const removeAll = store.removeAll.bind(store);
    let removals = 0;
    store.removeAll = async (prefix) => {
      removals++;
      await removeAll(prefix);
      throw Object.assign(new Error('in use'), { name: 'NoModificationAllowedError' });
    };
    const { data } = service({ store });
    const port = await data.openRun('run-1');
    store.setCapacity(0);
    const writer = new RunOutputWriter(port);
    writer.write(6, encoder.encode('row'));
    writer.end();
    await expect(data.commitRun('run-1')).rejects.toThrow('Not enough temporary storage');
    expect(removals).toBeGreaterThan(0);
  });

  it('counts the HSP records as readHits reads them: blank lines and empty chunks are not records', async () => {
    const { data } = service();
    const port = await data.openRun('run-1');
    // Raw messages, as an engine that does not filter empty chunks would send them.
    const chunks = [new Uint8Array(0), encoder.encode('{"index":0}\n\n{"inde'), encoder.encode('x":1}\n  \n{"index":2}')];
    for (const bytes of chunks) port.postMessage({ type: 'chunk', stream: 1, bytes });
    port.postMessage({ type: 'end', chunks: chunks.length, bytes: chunks.reduce((sum, c) => sum + c.length, 0) });
    expect((await data.commitRun('run-1')).hitCount).toBe(3);
    expect((await data.readHits('run-1')).map((hit) => hit.index)).toEqual([0, 1, 2]);
  });

  it('commits a run once when the commit is asked for twice', async () => {
    const { data } = service();
    const writer = new RunOutputWriter(await data.openRun('run-1'));
    writer.write(6, encoder.encode('row\n'));
    const [first, second] = [data.commitRun('run-1'), data.commitRun('run-1')];
    writer.end();
    expect(await first).toEqual(await second);
    expect(decoder.decode(await data.readOutput('run-1', 6))).toBe('row\n');
  });

  it('frees the bytes of a run that runs out of storage at once', async () => {
    const { data, store } = service();
    const port = await data.openRun('run-1');
    store.setCapacity(8);
    const writer = new RunOutputWriter(port);
    writer.write(6, encoder.encode('1234'));
    writer.write(6, encoder.encode('56789'));
    writer.end();
    await expect(data.commitRun('run-1')).rejects.toThrow('Not enough temporary storage');
    expect(store.usage()).toBe(0);
  });

  it('reports the memory budget of results kept in memory with its own message', async () => {
    const store = new MemoryBlockStore({ capacityBytes: 8, fullMessage: MEMORY_FULL_MESSAGE });
    const { data } = service({ store });
    const port = await data.openRun('run-1');
    const writer = new RunOutputWriter(port);
    writer.write(0, encoder.encode('123456789'));
    writer.end();
    await expect(data.commitRun('run-1')).rejects.toMatchObject({ name: 'StorageFullError', message: MEMORY_FULL_MESSAGE });
  });
});
