import { createHash } from 'node:crypto';
import { describe, expect, it } from 'vitest';
import type { IndexedRecord } from '../../src/domain/dataset';
import { sha256Hex } from '../../src/infra/browser/platform';
import { DataService, type DataServiceDeps } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeScanner } from '../../src/infra/fake/fake-fasta';
import { RunOutputWriter } from '../../src/infra/run-output/writer';
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
    const revision = await data.indexSource(source.sourceId, 0);
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
    const revision = await data.indexSource(source.sourceId, 0);
    expect(revision.records.map((r) => [r.index, r.id, r.length])).toEqual([
      [0, 'dup', 8],
      [1, 'dup', 4],
    ]);
    expect(revision.records[0]!.sha256).not.toBe(revision.records[1]!.sha256);
    expect(revision.records[0]!.residue_counts).toEqual({ a: 1, c: 1, g: 1, t: 1, A: 1, C: 1, G: 1, T: 1 });
    const second = await data.reviseDataset(revision.revisionId, [0]);
    const input = await data.buildRunInput([second.revisionId]);
    expect(decoder.decode(input.bytes)).toBe('>dup second\nGGCC\n');
  });

  it('rejects a source that the parser cannot read', async () => {
    const { data } = service();
    const source = await data.addSource(new File(['ACGT\n'], 'plain.txt'));
    await expect(data.indexSource(source.sourceId, 0)).rejects.toThrow('Expected > at record start.');
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
    await expect(data.indexSource(source.sourceId, 0)).rejects.toThrow('invalid record table');
  });

  it('makes a new revision for a new selection and refuses unknown records', async () => {
    const { data } = service();
    const source = await data.addSource(new File([THREE], 'three.fa'));
    const base = await data.indexSource(source.sourceId, 0);
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
    const revision = await data.indexSource(source.sourceId, 0);
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
    const revision = await data.reviseDataset((await data.indexSource(source.sourceId, 0)).revisionId, [1]);
    const input = await data.buildRunInput([revision.revisionId]);
    expect(decoder.decode(input.bytes)).toBe('>a one\nACGT\nAC\n>c three\nTTT');
    expect(input.records).toEqual([
      { id: 'a', length: 6 },
      { id: 'c', length: 3 },
    ]);
  });

  it('adds one newline after a source without a final newline when another follows', async () => {
    const { data } = service();
    const first = await data.indexSource((await data.addSource(new File(['>x\nAC'], 'x.fa'))).sourceId, 0);
    const second = await data.indexSource((await data.addSource(new File(['>y\nGG\n'], 'y.fa'))).sourceId, 0);
    const input = await data.buildRunInput([first.revisionId, second.revisionId, first.revisionId]);
    expect(decoder.decode(input.bytes)).toBe('>x\nAC\n>y\nGG\n>x\nAC');
    expect(input.records.map((r) => r.id)).toEqual(['x', 'y', 'x']);
    // A source whose records are all left out adds nothing, not even the newline.
    const none = await data.reviseDataset(second.revisionId, [0]);
    expect(decoder.decode((await data.buildRunInput([first.revisionId, none.revisionId])).bytes)).toBe('>x\nAC');
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
});
