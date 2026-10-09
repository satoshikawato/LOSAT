// The data layer of one working session (plan §3.1, §5.4, §5.6). It runs inside the Data
// worker (src/infra/data-worker/data-worker.ts), and directly in unit tests. It keeps:
// - the sources: File references, read with `File.slice`, never copied into storage;
// - the record tables (dataset revisions), built with the index scan port;
// - the run registry: staged runs, whose outputs arrive over a MessagePort from the
//   engine, and committed runs. Staged and committed runs use the same blocks; the
//   registry, not a move of files, says which runs are committed (plan §2.3).
import {
  includedRecords,
  normalizeExclusion,
  recordKey,
  type DatasetRevision,
  type FastaParserKind,
  type IndexedRecord,
  type RecordKey,
} from '../../domain/dataset';
import type { OutputFormat } from '../../domain/output-format';
import type { InputRole, ProgramId } from '../../domain/programs';
import type {
  CleanupState,
  DataGateway,
  ResultSetRef,
  RunInput,
  SourceRef,
  StorageInfo,
} from '../../ports/data';
import type { HspRecord } from '../../ports/engine';
import type { InputCheck, InputChecker } from '../../ports/input-check';
import { DIAGNOSTICS_STREAM, HITS_STREAM, OUTPUT_STREAMS, type OutputStream } from '../../ports/run-output';
import type { RecordScanner } from '../../ports/scan';
import { RunOutputReceiver } from '../run-output/receiver';
import { concatBytes } from '../bytes';
import { asStorageFull, isStorageFull, type BlockStore, type BlockWriter } from './block-store';

export interface DataServiceDeps {
  readonly store: BlockStore;
  readonly scanner: RecordScanner;
  /** The engine's reading of an input (`register`), for checkInput. */
  readonly checker: InputChecker;
  /** Lower-case hex SHA-256. */
  readonly digest: (bytes: Uint8Array) => Promise<string>;
  /** A new random token; tokens name blocks, so they never contain input names. */
  readonly newToken: () => string;
  readonly fallbackReason?: string;
  /** Settles when the removal of abandoned sessions has finished (session.ts). */
  readonly cleanup: Promise<CleanupState>;
  readonly estimate?: () => Promise<{ readonly usage: number; readonly quota: number } | undefined>;
  /** Largest range read from a source at once; 8 MiB by default. */
  readonly readChunkBytes?: number;
}

/** Block names of a run; tmp/<session-token>/runs/<run-token>/<name> in OPFS. */
const BLOCK_NAMES: Readonly<Record<OutputStream, string>> = {
  0: 'out0',
  6: 'out6',
  7: 'out7',
  [HITS_STREAM]: 'hits',
  [DIAGNOSTICS_STREAM]: 'diagnostics',
};
const LF = 0x0a;
const NEWLINE = new Uint8Array([LF]);

type Lengths = Record<OutputStream, number>;

interface StagedRun {
  readonly state: 'staged';
  readonly token: string;
  readonly writers: ReadonlyMap<OutputStream, BlockWriter>;
  readonly receiver: RunOutputReceiver;
  readonly lengths: Lengths;
  /** HSP records so far: the lines of stream 1 that are not blank, as `readHits` reads them. */
  hitLines: number;
  /** The current line of stream 1 has a character other than white space. */
  hitLineOpen: boolean;
  failure: Error | undefined;
  /** The commit in progress; a second commitRun waits for the same one. */
  commit: Promise<ResultSetRef> | undefined;
}

interface CommittedRun {
  readonly state: 'committed';
  readonly token: string;
  readonly lengths: Readonly<Lengths>;
}

export class DataService implements DataGateway {
  private readonly sources = new Map<string, File>();
  private readonly revisions = new Map<string, DatasetRevision>();
  private readonly runs = new Map<string, StagedRun | CommittedRun>();
  private readonly chunkBytes: number;
  private cleanup: CleanupState = { state: 'pending' };

  constructor(private readonly deps: DataServiceDeps) {
    this.chunkBytes = deps.readChunkBytes ?? 8 * 1024 * 1024;
    void deps.cleanup.then((state) => {
      this.cleanup = state;
    });
  }

  // --- sources and record tables -------------------------------------------------------

  async addSource(file: File): Promise<SourceRef> {
    if (!(file instanceof Blob)) throw new TypeError('a source must be a File');
    const sourceId = this.deps.newToken();
    this.sources.set(sourceId, file);
    return { sourceId, name: file.name, size: file.size };
  }

  async indexSource(sourceId: string, parser: FastaParserKind): Promise<DatasetRevision> {
    const file = this.source(sourceId);
    const { records } = await this.deps.scanner.scan(parser, readChunks(file, this.chunkBytes));
    checkRecordTable(records, file.size);
    const hashes = await this.hashRecords(file, records);
    return this.keep({
      revisionId: this.deps.newToken(),
      sourceId,
      parser,
      records: Object.freeze(records.map((record, i) => Object.freeze({ ...record, sha256: hashes[i]! }))),
      excluded: Object.freeze([]),
    });
  }

  async reviseDataset(revisionId: string, excluded: readonly number[]): Promise<DatasetRevision> {
    const base = this.revision(revisionId);
    return this.keep({
      ...base,
      revisionId: this.deps.newToken(),
      excluded: normalizeExclusion(excluded, base.records.length),
    });
  }

  async buildRunInput(revisionIds: readonly string[]): Promise<RunInput> {
    const { bytes, records } = await this.runInputBytes(revisionIds);
    return { bytes, sha256: await this.deps.digest(bytes), records };
  }

  /** The bytes and record keys of a run input, without its SHA-256. */
  private async runInputBytes(revisionIds: readonly string[]): Promise<Omit<RunInput, 'sha256'>> {
    if (revisionIds.length === 0) throw new Error('a run input needs at least one dataset revision');
    const parts: Uint8Array[] = [];
    const records: RecordKey[] = [];
    for (const revisionId of revisionIds) {
      const revision = this.revision(revisionId);
      const file = this.source(revision.sourceId);
      const included = includedRecords(revision);
      const bytes =
        revision.excluded.length === 0
          ? await readRange(file, 0, file.size)
          : await readRanges(
              file,
              included.map((record) => [record.header_offset, record.end_offset] as const),
            );
      const previous = parts[parts.length - 1];
      if (bytes.length > 0 && previous !== undefined && previous[previous.length - 1] !== LF) {
        parts.push(NEWLINE);
      }
      if (bytes.length > 0) parts.push(bytes);
      // One by one: spreading a large table into push() overflows the call stack.
      for (const record of included) records.push(recordKey(record));
    }
    return { bytes: concatBytes(parts), records: Object.freeze(records) };
  }

  async checkInput(program: ProgramId, role: InputRole, revisionIds: readonly string[]): Promise<InputCheck> {
    const { bytes } = await this.runInputBytes(revisionIds);
    return this.deps.checker.check(program, role, bytes);
  }

  async previewSource(sourceId: string, maxBytes: number): Promise<Uint8Array> {
    const file = this.source(sourceId);
    return readRange(file, 0, Math.min(file.size, Math.max(0, maxBytes)));
  }

  // --- runs -----------------------------------------------------------------------------

  async openRun(runId: string): Promise<MessagePort> {
    if (this.runs.has(runId)) throw new Error(`run ${runId} already exists`);
    const token = this.deps.newToken();
    const writers = new Map<OutputStream, BlockWriter>();
    try {
      for (const stream of OUTPUT_STREAMS) writers.set(stream, await this.deps.store.create(blockPath(token, stream)));
    } catch (error) {
      await this.deps.store.removeAll(runPrefix(token)).catch(() => undefined);
      throw isStorageFull(error) ? asStorageFull(error) : error;
    }
    const channel = new MessageChannel();
    const lengths: Lengths = { 0: 0, 6: 0, 7: 0, [HITS_STREAM]: 0, [DIAGNOSTICS_STREAM]: 0 };
    const run: StagedRun = {
      state: 'staged',
      token,
      writers,
      receiver: new RunOutputReceiver(channel.port2, (stream, bytes) => this.append(run, stream, bytes)),
      lengths,
      hitLines: 0,
      hitLineOpen: false,
      failure: undefined,
      commit: undefined,
    };
    this.runs.set(runId, run);
    return channel.port1;
  }

  commitRun(runId: string): Promise<ResultSetRef> {
    const run = this.runs.get(runId);
    if (run?.state !== 'staged') return Promise.reject(new Error(`run ${runId} is not staged`));
    run.commit ??= this.commit(runId, run);
    return run.commit;
  }

  private async commit(runId: string, run: StagedRun): Promise<ResultSetRef> {
    const result = await run.receiver.finished;
    if (this.runs.get(runId) !== run) throw new Error(`run ${runId} was discarded`);
    const failure =
      run.failure ??
      (result.state === 'broken' ? new Error(`The output of the run arrived incomplete: ${result.detail}`) : undefined);
    // A failure to clean up must not hide why the run failed.
    if (failure !== undefined) {
      await this.discardRun(runId).catch(() => undefined);
      throw failure;
    }
    try {
      for (const writer of run.writers.values()) await writer.seal();
    } catch (error) {
      await this.discardRun(runId).catch(() => undefined);
      throw isStorageFull(error) ? asStorageFull(error) : error;
    }
    this.runs.set(runId, { state: 'committed', token: run.token, lengths: Object.freeze({ ...run.lengths }) });
    const hitCount = run.hitLines + (run.hitLineOpen ? 1 : 0);
    return { runId, byteLengths: { 0: run.lengths[0], 6: run.lengths[6], 7: run.lengths[7] }, hitCount };
  }

  async discardRun(runId: string): Promise<void> {
    const run = this.runs.get(runId);
    if (run?.state !== 'staged') return;
    await this.drop(runId, run);
  }

  async readOutput(runId: string, format: OutputFormat): Promise<Uint8Array> {
    return this.readStream(runId, format);
  }

  async readHits(runId: string): Promise<readonly HspRecord[]> {
    return parseHitLines(new TextDecoder().decode(await this.readStream(runId, HITS_STREAM)));
  }

  async readDiagnostics(runId: string): Promise<string> {
    return new TextDecoder().decode(await this.readStream(runId, DIAGNOSTICS_STREAM));
  }

  async deleteRun(runId: string): Promise<void> {
    const run = this.runs.get(runId);
    if (run !== undefined) await this.drop(runId, run);
  }

  async storageInfo(): Promise<StorageInfo> {
    const estimate = await this.deps.estimate?.().catch(() => undefined);
    return {
      backend: this.deps.store.backend,
      ...(this.deps.fallbackReason === undefined ? {} : { fallbackReason: this.deps.fallbackReason }),
      sessionBytes: this.deps.store.usage(),
      ...(estimate === undefined ? {} : { estimate }),
      cleanup: this.cleanup,
    };
  }

  // --- helpers --------------------------------------------------------------------------

  /** Stores one chunk of a staged run. After a failure, the run keeps nothing. */
  private append(run: StagedRun, stream: OutputStream, bytes: Uint8Array): void {
    if (run.failure !== undefined || bytes.length === 0) return;
    try {
      run.writers.get(stream)!.append(bytes);
    } catch (error) {
      run.failure = isStorageFull(error) ? asStorageFull(error) : toError(error);
      // Give the space back at once; the commit reports the failure.
      void this.deps.store.removeAll(runPrefix(run.token)).catch(() => undefined);
      return;
    }
    run.lengths[stream] += bytes.length;
    if (stream === HITS_STREAM) {
      for (const byte of bytes) {
        if (byte === LF) {
          if (run.hitLineOpen) run.hitLines++;
          run.hitLineOpen = false;
        } else if (byte !== 0x20 && byte !== 0x09 && byte !== 0x0d) {
          run.hitLineOpen = true;
        }
      }
    }
  }

  private async readStream(runId: string, stream: OutputStream): Promise<Uint8Array> {
    const run = this.runs.get(runId);
    if (run?.state !== 'committed') throw new Error(`run ${runId} has no committed result`);
    return this.deps.store.read(blockPath(run.token, stream), 0, run.lengths[stream]);
  }

  private async drop(runId: string, run: StagedRun | CommittedRun): Promise<void> {
    this.runs.delete(runId);
    if (run.state === 'staged') run.receiver.close();
    await this.deps.store.removeAll(runPrefix(run.token));
  }

  private async hashRecords(file: File, records: readonly IndexedRecord[]): Promise<string[]> {
    const hashes: string[] = [];
    let first = 0;
    while (first < records.length) {
      // Read the following records together while their span fits in one read.
      const start = records[first]!.header_offset;
      let last = first;
      while (last + 1 < records.length && records[last + 1]!.end_offset - start <= this.chunkBytes) last++;
      const span = await readRange(file, start, records[last]!.end_offset);
      for (let i = first; i <= last; i++) {
        const record = records[i]!;
        hashes.push(await this.deps.digest(span.subarray(record.header_offset - start, record.end_offset - start)));
      }
      first = last + 1;
    }
    return hashes;
  }

  private keep(revision: DatasetRevision): DatasetRevision {
    const frozen = Object.freeze(revision);
    this.revisions.set(frozen.revisionId, frozen);
    return frozen;
  }

  private source(sourceId: string): File {
    const file = this.sources.get(sourceId);
    if (file === undefined) throw new Error(`unknown source ${sourceId}`);
    return file;
  }

  private revision(revisionId: string): DatasetRevision {
    const revision = this.revisions.get(revisionId);
    if (revision === undefined) throw new Error(`unknown dataset revision ${revisionId}`);
    return revision;
  }
}

/** Decodes the HSP records of ABI stream 1 (JSON Lines, docs/web/abi_v2.md §8). */
export function parseHitLines(text: string): readonly HspRecord[] {
  const records: HspRecord[] = [];
  for (const line of text.split('\n')) {
    if (line.trim() !== '') records.push(JSON.parse(line) as HspRecord);
  }
  return Object.freeze(records);
}

function runPrefix(token: string): string {
  return `runs/${token}`;
}

function blockPath(token: string, stream: OutputStream): string {
  return `${runPrefix(token)}/${BLOCK_NAMES[stream]}`;
}

/** Rejects a record table whose offsets do not describe ordered records inside the source. */
function checkRecordTable(records: readonly IndexedRecord[], size: number): void {
  let previousEnd = 0;
  records.forEach((record, i) => {
    const ordered =
      record.index === i &&
      previousEnd <= record.header_offset &&
      record.header_offset < record.sequence_offset &&
      record.sequence_offset <= record.end_offset &&
      record.end_offset <= size;
    if (!ordered) throw new Error(`the index scan returned an invalid record table (record ${i + 1})`);
    previousEnd = record.end_offset;
  });
}

async function* readChunks(file: Blob, chunkBytes: number): AsyncGenerator<Uint8Array> {
  for (let offset = 0; offset < file.size; offset += chunkBytes) {
    yield readRange(file, offset, Math.min(file.size, offset + chunkBytes));
  }
}

async function readRange(file: Blob, start: number, end: number): Promise<Uint8Array> {
  return new Uint8Array(await file.slice(start, end).arrayBuffer());
}

/** Reads ranges in order, joining ranges that touch into one read. */
async function readRanges(file: Blob, ranges: ReadonlyArray<readonly [number, number]>): Promise<Uint8Array> {
  const joined: Array<[number, number]> = [];
  for (const [start, end] of ranges) {
    const last = joined[joined.length - 1];
    if (last !== undefined && last[1] === start) last[1] = end;
    else joined.push([start, end]);
  }
  const parts: Uint8Array[] = [];
  for (const [start, end] of joined) parts.push(await readRange(file, start, end));
  return concatBytes(parts);
}

function toError(error: unknown): Error {
  return error instanceof Error ? error : new Error(String(error));
}
