// The data layer of one working session (plan §3.1, §5.4, §5.6). It runs inside the Data
// worker (src/infra/data-worker/data-worker.ts), and directly in unit tests. It keeps:
// - the sources: File references, read with `File.slice`, never copied into storage;
// - the record tables (dataset revisions), built with the index scan port, and the reading
//   of a record's residues from its source for extraction (domain/sequence-layout.ts);
// - the run registry: staged runs, whose outputs arrive over a MessagePort from the
//   engine, and committed runs. Staged and committed runs use the same blocks; the
//   registry, not a move of files, says which runs are committed (plan §2.3).
import { intervalText, withinRecord, type Interval } from '../../domain/coordinates';
import {
  includedRecords,
  normalizeExclusion,
  recordKey,
  type DatasetRecord,
  type DatasetRevision,
  type FastaParserKind,
  type IndexedRecord,
  type RecordKey,
} from '../../domain/dataset';
import { hspTable, type HspTable } from '../../domain/hsp-table';
import { ForwardReader, readLimit, readerKind, readStart, residueCounts, type ReaderKind } from '../../domain/sequence-layout';
import type { OutputFormat } from '../../domain/output-format';
import type { InputRole, ProgramId } from '../../domain/programs';
import type {
  CleanupState,
  DataGateway,
  RecordResidues,
  ResultSetRef,
  RunInput,
  RunInputDescription,
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
import { recordOfMessageLine } from './message-line';

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

/** A run input without its SHA-256, and the [start, end) byte range of each of its records. */
interface RunInputBytes extends Omit<RunInput, 'sha256'> {
  readonly spans: ReadonlyArray<readonly [number, number]>;
}

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
  /** Bytes delivered by the port so far, and the `stagedBytes` calls that wait for more. */
  received: number;
  readonly waiters: Array<{ readonly atLeast: number; readonly resolve: (received: number) => void }>;
}

interface CommittedRun {
  readonly state: 'committed';
  readonly token: string;
  readonly lengths: Readonly<Lengths>;
  /** Where each HSP record's line is in stream 1, found on the first `readHspRecords`. */
  hitLines?: Promise<HitLines> | undefined;
}

/**
 * The lines of stream 1 that hold HSP records (the lines that are not blank, as `readHits`
 * reads them): line `k` is the bytes [starts[k], ends[k]). `lineOf` maps an HSP index to its
 * line where the lines are not in the order of the index (the adapter writes them in order).
 */
interface HitLines {
  readonly starts: Float64Array;
  readonly ends: Float64Array;
  lineOf?: ReadonlyMap<number, number>;
}

/** The parts of a run input in order: each revision's source and included records. */
interface RunInputPart {
  readonly revision: DatasetRevision;
  readonly file: File;
  readonly included: readonly DatasetRecord[];
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

  /**
   * The parts of the run input of the revisions, in order: the included records of each revision
   * in turn. `buildRunInput` joins them and `readResidues` finds a record's position in them, so
   * the positions are those of the engine's records (`q_idx`, `s_idx`).
   */
  private runInputParts(revisionIds: readonly string[]): RunInputPart[] {
    if (revisionIds.length === 0) throw new Error('a run input needs at least one dataset revision');
    return revisionIds.map((revisionId) => {
      const revision = this.revision(revisionId);
      return { revision, file: this.source(revision.sourceId), included: includedRecords(revision) };
    });
  }

  /** The bytes and record keys of a run input, without its SHA-256, and where its records lie. */
  private async runInputBytes(revisionIds: readonly string[]): Promise<RunInputBytes> {
    const parts: Uint8Array[] = [];
    const records: RecordKey[] = [];
    const spans: Array<readonly [number, number]> = [];
    let length = 0;
    for (const { revision, file, included } of this.runInputParts(revisionIds)) {
      const whole = revision.excluded.length === 0;
      const bytes = whole
        ? await readRange(file, 0, file.size)
        : await readRanges(
            file,
            included.map((record) => [record.header_offset, record.end_offset] as const),
          );
      const previous = parts[parts.length - 1];
      if (bytes.length > 0 && previous !== undefined && previous[previous.length - 1] !== LF) {
        parts.push(NEWLINE);
        length += NEWLINE.length;
      }
      const first = included[0];
      if (first !== undefined && first.header_offset === first.sequence_offset && spans.length > 0) {
        // The engine would read the residues of a first record without a defline as part of
        // the record before it: only whole original records may form a run input.
        throw new Error(
          `the first record of ${file.name} has no defline (a line that begins with ">"), so it cannot follow ` +
            'other records in one search input; search it separately or give it a defline',
        );
      }
      let at = length;
      // One by one: spreading a large table into push() overflows the call stack.
      for (const record of included) {
        records.push(recordKey(record));
        const start = whole ? length + record.header_offset : at;
        const end = start + record.end_offset - record.header_offset;
        spans.push([start, end]);
        at = end;
      }
      if (bytes.length > 0) parts.push(bytes);
      length += bytes.length;
    }
    return { bytes: concatBytes(parts), records: Object.freeze(records), spans };
  }

  async checkInput(program: ProgramId, role: InputRole, revisionIds: readonly string[]): Promise<InputCheck> {
    const { bytes, spans } = await this.runInputBytes(revisionIds);
    const check = await this.deps.checker.check(program, role, bytes);
    if (check.ok || check.recordPosition !== undefined) return check;
    // NCBI's reader names the line it refuses; the record that holds the line can be excluded.
    const recordPosition = recordOfMessageLine(check.message, bytes, spans);
    return recordPosition === undefined ? check : { ...check, recordPosition };
  }

  async previewSource(sourceId: string, maxBytes: number): Promise<Uint8Array> {
    const file = this.source(sourceId);
    return readRange(file, 0, Math.min(file.size, Math.max(0, maxBytes)));
  }

  async readResidues(revisionIds: readonly string[], position: number, intervals: readonly Interval[]): Promise<RecordResidues> {
    const { revision, file, record } = this.recordAt(revisionIds, position);
    const kind = readerKind(revision.parser);
    for (const each of intervals) {
      if (!withinRecord(each, record.length)) {
        throw new RangeError(`${each.from}-${each.to} is not an interval of record ${position + 1} ("${record.id}", ${record.length} letters)`);
      }
    }
    const residues: Uint8Array[] = [];
    for (const each of intervals) residues.push(await this.readInterval(file, record, kind, each));
    return {
      origin: {
        sourceId: revision.sourceId,
        sourceName: file.name,
        revisionId: revision.revisionId,
        recordIndex: record.index,
        id: record.id,
        length: record.length,
        sha256: record.sha256,
      },
      residues,
    };
  }

  async checkRecord(revisionIds: readonly string[], position: number): Promise<void> {
    const { revision, file, record } = this.recordAt(revisionIds, position);
    const kind = readerKind(revision.parser);
    // In parts of at most `readChunkBytes` residues, read as `readResidues` reads an interval.
    const counts: Record<string, number> = {};
    let found = 0;
    for (let from = 1; from <= record.length; from += this.chunkBytes) {
      const part = { from, to: Math.min(record.length, from + this.chunkBytes - 1) };
      const residues = await this.readForward(file, record, kind, part);
      found += residues.length;
      for (const [letter, n] of Object.entries(residueCounts(kind, residues))) counts[letter] = (counts[letter] ?? 0) + n;
      if (residues.length < part.to - part.from + 1) break;
    }
    const problem =
      found !== record.length
        ? `only ${found} of its ${record.length} residues were found`
        : !sameCounts(counts, record.residue_counts)
          ? 'its residues are not those counted in the record table'
          : undefined;
    if (problem !== undefined) throw sourceChanged(file, record, { from: 1, to: record.length }, problem);
  }

  /** The record at `position` of the run input of the revisions (see `runInputParts`). */
  private recordAt(revisionIds: readonly string[], position: number): RunInputPart & { readonly record: DatasetRecord } {
    if (!Number.isSafeInteger(position) || position < 0) throw new RangeError(`${position} is not a record position`);
    let rest = position;
    for (const part of this.runInputParts(revisionIds)) {
      if (rest < part.included.length) return { ...part, record: part.included[rest]! };
      rest -= part.included.length;
    }
    throw new RangeError(`the run input has ${position - rest} records, so it has no record ${position + 1}`);
  }

  /**
   * Reads the residues of an interval forward from where the record's layout places its first
   * residue, in ranges of at most `readChunkBytes`, and stops at its last residue. The bytes
   * read are bounded by the interval (`readLimit`), not by the record: a short interval of a
   * 100 Mbp record reads a few hundred bytes (or up to two checkpoint spans), not a chunk.
   */
  private async readInterval(file: File, record: DatasetRecord, kind: ReaderKind, wanted: Interval): Promise<Uint8Array> {
    const count = wanted.to - wanted.from + 1;
    const residues = await this.readForward(file, record, kind, wanted);
    const whole = wanted.from === 1 && wanted.to === record.length;
    const problem =
      residues.length !== count
        ? `only ${residues.length} of its ${count} residues were found`
        : whole && !sameCounts(residueCounts(kind, residues), record.residue_counts)
          ? 'its residues are not those counted in the record table'
          : undefined;
    if (problem !== undefined) throw sourceChanged(file, record, wanted, problem);
    return residues;
  }

  /** The residues of an interval as the source File has them now: fewer than asked where it ran out. */
  private async readForward(file: File, record: DatasetRecord, kind: ReaderKind, wanted: Interval): Promise<Uint8Array> {
    const start = readStart(record.line_layout, record.sequence_offset, wanted.from - 1);
    const reader = new ForwardReader(kind, start.skip, wanted.to - wanted.from + 1);
    const end = Math.min(readLimit(record, wanted.to - 1), file.size);
    for (let offset = start.offset; offset < end && !reader.done; offset += this.chunkBytes) {
      reader.feed(await readRange(file, offset, Math.min(end, offset + this.chunkBytes)));
    }
    return reader.residues();
  }

  async describeRunInput(revisionIds: readonly string[]): Promise<RunInputDescription> {
    const parts = this.runInputParts(revisionIds);
    const readers = new Set(parts.map(({ revision }) => revision.parser));
    if (readers.size !== 1) throw new Error('the revisions of one run input were indexed with different reader kinds');
    const id: string[] = [];
    const length: number[] = [];
    const sha256: string[] = [];
    // One by one: spreading a large table into push() overflows the call stack.
    for (const { included } of parts) {
      for (const record of included) {
        id.push(record.id);
        length.push(record.length);
        sha256.push(record.sha256);
      }
    }
    return {
      reader: parts[0]!.revision.parser,
      records: { id, length, sha256 },
      sources: parts.map(({ revision, file }) => ({
        name: file.name,
        size: file.size,
        records: revision.records.length,
        excluded: [...revision.excluded],
      })),
    };
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
      received: 0,
      waiters: [],
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

  async readOutputRange(runId: string, format: OutputFormat, start: number, end: number): Promise<Uint8Array> {
    const run = this.runs.get(runId);
    if (run?.state !== 'committed') throw new Error(`run ${runId} has no committed result`);
    if (!Number.isSafeInteger(start) || !Number.isSafeInteger(end) || start < 0 || end < start || end > run.lengths[format]) {
      throw new RangeError(`bytes [${start}, ${end}) are not in outfmt ${format} of run ${runId} (${run.lengths[format]} bytes)`);
    }
    return this.deps.store.read(blockPath(run.token, format), start, end - start);
  }

  async readHits(runId: string): Promise<readonly HspRecord[]> {
    return parseHitLines(new TextDecoder().decode(await this.readStream(runId, HITS_STREAM)));
  }

  async readHspRecords(runId: string, indices: readonly number[]): Promise<readonly HspRecord[]> {
    const run = this.runs.get(runId);
    if (run?.state !== 'committed') throw new Error(`run ${runId} has no committed result`);
    // A failed search for the lines is not kept, so the next call tries again.
    run.hitLines ??= this.findHitLines(run).catch((error: unknown) => {
      run.hitLines = undefined;
      throw error;
    });
    const lines = await run.hitLines;
    const count = lines.starts.length;
    for (const index of indices) {
      if (!Number.isSafeInteger(index) || index < 0 || index >= count) {
        throw new RangeError(`run ${runId} has no HSP ${index} (it has ${count} HSP records)`);
      }
    }
    const path = blockPath(run.token, HITS_STREAM);
    const decoder = new TextDecoder();
    const readLine = async (line: number): Promise<HspRecord> => {
      const start = lines.starts[line]!;
      return JSON.parse(decoder.decode(await this.deps.store.read(path, start, lines.ends[line]! - start))) as HspRecord;
    };
    const records: HspRecord[] = [];
    for (const index of indices) {
      // The adapter writes the records in the order of their index, so line `index` is the
      // record; where it is not, every line is read once to map the indices to lines.
      let record = await readLine(lines.lineOf?.get(index) ?? index);
      if (record.index !== index) {
        lines.lineOf ??= await this.mapHitLines(path, lines);
        const line = lines.lineOf.get(index);
        if (line === undefined) throw new RangeError(`run ${runId} has no HSP ${index}`);
        record = await readLine(line);
      }
      records.push(record);
    }
    return Object.freeze(records);
  }

  async readHitTable(runId: string): Promise<HspTable> {
    return hspTable(await this.readHits(runId));
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

  async runBlockLengths(runId: string): Promise<Readonly<Record<OutputStream, number>>> {
    const run = this.runs.get(runId);
    if (run?.state !== 'committed') throw new Error(`run ${runId} has no committed result`);
    return { ...run.lengths };
  }

  async readRunBlock(runId: string, stream: OutputStream, start: number, end: number): Promise<Uint8Array> {
    const run = this.runs.get(runId);
    if (run?.state !== 'committed') throw new Error(`run ${runId} has no committed result`);
    if (!OUTPUT_STREAMS.includes(stream)) throw new RangeError(`${String(stream)} is not a stream of a run`);
    if (!Number.isSafeInteger(start) || !Number.isSafeInteger(end) || start < 0 || end < start || end > run.lengths[stream]) {
      throw new RangeError(`bytes [${start}, ${end}) are not in stream ${stream} of run ${runId} (${run.lengths[stream]} bytes)`);
    }
    return this.deps.store.read(blockPath(run.token, stream), start, end - start);
  }

  stagedBytes(runId: string, atLeast: number): Promise<number> {
    const run = this.runs.get(runId);
    if (run?.state !== 'staged' || run.received >= atLeast || run.failure !== undefined) {
      return Promise.resolve(run?.state === 'staged' ? run.received : 0);
    }
    // A port that broke the protocol or closed delivers nothing more.
    const ended = run.receiver.finished.then(() => run.received);
    const reached = new Promise<number>((resolve) => run.waiters.push({ atLeast, resolve }));
    return Promise.race([reached, ended]);
  }

  // --- helpers --------------------------------------------------------------------------

  /** Stores one chunk of a staged run. After a failure, the run keeps nothing. */
  private append(run: StagedRun, stream: OutputStream, bytes: Uint8Array): void {
    run.received += bytes.length;
    try {
      this.store(run, stream, bytes);
    } finally {
      this.release(run);
    }
  }

  /** Resolves the `stagedBytes` calls that the delivered bytes (or a failure) satisfy. */
  private release(run: StagedRun, all = false): void {
    for (let i = run.waiters.length - 1; i >= 0; i--) {
      const waiter = run.waiters[i]!;
      if (all || run.failure !== undefined || run.received >= waiter.atLeast) {
        run.waiters.splice(i, 1);
        waiter.resolve(run.received);
      }
    }
  }

  private store(run: StagedRun, stream: OutputStream, bytes: Uint8Array): void {
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
    if (run.state === 'staged') {
      run.receiver.close();
      this.release(run, true);
    }
    await this.deps.store.removeAll(runPrefix(run.token));
  }

  /** Finds the lines of stream 1 that are not blank (as `append` counts them), reading the block in bounded ranges. */
  private async findHitLines(run: CommittedRun): Promise<HitLines> {
    const path = blockPath(run.token, HITS_STREAM);
    const length = run.lengths[HITS_STREAM];
    const starts: number[] = [];
    const ends: number[] = [];
    let lineStart = 0;
    let open = false;
    for (let offset = 0; offset < length; offset += this.chunkBytes) {
      const bytes = await this.deps.store.read(path, offset, Math.min(this.chunkBytes, length - offset));
      for (let i = 0; i < bytes.length; i++) {
        const byte = bytes[i]!;
        if (byte === LF) {
          if (open) {
            starts.push(lineStart);
            ends.push(offset + i);
          }
          open = false;
          lineStart = offset + i + 1;
        } else if (byte !== 0x20 && byte !== 0x09 && byte !== 0x0d) {
          open = true;
        }
      }
    }
    if (open) {
      starts.push(lineStart);
      ends.push(length);
    }
    return { starts: Float64Array.from(starts), ends: Float64Array.from(ends) };
  }

  /** The line of every HSP index, from parsing every line once (lines read together in bounded ranges). */
  private async mapHitLines(path: string, lines: HitLines): Promise<Map<number, number>> {
    const decoder = new TextDecoder();
    const lineOf = new Map<number, number>();
    let first = 0;
    while (first < lines.starts.length) {
      const start = lines.starts[first]!;
      let last = first;
      while (last + 1 < lines.starts.length && lines.ends[last + 1]! - start <= this.chunkBytes) last++;
      const bytes = await this.deps.store.read(path, start, lines.ends[last]! - start);
      for (let line = first; line <= last; line++) {
        const text = decoder.decode(bytes.subarray(lines.starts[line]! - start, lines.ends[line]! - start));
        lineOf.set((JSON.parse(text) as HspRecord).index, line);
      }
      first = last + 1;
    }
    return lineOf;
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

/** The refusal of a source that no longer matches its record table (`readResidues`, `checkRecord`). */
function sourceChanged(file: File, record: DatasetRecord, wanted: Interval, problem: string): Error {
  return new Error(
    `The source "${file.name}" no longer matches its record table: record ${record.index + 1} ("${record.id}"), ` +
      `read for ${intervalText(wanted)}: ${problem}. Add the file again.`,
  );
}

/** Whether residue counts agree, a count of 0 being the same as no count. */
function sameCounts(a: Readonly<Record<string, number>>, b: Readonly<Record<string, number>>): boolean {
  const keys = new Set([...Object.keys(a), ...Object.keys(b)]);
  for (const key of keys) if ((a[key] ?? 0) !== (b[key] ?? 0)) return false;
  return true;
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

/**
 * Rejects a record table whose offsets do not describe ordered records inside the source.
 * Every record has a defline before its residues (`header_offset < sequence_offset`), except
 * a first record without one, whose offsets are both 0. A table without records (an input
 * of white space, blank and comment lines) is valid.
 */
function checkRecordTable(records: readonly IndexedRecord[], size: number): void {
  let previousEnd = 0;
  records.forEach((record, i) => {
    const headerless = i === 0 && record.header_offset === 0 && record.sequence_offset === 0;
    const ordered =
      record.index === i &&
      previousEnd <= record.header_offset &&
      (headerless || record.header_offset < record.sequence_offset) &&
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
