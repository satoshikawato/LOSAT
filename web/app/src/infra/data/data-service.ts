// The data layer of one working session (plan §3.1, §5.4, §5.6). It runs inside the Data
// worker (src/infra/data-worker/data-worker.ts), and directly in unit tests. It keeps:
// - the sources: File references, read with `File.slice`, never copied into storage;
// - the record tables (dataset revisions), built with the index scan port;
//   both until `dropSource` (the application's SourceRefs says when no one uses a source);
// - the run registry: staged runs, whose outputs arrive over a MessagePort from the
//   engine, and committed runs. Staged and committed runs use the same blocks; the
//   registry, not a move of files, says which runs are committed (plan §2.3).
import { RunOutputReceiver } from 'losat';
import {
  CUSTOM_OUTPUT,
  OUTPUT_FORMATS,
  type FastaReader,
  type IndexedRecord,
  type OutputPlace,
  type RecordKey,
  type RunOutput,
} from 'losat/domain';
import {
  DIAGNOSTICS_STREAM,
  HITS_STREAM,
  type FileInflater,
  type HspRecord,
  type InputCheck,
  type InputChecker,
  type OutputStream,
  type RecordScanner,
} from 'losat/ports';
import { includedRecords, normalizeExclusion, recordKey, type DatasetRevision } from '../../domain/dataset';
import { fastaParser, programById } from '../../domain/programs';
import { hspTable, type HspTable } from '../../domain/hsp-table';
import type {
  CleanupState,
  DataGateway,
  KeptRun,
  ResultSetRef,
  RunInput,
  SessionFileRef,
  SourceRef,
  StorageInfo,
  StoredResultsCheck,
  StoredToolFile,
  ToolOutputCommand,
  ToolOutputReply,
  ToolOutputsSummary,
} from '../../ports/data';
import { concatBytes } from '../bytes';
import { asStorageFull, isStorageFull, isStorageRemoved, type BlockStore, type BlockWriter } from './block-store';
import { inflateFile, inflateToWriter } from './inflate';
import { hasEntryBytes } from './zip-entries';

export interface DataServiceDeps {
  readonly store: BlockStore;
  readonly scanner: RecordScanner;
  /** The engine's inflation of a file named `.gz`, `.bz2` or `.zst`, for addSource. */
  readonly inflater: FileInflater;
  /** The engine's reading of an input (`register`), for checkInput. */
  readonly checker: InputChecker;
  /** Lower-case hex SHA-256. */
  readonly digest: (bytes: Uint8Array) => Promise<string>;
  /** A new random token; tokens name blocks, so they never contain input names. */
  readonly newToken: () => string;
  readonly fallbackReason?: string;
  /** The error behind fallbackReason, "<name>: <message>", for developers (session.ts). */
  readonly fallbackError?: string;
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
  [CUSTOM_OUTPUT]: 'outcustom',
  [HITS_STREAM]: 'hits',
  [DIAGNOSTICS_STREAM]: 'diagnostics',
};
const LF = 0x0a;
const NEWLINE = new Uint8Array([LF]);

/** The bytes of each stream of a run: those of its outputs (`openRun`), the HSP records and the diagnostics. */
type Lengths = Partial<Record<OutputStream, number>>;

/** What `checkStoredResults` found for a run. */
type StoredOutcome = 'kept' | 'removed' | 'unreadable';

/** How many reads `checkStoredResults` makes at a time. */
export const STORED_CHECK_READS = 16;

interface StagedRun {
  readonly state: 'staged';
  readonly token: string;
  readonly writers: ReadonlyMap<OutputStream, BlockWriter>;
  readonly receiver: RunOutputReceiver;
  readonly lengths: Lengths;
  /** The outputs that the engine said failed, with where their errors are (its discards, `clear`): `ResultSetRef.failed`. */
  readonly failed: Map<RunOutput, OutputPlace>;
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

interface StagedToolFile {
  readonly key: string;
  readonly name: string;
  readonly blockName: string;
  readonly writer: BlockWriter;
  length: number;
  closed: boolean;
}

interface StagedToolRun {
  readonly state: 'staged';
  readonly token: string;
  readonly port: MessagePort;
  readonly stdout: BlockWriter;
  readonly stderr: BlockWriter;
  stdoutSize: number;
  stderrSize: number;
  readonly files: Map<string, StagedToolFile>;
  readonly fileOrder: string[];
  failure: Error | undefined;
}

interface CommittedToolFile extends StoredToolFile {
  readonly blockName: string;
}

interface CommittedToolRun {
  readonly state: 'committed';
  readonly token: string;
  readonly stdoutSize: number;
  readonly stderrSize: number;
  readonly files: ReadonlyMap<string, CommittedToolFile>;
}

const TOOL_FILE_KEY = /^file:(\d+)$/;

export class DataService implements DataGateway {
  private readonly sources = new Map<string, File>();
  private readonly revisions = new Map<string, DatasetRevision>();
  private readonly runs = new Map<string, StagedRun | CommittedRun>();
  private readonly toolRuns = new Map<string, StagedToolRun | CommittedToolRun>();
  private readonly sessionFiles = new Map<string, { readonly prefix: string; readonly path: string }>();
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
    file = await inflateFile(file, this.deps.inflater);
    const sourceId = this.deps.newToken();
    this.sources.set(sourceId, file);
    return { sourceId, name: file.name, size: file.size };
  }

  async indexSource(sourceId: string, reader: FastaReader): Promise<DatasetRevision> {
    const file = this.source(sourceId);
    const parser = fastaParser(programById(reader.program), reader.role);
    const { records, reader_error: stop } = await this.deps.scanner.scan(reader, readChunks(file, this.chunkBytes));
    checkRecordTable(records, file.size, stop?.offset);
    const hashes = await this.hashRecords(file, records);
    return this.keep({
      revisionId: this.deps.newToken(),
      sourceId,
      parser,
      records: Object.freeze(records.map((record, i) => Object.freeze({ ...record, sha256: hashes[i]! }))),
      excluded: Object.freeze([]),
      ...(stop === undefined ? {} : { readerError: Object.freeze({ offset: stop.offset, defline: stop.defline, message: stop.message }) }),
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
   * The bytes and record keys of a run input, without its SHA-256. A source whose reading ends with
   * the reader's error (`DatasetRevision.readerError`) keeps the bytes from the error on after its
   * included records, so that the run meets the error as the command line meets it in the file, and
   * ends the input: the engine reads no record after it (line WQA's D10, 2026-11-19).
   */
  private async runInputBytes(revisionIds: readonly string[]): Promise<Omit<RunInput, 'sha256'>> {
    if (revisionIds.length === 0) throw new Error('a run input needs at least one dataset revision');
    const parts: Uint8Array[] = [];
    const records: RecordKey[] = [];
    for (const revisionId of revisionIds) {
      const revision = this.revision(revisionId);
      const file = this.source(revision.sourceId);
      const included = includedRecords(revision);
      const stop = revision.readerError?.offset;
      const bytes =
        revision.excluded.length === 0
          ? await readRange(file, 0, file.size)
          : await readRanges(file, [
              ...included.map((record) => [record.header_offset, record.end_offset] as const),
              ...(stop === undefined ? [] : [[stop, file.size] as const]),
            ]);
      const previous = parts[parts.length - 1];
      if (bytes.length > 0 && previous !== undefined && previous[previous.length - 1] !== LF) {
        parts.push(NEWLINE);
      }
      if (bytes.length > 0) parts.push(bytes);
      // One by one: spreading a large table into push() overflows the call stack.
      for (const record of included) records.push(recordKey(record));
      if (stop !== undefined) break;
    }
    return { bytes: concatBytes(parts), records: Object.freeze(records) };
  }

  async checkInput(reader: FastaReader, revisionIds: readonly string[]): Promise<InputCheck> {
    const { bytes } = await this.runInputBytes(revisionIds);
    return this.deps.checker.check(reader, bytes);
  }

  async dropSource(sourceId: string): Promise<void> {
    this.sources.delete(sourceId);
    for (const [revisionId, revision] of this.revisions) {
      if (revision.sourceId === sourceId) this.revisions.delete(revisionId);
    }
  }

  async previewSource(sourceId: string, maxBytes: number): Promise<Uint8Array> {
    const file = this.source(sourceId);
    return readRange(file, 0, Math.min(file.size, Math.max(0, maxBytes)));
  }

  // --- runs -----------------------------------------------------------------------------

  async openRun(runId: string, outputs: readonly RunOutput[] = OUTPUT_FORMATS): Promise<MessagePort> {
    if (this.runs.has(runId)) throw new Error(`run ${runId} already exists`);
    const token = this.deps.newToken();
    const writers = new Map<OutputStream, BlockWriter>();
    const streams: readonly OutputStream[] = [...new Set([...OUTPUT_FORMATS, ...outputs]), HITS_STREAM, DIAGNOSTICS_STREAM];
    try {
      for (const stream of streams) writers.set(stream, await this.deps.store.create(blockPath(token, stream)));
    } catch (error) {
      await this.deps.store.removeAll(runPrefix(token)).catch(() => undefined);
      throw isStorageFull(error) ? asStorageFull(error) : error;
    }
    const channel = new MessageChannel();
    const lengths: Lengths = Object.fromEntries(streams.map((stream) => [stream, 0]));
    const run: StagedRun = {
      state: 'staged',
      token,
      writers,
      receiver: new RunOutputReceiver(
        channel.port2,
        (stream, bytes) => this.append(run, stream, bytes),
        (stream, error) => this.clear(run, stream, error),
      ),
      lengths,
      failed: new Map(),
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
      (result.state === 'broken'
        ? new Error(`The output of the run arrived incomplete: ${result.detail}`)
        : result.state === 'failed'
          ? new Error('The engine failed the run')
          : undefined);
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
    return this.resultOf(runId, run);
  }

  /** The result of a run that was committed: the sizes of its outputs, its HSP records and the outputs that it could not write. */
  private resultOf(runId: string, run: StagedRun): ResultSetRef {
    const hitCount = run.hitLines + (run.hitLineOpen ? 1 : 0);
    const custom = run.lengths[CUSTOM_OUTPUT];
    const byteLengths = { 0: run.lengths[0]!, 6: run.lengths[6]!, 7: run.lengths[7]!, ...(custom === undefined ? {} : { [CUSTOM_OUTPUT]: custom }) };
    const failed = [...run.failed].sort(([a], [b]) => a - b).map(([output, error]) => ({ output, error }));
    return { runId, byteLengths, hitCount, failed };
  }

  async discardRun(runId: string): Promise<void> {
    const run = this.runs.get(runId);
    if (run?.state !== 'staged') return;
    await this.drop(runId, run);
  }

  async keepFailedRun(runId: string): Promise<KeptRun | undefined> {
    const run = this.runs.get(runId);
    if (run?.state !== 'staged') return undefined;
    // The chunks of a failed run, then its `fail`: the outputs and stream 3 are whole once `fail` arrived.
    const result = await run.receiver.finished;
    if (result.state === 'failed' && run.failure === undefined && this.runs.get(runId) === run) {
      try {
        for (const writer of run.writers.values()) await writer.seal();
        const length = run.lengths[DIAGNOSTICS_STREAM]!;
        const diagnostics =
          length === 0 ? '' : new TextDecoder().decode(await this.deps.store.read(blockPath(run.token, DIAGNOSTICS_STREAM), 0, length));
        if (this.runs.get(runId) === run) {
          if (!(await this.keptBytes(run))) {
            // Nothing was written to a kept output (an empty file, an archive of empty entries): the run keeps only
            // its warnings, as a failed run did before version 7.
            await this.drop(runId, run).catch(() => undefined);
            return { diagnostics };
          }
          // The commit of a completed run: its outputs are read, and exported, as a completed run's are.
          this.runs.set(runId, { state: 'committed', token: run.token, lengths: Object.freeze({ ...run.lengths }) });
          return { diagnostics, result: this.resultOf(runId, run) };
        }
      } catch {
        // The outputs could not be kept: the run is dropped, as a failed run was.
      }
    }
    // As the callers of `discardRun` do, a failure to remove the blocks does not hide why the run failed.
    const current = this.runs.get(runId);
    if (current !== undefined) await this.drop(runId, current).catch(() => undefined);
    return undefined;
  }

  /**
   * Whether a failed run wrote something to an output that it keeps: bytes of its own, or, for an archive (outfmt 13
   * and 14), bytes in its entries, not in the archive's own records (the empty master file of a run that failed
   * before writing it is an archive of 22 + 46 + 30 bytes of records and no data).
   */
  private async keptBytes(run: StagedRun): Promise<boolean> {
    for (const output of [0, 6, 7, CUSTOM_OUTPUT] as const) {
      const length = run.lengths[output] ?? 0;
      if (length === 0 || run.failed.has(output)) continue;
      const read = (start: number, end: number) => this.deps.store.read(blockPath(run.token, output), start, end - start);
      if (await hasEntryBytes(read, length)) return true;
    }
    return false;
  }

  async readOutput(runId: string, output: RunOutput): Promise<Uint8Array> {
    return this.readStream(runId, output);
  }

  async readOutputRange(runId: string, format: RunOutput, start: number, end: number): Promise<Uint8Array> {
    const length = this.committed(runId, format);
    const run = this.runs.get(runId)!;
    if (!Number.isSafeInteger(start) || !Number.isSafeInteger(end) || start < 0 || end < start || end > length) {
      throw new RangeError(`bytes [${start}, ${end}) are not in outfmt ${format} of run ${runId} (${length} bytes)`);
    }
    return this.deps.store.read(blockPath(run.token, format), start, end - start);
  }

  async readHits(runId: string): Promise<readonly HspRecord[]> {
    return parseHitLines(new TextDecoder().decode(await this.readStream(runId, HITS_STREAM)));
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
    const tool = this.toolRuns.get(runId);
    if (tool !== undefined) await this.dropTool(runId, tool);
  }

  // --- tool runs and OPFS-backed session files (WW3b) -----------------------------------

  async openToolRun(runId: string): Promise<MessagePort> {
    if (this.runs.has(runId) || this.toolRuns.has(runId)) throw new Error(`run ${runId} already exists`);
    const token = this.deps.newToken();
    let stdout: BlockWriter;
    let stderr: BlockWriter;
    try {
      stdout = await this.deps.store.create(toolBlockPath(token, 'stdout'));
      stderr = await this.deps.store.create(toolBlockPath(token, 'stderr'));
    } catch (error) {
      await this.deps.store.removeAll(toolRunPrefix(token)).catch(() => undefined);
      throw isStorageFull(error) ? asStorageFull(error) : error;
    }
    const channel = new MessageChannel();
    const run: StagedToolRun = {
      state: 'staged',
      token,
      port: channel.port2,
      stdout,
      stderr,
      stdoutSize: 0,
      stderrSize: 0,
      files: new Map(),
      fileOrder: [],
      failure: undefined,
    };
    this.toolRuns.set(runId, run);
    let chain: Promise<void> = Promise.resolve();
    channel.port2.onmessage = (event: MessageEvent<ToolOutputCommand>) => {
      const cmd = event.data;
      chain = chain.then(() => this.handleToolCommand(runId, run, cmd));
    };
    return channel.port1;
  }

  async readToolOutputRange(runId: string, key: string, start: number, end: number): Promise<Uint8Array> {
    const { path, size } = this.committedToolOutput(runId, key);
    if (!Number.isSafeInteger(start) || !Number.isSafeInteger(end) || start < 0 || end < start || end > size) {
      throw new RangeError(`bytes [${start}, ${end}) are not in output ${key} of run ${runId} (${size} bytes)`);
    }
    return this.deps.store.read(path, start, end - start);
  }

  async getToolOutputFile(runId: string, key: string, fileName: string): Promise<File> {
    const { path } = this.committedToolOutput(runId, key);
    return this.deps.store.readFile(path, fileName, 'application/octet-stream');
  }

  async inflateToSessionFile(file: File): Promise<SessionFileRef> {
    if (!(file instanceof Blob)) throw new TypeError('an input must be a File');
    const token = this.deps.newToken();
    const prefix = sessionFilePrefix(token);
    const path = `${prefix}/data`;
    const inflated = await inflateToWriter(file, this.deps.inflater, () => this.deps.store.create(path));
    if (inflated === undefined) {
      return { file, name: file.name, size: file.size };
    }
    this.sessionFiles.set(token, { prefix, path });
    const stored = await this.deps.store.readFile(path, inflated.plainName, 'text/plain');
    return { token, file: stored, name: inflated.plainName, size: inflated.size, inflatedFrom: file.name };
  }

  async dropSessionFile(token: string): Promise<void> {
    const entry = this.sessionFiles.get(token);
    if (entry === undefined) return;
    this.sessionFiles.delete(token);
    await this.deps.store.removeAll(entry.prefix).catch(() => undefined);
  }

  async checkStoredResults(): Promise<StoredResultsCheck> {
    type Entry =
      | { readonly kind: 'search'; readonly runId: string; readonly run: StagedRun | CommittedRun }
      | { readonly kind: 'tool'; readonly runId: string; readonly run: StagedToolRun | CommittedToolRun };
    const entries: Entry[] = [
      ...[...this.runs].map(([runId, run]): Entry => ({ kind: 'search', runId, run })),
      ...[...this.toolRuns].map(([runId, run]): Entry => ({ kind: 'tool', runId, run })),
    ];
    const outcomes: StoredOutcome[] = entries.map(() => 'kept');
    let next = 0;
    // At most STORED_CHECK_READS reads at a time, so that hundreds of runs do not flood the browser with requests.
    const reader = async () => {
      for (let i = next++; i < entries.length; i = next++) {
        const entry = entries[i]!;
        outcomes[i] =
          entry.kind === 'search'
            ? await this.checkStored(entry.runId, entry.run)
            : await this.checkStoredTool(entry.runId, entry.run);
      }
    };
    await Promise.all(Array.from({ length: Math.min(STORED_CHECK_READS, entries.length) }, reader));
    const of = (outcome: StoredOutcome) => entries.filter((_, i) => outcomes[i] === outcome).map((entry) => entry.runId);
    return { removed: of('removed'), unreadable: of('unreadable') };
  }

  /**
   * Opens the HSP records of a committed run, the first block that opening its results asks for, and reads its
   * first byte (none of an empty block). The other blocks of the run are in the same directory, which the
   * application removes only as a whole.
   */
  private async checkStored(runId: string, run: StagedRun | CommittedRun): Promise<StoredOutcome> {
    if (run.state !== 'committed') return 'kept';
    let failure: { readonly error: unknown } | undefined;
    try {
      await this.deps.store.read(blockPath(run.token, HITS_STREAM), 0, Math.min(1, run.lengths[HITS_STREAM]!));
    } catch (error) {
      failure = { error };
    }
    // A run deleted while it was checked is not reported.
    if (failure === undefined || this.runs.get(runId) !== run) return 'kept';
    return isStorageRemoved(failure.error) ? 'removed' : 'unreadable';
  }

  private async checkStoredTool(runId: string, run: StagedToolRun | CommittedToolRun): Promise<StoredOutcome> {
    if (run.state !== 'committed') return 'kept';
    let failure: { readonly error: unknown } | undefined;
    try {
      await this.deps.store.read(toolBlockPath(run.token, 'stdout'), 0, Math.min(1, run.stdoutSize));
    } catch (error) {
      failure = { error };
    }
    if (failure === undefined || this.toolRuns.get(runId) !== run) return 'kept';
    return isStorageRemoved(failure.error) ? 'removed' : 'unreadable';
  }

  private async handleToolCommand(runId: string, run: StagedToolRun, cmd: ToolOutputCommand): Promise<void> {
    const reply = (msg: ToolOutputReply, close = false) => {
      try {
        run.port.postMessage(msg);
      } finally {
        if (close) run.port.close();
      }
    };
    const removeToolBlocks = () => this.deps.store.removeAll(toolRunPrefix(run.token)).catch(() => undefined);
    const requireOpenFile = (key: string) => {
      const file = run.files.get(key);
      if (file === undefined || file.closed) throw new Error(`tool output ${key} is not open`);
      return file;
    };
    const fail = (error: unknown, close = false) => {
      const normalized = isStorageFull(error) ? asStorageFull(error) : toError(error);
      run.failure ??= normalized;
      void removeToolBlocks();
      reply({ id: cmd.id, ok: false, error: { name: normalized.name, message: normalized.message } }, close);
    };
    if (cmd.op === 'discard') {
      if (this.toolRuns.get(runId) === run) {
        this.toolRuns.delete(runId);
        await removeToolBlocks();
      }
      reply({ id: cmd.id, ok: true }, true);
      return;
    }
    if (this.toolRuns.get(runId) !== run) {
      reply({ id: cmd.id, ok: false, error: { name: 'Error', message: `run ${runId} was discarded` } }, cmd.op === 'commit');
      return;
    }
    if (run.failure !== undefined) {
      if (cmd.op === 'commit') {
        this.toolRuns.delete(runId);
        await removeToolBlocks();
      }
      reply({ id: cmd.id, ok: false, error: { name: run.failure.name, message: run.failure.message } }, cmd.op === 'commit');
      return;
    }
    try {
      if (cmd.op === 'stdout') {
        run.stdout.append(cmd.bytes);
        run.stdoutSize += cmd.bytes.length;
        reply({ id: cmd.id, ok: true });
      } else if (cmd.op === 'stderr') {
        run.stderr.append(cmd.bytes);
        run.stderrSize += cmd.bytes.length;
        reply({ id: cmd.id, ok: true });
      } else if (cmd.op === 'open') {
        const match = TOOL_FILE_KEY.exec(cmd.key);
        if (match === null) throw new Error(`invalid tool output key ${cmd.key}`);
        if (run.files.has(cmd.key)) throw new Error(`tool output ${cmd.key} is already open`);
        const blockName = `f${match[1]}`;
        const writer = await this.deps.store.create(toolBlockPath(run.token, blockName));
        run.files.set(cmd.key, { key: cmd.key, name: cmd.name, blockName, writer, length: 0, closed: false });
        run.fileOrder.push(cmd.key);
        reply({ id: cmd.id, ok: true });
      } else if (cmd.op === 'write') {
        const file = requireOpenFile(cmd.key);
        file.writer.writeAt(cmd.offset, cmd.bytes);
        if (cmd.bytes.length > 0) file.length = Math.max(file.length, cmd.offset + cmd.bytes.length);
        reply({ id: cmd.id, ok: true });
      } else if (cmd.op === 'truncate') {
        const file = requireOpenFile(cmd.key);
        file.writer.truncate(cmd.length);
        file.length = cmd.length;
        reply({ id: cmd.id, ok: true });
      } else if (cmd.op === 'close') {
        const file = run.files.get(cmd.key);
        if (file !== undefined && !file.closed) {
          file.length = await file.writer.seal();
          file.closed = true;
        }
        reply({ id: cmd.id, ok: true });
      } else if (cmd.op === 'commit') {
        const stdoutSize = await run.stdout.seal();
        const stderrSize = await run.stderr.seal();
        const committedFiles = new Map<string, CommittedToolFile>();
        const summaryFiles: StoredToolFile[] = [];
        for (const key of run.fileOrder) {
          const file = run.files.get(key)!;
          if (!file.closed) {
            file.length = await file.writer.seal();
            file.closed = true;
          }
          committedFiles.set(key, { key: file.key, name: file.name, size: file.length, blockName: file.blockName });
          summaryFiles.push({ key: file.key, name: file.name, size: file.length });
        }
        this.toolRuns.set(runId, {
          state: 'committed',
          token: run.token,
          stdoutSize,
          stderrSize,
          files: committedFiles,
        });
        const summary: ToolOutputsSummary = { stdoutSize, stderrSize, files: summaryFiles };
        reply({ id: cmd.id, ok: true, summary }, true);
      }
    } catch (error) {
      if (cmd.op === 'commit') this.toolRuns.delete(runId);
      fail(error, cmd.op === 'commit');
    }
  }

  private committedToolOutput(runId: string, key: string): { readonly path: string; readonly size: number } {
    const run = this.toolRuns.get(runId);
    if (run?.state !== 'committed') throw new Error(`tool run ${runId} has no committed output`);
    if (key === 'stdout') return { path: toolBlockPath(run.token, 'stdout'), size: run.stdoutSize };
    if (key === 'stderr') return { path: toolBlockPath(run.token, 'stderr'), size: run.stderrSize };
    const file = run.files.get(key);
    if (file === undefined) throw new Error(`tool run ${runId} has no output ${key}`);
    return { path: toolBlockPath(run.token, file.blockName), size: file.size };
  }

  private async dropTool(runId: string, run: StagedToolRun | CommittedToolRun): Promise<void> {
    this.toolRuns.delete(runId);
    if (run.state === 'staged') run.port.close();
    await this.deps.store.removeAll(toolRunPrefix(run.token));
  }

  async storageInfo(): Promise<StorageInfo> {
    const estimate = await this.deps.estimate?.().catch(() => undefined);
    return {
      backend: this.deps.store.backend,
      ...(this.deps.fallbackReason === undefined ? {} : { fallbackReason: this.deps.fallbackReason }),
      ...(this.deps.fallbackError === undefined ? {} : { fallbackError: this.deps.fallbackError }),
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
      const writer = run.writers.get(stream);
      // The engine writes the outputs of the run's argv (`openRun`); another output is not this run's.
      if (writer === undefined) throw new Error(`The engine wrote output stream ${stream}, which this run does not write`);
      writer.append(bytes);
    } catch (error) {
      run.failure = isStorageFull(error) ? asStorageFull(error) : toError(error);
      // Give the space back at once; the commit reports the failure.
      void this.deps.store.removeAll(runPrefix(run.token)).catch(() => undefined);
      return;
    }
    run.lengths[stream] = run.lengths[stream]! + bytes.length;
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

  /**
   * Records that one output of a staged run failed, and where its error is, as the engine said (docs/web/abi_v2.md
   * §5: it failed alone, or its search failed in a run of two searches), and drops what the run received of it: the
   * output stays empty. After a failure, the run keeps nothing.
   */
  private clear(run: StagedRun, stream: RunOutput, error: OutputPlace): void {
    if (run.failure !== undefined) return;
    try {
      const writer = run.writers.get(stream);
      if (writer === undefined) throw new Error(`The engine discarded output stream ${stream}, which this run does not write`);
      writer.clear();
    } catch (error) {
      run.failure = toError(error);
      void this.deps.store.removeAll(runPrefix(run.token)).catch(() => undefined);
      return;
    }
    run.lengths[stream] = 0;
    run.failed.set(stream, error);
  }

  private async readStream(runId: string, stream: OutputStream): Promise<Uint8Array> {
    const length = this.committed(runId, stream);
    return this.deps.store.read(blockPath(this.runs.get(runId)!.token, stream), 0, length);
  }

  /** The bytes of a stream of a committed run; throws for a run that is not committed, and for an output that it does not write. */
  private committed(runId: string, stream: OutputStream): number {
    const run = this.runs.get(runId);
    if (run?.state !== 'committed') throw new Error(`run ${runId} has no committed result`);
    const length = run.lengths[stream];
    if (length === undefined) throw new Error(`run ${runId} does not write output ${stream}`);
    return length;
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

function toolRunPrefix(token: string): string {
  return `tools/${token}`;
}

function toolBlockPath(token: string, name: string): string {
  return `${toolRunPrefix(token)}/${name}`;
}

function sessionFilePrefix(token: string): string {
  return `tool-inputs/${token}`;
}

/**
 * Rejects a record table whose offsets do not describe ordered records inside the source, and the
 * offset of a reader's error (`stop`) that is not after them inside the source. Only a first record
 * without a defline starts its sequence at its header offset (both 0).
 */
function checkRecordTable(records: readonly IndexedRecord[], size: number, stop: number | undefined): void {
  let previousEnd = 0;
  records.forEach((record, i) => {
    const ordered =
      record.index === i &&
      previousEnd <= record.header_offset &&
      (record.header_offset < record.sequence_offset || (i === 0 && record.header_offset === 0 && record.sequence_offset === 0)) &&
      record.sequence_offset <= record.end_offset &&
      record.end_offset <= size;
    if (!ordered) throw new Error(`the index scan returned an invalid record table (record ${i + 1})`);
    previousEnd = record.end_offset;
  });
  if (stop !== undefined && !(Number.isSafeInteger(stop) && previousEnd <= stop && stop <= size)) {
    throw new Error('the index scan returned an invalid offset of the reader\'s error');
  }
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
