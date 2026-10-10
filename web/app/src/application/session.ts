// Session files (plan §5.8, design §12.2, S15 items 4 and 6; the format is docs/web/session_file.md
// and domain/session-file.ts): saving the completed runs of the working session, with the tray's
// candidates and notes when chosen (S15 decision 1); opening a file in another working session
// without searching, its runs added as completed runs that are never queued and never given to
// the engine; and the explicit re-attachment of a loaded run's original FASTA (REQ-23), which only
// the extraction of sequences (and the download of a run's input FASTA) needs.
//
// Saving writes through the Writer contract (design §12.1): the container's lines and the run
// blocks, read from the Data worker in bounded ranges, go through a gzip sink into the download,
// block by block. Opening reads the file incrementally: the decompressed chunks go through the
// container reader, and each run's blocks go straight to the Data worker over the run output
// channel (ports/run-output.ts), as an engine's would, so a loaded run is in the RunStore like a
// searched one. A file is refused with a message that says what is wrong and where, and then
// nothing stays: every run staged or committed for it is deleted.
import { includedRecords, type DatasetRevision } from '../domain/dataset';
import { hspLabel } from '../domain/extraction';
import { splitOutfmt6Row, type Outfmt6Row } from '../domain/outfmt6';
import { programById, residueUnit, type InputRole } from '../domain/programs';
import type { InputSnapshot, RunRecord } from '../domain/run';
import {
  blockHeader,
  CANDIDATES_BLOCK,
  checkCandidates,
  checkManifest,
  containerEnd,
  containerHeader,
  hspRecordBounds,
  isGzip,
  MANIFEST_BLOCK,
  matchSources,
  recordsMismatch,
  runBlockName,
  SESSION_FORMAT,
  SESSION_LIMITS,
  SESSION_MIME,
  SESSION_SCHEMA,
  SESSION_STREAMS,
  SessionFileError,
  SessionFileReader,
  sessionFileName,
  type RecordIdentity,
  type SessionApp,
  type SessionCandidate,
  type SessionGroup,
  type SessionInput,
  type SessionLimits,
  type SessionManifest,
  type SessionRun,
  type SessionRunRecord,
  type SessionStream,
} from '../domain/session-file';
import type { Compression } from '../ports/compression';
import type { DataGateway, ResultSetRef, RunInput, RunStore } from '../ports/data';
import type { Downloader } from '../ports/download';
import type { HspRecord } from '../ports/engine';
import { DIAGNOSTICS_STREAM, HITS_STREAM, type OutputStream, type RunOutputMessage } from '../ports/run-output';
import type { CandidateRecord, CandidateRun, CandidateSource, CandidateTray } from './candidates';
import type { AppState, Coordinator, RunView, SessionRunInit } from './coordinator';
import { ExportWriter, writeFile } from './export-writer';
import { Store } from './store';

/** The ABI stream of each block of a run (ports/run-output.ts). */
const STREAM_OF: Readonly<Record<SessionStream, OutputStream>> = { out0: 0, out6: 6, out7: 7, hits: HITS_STREAM, diagnostics: DIAGNOSTICS_STREAM };
/** The RunRecord fields that a session file keeps (a completed run has no error). */
const RECORD_FIELDS = [
  'runtimePath',
  'threads',
  'fallbackReason',
  'engineBuild',
  'runtimeGeneration',
  'memory',
  'subjectRetained',
  'startedAt',
  'phaseTimes',
  'endedAt',
] as const satisfies readonly (keyof SessionRunRecord)[];
/** Bytes of a run block read from the Data worker at once when saving. */
const READ_BYTES = 8 * 1024 * 1024;
/** Bytes of the manifest and the candidates written to the gzip sink at once. */
const WRITE_BYTES = 1024 * 1024;
/** Bytes of a loaded run sent to the Data worker per message. */
const SEND_BYTES = 1024 * 1024;
/** Messages of a loaded run sent and not yet stored by the Data worker, at most (16 MiB). */
const IN_FLIGHT_MESSAGES = 16;
/** HSP records of a run read at once to restore candidates; outfmt 6 rows read together while their span is this small. */
const RECORD_BATCH = 1000;
const ROW_SPAN_BYTES = 4 * 1024 * 1024;

export interface SessionDeps {
  readonly coordinator: Pick<Coordinator, 'state' | 'addSessionRuns' | 'attach'>;
  readonly tray: Pick<CandidateTray, 'state' | 'restore'>;
  readonly data: Pick<
    DataGateway,
    | 'openRun'
    | 'commitRun'
    | 'deleteRun'
    | 'stagedBytes'
    | 'checkHspRecords'
    | 'readHspRecords'
    | 'readOutputRange'
    | 'runBlockLengths'
    | 'readRunBlock'
    | 'describeRunInput'
    | 'addSource'
    | 'indexSource'
    | 'reviseDataset'
    | 'releaseSources'
    | 'buildRunInput'
  >;
  readonly compression: Compression;
  readonly downloader: Pick<Downloader, 'open'>;
  /** This LOSAT Web's version and build, written into the manifest. */
  readonly app: SessionApp;
  readonly now: () => number;
  /** New run IDs (and group IDs) for loaded runs; the file's are never used. */
  readonly newRunId: () => string;
  /** Lower the sizes in tests. */
  readonly limits?: SessionLimits;
  readonly readBytes?: number;
  readonly sendBytes?: number;
}

export interface SessionState {
  readonly busy?: 'saving' | 'loading';
  /** What the last save or open did, or why it was refused (English, for the screen). */
  readonly message?: { readonly text: string; readonly error: boolean };
  /** Re-attachments being checked, and the last refusal of each, by `attachKey(runId, role)`. */
  readonly attaching: ReadonlyMap<string, { readonly busy: boolean; readonly message?: string }>;
}

export type SaveResult =
  | { readonly ok: true; readonly fileName: string; readonly runs: number; readonly candidates: number; readonly bytes: number }
  | { readonly ok: false; readonly message: string };
export type LoadResult = { readonly ok: true; readonly runs: number; readonly candidates: number } | { readonly ok: false; readonly message: string };
export type AttachResult = { readonly ok: true } | { readonly ok: false; readonly message: string };

export const attachKey = (runId: string, role: InputRole): string => `${runId}/${role}`;

/** A loaded run once its blocks are committed in the Data worker. */
interface CommittedRun {
  readonly runId: string;
  readonly result: ResultSetRef;
}

/** What the tray keeps of a candidate's HSP record: not its aligned rows, which an export reads again. */
type CandidateHsp = Pick<HspRecord, 'index' | 'q_idx' | 's_idx' | 'rank' | 'q_start' | 'q_end' | 's_start' | 's_end' | 'query_frame' | 'subject_frame'> & {
  readonly out6: readonly [number, number];
};

/** A tray entry of a candidate of the file, made before its run joins the working session: all but the run. */
interface PendingCandidate {
  /** 1-based position of the run in the file. */
  readonly run: number;
  readonly source: Omit<CandidateSource, 'run'>;
  readonly note: string;
  readonly addedAt: number;
}

export class Session {
  readonly state = new Store<SessionState>({ attaching: new Map() });
  private readonly limits: SessionLimits;
  /** The Data worker's sources of each attached original, by `attachKey`: released when another replaces it. */
  private readonly attachedSources = new Map<string, readonly string[]>();

  constructor(private readonly deps: SessionDeps) {
    this.limits = deps.limits ?? SESSION_LIMITS;
  }

  /** The coordinator's state: the runs, the loaded ones with their origin and their attached originals. */
  get runs(): Store<AppState> {
    return this.deps.coordinator.state;
  }

  // --- saving ----------------------------------------------------------------------------------------

  /**
   * Saves the completed runs, in run order, and the candidates of those runs with their notes
   * (unless `includeCandidates` is false), as one gzip session file. Queued, running, cancelled
   * and failed runs are never saved. A failure saves nothing and says why.
   */
  async save(options: { readonly includeCandidates?: boolean } = {}): Promise<SaveResult> {
    if (this.state.get().busy !== undefined) return this.refuse('A session file is being saved or opened; wait until that ends.');
    const views = this.deps.coordinator.state.get().runs.filter((run) => run.status === 'completed');
    if (views.length === 0) return this.refuse('There are no completed runs to save. Queued, running, cancelled and failed runs are never saved.');
    const include = options.includeCandidates ?? true;
    this.set({ busy: 'saving', message: undefined });
    try {
      const groups = new Map<string, number>();
      const runs: SessionRun[] = [];
      for (const view of views) runs.push(await this.sessionRun(view, groups));
      const positions = new Map(views.map((view, k) => [view.snapshot.runId, k + 1]));
      const candidates: SessionCandidate[] = include
        ? this.deps.tray.state
            .get()
            .candidates.filter((candidate) => positions.has(candidate.run.runId))
            .map((candidate) => ({
              run: positions.get(candidate.run.runId)!,
              index: candidate.index,
              qIdx: candidate.id.qIdx,
              rank: candidate.id.rank,
              note: candidate.note,
              addedAt: candidate.addedAt,
            }))
        : [];
      // The checks of opening, so that a saved file always opens.
      const manifest = checkManifest(
        { format: SESSION_FORMAT, schema: SESSION_SCHEMA, app: this.deps.app, savedAt: this.deps.now(), candidates: include, runs },
        this.limits,
      );
      const manifestBytes = this.encodeBlock(MANIFEST_BLOCK, manifest, this.limits.manifestBytes);
      let candidatesBytes: Uint8Array | undefined;
      if (include) {
        checkCandidates({ candidates }, manifest, this.limits);
        candidatesBytes = this.encodeBlock(CANDIDATES_BLOCK, { candidates }, this.limits.candidatesBytes);
      }
      const fileName = sessionFileName(manifest.savedAt);
      const bytes = await writeFile(this.deps.downloader, fileName, SESSION_MIME, (out) =>
        this.writeContainer(out, views, manifest, manifestBytes, candidatesBytes),
      );
      const what = `${count(views.length, 'run')}${include ? ` and ${count(candidates.length, 'candidate')} with their notes` : ''}`;
      this.set({ busy: undefined, message: { text: `Saved ${fileName}: ${what}.`, error: false } });
      return { ok: true, fileName, runs: views.length, candidates: candidates.length, bytes };
    } catch (error) {
      const reason = error instanceof SessionFileError ? `a session file cannot hold it (${error.message})` : errorMessage(error);
      return this.fail(`The session could not be saved: ${reason}`);
    }
  }

  /** What the file records of a completed run. */
  private async sessionRun(view: RunView, groups: Map<string, number>): Promise<SessionRun> {
    const { snapshot } = view;
    if (view.result === undefined) throw new Error(`run ${snapshot.number} has no result`);
    const lengths = await this.deps.data.runBlockLengths(snapshot.runId);
    let group: SessionGroup | undefined;
    if (snapshot.group !== undefined) {
      const index = groups.get(snapshot.group.groupId) ?? groups.size + 1;
      groups.set(snapshot.group.groupId, index);
      group = { index, position: snapshot.group.position, size: snapshot.group.size };
    }
    const record: Partial<Record<keyof SessionRunRecord, unknown>> = {};
    for (const field of RECORD_FIELDS) if (view.record[field] !== undefined) record[field] = view.record[field];
    return {
      number: snapshot.number,
      ...(snapshot.title === undefined ? {} : { title: snapshot.title }),
      program: snapshot.program,
      argv: snapshot.argv,
      requestedThreads: snapshot.requestedThreads,
      ...(group === undefined ? {} : { group }),
      queuedAt: snapshot.queuedAt,
      record: record as SessionRunRecord,
      hitCount: view.result.hitCount,
      blocks: Object.fromEntries(SESSION_STREAMS.map((stream) => [stream, lengths[STREAM_OF[stream]]])) as Record<SessionStream, number>,
      query: await this.inputIdentity(view, 'query'),
      subject: await this.inputIdentity(view, 'subject'),
    };
  }

  /** The identity of a run's input: the Data worker's record of it, or for a loaded run, the file's. */
  private async inputIdentity(view: RunView, role: InputRole): Promise<SessionInput> {
    if (view.fromSession !== undefined) return view.fromSession.inputs[role];
    const input = view.snapshot[role];
    const description = await this.deps.data.describeRunInput(input.revisionIds);
    if (description.records.id.length !== input.records.length) {
      throw new Error(`the ${role} record table of run ${view.snapshot.number} no longer matches its input`);
    }
    return {
      name: input.name,
      sha256: input.sha256,
      length: input.bytes?.length ?? 0,
      reader: description.reader,
      records: description.records,
      sources: description.sources,
    };
  }

  private encodeBlock(name: string, value: unknown, max: number): Uint8Array {
    const bytes = new TextEncoder().encode(JSON.stringify(value));
    if (bytes.length > max) throw new SessionFileError(`its ${name} would be ${bytes.length} bytes, more than the ${max} that a session file may have`);
    return bytes;
  }

  private async writeContainer(
    out: ExportWriter,
    views: readonly RunView[],
    manifest: SessionManifest,
    manifestBytes: Uint8Array,
    candidatesBytes: Uint8Array | undefined,
  ): Promise<void> {
    const gzip = this.deps.compression.gzip((bytes) => out.bytes(bytes));
    try {
      const writer = new ExportWriter(gzip);
      await writer.text(containerHeader());
      await writeJsonBlock(writer, MANIFEST_BLOCK, manifestBytes);
      const readBytes = this.deps.readBytes ?? READ_BYTES;
      for (const [k, view] of views.entries()) {
        for (const stream of SESSION_STREAMS) {
          const length = manifest.runs[k]!.blocks[stream];
          await writer.text(blockHeader(runBlockName(k + 1, stream), length));
          for (let start = 0; start < length; start += readBytes) {
            const end = Math.min(length, start + readBytes);
            await writer.bytes(await this.deps.data.readRunBlock(view.snapshot.runId, STREAM_OF[stream], start, end));
          }
          await writer.text('\n');
        }
      }
      if (candidatesBytes !== undefined) await writeJsonBlock(writer, CANDIDATES_BLOCK, candidatesBytes);
      await writer.text(containerEnd());
      await writer.flush();
      await gzip.close();
    } catch (error) {
      gzip.abort();
      throw error;
    }
  }

  // --- opening ---------------------------------------------------------------------------------------

  /**
   * Opens a session file: its runs become completed runs of this working session, numbered after
   * its runs, and its candidates (if it has them) are added to the tray with their notes. Nothing
   * is searched, and the engine is never asked. A file that is refused loads nothing.
   */
  async load(file: File): Promise<LoadResult> {
    if (this.state.get().busy !== undefined) return this.refuse('A session file is being saved or opened; wait until that ends.');
    this.set({ busy: 'loading', message: undefined });
    const runIds: string[] = [];
    try {
      const { runs, candidates } = await this.read(file, runIds);
      const numbers = runs.length === 1 ? `run ${runs[0]!.snapshot.number}` : `runs ${runs[0]!.snapshot.number} to ${runs.at(-1)!.snapshot.number}`;
      const restored = candidates === undefined ? '' : ` and ${count(candidates, 'candidate')} with their notes`;
      const text = `Opened ${file.name}: ${count(runs.length, 'run')} (${numbers} here)${restored}. Nothing was searched again.`;
      this.set({ busy: undefined, message: { text, error: false } });
      return { ok: true, runs: runs.length, candidates: candidates ?? 0 };
    } catch (error) {
      for (const runId of runIds) await this.deps.data.deleteRun(runId).catch(() => undefined);
      const reason = error instanceof SessionFileError ? error.message : `It could not be read: ${errorMessage(error)}`;
      return this.fail(`${file.name} was not opened, and nothing was loaded. ${reason}`);
    }
  }

  /** Reads the file into the Data worker, checks it, and only then adds its runs and candidates. */
  private async read(file: File, runIds: string[]): Promise<{ runs: readonly RunView[]; candidates?: number }> {
    if (!isGzip(new Uint8Array(await file.slice(0, 2).arrayBuffer()))) {
      throw new SessionFileError('This is not a LOSAT Web session file: it is not gzip data (a session file is saved as a .losat-session.gz file).');
    }
    const reader = new SessionFileReader(this.limits);
    const chunks = this.deps.compression.gunzip(file)[Symbol.asyncIterator]();
    let manifest: SessionManifest | undefined;
    let candidates: readonly SessionCandidate[] = [];
    const committed: CommittedRun[] = [];
    let sender: RunSender | undefined;
    try {
      for (;;) {
        let next: IteratorResult<Uint8Array>;
        try {
          next = await chunks.next();
        } catch (error) {
          throw new SessionFileError(`The session file is damaged: its gzip data is not valid or is cut short${detail(error)}.`);
        }
        if (next.done === true) break;
        for (const event of reader.push(next.value)) {
          switch (event.type) {
            case 'manifest':
              manifest = event.manifest;
              break;
            case 'run-start': {
              const runId = this.deps.newRunId();
              runIds.push(runId);
              sender = new RunSender(runId, await this.deps.data.openRun(runId), this.deps.data, this.deps.sendBytes ?? SEND_BYTES);
              break;
            }
            case 'data':
              await sender!.write(STREAM_OF[event.stream], event.bytes);
              break;
            case 'run-end':
              committed.push(await this.commit(sender!, manifest!.runs[event.run - 1]!, event.run));
              sender = undefined;
              break;
            case 'candidates':
              candidates = event.candidates;
              break;
            case 'end':
              break;
          }
        }
      }
      reader.finish();
    } finally {
      await chunks.return?.().catch(() => undefined);
    }
    // Everything that can fail comes before the runs join this working session, so that a refused
    // file leaves nothing in it (the caller deletes the runs from the Data worker).
    const groupIds = new Map<number, string>();
    const inits = manifest!.runs.map((run, k) => this.runInit(run, committed[k]!, file.name, manifest!, groupIds));
    const pending = await this.candidateEntries(manifest!, committed, candidates);
    // Every check has passed: the runs join, then their candidates (the tray takes HSPs of completed runs only).
    const runs = this.deps.coordinator.addSessionRuns(inits);
    if (!manifest!.candidates) return { runs };
    const restored = this.deps.tray.restore(
      pending.map(({ run, source, note, addedAt }) => ({ source: { ...source, run: candidateRun(runs[run - 1]!) }, note, addedAt })),
    );
    return { runs, candidates: restored.ok ? restored.added : 0 };
  }

  /** Ends a run's output, commits it, and checks what the Data worker stored against the manifest. */
  private async commit(sender: RunSender, run: SessionRun, position: number): Promise<CommittedRun> {
    sender.end();
    const result = await this.deps.data.commitRun(sender.runId);
    const where = `run ${position} in the file (run ${run.number} there)`;
    if (result.hitCount !== run.hitCount) {
      throw new SessionFileError(`The session file is damaged: ${where} has ${result.hitCount} HSP records, but the manifest gives ${run.hitCount}.`);
    }
    // The records as the file holds them, before anything coerces them: every later reader (the
    // results' table, the restored candidates, the exports) takes them as an engine wrote them.
    let problem: string | undefined;
    try {
      problem = await this.deps.data.checkHspRecords(sender.runId, hspRecordBounds(run));
    } catch (error) {
      throw new SessionFileError(`The session file is damaged: the HSP records of ${where} cannot be read${detail(error)}.`);
    }
    if (problem !== undefined) throw new SessionFileError(`The session file is damaged: in the HSP records of ${where}, ${problem}.`);
    return { runId: sender.runId, result };
  }

  /**
   * The tray's entries of the file's candidates, in the file's order, made from the loaded runs
   * (only the reference comes from the file): the HSP record at a candidate's index must be the HSP
   * that it names (its query record and rank), and its outfmt 6 row is read from the run's
   * outfmt 6. The records are read RECORD_BATCH at a time, as the exports read them, and only
   * the fields that the tray keeps stay (no aligned rows).
   */
  private async candidateEntries(
    manifest: SessionManifest,
    committed: readonly CommittedRun[],
    candidates: readonly SessionCandidate[],
  ): Promise<readonly PendingCandidate[]> {
    const entries = new Array<PendingCandidate>(candidates.length);
    const byRun = new Map<number, number[]>();
    candidates.forEach((candidate, i) => {
      const list = byRun.get(candidate.run);
      if (list === undefined) byRun.set(candidate.run, [i]);
      else list.push(i);
    });
    for (const [position, items] of byRun) {
      const run = manifest.runs[position - 1]!;
      const { runId } = committed[position - 1]!;
      const recordOf = candidateRecords(run);
      for (let first = 0; first < items.length; first += RECORD_BATCH) {
        const batch = items.slice(first, first + RECORD_BATCH);
        const read = await this.deps.data.readHspRecords(runId, batch.map((i) => candidates[i]!.index));
        const hsps = batch.map((i, j) => candidateHsp(candidates[i]!, i, position, read[j]!));
        const rows = await this.readRows(runId, hsps, run.blocks.out6);
        batch.forEach((i, j) => {
          const hsp = hsps[j]!;
          const { note, addedAt } = candidates[i]!;
          const source: Omit<CandidateSource, 'run'> = {
            id: { runId, qIdx: hsp.q_idx, rank: hsp.rank },
            index: hsp.index,
            query: recordOf('query', hsp.q_idx),
            subject: recordOf('subject', hsp.s_idx),
            coordinates: {
              q_start: hsp.q_start,
              q_end: hsp.q_end,
              s_start: hsp.s_start,
              s_end: hsp.s_end,
              query_frame: hsp.query_frame,
              subject_frame: hsp.subject_frame,
            },
            row: rows[j]!,
          };
          entries[i] = { run: position, source, note, addedAt };
        });
      }
    }
    return entries;
  }

  /** The outfmt 6 rows of HSP records of a run, read in spans of nearby rows. */
  private async readRows(runId: string, records: readonly { readonly out6: readonly [number, number] }[], out6Length: number): Promise<Outfmt6Row[]> {
    const order = records.map((_, j) => j).sort((a, b) => records[a]!.out6[0] - records[b]!.out6[0]);
    const rows = new Array<Outfmt6Row>(records.length);
    const decoder = new TextDecoder();
    let first = 0;
    while (first < order.length) {
      const start = records[order[first]!]!.out6[0];
      let last = first;
      let end = records[order[first]!]!.out6[1];
      while (last + 1 < order.length && records[order[last + 1]!]!.out6[1] - start <= ROW_SPAN_BYTES) {
        last++;
        end = Math.max(end, records[order[last]!]!.out6[1]);
      }
      const bytes = await this.deps.data.readOutputRange(runId, 6, start, Math.min(end, out6Length));
      for (let k = first; k <= last; k++) {
        const [from, to] = records[order[k]!]!.out6;
        try {
          rows[order[k]!] = splitOutfmt6Row(decoder.decode(bytes.subarray(from - start, to - start)));
        } catch (error) {
          throw new SessionFileError(`The session file is damaged: the outfmt 6 row of a candidate is not a row${detail(error)}.`);
        }
      }
      first = last + 1;
    }
    return rows;
  }

  private runInit(run: SessionRun, committed: CommittedRun, fileName: string, manifest: SessionManifest, groupIds: Map<number, string>): SessionRunInit {
    let group;
    if (run.group !== undefined) {
      const groupId = groupIds.get(run.group.index) ?? this.deps.newRunId();
      groupIds.set(run.group.index, groupId);
      group = Object.freeze({ groupId, position: run.group.position, size: run.group.size });
    }
    return {
      snapshot: Object.freeze({
        runId: committed.runId,
        program: run.program,
        ...(run.title === undefined ? {} : { title: run.title }),
        argv: Object.freeze([...run.argv]),
        query: inputSnapshot(run.query),
        subject: inputSnapshot(run.subject),
        requestedThreads: run.requestedThreads,
        queuedAt: run.queuedAt,
        ...(group === undefined ? {} : { group }),
      }),
      record: { ...run.record } as RunRecord,
      result: committed.result,
      fromSession: {
        fileName,
        number: run.number,
        savedAt: manifest.savedAt,
        app: manifest.app,
        ...(run.record.engineBuild === undefined ? {} : { engineBuild: run.record.engineBuild }),
        inputs: { query: run.query, subject: run.subject },
      },
    };
  }

  // --- re-attaching the original FASTA ---------------------------------------------------------------

  /**
   * Attaches the original FASTA of one role of a run loaded from a session file (REQ-23), only on
   * this explicit choice: the files (several, in any order, for a joined input) are indexed with
   * the reader kind that the file recorded and matched to the recorded sources by their records
   * (domain/session-file.ts `matchSources`); the run input is then made in the recorded order
   * with the recorded exclusions, and attached only if its records and its SHA-256 are those
   * that the run searched. Otherwise nothing changes, the message names the first record or file
   * that differs, and the attempt's sources and record tables are released from the Data worker,
   * as are those of an original that a new attachment replaces.
   */
  async attach(runId: string, role: InputRole, files: readonly File[]): Promise<AttachResult> {
    const view = this.deps.coordinator.state.get().runs.find((run) => run.snapshot.runId === runId);
    if (view?.fromSession === undefined) return { ok: false, message: 'Only a run loaded from a session file takes its original FASTA again.' };
    const key = attachKey(runId, role);
    if (this.state.get().attaching.get(key)?.busy === true) return { ok: false, message: 'These files are being checked; wait until that ends.' };
    this.setAttaching(key, { busy: true });
    const number = view.snapshot.number;
    const sourceIds: string[] = [];
    try {
      const saved = view.fromSession.inputs[role];
      if (files.length !== saved.sources.length) {
        const names = saved.sources.map((source) => JSON.stringify(source.name)).join(', ');
        const chosen = `${files.length} ${files.length === 1 ? 'was' : 'were'} chosen`;
        throw new Error(
          saved.sources.length === 1
            ? `choose one file (the run's ${role} was ${names}); ${chosen}`
            : `choose the ${saved.sources.length} files that the run joined (${names}), together and in any order; ${chosen}`,
        );
      }
      const indexed: DatasetRevision[] = [];
      for (const file of files) {
        const source = await this.deps.data.addSource(file);
        sourceIds.push(source.sourceId);
        try {
          indexed.push(await this.deps.data.indexSource(source.sourceId, saved.reader));
        } catch (error) {
          throw new Error(`${JSON.stringify(file.name)} could not be read as the run read its ${role}: ${errorMessage(error)}`);
        }
      }
      const match = matchSources(saved, indexed.map((revision, j) => ({ name: files[j]!.name, records: revision.records })));
      if (!match.ok) throw new Error(match.message);
      const revisionIds: string[] = [];
      const chosen: RecordIdentity[] = [];
      for (const [i, j] of match.files.entries()) {
        const excluded = saved.sources[i]!.excluded;
        const revised = excluded.length === 0 ? indexed[j]! : await this.deps.data.reviseDataset(indexed[j]!.revisionId, excluded);
        revisionIds.push(revised.revisionId);
        for (const record of includedRecords(revised)) chosen.push({ id: record.id, length: record.length, sha256: record.sha256 });
      }
      const mismatch = recordsMismatch(saved, chosen);
      if (mismatch !== undefined) throw new Error(mismatch);
      const input = await this.deps.data.buildRunInput(revisionIds);
      if (input.sha256 !== saved.sha256) {
        throw new Error(
          `the chosen files have the run's records, but the input that they make (SHA-256 ${input.sha256}) is not the one that run ${number} searched ` +
            `(SHA-256 ${saved.sha256}): a file has other lines before or between its records`,
        );
      }
      this.deps.coordinator.attach(runId, role, { revisionIds, fileNames: match.files.map((j) => files[j]!.name) });
      const replaced = this.attachedSources.get(key);
      this.attachedSources.set(key, sourceIds);
      if (replaced !== undefined) await this.deps.data.releaseSources(replaced).catch(() => undefined);
      this.setAttaching(key, undefined);
      return { ok: true };
    } catch (error) {
      await this.deps.data.releaseSources(sourceIds).catch(() => undefined);
      const message = `The ${role} FASTA was not attached to run ${number}: ${errorMessage(error)}.`;
      this.setAttaching(key, { busy: false, message });
      return { ok: false, message };
    }
  }

  /**
   * The run input (the bytes that the engine searched) of a loaded run's attached original FASTA,
   * or undefined while none is attached: for the download of a run's input FASTA.
   */
  async attachedInput(runId: string, role: InputRole): Promise<RunInput | undefined> {
    const attached = this.deps.coordinator.state.get().runs.find((run) => run.snapshot.runId === runId)?.attached?.[role];
    return attached === undefined ? undefined : this.deps.data.buildRunInput(attached.revisionIds);
  }

  // --- state -----------------------------------------------------------------------------------------

  private setAttaching(key: string, value: { readonly busy: boolean; readonly message?: string } | undefined): void {
    const attaching = new Map(this.state.get().attaching);
    if (value === undefined) attaching.delete(key);
    else attaching.set(key, value);
    this.set({ attaching });
  }

  private refuse(message: string): { readonly ok: false; readonly message: string } {
    this.set({ message: { text: message, error: true } });
    return { ok: false, message };
  }

  private fail(message: string): { readonly ok: false; readonly message: string } {
    this.set({ busy: undefined, message: { text: message, error: true } });
    return { ok: false, message };
  }

  private set(change: { [K in keyof SessionState]?: SessionState[K] | undefined }): void {
    const next = { ...this.state.get() } as Record<string, unknown>;
    for (const [key, value] of Object.entries(change)) {
      if (value === undefined) delete next[key];
      else next[key] = value;
    }
    this.state.set(next as unknown as SessionState);
  }
}

/**
 * Sends a loaded run's blocks to the Data worker over its run output port (ports/run-output.ts),
 * in messages of about `sendBytes`, and waits while more than IN_FLIGHT_MESSAGES of them are not
 * yet stored, so that a large file does not pile up in the worker's queue. Once the worker cannot
 * store the run (the storage ran out), that wait rejects with the reason, and the load stops
 * there instead of decompressing and sending the rest of the file (code review L4).
 */
class RunSender {
  private buffer: Uint8Array;
  private filled = 0;
  private stream: OutputStream | undefined;
  private chunks = 0;
  private sent = 0;
  private stored = 0;

  constructor(
    readonly runId: string,
    private readonly port: MessagePort,
    private readonly data: Pick<RunStore, 'stagedBytes'>,
    private readonly sendBytes: number,
  ) {
    this.buffer = new Uint8Array(sendBytes);
  }

  async write(stream: OutputStream, bytes: Uint8Array): Promise<void> {
    if (stream !== this.stream) this.post();
    this.stream = stream;
    let at = 0;
    while (at < bytes.length) {
      const take = Math.min(this.buffer.length - this.filled, bytes.length - at);
      this.buffer.set(bytes.subarray(at, at + take), this.filled);
      this.filled += take;
      at += take;
      if (this.filled === this.buffer.length) {
        this.post();
        const window = IN_FLIGHT_MESSAGES * this.sendBytes;
        if (this.sent - this.stored > window) this.stored = await this.data.stagedBytes(this.runId, this.sent - window / 2);
      }
    }
  }

  /** Sends what is left and the `end` with the totals. */
  end(): void {
    this.post();
    const message: RunOutputMessage = { type: 'end', chunks: this.chunks, bytes: this.sent };
    this.port.postMessage(message);
  }

  private post(): void {
    if (this.filled === 0 || this.stream === undefined) return;
    const full = this.filled === this.buffer.length;
    const bytes = full ? this.buffer : this.buffer.slice(0, this.filled);
    const message: RunOutputMessage = { type: 'chunk', stream: this.stream, bytes };
    // Count before the transfer detaches the bytes.
    this.chunks++;
    this.sent += bytes.length;
    this.port.postMessage(message, [bytes.buffer]);
    this.filled = 0;
    // A full buffer was transferred; a part was copied, and the buffer is used again.
    if (full) this.buffer = new Uint8Array(this.sendBytes);
  }
}

/** The InputSnapshot of a loaded run: the file's identity of the input, without bytes or revisions. */
function inputSnapshot(input: SessionInput): InputSnapshot {
  const { id, length } = input.records;
  const records = new Array<{ readonly id: string; readonly length: number }>(id.length);
  for (let k = 0; k < id.length; k++) records[k] = { id: id[k]!, length: length[k]! };
  return Object.freeze({ name: input.name, sha256: input.sha256, revisionIds: Object.freeze([]), records: Object.freeze(records) });
}

/**
 * The HSP record of candidate `i` of the file (run `position` there), read from the loaded run,
 * with only the fields that the tray keeps; refused unless it is the HSP that the candidate names
 * and has an outfmt 6 row. The run's records were checked when it was committed (`checkHspRecords`).
 */
function candidateHsp(candidate: SessionCandidate, i: number, position: number, record: HspRecord): CandidateHsp {
  if (record.q_idx !== candidate.qIdx || record.rank !== candidate.rank) {
    throw new SessionFileError(
      `The session file's candidates block is not valid: candidates[${i}] is HSP ${hspLabel(candidate.qIdx, candidate.rank)} at index ` +
        `${candidate.index} of run ${position} in the file, but the HSP record there is HSP ${hspLabel(record.q_idx, record.rank)}.`,
    );
  }
  if (record.out6 === null) {
    throw new SessionFileError(`The session file is damaged: the HSP record of candidates[${i}] (run ${position} in the file) has no outfmt 6 row.`);
  }
  const { index, q_idx, s_idx, rank, q_start, q_end, s_start, s_end, query_frame, subject_frame, out6 } = record;
  return { index, q_idx, s_idx, rank, q_start, q_end, s_start, s_end, query_frame, subject_frame, out6: [out6[0], out6[1]] };
}

/** The tray's records of a loaded run's inputs, made once each from the file's record tables. */
function candidateRecords(run: SessionRun): (role: InputRole, position: number) => CandidateRecord {
  const program = programById(run.program);
  const made = new Map<string, CandidateRecord>();
  return (role, position) => {
    const key = `${role}/${position}`;
    let record = made.get(key);
    if (record === undefined) {
      const kind = role === 'query' ? program.query : program.subject;
      record = { position, id: run[role].records.id[position]!, length: run[role].records.length[position]!, kind, unit: residueUnit(kind) };
      made.set(key, record);
    }
    return record;
  };
}

/** How the tray names a loaded run. */
function candidateRun(view: RunView): CandidateRun {
  const { snapshot } = view;
  return { runId: snapshot.runId, number: snapshot.number, ...(snapshot.title === undefined ? {} : { title: snapshot.title }), program: snapshot.program };
}

async function writeJsonBlock(writer: ExportWriter, name: string, bytes: Uint8Array): Promise<void> {
  await writer.text(blockHeader(name, bytes.length));
  for (let start = 0; start < bytes.length; start += WRITE_BYTES) await writer.bytes(bytes.subarray(start, start + WRITE_BYTES));
  await writer.text('\n');
}

const count = (n: number, one: string) => `${n} ${n === 1 ? one : `${one}s`}`;

function detail(error: unknown): string {
  // The browser's own text, joined cleanly: its final full stop would double the sentence's ("…stream.).").
  const message = (error instanceof Error ? error.message : String(error)).replace(/[.\s]+$/, '');
  return message === '' ? '' : ` (${message})`;
}

function errorMessage(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
