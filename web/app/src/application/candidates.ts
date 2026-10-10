// The candidate tray (plan §5.7 "候補トレイは application 層の状態", design §11.4, REQ-14): HSPs
// collected from the results of several runs, each held by its identity (run, query record,
// rank) with the values that the tray shows, a note and when it was added; the user's order; a
// selection; and the two outputs made from the selected candidates, kept apart (design §11.4):
// the original residues of their records (domain/extraction.ts plans the intervals, the Data
// worker reads them from the source Files) and their gapped alignments, from the aligned rows of
// their HSP records. Candidates come only from completed runs (REQ-10): a cancelled or failed run
// keeps no results and a run that has not ended has none yet, so the tray refuses them here,
// whatever a screen offers. The tray holds plain data, never a run's tables, and each operation
// is linear in the number of candidates (sorting n log n), so that adding thousands of HSPs at
// once, moving and sorting stay fast.
import type { HspCoordinates, Interval, Unit } from '../domain/coordinates';
import {
  alignmentFasta,
  extractionTarget,
  hspLabel,
  missingAlignmentNote,
  planExtraction,
  readRequests,
  recordLabel,
  sequenceFasta,
  type ExtractedHsp,
  type ExtractionJoin,
  type ExtractionRegion,
  type RunRef,
} from '../domain/extraction';
import type { Outfmt6Row } from '../domain/outfmt6';
import type { InputRole, ProgramId, SequenceKind } from '../domain/programs';
import type { DataGateway } from '../ports/data';
import type { Downloader } from '../ports/download';
import type { HspRecord } from '../ports/engine';
import type { AppState, RunView } from './coordinator';
import type { HspId } from './results';
import { Store } from './store';

/** The file names of the two outputs; they never name an input (an application choice). */
export const SEQUENCES_FILE = 'losat-candidates.fa';
export const ALIGNMENTS_FILE = 'losat-candidates-aligned.fa';
/** The outputs are text (FASTA); the other exports of the application are saved as text/plain too. */
export const FASTA_MIME = 'text/plain';

/** A record of a run's input that an HSP lies on. */
export interface CandidateRecord {
  /** 0-based position of the record in the run's input of its role: the HSP's `q_idx` or `s_idx`. */
  readonly position: number;
  /** The record's ID in the run's record table; empty for a record without one (`recordLabel` names it). */
  readonly id: string;
  /** The record's length in `unit` (the record table's length). */
  readonly length: number;
  readonly kind: SequenceKind;
  readonly unit: Unit;
}

/** The run that a candidate comes from, as the tray names it. */
export interface CandidateRun {
  readonly runId: string;
  /** The run's number in the working session (`RunSnapshot.number`). */
  readonly number: number;
  /** The run's name (NCBI's "Job Title"), if it has one. */
  readonly title?: string;
  readonly program: ProgramId;
}

/** What the tray takes of an HSP (the results browser's `candidateSources`). */
export interface CandidateSource {
  readonly id: HspId;
  /** The engine's HSP index (docs/web/abi_v2.md §8): its position in the run's final hit list. */
  readonly index: number;
  readonly run: CandidateRun;
  readonly query: CandidateRecord;
  readonly subject: CandidateRecord;
  /** The HSP record's coordinates and frames, as written (domain/coordinates.ts). */
  readonly coordinates: HspCoordinates;
  /** The HSP's outfmt 6 row, split into its fields as the engine wrote them. */
  readonly row: Outfmt6Row;
}

export interface Candidate extends CandidateSource {
  /** `candidateKey(id)`: stable, the same HSP has the same key. */
  readonly key: string;
  /** The user's note (free text). */
  readonly note: string;
  /** When the candidate was added (ms since the epoch). */
  readonly addedAt: number;
  /** The order of addition over the tray's life: the `added` order, also among candidates added at once. */
  readonly serial: number;
}

/** The orders that `sortBy` writes into the tray's own order. */
export type TrayOrder =
  /** The order of addition. */
  | 'added'
  /** The run's number, then the engine's order of the HSPs (their index). */
  | 'run'
  /** The subject's ID, then the HSP's position on it (its smaller, then larger coordinate), then the run. */
  | 'subject';

/** The intervals of an extracted sequence that a record's end cut. */
export interface ClippedSequence {
  /** The record's name in the header (`recordLabel`). */
  readonly name: string;
  readonly runNumber: number;
  readonly role: InputRole;
  /** 0-based position of the record in the run's input (the header's `{role}_record` is this plus 1). */
  readonly position: number;
  /** The HSP labels (`Q.R`) of the sequence. */
  readonly hsps: readonly string[];
  readonly requested: Interval;
  readonly actual: Interval;
  /** The record's length, in `unit` (the header's `length=` and `unit=`). */
  readonly recordLength: number;
  readonly unit: Unit;
}

/** What `extract` wrote. */
export interface SequencesSummary {
  readonly output: 'sequences';
  readonly fileName: string;
  readonly role: InputRole;
  /** Candidates extracted. */
  readonly candidates: number;
  /** FASTA records written. */
  readonly sequences: number;
  readonly bytes: number;
  /** The sequences that a record's end cut, with the interval asked for and the one written. */
  readonly clipped: readonly ClippedSequence[];
  /**
   * Why the strand of some HSPs is unknown (BLASTN HSPs of one letter, whose record does not give
   * it; their headers say `hit_strand=unknown`): one sentence each, naming the HSP.
   */
  readonly unknownStrand: readonly string[];
}

/** What `exportAlignments` wrote. */
export interface AlignmentsSummary {
  readonly output: 'alignments';
  readonly fileName: string;
  /** Candidates exported. */
  readonly candidates: number;
  /** HSPs written (a query and a subject record each). */
  readonly alignments: number;
  readonly bytes: number;
  /** The HSPs whose record has no aligned rows, so nothing was written for them: one sentence each. */
  readonly missing: readonly string[];
}

export type OutputResult<S> = { readonly ok: true; readonly summary: S } | { readonly ok: false; readonly message: string };
export type AddResult =
  | { readonly ok: true; readonly added: number; /** Already in the tray (or given twice). */ readonly already: number }
  | { readonly ok: false; readonly message: string };

export interface TrayState {
  /** In the tray's own order (added at the end, moved, or rewritten by `sortBy`). */
  readonly candidates: readonly Candidate[];
  /** Keys of the candidates that `extract` and `exportAlignments` use. */
  readonly selected: ReadonlySet<string>;
  /** The output being written; another waits until it ends. */
  readonly busy?: 'sequences' | 'alignments';
  /** Why the last addition or output was refused or failed (English, for the screen). */
  readonly message?: string;
  /** The last output written. */
  readonly last?: SequencesSummary | AlignmentsSummary;
}

/** Where the candidates of one run come from (S14 item 5, the tray's list of origins). */
export interface CandidateOrigin {
  readonly runId: string;
  readonly number: number;
  readonly title?: string;
  readonly program: ProgramId;
  /** The options of the run's argv: the words after the program and the two inputs (plan §5.3). */
  readonly options: readonly string[];
  /** The -query and -subject names and the SHA-256 of the bytes that the engine searched. */
  readonly query: { readonly name: string; readonly sha256: string };
  readonly subject: { readonly name: string; readonly sha256: string };
  readonly engineBuild?: string;
  /** When the run ended (ms since the epoch). */
  readonly endedAt?: number;
  /** Candidates of this run in the tray. */
  readonly candidates: number;
}

export interface ExtractOptions {
  /** The records to read: the subjects' (default) or the queries'. */
  readonly role?: InputRole;
  /** The HSP's interval (default), with flanks, or the whole record. */
  readonly region?: ExtractionRegion;
  /** One sequence per HSP (default), or one spanning interval per record. */
  readonly join?: ExtractionJoin;
}

export interface CandidateTrayDeps {
  /** The coordinator's state: which runs completed, and their snapshots. */
  readonly runs: Store<AppState>;
  readonly data: Pick<DataGateway, 'readResidues' | 'readHspRecords'>;
  readonly downloader: Downloader;
  readonly now: () => number;
}

/** The key of an HSP in the tray: `runId/qIdx/rank`. */
export const candidateKey = (id: HspId): string => `${id.runId}/${id.qIdx}/${id.rank}`;

const INITIAL: TrayState = Object.freeze<TrayState>({ candidates: Object.freeze([]), selected: new Set<string>() });

type Action = 'add' | 'extract' | 'export';
const ONLY_COMPLETED: Readonly<Record<Action, string>> = {
  add: 'Only HSPs of completed runs can be added to the candidates.',
  extract: 'Only candidates of completed runs can be extracted.',
  export: 'Only candidates of completed runs can be exported.',
};

export class CandidateTray {
  readonly state = new Store<TrayState>(INITIAL);
  private serial = 0;

  constructor(private readonly deps: CandidateTrayDeps) {}

  /**
   * Adds HSPs at the end of the tray, selected, in the order given; an HSP already in the tray
   * stays where it is, with its note. Refused, with nothing added, when an HSP's run has not
   * completed (REQ-10).
   */
  add(sources: readonly CandidateSource[]): AddResult {
    const refused = this.refusal(sources.map((source) => source.run), 'add');
    if (refused !== undefined) return this.refuse(refused);
    const state = this.state.get();
    const present = new Set(state.candidates.map((candidate) => candidate.key));
    const selected = new Set(state.selected);
    const addedAt = this.deps.now();
    const added: Candidate[] = [];
    for (const source of sources) {
      const key = candidateKey(source.id);
      if (present.has(key)) continue;
      present.add(key);
      selected.add(key);
      const { id, index, run, query, subject, coordinates, row } = source;
      added.push({ id, index, run, query, subject, coordinates, row, key, note: '', addedAt, serial: this.serial++ });
    }
    this.set({ candidates: [...state.candidates, ...added], selected, message: undefined });
    return { ok: true, added: added.length, already: sources.length - added.length };
  }

  remove(keys: readonly string[]): void {
    const gone = new Set(keys);
    const state = this.state.get();
    this.set({
      candidates: state.candidates.filter((candidate) => !gone.has(candidate.key)),
      selected: new Set([...state.selected].filter((key) => !gone.has(key))),
    });
  }

  /** Moves a candidate to position `toIndex` (0-based, kept within the tray) of the tray's order. */
  move(key: string, toIndex: number): void {
    const candidates = [...this.state.get().candidates];
    const from = candidates.findIndex((candidate) => candidate.key === key);
    if (from < 0) return;
    const to = Math.max(0, Math.min(candidates.length - 1, Math.trunc(toIndex)));
    if (to === from) return;
    const [moved] = candidates.splice(from, 1);
    candidates.splice(to, 0, moved!);
    this.set({ candidates });
  }

  /** Rewrites the tray's order by `order` (see `TrayOrder`); later moves start from it. */
  sortBy(order: TrayOrder): void {
    const candidates = [...this.state.get().candidates];
    candidates.sort(COMPARE[order]);
    this.set({ candidates });
  }

  setNote(key: string, text: string): void {
    let found = false;
    const candidates = this.state.get().candidates.map((candidate) => {
      if (candidate.key !== key) return candidate;
      found = true;
      return { ...candidate, note: text };
    });
    if (found) this.set({ candidates });
  }

  /** Selects or deselects candidates for the outputs. */
  select(keys: readonly string[], on: boolean): void {
    const state = this.state.get();
    const present = new Set(state.candidates.map((candidate) => candidate.key));
    const selected = new Set(state.selected);
    for (const key of keys) {
      if (!present.has(key)) continue;
      if (on) selected.add(key);
      else selected.delete(key);
    }
    this.set({ selected });
  }

  selectAll(on: boolean): void {
    this.set({ selected: new Set(on ? this.state.get().candidates.map((candidate) => candidate.key) : []) });
  }

  /** Empties the tray. The runs and their results are not touched. */
  clear(): void {
    this.set({ candidates: [], selected: new Set(), message: undefined, last: undefined });
  }

  /**
   * One entry per run that has candidates, in run order: what the run searched (its title,
   * program, options, the names and SHA-256 of its inputs), the engine build and when it ended,
   * from the coordinator's record of the run.
   */
  origins(): readonly CandidateOrigin[] {
    const counts = new Map<string, number>();
    for (const candidate of this.state.get().candidates) counts.set(candidate.run.runId, (counts.get(candidate.run.runId) ?? 0) + 1);
    return this.deps.runs
      .get()
      .runs.filter((run) => counts.has(run.snapshot.runId))
      .map(({ snapshot, record }) => ({
        runId: snapshot.runId,
        number: snapshot.number,
        ...(snapshot.title === undefined ? {} : { title: snapshot.title }),
        program: snapshot.program,
        options: snapshot.argv.slice(5),
        query: { name: snapshot.query.name, sha256: snapshot.query.sha256 },
        subject: { name: snapshot.subject.name, sha256: snapshot.subject.sha256 },
        ...(record.engineBuild === undefined ? {} : { engineBuild: record.engineBuild }),
        ...(record.endedAt === undefined ? {} : { endedAt: record.endedAt }),
        candidates: counts.get(snapshot.runId)!,
      }))
      .sort((a, b) => a.number - b.number);
  }

  /**
   * Saves the original residues of the selected candidates' records as one multi-FASTA
   * (`SEQUENCES_FILE`), in tray order: the intervals that domain/extraction.ts plans for the
   * options, read with one `readResidues` call per record. A failure saves nothing, leaves the
   * tray as it was, and says why.
   */
  async extract(options: ExtractOptions = {}): Promise<OutputResult<SequencesSummary>> {
    const role = options.role ?? 'subject';
    const chosen = this.chosen('extract');
    if (typeof chosen === 'string') return this.refuse(chosen);
    const views = this.runViews();
    this.set({ busy: 'sequences', message: undefined });
    try {
      const targets = chosen.map((candidate) => extractionTarget(runRef(candidate), extractedHsp(candidate), role, candidate[role]));
      const plan = planExtraction(targets, { region: options.region ?? { kind: 'hit' }, join: options.join ?? 'separate' });
      const parts: Uint8Array[] = new Array<Uint8Array>(plan.pieces.length);
      for (const request of readRequests(plan)) {
        const view = views.get(request.runId)!;
        const read = await this.deps.data.readResidues(view.snapshot[role].revisionIds, request.position, request.intervals);
        request.pieces.forEach((p, i) => {
          const piece = plan.pieces[p]!;
          if (read.origin.id !== piece.recordId || read.origin.length !== piece.recordLength) {
            throw new Error(
              `${role} record ${piece.position + 1} of run ${piece.runNumber} is "${piece.recordId}" (${piece.recordLength} ${piece.unit}) ` +
                `in the run, but its source now gives "${read.origin.id}" (${read.origin.length})`,
            );
          }
          parts[p] = sequenceFasta(piece, read.residues[i]!);
        });
      }
      const bytes = concat(parts);
      this.deps.downloader.save(SEQUENCES_FILE, bytes, FASTA_MIME);
      const summary: SequencesSummary = {
        output: 'sequences',
        fileName: SEQUENCES_FILE,
        role,
        candidates: chosen.length,
        sequences: plan.pieces.length,
        bytes: bytes.length,
        clipped: plan.pieces
          .filter((piece) => piece.interval.clippedLeft || piece.interval.clippedRight)
          .map((piece) => ({
            name: recordLabel(piece.recordId, piece.role, piece.position),
            runNumber: piece.runNumber,
            role: piece.role,
            position: piece.position,
            hsps: piece.hsps,
            requested: piece.interval.requested,
            actual: piece.interval.actual,
            recordLength: piece.recordLength,
            unit: piece.unit,
          })),
        unknownStrand: plan.notes,
      };
      this.set({ busy: undefined, last: summary });
      return { ok: true, summary };
    } catch (error) {
      return this.fail(`The sequences could not be extracted: ${errorMessage(error)}`);
    }
  }

  /**
   * Saves the gapped alignments of the selected candidates (`ALIGNMENTS_FILE`), in tray order:
   * the query and subject rows of each HSP record exactly as the engine wrote them, read with one
   * `readHspRecords` call per run; never mixed with extracted residues (design §11.4). A failure
   * saves nothing, leaves the tray as it was, and says why.
   */
  async exportAlignments(): Promise<OutputResult<AlignmentsSummary>> {
    const chosen = this.chosen('export');
    if (typeof chosen === 'string') return this.refuse(chosen);
    this.set({ busy: 'alignments', message: undefined });
    try {
      const byRun = new Map<string, Candidate[]>();
      for (const candidate of chosen) {
        const list = byRun.get(candidate.run.runId);
        if (list === undefined) byRun.set(candidate.run.runId, [candidate]);
        else list.push(candidate);
      }
      const records = new Map<string, HspRecord>();
      for (const [runId, list] of byRun) {
        const read = await this.deps.data.readHspRecords(runId, list.map((candidate) => candidate.index));
        list.forEach((candidate, i) => {
          const record = read[i];
          if (record === undefined || record.q_idx !== candidate.query.position || record.s_idx !== candidate.subject.position || record.rank !== candidate.id.rank) {
            throw new Error(`HSP record ${candidate.index} of run ${candidate.run.number} is not HSP ${hspLabel(candidate.id.qIdx, candidate.id.rank)}`);
          }
          records.set(candidate.key, record);
        });
      }
      const parts: Uint8Array[] = [];
      const missing: string[] = [];
      for (const candidate of chosen) {
        const record = records.get(candidate.key)!;
        const fasta = alignmentFasta(runRef(candidate), record, { query: candidate.query.id, subject: candidate.subject.id });
        if (fasta === undefined) missing.push(missingAlignmentNote(runRef(candidate), record));
        else parts.push(fasta);
      }
      if (parts.length === 0) throw new Error(`no selected candidate has aligned sequences in its HSP record. ${missing.join(' ')}`);
      const bytes = concat(parts);
      this.deps.downloader.save(ALIGNMENTS_FILE, bytes, FASTA_MIME);
      const summary: AlignmentsSummary = {
        output: 'alignments',
        fileName: ALIGNMENTS_FILE,
        candidates: chosen.length,
        alignments: parts.length,
        bytes: bytes.length,
        missing,
      };
      this.set({ busy: undefined, last: summary });
      return { ok: true, summary };
    } catch (error) {
      return this.fail(`The alignments could not be exported: ${errorMessage(error)}`);
    }
  }

  /** The selected candidates in tray order, or why an output cannot be written now. */
  private chosen(action: 'extract' | 'export'): readonly Candidate[] | string {
    const state = this.state.get();
    if (state.busy !== undefined) return 'The candidates are being written to a file; wait until that ends.';
    const chosen = state.candidates.filter((candidate) => state.selected.has(candidate.key));
    if (chosen.length === 0) return action === 'extract' ? 'Select the candidates to extract.' : 'Select the candidates to export.';
    return this.refusal(chosen.map((candidate) => candidate.run), action) ?? chosen;
  }

  /** Why candidates of these runs are refused (REQ-10), or undefined when every run has completed. */
  private refusal(runs: readonly CandidateRun[], action: Action): string | undefined {
    const views = this.runViews();
    const seen = new Set<string>();
    for (const run of runs) {
      if (seen.has(run.runId)) continue;
      seen.add(run.runId);
      const view = views.get(run.runId);
      if (view?.status !== 'completed') return `${ONLY_COMPLETED[action]} ${whyNotCompleted(view, run.number)}`;
    }
    return undefined;
  }

  private runViews(): ReadonlyMap<string, RunView> {
    return new Map(this.deps.runs.get().runs.map((run) => [run.snapshot.runId, run]));
  }

  private refuse(message: string): { readonly ok: false; readonly message: string } {
    this.set({ message });
    return { ok: false, message };
  }

  private fail(message: string): { readonly ok: false; readonly message: string } {
    this.set({ busy: undefined, message });
    return { ok: false, message };
  }

  private set(change: { [K in keyof TrayState]?: TrayState[K] | undefined }): void {
    const next = { ...this.state.get() } as Record<string, unknown>;
    for (const [key, value] of Object.entries(change)) {
      if (value === undefined) delete next[key];
      else next[key] = value;
    }
    this.state.set(next as unknown as TrayState);
  }
}

function whyNotCompleted(view: RunView | undefined, number: number): string {
  if (view === undefined) return `Run ${number} is not in this working session.`;
  switch (view.status) {
    case 'cancelled':
      return `Run ${number} was cancelled, and a cancelled run keeps no results to take candidates from.`;
    case 'failed':
      return `Run ${number} failed, and a failed run keeps no results to take candidates from.`;
    default:
      return `Run ${number} has not completed yet (it is ${view.status}).`;
  }
}

const runRef = (candidate: Candidate): RunRef => ({ runId: candidate.run.runId, number: candidate.run.number, program: candidate.run.program });

/** The HSP record fields that extraction plans with, from the candidate's copy of them. */
const extractedHsp = (candidate: Candidate): ExtractedHsp => ({
  ...candidate.coordinates,
  q_idx: candidate.query.position,
  s_idx: candidate.subject.position,
  rank: candidate.id.rank,
});

const compareText = (a: string, b: string): number => (a < b ? -1 : a > b ? 1 : 0);
const subjectFrom = (candidate: Candidate) => Math.min(candidate.coordinates.s_start, candidate.coordinates.s_end);
const subjectTo = (candidate: Candidate) => Math.max(candidate.coordinates.s_start, candidate.coordinates.s_end);
const byRun = (a: Candidate, b: Candidate) => a.run.number - b.run.number || a.index - b.index;

const COMPARE: Readonly<Record<TrayOrder, (a: Candidate, b: Candidate) => number>> = {
  added: (a, b) => a.serial - b.serial,
  run: byRun,
  subject: (a, b) =>
    compareText(recordLabel(a.subject.id, 'subject', a.subject.position), recordLabel(b.subject.id, 'subject', b.subject.position)) ||
    subjectFrom(a) - subjectFrom(b) ||
    subjectTo(a) - subjectTo(b) ||
    a.run.number - b.run.number ||
    a.subject.position - b.subject.position ||
    a.index - b.index,
};

function concat(parts: readonly Uint8Array[]): Uint8Array {
  const out = new Uint8Array(parts.reduce((total, part) => total + part.length, 0));
  let at = 0;
  for (const part of parts) {
    out.set(part, at);
    at += part.length;
  }
  return out;
}

function errorMessage(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
