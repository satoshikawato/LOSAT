// The search draft (design §3.1, §3.2; plan §5.2-§5.4): the job being edited on the search
// screen. It is independent of the runs, so the next job can be edited while a search
// runs, and "Run LOSAT" freezes the draft's inputs and conditions into run snapshots
// (Coordinator.enqueueAll; runs search one at a time); the draft itself stays, so a subject
// can be kept while the query changes.
//
// Inputs: every pasted text and every file becomes a source of the data layer and is
// indexed when it is added (DatasetStore). The record table of a source is its
// `base` revision; the records that the user leaves out make other revisions, created on
// demand and reused for the same selection, so that a subject whose selection does not
// change is the same revision from search to search (the engine then keeps it, R1).
// Whether the engine reads a source is the engine's verdict (`checkInput`, the program's
// `register`), shown with its message. Several sources of a role are one combined search
// input, or separate searches (one run per source, started as a group). A box of accessions
// only is fetched when the user asks (fetchAccessions, line WF-1) and becomes a pasted text;
// with "Fastest" (line WF-3), the first fetch of a kind races the sources (source-race.ts); a
// box that no one source serves is fetched by groups, one source a press (fetchGroup, WF-5).
// An Example (useExample, line WF-2) replaces both roles: its query becomes the query box's
// text and its subject set the only subject source.
import {
  accessionLines,
  FETCH_ACCESSION_LIMIT,
  FETCH_BYTE_LIMIT,
  FETCH_SOURCES,
  fetchGroups,
  servingSources,
  sourceLabel,
  sourceList,
  unserved,
  type FetchChoice,
  type FetchGroup,
  type FetchSourceId,
} from '../domain/accession';
import { buildArgv, COMBINED_NAMES, PASTED_NAMES, splitOptions } from '../domain/argv';
import { includedRecords, type DatasetRecord, type DatasetRevision, type FastaParserKind } from '../domain/dataset';
import type { ExampleId } from '../domain/examples';
import { formControlOf, formParameters, setFields, type FieldValue, type FormValues } from '../domain/parameters';
import {
  fastaParser,
  PROGRAMS,
  programById,
  sequenceKind,
  type InputRole,
  type ProgramId,
  type SequenceKind,
} from '../domain/programs';
import { REGION_FLAG, regionProblem, regionValue, type RegionText } from '../domain/region';
import type { DataGateway } from '../ports/data';
import type { EngineGateway, ProgramDescription } from '../ports/engine';
import type { AccessionFetcher, FetchResult } from '../ports/fetch';
import type { InputCheck } from '../ports/input-check';
import type { EnqueueAllResult, SearchRequest } from './coordinator';
import { SourceRace, type RaceStanding } from './source-race';
import { SourceRefs } from './source-refs';
import { Store } from './store';

export type InputMode = 'combined' | 'separate';

/** The engine's verdict on a source (ports/input-check.ts). */
export type SourceCheck =
  | { readonly state: 'pending' }
  | { readonly state: 'ok'; readonly records: number }
  /** `record` is the index, in the record table, of the record that the message names. */
  | { readonly state: 'refused'; readonly message: string; readonly record?: number }
  /** The check could not run (not a verdict on the input). */
  | { readonly state: 'error'; readonly message: string };

/** A record table of a source and the records that runs leave out of it. */
export interface SourceTable {
  readonly base: DatasetRevision;
  readonly excluded: readonly number[];
  readonly indexMs?: number;
}

export interface DraftSource {
  /** Identifies the source within the draft. */
  readonly key: string;
  readonly origin: 'file' | 'paste';
  readonly name: string;
  readonly size: number;
  /** The File of the source (pasted text is a File too); a source is indexed again from it. */
  readonly file: File;
  /** The Example set whose copy the file is (useExample). */
  readonly example?: ExampleId;
  readonly status: 'indexing' | 'ready' | 'failed';
  /** The FASTA reader kind of the index (running, failed or its record table), the program's for the role (TD-8). */
  readonly parser: FastaParserKind;
  /** Why the source could not be indexed (the index scan's error). */
  readonly error?: string;
  /** The record table: the revision of the data layer that includes every record. */
  readonly base?: DatasetRevision;
  /** Indices of the records that runs leave out, ascending. */
  readonly excluded: readonly number[];
  /** The data layer's id of the source, once it is added (the other kind's index reads the same source). */
  readonly sourceId?: string;
  /**
   * The record tables of the other FASTA reader kind (the one that `parser` is not), with the
   * records left out of each: a program that reads the role with another kind puts the current
   * table here, and another program that reads it with this kind takes it back.
   */
  readonly tables?: Readonly<Partial<Record<FastaParserKind, SourceTable>>>;
  /** The first lines of the source, as text, for its preview. */
  readonly head?: string;
  readonly check?: SourceCheck;
  /** Milliseconds that the index (scan and record hashes) took. */
  readonly indexMs?: number;
  /** Milliseconds that the latest engine check took (the run input and `register`). */
  readonly checkMs?: number;
}

/** The fetch of the paste box's accessions (fetchAccessions). */
export type FetchState =
  /** Waiting for the race of the sources ("Fastest"): no accession was sent yet. */
  | { readonly state: 'racing' }
  /** `done` of `total` accessions have come. */
  | { readonly state: 'fetching'; readonly source: FetchSourceId; readonly done: number; readonly total: number }
  | { readonly state: 'failed'; readonly message: string };

export interface RoleDraft {
  /**
   * The text of the paste box; its source is the role's paste source. A text of accessions only
   * has no source until it is fetched.
   */
  readonly paste: string;
  readonly fetch?: FetchState;
  /**
   * The FASTA that the groups of a box that no one source serves have brought (fetchGroup), by
   * line of the box, with the source it came from; the box is written once every line has come.
   */
  readonly fetched?: Readonly<Record<number, { readonly fasta: string; readonly source: FetchSourceId }>>;
  /** The paste source first, then the files in the order they were added. */
  readonly sources: readonly DraftSource[];
  readonly mode: InputMode;
  /**
   * The regions chosen on the role's only record (DW-9), by that record (`regionIdentity`). One
   * applies only while its record is still the role's only included record; another record, or
   * another count of records, leaves it aside, and it comes back with its record (a mode that
   * reads the source with the other kind has other records, and has its own regions).
   */
  readonly regions?: Readonly<Record<string, RegionText>>;
}

export type Validation =
  | { readonly state: 'idle' | 'checking' | 'ok' }
  | { readonly state: 'invalid'; readonly message: string };

export interface DraftState {
  readonly program: ProgramId;
  readonly query: RoleDraft;
  readonly subject: RoleDraft;
  /** The form's values of each program: what the user set, nothing else. */
  readonly values: Readonly<Record<ProgramId, FormValues>>;
  /** The "Additional options" of each program, as typed: other options, as on a command line (`splitOptions`). */
  readonly extra: Readonly<Record<ProgramId, string>>;
  /** The Job Title of the next runs, as typed ('' for none); never part of the argv. */
  readonly title: string;
  readonly threads: number | 'auto';
  /** Where "Fetch" asks for accessions; the choice stays for the session. */
  readonly fetchSource: FetchChoice;
  /** What the race of the sources found for each kind (source-race.ts), kept for the tab. */
  readonly race: Readonly<Partial<Record<SequenceKind, RaceStanding>>>;
  /** The engine's `describe` of the program, once it has arrived. */
  readonly description?: ProgramDescription;
  readonly descriptionError?: string;
  /** The engine's `validate` of the argv that the draft makes now. */
  readonly validation: Validation;
  readonly submitting: boolean;
  readonly message?: { readonly kind: 'info' | 'error'; readonly text: string };
}

export interface DraftDeps {
  readonly engine: Pick<EngineGateway, 'describe' | 'validate' | 'prepare'>;
  readonly data: Pick<
    DataGateway,
    'addSource' | 'indexSource' | 'reviseDataset' | 'checkInput' | 'previewSource' | 'dropSource'
  >;
  /**
   * Who uses the data layer's sources, shared with the coordinator, whose runs hold theirs; a
   * draft without it drops a source once the draft is done with it.
   */
  readonly refs?: SourceRefs;
  readonly enqueueAll: (requests: readonly SearchRequest[]) => Promise<EnqueueAllResult>;
  readonly fetcher: AccessionFetcher;
  readonly now: () => number;
  /** Milliseconds of quiet before a pasted text is indexed and the argv is validated. */
  readonly debounceMs?: number;
}

const ROLES: readonly InputRole[] = ['query', 'subject'];
const PASTE_KEY = 'paste';
/** Bytes of a source read for its preview, and the lines shown. */
const HEAD_BYTES = 2048;
const HEAD_LINES = 6;
const HEAD_LINE_CHARS = 120;
/** "query record 3 (id) …": the record that an engine message names (1-based, among the input's records). */
const RECORD_IN_MESSAGE = /\b(?:query|subject) record (\d+)\b/;

const emptyRole = (): RoleDraft => ({ paste: '', sources: [], mode: 'combined' });

export class SearchDraft {
  readonly state: Store<DraftState>;
  private readonly descriptions = new Map<ProgramId, Promise<ProgramDescription>>();
  /** The descriptions that have arrived, so that a program's form is complete when it is chosen again. */
  private readonly described = new Map<ProgramId, ProgramDescription>();
  /** The engine's verdicts, by program, role and revision: a selection checked once is not checked again. */
  private readonly verdicts = new Map<string, Promise<InputCheck>>();
  /** Revisions of a selection, by base revision and excluded records (revisionFor). */
  private readonly revisions = new Map<string, Promise<string>>();
  /** The data layer's sources that the form shows, which the draft holds (SourceRefs). */
  private formSources: ReadonlySet<string> = new Set();
  /** The revisions of each of those sources, whose cached verdicts and selections go with the source. */
  private readonly sourceRevisions = new Map<string, Set<string>>();
  private readonly refs: SourceRefs;
  /** The latest request of each asynchronous task; older results are dropped. */
  private readonly generations = new Map<string, number>();
  private readonly timers = new Map<string, ReturnType<typeof setTimeout>>();
  private readonly workers = new Map<string, () => void>();
  private readonly pending = new Set<Promise<unknown>>();
  /** The running fetch of each role's accessions. */
  private readonly fetches = new Map<InputRole, AbortController>();
  private readonly race: SourceRace;
  /** The argv of the latest validation that the engine accepted (the prepare hint's argv). */
  private validArgv: readonly string[] | undefined;
  private nextKey = 1;

  constructor(private readonly deps: DraftDeps) {
    this.refs = deps.refs ?? new SourceRefs((sourceId) => void deps.data.dropSource(sourceId).catch(() => undefined));
    this.race = new SourceRace({
      fetcher: deps.fetcher,
      now: deps.now,
      onChange: (kind, standing) => {
        const race = { ...this.state.get().race };
        if (standing === undefined) delete race[kind];
        else race[kind] = standing;
        this.set({ race });
      },
    });
    const values = Object.fromEntries(PROGRAMS.map((program) => [program.id, Object.freeze({})])) as Record<
      ProgramId,
      FormValues
    >;
    this.state = new Store<DraftState>({
      program: 'blastn',
      query: emptyRole(),
      subject: emptyRole(),
      values,
      extra: Object.fromEntries(PROGRAMS.map((program) => [program.id, ''])) as Record<ProgramId, string>,
      title: '',
      threads: 'auto',
      fetchSource: 'fastest',
      race: {},
      validation: { state: 'idle' },
      submitting: false,
    });
    this.loadDescription('blastn');
    this.scheduleValidation();
  }

  // --- program and form -------------------------------------------------------------------

  setProgram(program: ProgramId): void {
    if (this.state.get().program === program) return;
    const previous = programById(this.state.get().program);
    const next = programById(program);
    for (const role of ROLES) {
      // A fetch of the other kind's record would write it into the box.
      if (sequenceKind(previous, role) !== sequenceKind(next, role)) this.cancelFetch(role);
    }
    const description = this.described.get(program);
    this.set({ program, description, descriptionError: undefined, message: undefined });
    this.loadDescription(program);
    for (const role of ROLES) {
      const parser = fastaParser(next, role);
      for (const source of this.role(role).sources) {
        if (source.parser !== parser) {
          // A program that reads the role with another FASTA reader kind needs another record table (plan TD-8).
          this.reindex(role, source.key, parser);
        } else {
          this.scheduleCheck(role, source.key, 0);
        }
      }
    }
    this.scheduleValidation();
  }

  setField(flag: string, value: FieldValue): void {
    this.setFields([[flag, value]]);
  }

  /** Sets several fields at once (NCBI's select of a pair: -reward and -penalty, or the gap costs). */
  setFields(entries: ReadonlyArray<readonly [string, FieldValue]>): void {
    const { program, values, description } = this.state.get();
    const next = setFields(programById(program), values[program], entries, description?.parameters);
    this.set({ values: { ...values, [program]: next }, message: undefined });
    this.scheduleValidation();
  }

  /**
   * Clears the values of some fields of the program ("Restore default search parameters"):
   * the engine's defaults apply to them again.
   */
  resetFields(flags: readonly string[]): void {
    const { program, values } = this.state.get();
    const next: Record<string, FieldValue> = { ...values[program] };
    for (const flag of flags) delete next[flag];
    this.set({ values: { ...values, [program]: Object.freeze(next) }, message: undefined });
    this.scheduleValidation();
  }

  /** The program's "Additional options": words of other options, after the form's in the argv. */
  setExtra(text: string): void {
    const { program, extra } = this.state.get();
    this.set({ extra: { ...extra, [program]: text }, message: undefined });
    this.scheduleValidation();
  }

  /** The Job Title of the runs that "Run LOSAT" makes next. */
  setTitle(title: string): void {
    this.set({ title, message: undefined });
  }

  setThreads(threads: number | 'auto'): void {
    this.set({ threads, message: undefined });
    // After a pause, as for typing: stepping through the choices hints only the last one.
    this.debounce('prepare', () => this.prepareEngine());
  }

  // --- inputs -----------------------------------------------------------------------------

  /**
   * The paste box's text. Its source shows that it is being read at once, and is indexed
   * after a pause in typing.
   */
  setPaste(role: InputRole, text: string): void {
    // A newer text supersedes a fetch that is still running.
    this.cancelFetch(role);
    this.updateRole(role, (draft) => ({ ...draft, paste: text }));
    this.set({ message: undefined });
    if (text === '' || accessionLines(text) !== undefined) {
      this.dropSource(role, PASTE_KEY);
      return;
    }
    const file = new File([text], PASTED_NAMES[role], { type: 'text/plain' });
    const source: DraftSource = {
      key: PASTE_KEY,
      origin: 'paste',
      name: PASTED_NAMES[role],
      size: file.size,
      file,
      status: 'indexing',
      parser: this.parserOf(role),
      excluded: [],
    };
    // A newer text supersedes an index of the older one that is still running.
    this.nextGeneration(`${role}:${PASTE_KEY}`);
    this.cancelTask(`check:${role}:${PASTE_KEY}`);
    this.updateRole(role, (draft) => ({
      ...draft,
      sources: draft.sources.some((s) => s.key === PASTE_KEY)
        ? draft.sources.map((s) => (s.key === PASTE_KEY ? source : s))
        : [source, ...draft.sources],
    }));
    // The program may change before the pause ends (setProgram then sets the source's kind).
    this.debounce(`paste:${role}`, () => this.track(this.index(role, PASTE_KEY, file, this.parserOf(role))));
    this.scheduleValidation();
  }

  /**
   * Puts the line `>pasted_<role>` before a pasted sequence that has no defline, so that it
   * is FASTA. The user asks for it; the text in the box changes to what is searched.
   */
  addDefline(role: InputRole): void {
    this.setPaste(role, `>pasted_${role}\n${this.role(role).paste.trimStart()}`);
  }

  /** Adds files; each is indexed at once (S12 instructions). */
  addFiles(role: InputRole, files: readonly File[]): void {
    for (const file of files) this.addFile(role, file, undefined);
    this.set({ message: undefined });
    this.scheduleValidation();
  }

  private addFile(role: InputRole, file: File, example: ExampleId | undefined): void {
    const key = `file-${this.nextKey++}`;
    const parser = this.parserOf(role);
    const source: DraftSource = {
      key,
      origin: 'file',
      name: file.name,
      size: file.size,
      file,
      ...(example === undefined ? {} : { example }),
      status: 'indexing',
      parser,
      excluded: [],
    };
    this.updateRole(role, (draft) => ({ ...draft, sources: [...draft.sources, source] }));
    this.track(this.index(role, key, file, parser));
  }

  removeSource(role: InputRole, key: string): void {
    this.dropSource(role, key);
    if (key === PASTE_KEY) this.updateRole(role, (draft) => ({ ...draft, paste: '' }));
  }

  // --- Example sets (line WF-2) ------------------------------------------------------------

  /**
   * Replaces both roles with an Example (as "Use example sequence" does at EBI): `query`, its
   * FASTA, becomes the query box's text, and `subject`, the copy of its set, the only subject
   * source. The other sources and the regions go.
   */
  useExample(example: ExampleId, query: string, subject: File): void {
    for (const role of ROLES) {
      for (const source of this.role(role).sources) {
        if (source.key !== PASTE_KEY) this.dropSource(role, source.key);
      }
      this.updateRole(role, (draft) => {
        const next = { ...draft };
        delete next.regions;
        return next;
      });
    }
    this.setPaste('subject', '');
    this.setPaste('query', query);
    this.addFile('subject', subject, example);
    this.scheduleValidation();
  }

  /** Removes the sources made from an Example's copy, which was deleted and can no longer be read. */
  dropExample(example: ExampleId): void {
    for (const role of ROLES) {
      for (const source of this.role(role).sources) {
        if (source.example === example) this.dropSource(role, source.key);
      }
    }
  }

  // --- accessions (line WF-1) ---------------------------------------------------------------

  /** The accessions of the role's paste box, when it holds accessions only. */
  accessions(role: InputRole): readonly string[] | undefined {
    return accessionLines(this.role(role).paste);
  }

  setFetchSource(source: FetchChoice): void {
    this.set({ fetchSource: source });
    for (const role of ROLES) this.clearFetchFailure(role);
  }

  /** Drops the message of the role's failed fetch (another source is chosen for it). */
  clearFetchFailure(role: InputRole): void {
    if (this.role(role).fetch?.state === 'failed') this.updateRole(role, withoutFetch);
  }

  /**
   * Fetches the FASTA of the paste box's accessions from the chosen source, one after another,
   * and writes it into the box, which is then read as any pasted text. The user asks for it
   * ("Fetch"); only the accessions are sent. With "Fastest", they go to the source that the
   * race of the kind found first among those that serve them. An edit of the box drops the
   * fetch; a failure leaves the box as it was and says why.
   */
  fetchAccessions(role: InputRole): void {
    const accessions = this.accessions(role);
    if (accessions === undefined) return;
    const { fetchSource, program } = this.state.get();
    this.cancelFetch(role);
    const controller = new AbortController();
    this.fetches.set(role, controller);
    this.track(this.fetchAll(role, accessions, fetchSource, sequenceKind(programById(program), role), controller.signal));
  }

  private async fetchAll(
    role: InputRole,
    accessions: readonly string[],
    choice: FetchChoice,
    kind: SequenceKind,
    signal: AbortSignal,
  ): Promise<void> {
    const fail = (message: string) => {
      if (!signal.aborted) this.updateRole(role, (draft) => ({ ...draft, fetch: { state: 'failed', message } }));
    };
    if (accessions.length > FETCH_ACCESSION_LIMIT) {
      fail(`Fetch takes up to ${FETCH_ACCESSION_LIMIT} accessions at a time.`);
      return;
    }
    // Before any request: an accession the source cannot serve stops the fetch without sending the others.
    for (const accession of accessions) {
      const reason = choice === 'fastest' ? undefined : unserved(choice, kind, accession);
      if (reason !== undefined) {
        fail(`${reason}: fetch ${accession} from ${sourceList(servingSources(kind, [accession]))}.`);
        return;
      }
    }
    let source: FetchSourceId;
    if (choice === 'fastest') {
      const eligible = servingSources(kind, accessions);
      if (eligible.length === 0) {
        const only = accessions
          .map((accession) => ({ accession, sources: servingSources(kind, [accession]) }))
          .filter((entry) => entry.sources.length < FETCH_SOURCES.length)
          .map((entry) => `${entry.accession} only from ${sourceList(entry.sources)}`);
        fail(`No one site has all of these accessions: ${only.join(', ')}.`);
        return;
      }
      if (eligible.length === 1) {
        // One site has them: nothing to race.
        source = eligible[0]!;
      } else {
        this.updateRole(role, (draft) => ({ ...draft, fetch: { state: 'racing' } }));
        try {
          source = await this.race.fastest(kind, eligible, signal);
        } catch (error) {
          fail(`The race of the sites failed: ${messageOf(error)}`);
          return;
        }
      }
    } else {
      source = choice;
    }
    const texts: string[] = [];
    let bytes = 0;
    for (const [done, accession] of accessions.entries()) {
      this.updateRole(role, (draft) => ({ ...draft, fetch: { state: 'fetching', source, done, total: accessions.length } }));
      let result: FetchResult;
      try {
        result = await this.deps.fetcher.fetch({ source, kind, accession, maxBytes: FETCH_BYTE_LIMIT - bytes }, signal);
      } catch (error) {
        fail(`${accession} could not be fetched from ${sourceLabel(source)}: ${messageOf(error)}`);
        return;
      }
      if (signal.aborted) return;
      if (!result.ok) {
        fail(fetchProblem(result, source, kind, accession));
        return;
      }
      texts.push(result.fasta);
      bytes += new TextEncoder().encode(result.fasta).byteLength;
    }
    this.fetches.delete(role);
    this.setPaste(role, texts.join(''));
  }

  /**
   * Fetches one group of a box that no one source serves (fetchGroups, line WF-5) from `source`,
   * one of the group's: only that group's accessions are sent, only there, and no race runs. The
   * FASTA is held until every group has come, then written into the box in the box's order. A
   * failure keeps the groups that have come, so that this group can be fetched again alone.
   */
  fetchGroup(role: InputRole, group: number, source: FetchSourceId): void {
    const accessions = this.accessions(role);
    if (accessions === undefined) return;
    const kind = sequenceKind(programById(this.state.get().program), role);
    if (servingSources(kind, accessions).length > 0) return;
    const target = fetchGroups(kind, accessions)[group];
    if (target === undefined || !target.sources.includes(source)) return;
    this.fetches.get(role)?.abort();
    const controller = new AbortController();
    this.fetches.set(role, controller);
    this.track(this.fetchPart(role, accessions, target, source, kind, controller.signal));
  }

  private async fetchPart(
    role: InputRole,
    accessions: readonly string[],
    group: FetchGroup,
    source: FetchSourceId,
    kind: SequenceKind,
    signal: AbortSignal,
  ): Promise<void> {
    const fail = (message: string) => {
      if (!signal.aborted) this.updateRole(role, (draft) => ({ ...draft, fetch: { state: 'failed', message } }));
    };
    if (accessions.length > FETCH_ACCESSION_LIMIT) {
      fail(`Fetch takes up to ${FETCH_ACCESSION_LIMIT} accessions at a time.`);
      return;
    }
    // The limit counts the whole box: the groups that have come, then this one.
    let bytes = 0;
    for (const { fasta } of Object.values(this.role(role).fetched ?? {})) bytes += new TextEncoder().encode(fasta).byteLength;
    const texts: string[] = [];
    for (const [done, accession] of group.accessions.entries()) {
      this.updateRole(role, (draft) => ({ ...draft, fetch: { state: 'fetching', source, done, total: group.accessions.length } }));
      let result: FetchResult;
      try {
        result = await this.deps.fetcher.fetch({ source, kind, accession, maxBytes: FETCH_BYTE_LIMIT - bytes }, signal);
      } catch (error) {
        fail(`${accession} could not be fetched from ${sourceLabel(source)}: ${messageOf(error)}`);
        return;
      }
      if (signal.aborted) return;
      if (!result.ok) {
        const problem = fetchProblem(result, source, kind, accession);
        fail(problem.includes(accession) ? problem : `${accession}: ${problem}`);
        return;
      }
      texts.push(result.fasta);
      bytes += new TextEncoder().encode(result.fasta).byteLength;
    }
    this.fetches.delete(role);
    const fetched = { ...this.role(role).fetched };
    for (const [index, line] of group.lines.entries()) fetched[line] = { fasta: texts[index]!, source };
    if (accessions.every((_, line) => fetched[line] !== undefined)) {
      this.setPaste(role, accessions.map((_, line) => fetched[line]!.fasta).join(''));
    } else {
      this.updateRole(role, (draft) => ({ ...withoutFetch(draft), fetched }));
    }
  }

  /** Stops the role's fetch, and drops the groups that have come (fetchGroup). */
  private cancelFetch(role: InputRole): void {
    this.fetches.get(role)?.abort();
    this.fetches.delete(role);
    const draft = this.role(role);
    if (draft.fetch !== undefined || draft.fetched !== undefined) this.updateRole(role, (next) => withoutFetch(next, true));
  }

  /** Removes a source; the paste box keeps its text. */
  private dropSource(role: InputRole, key: string): void {
    this.cancelTask(`${role}:${key}`);
    if (key === PASTE_KEY) this.cancelTask(`paste:${role}`);
    this.updateRole(role, (draft) => ({ ...draft, sources: draft.sources.filter((source) => source.key !== key) }));
    this.set({ message: undefined });
    this.scheduleValidation();
  }

  setMode(role: InputRole, mode: InputMode): void {
    this.updateRole(role, (draft) => ({ ...draft, mode }));
    this.set({ message: undefined });
  }

  /** Includes or leaves out records of a source (by their index in its record table). */
  setIncluded(role: InputRole, key: string, indices: readonly number[], included: boolean): void {
    this.updateSource(role, key, (source) => {
      const excluded = new Set(source.excluded);
      for (const index of indices) {
        if (included) excluded.delete(index);
        else excluded.add(index);
      }
      return { ...source, excluded: [...excluded].sort((a, b) => a - b), check: { state: 'pending' } };
    });
    this.set({ message: undefined });
    this.scheduleCheck(role, key);
    this.scheduleValidation();
  }

  setRegion(role: InputRole, region: RegionText | undefined): void {
    const record = this.regionIdentity(role);
    if (record === undefined) return;
    this.updateRole(role, (draft) => {
      const regions: Record<string, RegionText> = { ...draft.regions };
      if (region === undefined) delete regions[record];
      else regions[record] = region;
      const { regions: _before, ...rest } = draft;
      void _before;
      return Object.keys(regions).length === 0 ? rest : { ...rest, regions };
    });
    this.set({ message: undefined });
    this.scheduleValidation();
  }

  // --- derived views ----------------------------------------------------------------------

  /** The record that the role's region applies to: the role's only included record. */
  regionRecord(role: InputRole): DatasetRecord | undefined {
    return this.onlyRecord(role)?.record;
  }

  /** The region of the role, if it was chosen on the record that `regionRecord` gives now. */
  region(role: InputRole): RegionText | undefined {
    const record = this.regionIdentity(role);
    return record === undefined ? undefined : this.role(role).regions?.[record];
  }

  private onlyRecord(role: InputRole): { readonly source: DraftSource; readonly record: DatasetRecord } | undefined {
    const included = this.role(role).sources.flatMap((source) =>
      source.status === 'ready' && source.base !== undefined ? includedOf(source).map((record) => ({ source, record })) : [],
    );
    return included.length === 1 ? included[0] : undefined;
  }

  /** The source, record table and record of the role's only record. */
  private regionIdentity(role: InputRole): string | undefined {
    const only = this.onlyRecord(role);
    return only === undefined ? undefined : `${only.source.key}|${only.source.base!.revisionId}|${only.record.index}`;
  }

  /** The parameters of the argv: the form's and the regions'. */
  parameters(): Array<readonly [string, string | true]> {
    const { program, values, description } = this.state.get();
    const parameters = formParameters(programById(program), values[program], description?.parameters);
    for (const role of ROLES) {
      const region = this.region(role);
      const record = this.regionRecord(role);
      if (region === undefined || record === undefined) continue;
      if (regionProblem(region, record.length) === undefined) parameters.push([REGION_FLAG[role], regionValue(region)]);
    }
    return parameters;
  }

  /**
   * The words of the program's "Additional options". Throws on an unclosed quote, and on an
   * option that the form sets (`formControlOf`) or a role's subrange, which are set there.
   */
  extraWords(): readonly string[] {
    const { program, extra, description } = this.state.get();
    const words = splitOptions(extra[program]);
    for (const word of words) {
      const flag = word.split('=', 1)[0]!;
      const role = ROLES.find((r) => REGION_FLAG[r] === flag);
      const control =
        role === undefined ? formControlOf(programById(program), description?.parameters, flag) : `${role === 'query' ? 'Query' : 'Subject'} subrange`;
      if (control !== undefined) throw new Error(`${flag} is set in the form: ${control}`);
    }
    return words;
  }

  /** The name of the role's input in the argv of a combined search (plan §5.3). */
  inputName(role: InputRole): string {
    const sources = this.role(role).sources.filter((source) => source.status !== 'failed');
    if (sources.length === 1) return nameOf(sources[0]!, role);
    return sources.length === 0 ? PASTED_NAMES[role] : COMBINED_NAMES[role];
  }

  /** The runs that "Run LOSAT" would make now. */
  runCount(): number {
    return ROLES.map((role) => this.inputGroups(role).length).reduce((a, b) => a * Math.max(1, b), 1);
  }

  // --- submit -----------------------------------------------------------------------------

  /** Freezes the draft into runs (one, or a group of separate searches) and starts them. */
  async submit(): Promise<EnqueueAllResult> {
    if (this.state.get().submitting) return { ok: false, message: 'The search is being started already.' };
    this.set({ submitting: true, message: undefined });
    try {
      let result: EnqueueAllResult;
      try {
        result = await this.prepareAndEnqueue();
      } catch (error) {
        result = { ok: false, message: `The search could not be started: ${messageOf(error)}` };
      }
      // A started search shows itself: the page moves to the run (stop-describing-ui rule 6).
      if (!result.ok) this.set({ message: { kind: 'error', text: result.message } });
      return result;
    } finally {
      this.refs.release('submit');
      this.set({ submitting: false });
    }
  }

  /**
   * Runs the waiting work at once and settles when no indexing, check, validation or
   * revision is in progress (submit, and tests).
   */
  async idle(): Promise<void> {
    for (;;) {
      if (this.timers.size > 0) {
        for (const name of [...this.timers.keys()]) {
          clearTimeout(this.timers.get(name));
          this.timers.delete(name);
          this.runTimer(name);
        }
        continue;
      }
      if (this.pending.size === 0) return;
      await Promise.allSettled([...this.pending]);
    }
  }

  /** Resolves when the race requests have ended; idle does not wait for them, nor does Run (tests). */
  async raceIdle(): Promise<void> {
    await this.race.idle();
  }

  private async prepareAndEnqueue(): Promise<EnqueueAllResult> {
    // Run never fetches, nor waits for a fetch that is running (fetchAccessions).
    for (const role of ROLES) {
      const problem = this.accessionProblem(role);
      if (problem !== undefined) return { ok: false, message: problem };
    }
    await this.idle();
    const { program, threads, title } = this.state.get();
    for (const role of ROLES) {
      const problem = this.inputProblem(role);
      if (problem !== undefined) return { ok: false, message: problem };
    }
    // Everything from one state, before any wait.
    const groups = { query: this.inputGroups('query'), subject: this.inputGroups('subject') };
    // Their sources stay until the runs hold them, also when the form lets them go meanwhile.
    this.refs.set(
      'submit',
      [...groups.query, ...groups.subject].flatMap((group) => group.sources.map((source) => source.sourceId!)),
    );
    const parameters = this.parameters();
    let extra: readonly string[];
    try {
      extra = this.extraWords();
    } catch (error) {
      return { ok: false, message: messageOf(error) };
    }
    const resolved = {
      query: await Promise.all(groups.query.map((group) => this.resolveGroup(group))),
      subject: await Promise.all(groups.subject.map((group) => this.resolveGroup(group))),
    };
    const requests: SearchRequest[] = [];
    for (const query of resolved.query) {
      for (const subject of resolved.subject) {
        requests.push({
          program,
          query: { dataset: query },
          subject: { dataset: subject },
          parameters,
          ...(extra.length === 0 ? {} : { extra }),
          requestedThreads: threads,
          ...(title.trim() === '' ? {} : { title }),
        });
      }
    }
    return this.deps.enqueueAll(requests);
  }

  /** Why the role cannot be searched as it is (the full message), or undefined. */
  private inputProblem(role: InputRole): string | undefined {
    return this.problemOf(role)?.message;
  }

  /**
   * What keeps the draft from being searched now, in short, for the line above the button;
   * undefined when nothing does. "Run LOSAT" reports the full message.
   */
  readiness(): string | undefined {
    for (const role of ROLES) {
      const problem = this.problemOf(role);
      if (problem !== undefined) return problem.short;
    }
    return undefined;
  }

  /** A box of accessions only, fetched or not yet: it is not searched. */
  private accessionProblem(role: InputRole): string | undefined {
    const accessions = accessionLines(this.role(role).paste);
    if (accessions === undefined) return undefined;
    const title = role === 'query' ? 'query' : 'subject';
    const noun = accessions.length === 1 ? 'accession' : 'accessions';
    const fetch = this.role(role).fetch?.state;
    return fetch === 'fetching' || fetch === 'racing'
      ? `The ${title} ${noun} ${accessions.length === 1 ? 'is' : 'are'} still being fetched.`
      : `Fetch the ${title} ${noun} first.`;
  }

  private problemOf(role: InputRole): { readonly message: string; readonly short: string } | undefined {
    const draft = this.role(role);
    const title = role === 'query' ? 'query' : 'subject';
    const label = role === 'query' ? 'Query' : 'Subject';
    const accessions = this.accessionProblem(role);
    if (accessions !== undefined) return { message: accessions, short: accessions };
    if (draft.sources.length === 0) {
      const message = `Add the ${title} sequences: paste them or open FASTA files.`;
      return { message, short: message };
    }
    for (const source of draft.sources) {
      const where = `${label} (${describeSource(source)})`;
      if (source.status === 'failed') {
        return { message: `${where}: ${source.error ?? 'the input cannot be read'}`, short: `${where} cannot be read.` };
      }
      if (source.status !== 'ready') {
        const message = `${where} is still being read.`;
        return { message, short: message };
      }
    }
    if (draft.sources.every((source) => includedOf(source).length === 0 && (source.base?.records.length ?? 0) > 0)) {
      const message = `Every ${title} record is excluded.`;
      return { message, short: message };
    }
    if (draft.mode === 'combined') {
      // The sources follow one another in one input, where a sequence without a defline
      // would join the record before it (NCBI's reader): only the first can start so.
      const late = draft.sources
        .filter((source) => includedOf(source).length > 0)
        .slice(1)
        .find((source) => includedOf(source)[0]!.sequence_offset === 0);
      if (late !== undefined) {
        const message = `${label} (${describeSource(late)}) starts with a sequence without a defline (">"), which cannot follow another input.`;
        return { message, short: message };
      }
    }
    for (const source of draft.sources) {
      const where = `${label} (${describeSource(source)})`;
      const check = source.check;
      if (check?.state === 'refused') {
        return {
          message: `${where}: ${check.message}`,
          short: `${where}: the engine refuses ${check.record === undefined ? 'this input' : 'a record'}. Exclude it or correct the input to search.`,
        };
      }
      // A check that could not run is not the engine's verdict: the search reads the input
      // with the engine again, and reports what it finds.
    }
    const record = this.regionRecord(role);
    const region = this.region(role);
    if (region !== undefined && record !== undefined) {
      const problem = regionProblem(region, record.length);
      if (problem !== undefined) {
        const message = `${label} region: ${problem}`;
        return { message, short: message };
      }
    }
    return undefined;
  }

  /** The inputs of the role's runs: all sources together, or one per source (Separate). */
  private inputGroups(role: InputRole): ReadonlyArray<{ readonly name: string; readonly sources: readonly DraftSource[] }> {
    const draft = this.role(role);
    const usable = draft.sources.filter((source) => source.status === 'ready');
    if (draft.mode === 'separate' && usable.length > 1) {
      // A source whose records are all excluded has nothing to search alone.
      return usable
        .filter((source) => includedOf(source).length > 0 || (source.base?.records.length ?? 0) === 0)
        .map((source) => ({ name: nameOf(source, role), sources: [source] }));
    }
    return usable.length === 0 ? [] : [{ name: this.inputName(role), sources: usable }];
  }

  private async resolveGroup(group: { readonly name: string; readonly sources: readonly DraftSource[] }) {
    const revisionIds = await Promise.all(group.sources.map((source) => this.revisionFor(source)));
    return { name: group.name, revisionIds, sourceIds: group.sources.map((source) => source.base!.sourceId) };
  }

  // --- asynchronous work ------------------------------------------------------------------

  private loadDescription(program: ProgramId): void {
    let description = this.descriptions.get(program);
    if (description === undefined) {
      description = this.deps.engine.describe(program);
      this.descriptions.set(program, description);
      description.catch(() => this.descriptions.delete(program));
    }
    this.track(
      description.then(
        (value) => {
          this.described.set(program, value);
          if (this.state.get().program === program) this.set({ description: value });
          this.scheduleValidation();
        },
        (error: unknown) => {
          if (this.state.get().program === program) this.set({ descriptionError: messageOf(error) });
        },
      ),
    );
  }

  private parserOf(role: InputRole): FastaParserKind {
    return fastaParser(programById(this.state.get().program), role);
  }

  private reindex(role: InputRole, key: string, parser: FastaParserKind): void {
    const source = this.role(role).sources.find((s) => s.key === key);
    if (source === undefined) return;
    // A pending index of the pasted text would read it again after this one (and drop the exclusions made since).
    if (key === PASTE_KEY) this.cancelTask(`paste:${role}`);
    // The table of the kind that the source leaves is kept, with its exclusions, for a program that reads it again.
    const tables: Partial<Record<FastaParserKind, SourceTable>> = { ...source.tables };
    if (source.base !== undefined) {
      tables[source.parser] = { base: source.base, excluded: source.excluded, ...(source.indexMs === undefined ? {} : { indexMs: source.indexMs }) };
    }
    const kept = tables[parser];
    delete tables[parser];
    this.updateSource(role, key, (s) => {
      const next: DraftSource = { ...s, parser, tables, ...(kept === undefined ? { status: 'indexing', excluded: [] } : { status: 'ready', base: kept.base, excluded: kept.excluded }) };
      if (kept === undefined) delete (next as { base?: unknown }).base;
      delete (next as { check?: unknown }).check;
      delete (next as { error?: unknown }).error;
      delete (next as { indexMs?: unknown }).indexMs;
      if (kept?.indexMs !== undefined) return { ...next, indexMs: kept.indexMs };
      return next;
    });
    if (kept !== undefined) {
      // An index of the kind that was left, still running, must not replace this table.
      this.nextGeneration(`${role}:${key}`);
      this.scheduleCheck(role, key, 0);
      return;
    }
    this.track(this.index(role, key, source.file, parser, source.sourceId));
  }

  private async index(role: InputRole, key: string, file: File, parser: FastaParserKind, known?: string): Promise<void> {
    const task = `${role}:${key}`;
    const generation = this.nextGeneration(task);
    const current = () => this.generations.get(task) === generation;
    // The index holds its source until the form shows it; a superseded or failed index lets it go.
    const holder = `index:${task}:${generation}`;
    const started = this.deps.now();
    try {
      // The source is added once: its other kind's index reads the same data (no second inflation of a .gz).
      const ref = known === undefined ? await this.deps.data.addSource(file) : undefined;
      const sourceId = known ?? ref!.sourceId;
      this.refs.set(holder, [sourceId]);
      const head =
        ref === undefined
          ? undefined
          : await this.deps.data.previewSource(sourceId, HEAD_BYTES).then(headText, () => undefined);
      const base = await this.deps.data.indexSource(sourceId, parser);
      if (!current()) return;
      this.updateSource(role, key, (source) => ({
        ...source,
        // A file named .gz is inflated by addSource: the card shows its plain name and inflated size.
        ...(ref === undefined ? {} : { name: ref.name, size: ref.size }),
        sourceId,
        status: 'ready',
        base,
        excluded: [],
        check: { state: 'pending' },
        indexMs: this.deps.now() - started,
        ...(head === undefined ? {} : { head }),
      }));
      this.noteRevision(sourceId, base.revisionId);
      this.scheduleCheck(role, key, 0);
    } catch (error) {
      if (!current()) return;
      this.updateSource(role, key, (source) => ({ ...source, status: 'failed', error: messageOf(error) }));
    } finally {
      this.refs.release(holder);
      this.scheduleValidation();
    }
  }

  private scheduleCheck(role: InputRole, key: string, delay?: number): void {
    this.updateSource(role, key, (source) => (source.status === 'ready' ? { ...source, check: { state: 'pending' } } : source));
    this.debounce(`check:${role}:${key}`, () => this.track(this.check(role, key)), delay);
  }

  private async check(role: InputRole, key: string): Promise<void> {
    const task = `check:${role}:${key}`;
    const generation = this.nextGeneration(task);
    const source = this.role(role).sources.find((s) => s.key === key);
    if (source?.status !== 'ready' || source.base === undefined) return;
    if (includedOf(source).length === 0) {
      this.updateSource(role, key, (s) => ({ ...s, check: { state: 'ok', records: 0 } }));
      return;
    }
    const program = this.state.get().program;
    const started = this.deps.now();
    let check: SourceCheck;
    try {
      const revisionId = await this.revisionFor(source);
      // A source removed meanwhile is not checked (nor cached).
      if (source.sourceId === undefined || !this.formSources.has(source.sourceId)) return;
      const key = `${program}|${role}|${revisionId}`;
      let verdict = this.verdicts.get(key);
      if (verdict === undefined) {
        verdict = this.deps.data.checkInput(program, role, [revisionId]);
        this.verdicts.set(key, verdict);
        // A check that could not run is tried again next time.
        verdict.catch(() => this.verdicts.delete(key));
      }
      check = toSourceCheck(await verdict, source);
    } catch (error) {
      check = { state: 'error', message: messageOf(error) };
    }
    if (this.generations.get(task) !== generation || this.state.get().program !== program) return;
    const checkMs = this.deps.now() - started;
    this.updateSource(role, key, (s) =>
      s.excluded === source.excluded && s.base === source.base ? { ...s, check, checkMs } : s,
    );
  }

  /** The revision of a source's selection: its record table, or a revision without the excluded records. */
  private revisionFor(source: DraftSource): Promise<string> {
    const base = source.base;
    if (base === undefined) return Promise.reject(new Error(`${source.name} has no record table`));
    if (source.excluded.length === 0) return Promise.resolve(base.revisionId);
    const key = `${base.revisionId}|${source.excluded.join(',')}`;
    let revision = this.revisions.get(key);
    if (revision === undefined) {
      revision = this.deps.data.reviseDataset(base.revisionId, source.excluded).then((r) => {
        this.noteRevision(base.sourceId, r.revisionId);
        return r.revisionId;
      });
      this.revisions.set(key, revision);
      revision.catch(() => this.revisions.delete(key));
    }
    return revision;
  }

  private scheduleValidation(): void {
    // A validation of an older argv that is still running is not shown.
    this.nextGeneration('validate');
    this.debounce('validate', () => this.track(this.validate()));
  }

  private async validate(): Promise<void> {
    const generation = this.generations.get('validate');
    const { program } = this.state.get();
    let argv: readonly string[];
    try {
      argv = buildArgv({
        program,
        queryName: this.inputName('query'),
        subjectName: this.inputName('subject'),
        parameters: this.parameters(),
        extra: this.extraWords(),
      });
    } catch (error) {
      this.set({ validation: { state: 'invalid', message: messageOf(error) } });
      return;
    }
    this.set({ validation: { state: 'checking' } });
    let validation: Validation;
    try {
      const result = await this.deps.engine.validate(argv);
      validation = result.ok ? { state: 'ok' } : { state: 'invalid', message: result.message };
    } catch (error) {
      validation = { state: 'invalid', message: messageOf(error) };
    }
    if (this.generations.get('validate') !== generation) return;
    this.validArgv = validation.state === 'ok' ? argv : undefined;
    this.set({ validation });
    this.prepareEngine();
  }

  /**
   * Tells the engine that a search like the draft's may come (EngineGateway.prepare), once
   * its conditions are valid and both roles have indexed inputs, so that the engine can
   * start its runtime and compile the code of the search before "Run LOSAT". The size is
   * that of the largest run.
   */
  private prepareEngine(): void {
    const { threads, validation } = this.state.get();
    if (this.deps.engine.prepare === undefined || validation.state !== 'ok') return;
    let inputBytes = 0;
    for (const role of ROLES) {
      const groups = this.inputGroups(role);
      if (groups.length === 0) return;
      inputBytes += Math.max(...groups.map((group) => group.sources.reduce((sum, source) => sum + source.size, 0)));
    }
    this.deps.engine.prepare({ requestedThreads: threads, inputBytes, ...(this.validArgv === undefined ? {} : { argv: this.validArgv }) });
  }

  // --- helpers ----------------------------------------------------------------------------

  private debounce(name: string, work: () => void, delay = this.deps.debounceMs ?? 300): void {
    const timer = this.timers.get(name);
    if (timer !== undefined) clearTimeout(timer);
    this.workers.set(name, work);
    this.timers.set(
      name,
      setTimeout(() => {
        this.timers.delete(name);
        this.runTimer(name);
      }, delay),
    );
  }

  private runTimer(name: string): void {
    const work = this.workers.get(name);
    this.workers.delete(name);
    work?.();
  }

  private cancelTask(name: string): void {
    const timer = this.timers.get(name);
    if (timer !== undefined) clearTimeout(timer);
    this.timers.delete(name);
    this.workers.delete(name);
    this.nextGeneration(name);
  }

  private nextGeneration(task: string): number {
    const generation = (this.generations.get(task) ?? 0) + 1;
    this.generations.set(task, generation);
    return generation;
  }

  private track(work: Promise<unknown>): void {
    this.pending.add(work);
    void work.finally(() => this.pending.delete(work));
  }

  private role(role: InputRole): RoleDraft {
    return this.state.get()[role];
  }

  private updateRole(role: InputRole, change: (draft: RoleDraft) => RoleDraft): void {
    this.set(role === 'query' ? { query: change(this.role(role)) } : { subject: change(this.role(role)) });
  }

  private updateSource(role: InputRole, key: string, change: (source: DraftSource) => DraftSource): void {
    this.updateRole(role, (draft) => ({
      ...draft,
      sources: draft.sources.map((source) => (source.key === key ? change(source) : source)),
    }));
  }

  /** Changes the state; a key given as undefined is removed. */
  private set(change: { readonly [K in keyof DraftState]?: DraftState[K] | undefined }): void {
    const next = { ...this.state.get(), ...change };
    for (const key of Object.keys(change) as Array<keyof DraftState>) {
      if (change[key] === undefined) delete next[key];
    }
    this.state.set(next as DraftState);
    this.holdFormSources(next as DraftState);
  }

  /** Holds the sources that the form shows; one it no longer shows goes, with the draft's caches of it. */
  private holdFormSources(state: DraftState): void {
    const shown = new Set<string>();
    for (const role of ROLES) {
      for (const source of state[role].sources) if (source.sourceId !== undefined) shown.add(source.sourceId);
    }
    if (shown.size === this.formSources.size && [...shown].every((sourceId) => this.formSources.has(sourceId))) return;
    for (const sourceId of this.formSources) if (!shown.has(sourceId)) this.forget(sourceId);
    this.formSources = shown;
    this.refs.set('draft', shown);
  }

  private noteRevision(sourceId: string, revisionId: string): void {
    // A source that the form no longer shows keeps nothing in the caches.
    if (!this.formSources.has(sourceId)) return;
    let revisions = this.sourceRevisions.get(sourceId);
    if (revisions === undefined) this.sourceRevisions.set(sourceId, (revisions = new Set()));
    revisions.add(revisionId);
  }

  /** Drops the verdicts and the selections of a source's revisions. */
  private forget(sourceId: string): void {
    const revisions = this.sourceRevisions.get(sourceId);
    this.sourceRevisions.delete(sourceId);
    if (revisions === undefined) return;
    for (const key of [...this.verdicts.keys()]) {
      if (revisions.has(key.slice(key.lastIndexOf('|') + 1))) this.verdicts.delete(key);
    }
    for (const key of [...this.revisions.keys()]) {
      if (revisions.has(key.slice(0, key.indexOf('|')))) this.revisions.delete(key);
    }
  }
}

/** The draft without its fetch's state, and with `groups` also without the groups that have come. */
function withoutFetch(draft: RoleDraft, groups = false): RoleDraft {
  const next = { ...draft };
  delete (next as { fetch?: unknown }).fetch;
  if (groups) delete (next as { fetched?: unknown }).fetched;
  return next;
}

/** Why a fetch brought no record, for the line under the box. */
function fetchProblem(result: Exclude<FetchResult, { ok: true }>, source: FetchSourceId, kind: SequenceKind, accession: string): string {
  const label = sourceLabel(source);
  switch (result.reason) {
    case 'not-found':
      return `${label} has no ${kind} record ${accession}.`;
    case 'too-large':
      return `${accession} makes the fetched text larger than ${FETCH_BYTE_LIMIT / 1_000_000} MB: download it from ${label} and open it with "Open FASTA files…".`;
    case 'unreachable':
      return `${label} could not be reached (${result.detail}).`;
    case 'refused':
      return result.status === 429
        ? `${label} asks to wait (too many requests): fetch again in a moment.`
        : `${label} answered HTTP ${result.status} for ${accession}.`;
  }
}

/** The included records of a ready source. */
export function includedOf(source: DraftSource): readonly DatasetRecord[] {
  return source.base === undefined ? [] : includedRecords({ ...source.base, excluded: source.excluded });
}

function nameOf(source: DraftSource, role: InputRole): string {
  return source.origin === 'paste' ? PASTED_NAMES[role] : source.name;
}

function describeSource(source: DraftSource): string {
  return source.origin === 'paste' ? 'pasted' : `file ${source.name}`;
}

/** The engine's verdict, with the record its message names mapped to the record table. */
function toSourceCheck(check: InputCheck, source: DraftSource): SourceCheck {
  if (check.ok) return { state: 'ok', records: check.records.length };
  const match = RECORD_IN_MESSAGE.exec(check.message);
  const record = match === null ? undefined : includedOf(source)[Number(match[1]) - 1];
  return { state: 'refused', message: check.message, ...(record === undefined ? {} : { record: record.index }) };
}

/** The first lines of a source, decoded leniently and shortened, for its preview. */
function headText(bytes: Uint8Array): string {
  const text = new TextDecoder('utf-8', { fatal: false }).decode(bytes);
  const lines = text.split(/\r?\n/).slice(0, HEAD_LINES);
  // The last line may be cut by the read, unless the source is shorter than the read.
  return lines.map((line) => (line.length > HEAD_LINE_CHARS ? `${line.slice(0, HEAD_LINE_CHARS)}…` : line)).join('\n');
}

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
