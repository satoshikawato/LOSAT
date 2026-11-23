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
  REGION_FLAG,
  WINDOW_MASKER_FLAG,
  WINDOW_MASKER_TAXID_FLAG,
  type FastaParserKind,
  type FastaReader,
  type InputRole,
  type PrintedKind,
  type ProgramId,
  type RegistryFile,
  type SequenceKind,
} from 'losat/domain';
import { isFileReadError, isUnavailable, type EngineGateway, type InputCheck, type ProgramDescription, type StartInfo } from 'losat/ports';
import {
  accessionLines,
  accessionsOnALine,
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
import { MOST_LEFT_OUT_IDS } from '../domain/input-parts';
import type { InputPart } from '../domain/run';
import { buildArgv, COMBINED_NAMES, INCLUDED_PREFIX, PASTED_NAMES, splitOptions } from '../domain/argv';
import { endsInput, includedRecords, type DatasetRecord, type DatasetRevision } from '../domain/dataset';
import { parseTaxid, type CountsOrigin } from '../domain/counts';
import type { ExampleId } from '../domain/examples';
import { inputLead, type InputLead } from '../domain/input-format';
import { describedSections, formControlOf, formParameters, setFields, type FieldValue, type FormValues } from '../domain/parameters';
import { fastaParser, programById, PROGRAMS, sequenceKind } from '../domain/programs';
import { regionProblem, regionValue, type RegionText } from '../domain/region';
import { withUnit } from '../domain/units';
import type { DataGateway, SourceRef } from '../ports/data';
import { failureText, withRemedy } from './availability';
import type { AccessionFetcher, FetchResult } from '../ports/fetch';
import { registryName, unreadableConffile } from './conffile';
import type { EnqueueAllResult, SearchRequest } from './coordinator';
import { SourceRace, type RaceStanding } from './source-race';
import { SourceRefs } from './source-refs';
import { Store } from './store';

export type InputMode = 'combined' | 'separate';

/** The engine's verdict on a source (losat/ports: input-check.ts). */
export type SourceCheck =
  | { readonly state: 'pending' }
  /**
   * `readerError`: the reader's error that ends the query's records, the command line's text; a run
   * reads the records before it and ends with it, as the command line does, so it is a warning, not
   * a refusal (line WQA's D10, 2026-11-19).
   */
  | { readonly state: 'ok'; readonly records: number; readonly readerError?: string }
  /** `record` is the index, in the record table, of the record that the message names. */
  | { readonly state: 'refused'; readonly message: string; readonly record?: number }
  /**
   * The check could not run (not a verdict on the input); `unavailable` when a part of LOSAT was (UnavailableError,
   * line WQA's sweep, D1).
   */
  | { readonly state: 'error'; readonly message: string; readonly unavailable?: true };

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
  /**
   * The index failed because a part of LOSAT could not be loaded or stopped (UnavailableError), not because of the
   * input: the card does not blame the input, and the page says what failed (line WQA's sweep, D1).
   */
  readonly unavailable?: true;
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
  /** What the first bytes say besides FASTA: a byte order mark, or another file format (WA-IN-04, WA-IN-05). */
  readonly lead?: InputLead;
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
  | { readonly state: 'invalid'; readonly message: string }
  /** The engine could not check the options: a part of LOSAT could not be loaded or stopped (line WQA's sweep, D1). */
  | { readonly state: 'unavailable'; readonly message: string };

/** The counts of a taxonomy ID that the draft can search with. */
export interface SavedCountsFile {
  readonly file: File;
  readonly origin: CountsOrigin;
}

export interface DraftState {
  readonly program: ProgramId;
  readonly query: RoleDraft;
  readonly subject: RoleDraft;
  /** The form's values of each program: what the user set, nothing else. */
  readonly values: Readonly<Record<ProgramId, FormValues>>;
  /**
   * The files chosen for the options of the kind `file` (BLASTN's WindowMasker counts), by program
   * and flag; the field's value holds the file's name. They stay in the working session, as the
   * sources do.
   */
  readonly files: Readonly<Record<ProgramId, Readonly<Record<string, File>>>>;
  /**
   * The WindowMasker counts that the browser keeps for the taxonomy IDs (line WW2), by the ID in
   * decimal: the file, and where it came from. A search with `-window_masker_taxid` reads the file
   * of its ID (CountsLoader hands them over: setCounts).
   */
  readonly counts: Readonly<Record<string, SavedCountsFile>>;
  /** The "Additional options" of each program, as typed: other options, as on a command line (`splitOptions`). */
  readonly extra: Readonly<Record<ProgramId, string>>;
  /**
   * The file chosen for the `-conffile` of each program's Additional options (line RG E1w-3b): the
   * validation and the run read its bytes under the name that the argv gives. It stays in the
   * working session, as the sources do.
   */
  readonly conffiles: Readonly<Partial<Record<ProgramId, File>>>;
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
  /** The describe failed because a part of LOSAT could not be loaded or stopped (UnavailableError): the page says so. */
  readonly descriptionUnavailable?: true;
  /** The engine's `validate` of the argv that the draft makes now. */
  readonly validation: Validation;
  /**
   * Whether the run of the draft's argv reads its inputs with NCBI's data loaders on, as the latest
   * validation that the engine answered says (`[BLAST] DATA_LOADERS` of the argv's registry; on before
   * one): the engine's readers of the sources' index and check have them so (line WQA's D16, 2026-11-20).
   */
  readonly dataLoaders: boolean;
  /**
   * What the engine's start read of the argv of the latest validation (`EngineGateway.start`): the
   * name of its `-conffile`, for which the form asks for a file; undefined before it, or when the
   * engine could not start the argv.
   */
  readonly start?: StartInfo;
  readonly submitting: boolean;
  readonly message?: { readonly kind: 'info' | 'error'; readonly text: string };
}

export interface DraftDeps {
  readonly engine: Pick<EngineGateway, 'describe' | 'start' | 'validate' | 'prepare'>;
  readonly data: Pick<
    DataGateway,
    'addSource' | 'indexSource' | 'reviseDataset' | 'checkInput' | 'previewSource' | 'dropSource'
  >;
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
  /** Who uses the data layer's sources: one that the form, an index and "Run LOSAT" no longer use is dropped. */
  private readonly refs: SourceRefs;
  /** The data layer's sources that the form shows. */
  private formSources: ReadonlySet<string> = new Set();
  /** The revisions of each of those sources, whose cached verdicts and selections go with the source. */
  private readonly sourceRevisions = new Map<string, Set<string>>();
  /** The latest request of each asynchronous task; older results are dropped. */
  private readonly generations = new Map<string, number>();
  private readonly timers = new Map<string, ReturnType<typeof setTimeout>>();
  private readonly workers = new Map<string, () => void>();
  private readonly pending = new Set<Promise<unknown>>();
  /**
   * What askAgain asked the engine again that has not answered yet. While it is on its way, the question reads as
   * asked (the options are described, the input indexing, the validation checking), so a later askAgain, which asks
   * only what failed, waits for it: an Example pressed again then asked nothing and was refused with "Reload the page"
   * also once the answers came (the fifth review of line WQA's D1 and D2, M1, 2026-11-19).
   */
  private readonly asking = new Set<Promise<unknown>>();
  /** The running fetch of each role's accessions. */
  private readonly fetches = new Map<InputRole, AbortController>();
  private readonly race: SourceRace;
  /** The argv of the latest validation that the engine accepted (the prepare hint's argv). */
  private validArgv: readonly string[] | undefined;
  /** The validation that runs, or the last one: the next waits for it (validateInTurn). */
  private validationTurn: Promise<void> = Promise.resolve();
  private nextKey = 1;

  constructor(private readonly deps: DraftDeps) {
    this.refs = new SourceRefs((sourceId) => void deps.data.dropSource(sourceId).catch(() => undefined));
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
      files: Object.fromEntries(PROGRAMS.map((program) => [program.id, Object.freeze({})])) as Record<ProgramId, Readonly<Record<string, File>>>,
      counts: Object.freeze({}),
      extra: Object.fromEntries(PROGRAMS.map((program) => [program.id, ''])) as Record<ProgramId, string>,
      conffiles: Object.freeze({}),
      title: '',
      threads: 'auto',
      fetchSource: 'fastest',
      race: {},
      validation: { state: 'idle' },
      dataLoaders: true,
      submitting: false,
    });
    void this.loadDescription('blastn');
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
    // The engine's start of the other program's argv is not this program's (its -conffile row goes until the validation).
    this.set({ program, description, descriptionError: undefined, descriptionUnavailable: undefined, message: undefined, start: undefined });
    void this.loadDescription(program);
    for (const role of ROLES) {
      const parser = fastaParser(next, role);
      for (const source of this.role(role).sources) {
        if (source.parser !== parser || source.status === 'failed') {
          // A program that reads the role with another FASTA reader kind needs another record table (plan TD-8); a
          // refusal of the reader names its program (line WQA's D13, 2026-11-20), so a source that failed is read again.
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
   * Chooses the file of an option of the kind `file`: the option's value becomes the file's name,
   * which the argv gives the option, and the run reads the file under it. Without a file the
   * option is cleared (the box is unchecked).
   */
  setFile(flag: string, file: File | undefined): void {
    const { program, files } = this.state.get();
    const next = { ...files[program] };
    if (file === undefined) delete next[flag];
    else next[flag] = file;
    this.set({ files: { ...files, [program]: Object.freeze(next) } });
    this.setField(flag, file === undefined ? '' : file.name);
  }

  /**
   * Chooses how the row "Species-specific repeats filter for:" names the counts (BLASTN): not at
   * all (the box is unchecked), by a taxonomy ID, or by a file. The box checked waits for the ID or
   * the file (the search cannot start without it), and a chosen file is let go whichever way the
   * row changes.
   */
  setRepeats(mode: 'off' | 'taxid' | 'file'): void {
    const { program, files } = this.state.get();
    const next = { ...files[program] };
    delete next[WINDOW_MASKER_FLAG];
    this.set({ files: { ...files, [program]: Object.freeze(next) } });
    this.setFields([
      [WINDOW_MASKER_FLAG, mode === 'file' ? true : ''],
      [WINDOW_MASKER_TAXID_FLAG, mode === 'taxid' ? true : ''],
    ]);
  }

  /** Types the taxonomy ID of the row (an empty box stays checked and waits for it). */
  setTaxid(text: string): void {
    this.setField(WINDOW_MASKER_TAXID_FLAG, text.trim() === '' ? true : text);
  }

  /**
   * Hands the draft the whole counts file that the browser keeps for a taxonomy ID, or takes it
   * back (undefined: it was deleted, or is not whole). A search with that ID reads this file.
   */
  setCounts(taxid: number, saved: SavedCountsFile | undefined): void {
    const counts = { ...this.state.get().counts };
    if (saved === undefined) delete counts[taxid];
    else counts[taxid] = saved;
    this.set({ counts: Object.freeze(counts) });
  }

  /**
   * Clears the values of some fields of the program ("Restore default search parameters"):
   * the engine's defaults apply to them again, and a chosen file is let go.
   */
  resetFields(flags: readonly string[]): void {
    const { program, values, files } = this.state.get();
    const next: Record<string, FieldValue> = { ...values[program] };
    const nextFiles: Record<string, File> = { ...files[program] };
    for (const flag of flags) {
      delete next[flag];
      delete nextFiles[flag];
    }
    this.set({
      values: { ...values, [program]: Object.freeze(next) },
      files: { ...files, [program]: Object.freeze(nextFiles) },
      message: undefined,
    });
    this.scheduleValidation();
  }

  /** The program's "Additional options": words of other options, after the form's in the argv. */
  setExtra(text: string): void {
    const { program, extra } = this.state.get();
    this.set({ extra: { ...extra, [program]: text }, message: undefined });
    this.scheduleValidation();
  }

  /**
   * Chooses the file of the `-conffile` that the program's Additional options give (undefined lets
   * it go): the validation and the run read it under the argv's name, whatever the file is called.
   * The page holds the File, not its bytes: a validation reads it on an instance of the engine of its
   * own, and "Run LOSAT" reads it once for the runs (Coordinator.enqueueAll).
   */
  setConffile(file: File | undefined): void {
    const { program, conffiles } = this.state.get();
    const next = { ...conffiles };
    if (file === undefined) delete next[program];
    else next[program] = file;
    this.set({ conffiles: Object.freeze(next), message: undefined });
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

  /**
   * Puts ">" before the first line of a pasted text that the engine reads as a Seq-id, which NCBI BLAST+
   * would fetch, so that the line is a defline (WA2-OP-04). The user asks for it; the first line is
   * the only change, as the reason says ("start the input with a '>' defline").
   */
  markFirstLineAsDefline(role: InputRole): void {
    this.setPaste(role, `>${this.role(role).paste.trimStart()}`);
  }

  /**
   * Puts each accession of a pasted text of accessions on its own line (WA2-OP-04), where Fetch
   * takes them; the user asks for it.
   */
  splitAccessions(role: InputRole): void {
    const accessions = accessionsOnALine(this.role(role).paste);
    if (accessions !== undefined) this.setPaste(role, accessions.join('\n'));
  }

  /** Removes the byte order mark (U+FEFF) that a pasted text starts with, which NCBI's reader cannot read (WA-IN-05). */
  removeByteOrderMark(role: InputRole): void {
    const paste = this.role(role).paste;
    if (paste.startsWith('\uFEFF')) this.setPaste(role, paste.slice(1));
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
    this.recheckReaderErrors(role);
    this.scheduleValidation();
  }

  setMode(role: InputRole, mode: InputMode): void {
    this.updateRole(role, (draft) => ({ ...draft, mode }));
    this.set({ message: undefined });
    this.recheckReaderErrors(role);
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
    this.recheckReaderErrors(role, key);
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
    const sources = (r: InputRole) => this.role(r).sources.filter((source) => source.status !== 'failed');
    const [query, subject] = [sources('query'), sources('subject')];
    return runNames(
      { name: nameOf(query, 'query'), made: !isFileAsGiven(query) },
      { name: nameOf(subject, 'subject'), made: !isFileAsGiven(subject) },
    )[role];
  }

  /**
   * The runs that "Run LOSAT" would make now: the runs that `submit` starts (one `inputGroups` for both). A draft whose argv
   * prints the command line's help or version starts the one run that prints it, whatever the inputs (line WQA's D20, S3-2,
   * 2026-11-23).
   */
  runCount(): number {
    if (this.printsKind() !== null) return 1;
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
        result = { ok: false, message: isUnavailable(error) ? failureText(error) : `The search could not be started: ${messageOf(error)}` };
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
    // Run never fetches, nor waits for a fetch that is running (fetchAccessions). A run that prints the command line's
    // help or version reads no input (the engine's start says so, `prints`; line WQA's D14, 2026-11-20): it is read once the
    // validation of the draft's argv has answered (not `idle`, which waits for a fetch too), so that a word typed just before
    // Run, or taken out, decides as the engine's answer for the argv says (line WQA's D20 round 4, 2026-11-23).
    await this.validationSettled();
    for (const role of ROLES) {
      const problem = this.printsText() ? undefined : this.accessionProblem(role);
      if (problem !== undefined) return { ok: false, message: problem };
    }
    // What the engine could not answer while a part of LOSAT could not be loaded is asked again first: Run LOSAT was
    // refused with "Reload the page" and asked nothing, also once the engine file downloaded again (the third review of
    // line WQA's D1 and D2, M1, 2026-11-19). A question still in flight at the press, which fails only after it, was not
    // among those, and the next press then succeeded (its fourth review, L1): such a question is asked once more.
    const failed = this.unavailableQuestions();
    void this.askAgain();
    await this.idle();
    if ([...this.unavailableQuestions()].some((question) => !failed.has(question))) {
      void this.askAgain();
      await this.idle();
    }
    const { program, threads, title } = this.state.get();
    if (this.printsText()) return this.enqueuePrinting();
    for (const role of ROLES) {
      const problem = this.inputProblem(role);
      if (problem !== undefined) return { ok: false, message: problem };
    }
    const fileProblem = this.fileProblem();
    if (fileProblem !== undefined) return { ok: false, message: fileProblem };
    // Everything from one state, before any wait.
    const groups = { query: this.inputGroups('query'), subject: this.inputGroups('subject') };
    // Their sources stay until the runs have their inputs, also when the form lets them go meanwhile.
    this.refs.set(
      'submit',
      [...groups.query, ...groups.subject].flatMap((sources) => sources.map((source) => source.sourceId!)),
    );
    const parameters = this.parameters();
    // The file of -window_masker_db, when the argv has it (a file chosen and then unchecked is not in the form).
    let windowMasker = parameters.some(([flag]) => flag === WINDOW_MASKER_FLAG) ? this.state.get().files[program][WINDOW_MASKER_FLAG] : undefined;
    // Or the file kept for the taxonomy ID of -window_masker_taxid (line WW2).
    const taxid = parseTaxid(String(parameters.find(([flag]) => flag === WINDOW_MASKER_TAXID_FLAG)?.[1] ?? ''));
    const kept = taxid === undefined ? undefined : this.state.get().counts[taxid];
    if (kept !== undefined) windowMasker = kept.file;
    let extra: readonly string[];
    try {
      extra = this.extraWords();
    } catch (error) {
      return { ok: false, message: messageOf(error) };
    }
    // The file of the argv's -conffile, which the engine's start named at the validation (idle).
    const state = this.state.get();
    const conffile = registryName(state.start) === undefined ? undefined : state.conffiles[program];
    const resolved = {
      query: await Promise.all(groups.query.map((sources) => this.resolveGroup('query', sources))),
      subject: await Promise.all(groups.subject.map((sources) => this.resolveGroup('subject', sources))),
    };
    const requests: SearchRequest[] = [];
    for (const query of resolved.query) {
      for (const subject of resolved.subject) {
        // The two names of one run differ, unless they are the user's own (runNames).
        const names = runNames({ name: query.name, made: query.from !== undefined }, { name: subject.name, made: subject.from !== undefined });
        requests.push({
          program,
          query: { dataset: { ...query, name: names.query } },
          subject: { dataset: { ...subject, name: names.subject } },
          parameters,
          ...(extra.length === 0 ? {} : { extra }),
          requestedThreads: threads,
          ...(windowMasker === undefined ? {} : { windowMasker }),
          ...(kept === undefined ? {} : { windowMaskerOrigin: kept.origin }),
          ...(conffile === undefined ? {} : { conffile }),
          ...(title.trim() === '' ? {} : { title }),
        });
      }
    }
    return this.deps.enqueueAll(requests);
  }

  /**
   * Whether the run of the draft's argv prints the command line's help or version in place of a search, as the engine's
   * start says (`StartInfo.prints`; line WQA's D14, 2026-11-20): it needs no inputs, files or fetches, as the command
   * line reads none for it.
   */
  private printsText(): boolean {
    return this.printsKind() !== null;
  }

  /**
   * What the run of the draft's argv prints in place of a search (`printsText`): help or version, or null for a search, as the
   * engine's latest answer says (`start.prints`, validation ok). The page reads it for the line under the button and the run
   * count, and Run LOSAT reads it once the validation is idle. Between an edit and its answer the older answer stands (the
   * app cannot know before the engine answers that a word makes the argv print or not); it is not set aside at every edit,
   * which would flip the line on each edit of an ordinary search (line WQA's D20 round 4, 2026-11-23).
   */
  printsKind(): PrintedKind | null {
    const { start, validation } = this.state.get();
    return validation.state === 'ok' ? (start?.prints ?? null) : null;
  }

  /**
   * The one run of an argv that prints the command line's help or version (`printsText`): without inputs (the
   * coordinator gives it none), with the registry file of its -conffile, which the command line reads before the help.
   */
  private enqueuePrinting(): Promise<EnqueueAllResult> {
    const { program, threads, title } = this.state.get();
    let extra: readonly string[];
    try {
      extra = this.extraWords();
    } catch (error) {
      return Promise.resolve({ ok: false, message: messageOf(error) });
    }
    const state = this.state.get();
    const conffile = registryName(state.start) === undefined ? undefined : state.conffiles[program];
    return this.deps.enqueueAll([
      {
        program,
        query: { text: '' },
        subject: { text: '' },
        parameters: this.parameters(),
        ...(extra.length === 0 ? {} : { extra }),
        requestedThreads: threads,
        ...(conffile === undefined ? {} : { conffile }),
        ...(title.trim() === '' ? {} : { title }),
      },
    ]);
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
    if (this.printsText()) return undefined;
    for (const role of ROLES) {
      const problem = this.problemOf(role);
      if (problem !== undefined) return problem.short;
    }
    return this.fileProblem();
  }

  /**
   * A box of an option of the kind `file` that is checked while no file is chosen (BLASTN's
   * WindowMasker counts), or while the taxonomy ID of its row is not typed or its counts are not
   * downloaded: the sentence that asks for it; undefined when there is none.
   */
  private fileProblem(): string | undefined {
    const { program, values, description, counts } = this.state.get();
    for (const section of describedSections(programById(program), description?.parameters)) {
      for (const field of section.fields) {
        if (field.kind !== 'file') continue;
        if (values[program][field.flag] === true) return field.file?.missing ?? `Choose the file for ${field.label}`;
        // The taxonomy ID of the same row: typed, and its counts downloaded.
        const taxid = field.file?.taxid;
        const typed = taxid === undefined ? undefined : values[program][taxid.flag];
        if (taxid === undefined || (typed !== true && (typeof typed !== 'string' || typed.trim() === ''))) continue;
        const id = typeof typed === 'string' ? parseTaxid(typed) : undefined;
        if (id === undefined) return taxid.missing;
        if (counts[id] === undefined) return `Download the counts for ${id}.`;
      }
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
      if (source.status === 'failed' && source.unavailable === true) {
        // LOSAT's own failure, not the input's (line WQA's sweep, D1): what the page says.
        const message = withRemedy(source.error ?? '');
        return { message, short: message };
      }
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
      // would join the record before it (NCBI's reader): only the first can start so. So can
      // neither a source whose reader's error is on its first line, read without a defline (a
      // byte order mark, a line of other text): the record before it would read that line
      // (review M2, 2026-11-20).
      const late = draft.sources.find((source) => cannotFollow(draft.sources, source));
      if (late !== undefined) {
        const message = `${label} (${describeSource(late)}) ${cannotFollowText(late)}`;
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
  private inputGroups(role: InputRole): ReadonlyArray<readonly DraftSource[]> {
    const draft = this.role(role);
    const usable = draft.sources.filter((source) => source.status === 'ready');
    if (draft.mode === 'separate' && usable.length > 1) {
      // A source whose records are all excluded has nothing to search alone.
      return usable
        .filter((source) => includedOf(source).length > 0 || (source.base?.records.length ?? 0) === 0)
        .map((source) => [source]);
    }
    return usable.length === 0 ? [] : [usable];
  }

  /**
   * An input of a run, from its sources: the name that its file has, the revisions and what it is made of. What it is made of
   * is made here, when the run is resolved, not when the runs are counted (runCount): it reads the record tables.
   */
  private async resolveGroup(role: InputRole, sources: readonly DraftSource[]) {
    const revisionIds = await Promise.all(sources.map((source) => this.revisionFor(source)));
    const from = fromOf(sources);
    return { name: nameOf(sources, role), revisionIds, ...(from === undefined ? {} : { from }) };
  }

  // --- asynchronous work ------------------------------------------------------------------

  /** Asks the engine for the program's options (once, while it answers); settles once the form has the answer. */
  private loadDescription(program: ProgramId): Promise<void> {
    let description = this.descriptions.get(program);
    if (description === undefined) {
      description = this.deps.engine.describe(program);
      this.descriptions.set(program, description);
      description.catch(() => this.descriptions.delete(program));
    }
    const answered = description.then(
      (value) => {
        this.described.set(program, value);
        if (this.state.get().program === program) this.set({ description: value, descriptionError: undefined, descriptionUnavailable: undefined });
        this.scheduleValidation();
        void this.askAgain();
      },
      (error: unknown) => {
        if (this.state.get().program === program) {
          this.set({ descriptionError: messageOf(error), descriptionUnavailable: isUnavailable(error) ? true : undefined });
        }
      },
    );
    this.track(answered);
    return answered;
  }

  /**
   * Asks the engine again what the form could not ask it while a part of LOSAT could not be loaded: the options, the
   * inputs that could not be read or checked, the validation. An engine file that failed to download is downloaded
   * again by the next call, so one network error at load left "Not read." on the options and the page's "Reload the
   * page" for the whole session, also after searches completed (the second review of line WQA's D1 and D2, M1,
   * 2026-11-19). It is asked when the engine answers the options, an index or a check, when Run LOSAT is pressed, and
   * before an Example is loaded (Availability.recheck). A question is pending until it answers; one that fails asks
   * nothing again. Settles once what it asked has answered (and a validation that it asked has run), and what an
   * earlier call asked that has not answered yet (asking), not the rest of the form's work: an Example waited for all
   * of it, a large index or a fetch too (the fourth review, M1).
   */
  askAgain(): Promise<void> {
    const state = this.state.get();
    const asked: Array<Promise<unknown>> = [];
    const options = this.askOptionsAgain();
    if (options !== undefined) asked.push(options);
    for (const role of ROLES) {
      for (const source of state[role].sources) {
        if (source.status === 'failed' && source.unavailable === true) asked.push(this.reindex(role, source.key, this.parserOf(role)));
        else if (source.status === 'ready' && source.check?.state === 'error' && source.check.unavailable === true) {
          asked.push(this.checkNow(role, source.key));
        }
      }
    }
    const validation = state.validation.state === 'unavailable';
    if (validation) this.scheduleValidation();
    // The answers schedule a validation of their own; one that was unavailable has answered once none is left.
    if (asked.length > 0 || validation) this.keepAsking(Promise.allSettled(asked).then(() => (validation ? this.validationSettled() : undefined)));
    return Promise.allSettled([...this.asking]).then(() => undefined);
  }

  /** Keeps a question asked again among those on their way (asking) until it answers. */
  private keepAsking<T>(work: Promise<T>): Promise<T> {
    this.asking.add(work);
    const done = () => this.asking.delete(work);
    work.then(done, done);
    return work;
  }

  /** What the engine could not answer because a part of LOSAT could not be loaded, by question (unavailableOf). */
  private unavailableQuestions(): Set<string> {
    const state = this.state.get();
    const questions = new Set<string>();
    if (state.descriptionUnavailable === true) questions.add('options');
    if (state.validation.state === 'unavailable') questions.add('validation');
    for (const role of ROLES) {
      for (const source of state[role].sources) {
        if (source.status === 'failed' && source.unavailable === true) questions.add(`index:${role}:${source.key}`);
        if (source.check?.state === 'error' && source.check.unavailable === true) questions.add(`check:${role}:${source.key}`);
      }
    }
    return questions;
  }

  /**
   * Asks the engine for the options again if it could not answer them before. A validation that the engine answered
   * asks only this: a validation follows every index, also one that failed, which would otherwise be read again after
   * each validation; the options that answer ask the rest again (askAgain). They are on their way (asking) until they
   * answer, also when a validation asked them.
   */
  private askOptionsAgain(): Promise<void> | undefined {
    const { descriptionUnavailable, program } = this.state.get();
    return descriptionUnavailable === true && !this.descriptions.has(program) ? this.keepAsking(this.loadDescription(program)) : undefined;
  }

  private parserOf(role: InputRole): FastaParserKind {
    return fastaParser(programById(this.state.get().program), role);
  }

  /** The engine's reader of the role's input now: the program's, with the data loaders of the latest validation. */
  private readerOf(role: InputRole): FastaReader {
    const { program, dataLoaders } = this.state.get();
    return { program, role, dataLoaders };
  }

  /**
   * The validation says that the run reads its inputs with the data loaders on or off (`dataLoaders`), as the CLI
   * reads them under the argv's registry (line WQA's D16, 2026-11-20): where that changes, the sources that the reader
   * refused are read again (a first line that may be a Seq-id is data with the loaders off), and the rest checked again.
   */
  private setDataLoaders(dataLoaders: boolean): void {
    if (this.state.get().dataLoaders === dataLoaders) return;
    this.set({ dataLoaders });
    for (const role of ROLES) {
      for (const source of this.role(role).sources) {
        if (source.status === 'failed') this.reindex(role, source.key, this.parserOf(role));
        else this.scheduleCheck(role, source.key, 0);
      }
    }
  }

  /** Reads a source again with `parser`; settles once it is read (or once its table of that kind is back). */
  private reindex(role: InputRole, key: string, parser: FastaParserKind): Promise<void> {
    const source = this.role(role).sources.find((s) => s.key === key);
    if (source === undefined) return Promise.resolve();
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
      delete (next as { unavailable?: unknown }).unavailable;
      delete (next as { indexMs?: unknown }).indexMs;
      if (kept?.indexMs !== undefined) return { ...next, indexMs: kept.indexMs };
      return next;
    });
    if (kept !== undefined) {
      // An index of the kind that was left, still running, must not replace this table.
      this.nextGeneration(`${role}:${key}`);
      this.scheduleCheck(role, key, 0);
      return Promise.resolve();
    }
    const read = this.index(role, key, source.file, parser, source.sourceId);
    this.track(read);
    return read;
  }

  private async index(role: InputRole, key: string, file: File, parser: FastaParserKind, known?: string): Promise<void> {
    const task = `${role}:${key}`;
    const generation = this.nextGeneration(task);
    // The engine's reader of the input now. What it finds (the records, or its refusal) stands only while it is the
    // reader in force: an index that ends after the program or the data loaders changed reads the input again
    // (`readAgain`), whatever the order of the change and the index (line WQA's D13 round 2, M1, 2026-11-20).
    const reader = this.readerOf(role);
    const current = () => this.generations.get(task) === generation;
    let ref: SourceRef | undefined;
    let sourceId: string | undefined;
    let head: string | undefined;
    // The index holds its source until the form shows it; a superseded or failed index lets it go.
    const holder = `index:${task}:${generation}`;
    const started = this.deps.now();
    // Kept for a source that cannot be read too: a byte order mark says why.
    let lead: InputLead | undefined;
    try {
      // The source is added once: its other kind's index reads the same data (no second inflation of a .gz).
      ref = known === undefined ? await this.deps.data.addSource(file) : undefined;
      const id = known ?? ref!.sourceId;
      sourceId = id;
      this.refs.set(holder, [id]);
      const preview =
        ref === undefined ? undefined : await this.deps.data.previewSource(id, HEAD_BYTES).catch(() => undefined);
      head = preview === undefined ? undefined : headText(preview);
      lead = preview === undefined ? undefined : inputLead(preview, preview.length >= ref!.size);
      const base = await this.deps.data.indexSource(id, reader);
      if (!current()) return;
      if (this.readAgain(role, key, file, parser, reader, { ref, sourceId, head, lead })) return;
      const added = ref;
      this.updateSource(role, key, (source) => ({
        ...source,
        // A file named .gz is inflated by addSource: the card shows its plain name and inflated size.
        ...(added === undefined ? {} : { name: added.name, size: added.size }),
        sourceId: id,
        status: 'ready',
        base,
        excluded: [],
        check: { state: 'pending' },
        indexMs: this.deps.now() - started,
        ...(head === undefined ? {} : { head }),
        ...(lead === undefined ? {} : { lead }),
      }));
      this.noteRevision(id, base.revisionId);
      this.scheduleCheck(role, key, 0);
      this.recheckReaderErrors(role, key);
      void this.askAgain();
    } catch (error) {
      if (!current()) return;
      if (this.readAgain(role, key, file, parser, reader, { ref, sourceId, head, lead })) return;
      this.updateSource(role, key, (source) => ({
        ...source,
        status: 'failed',
        error: messageOf(error),
        ...(isUnavailable(error) ? { unavailable: true } : {}),
        ...(lead === undefined ? {} : { lead }),
      }));
    } finally {
      this.refs.release(holder);
      this.scheduleValidation();
    }
  }

  /**
   * Reads the input of a source again when the reader that an index of it used (`reader`) is no longer the reader in
   * force: the program or the data loaders changed while it was read, and what that reader found (its records or its
   * refusal) is not what the run's reader finds. The name, the first lines and the lead of the first read stay with the
   * source; the data layer's source (`sourceId`) is read again, or added again if the first read could not add it.
   */
  private readAgain(
    role: InputRole,
    key: string,
    file: File,
    parser: FastaParserKind,
    reader: FastaReader,
    first: { ref: SourceRef | undefined; sourceId: string | undefined; head: string | undefined; lead: InputLead | undefined },
  ): boolean {
    const inForce = this.readerOf(role);
    if (inForce.program === reader.program && inForce.dataLoaders === reader.dataLoaders) return false;
    const { ref, sourceId, head, lead } = first;
    this.updateSource(role, key, (source) => ({
      ...source,
      ...(ref === undefined ? {} : { name: ref.name, size: ref.size }),
      ...(head === undefined ? {} : { head }),
      ...(lead === undefined ? {} : { lead }),
    }));
    // The new index holds the data layer's source before this one lets it go (its holder is set synchronously).
    this.track(this.index(role, key, file, parser, sourceId));
    return true;
  }

  private scheduleCheck(role: InputRole, key: string, delay?: number): void {
    this.updateSource(role, key, (source) => (source.status === 'ready' ? { ...source, check: { state: 'pending' } } : source));
    this.debounce(`check:${role}:${key}`, () => this.track(this.check(role, key)), delay);
  }

  /** Checks a source now, instead of a check of it that waits (askAgain); settles once it is checked. */
  private checkNow(role: InputRole, key: string): Promise<void> {
    this.updateSource(role, key, (source) => (source.status === 'ready' ? { ...source, check: { state: 'pending' } } : source));
    const name = `check:${role}:${key}`;
    clearTimeout(this.timers.get(name));
    this.timers.delete(name);
    this.workers.delete(name);
    const checked = this.check(role, key);
    this.track(checked);
    return checked;
  }

  private async check(role: InputRole, key: string): Promise<void> {
    const task = `check:${role}:${key}`;
    const generation = this.nextGeneration(task);
    const source = this.role(role).sources.find((s) => s.key === key);
    if (source?.status !== 'ready' || source.base === undefined) return;
    // An input whose reading ends with the reader's error is checked whatever it includes: its run
    // input keeps the error (line WQA's D10, 2026-11-19).
    if (includedOf(source).length === 0 && source.base.readerError === undefined) {
      this.updateSource(role, key, (s) => ({ ...s, check: { state: 'ok', records: 0 } }));
      return;
    }
    const reader = this.readerOf(role);
    const { program, dataLoaders } = reader;
    const started = this.deps.now();
    let check: SourceCheck;
    try {
      const revisionId = await this.revisionFor(source);
      // A source removed meanwhile is not checked (nor cached).
      if (source.sourceId === undefined || !this.formSources.has(source.sourceId)) return;
      const key = `${program}|${role}|${dataLoaders}|${revisionId}`;
      let verdict = this.verdicts.get(key);
      if (verdict === undefined) {
        verdict = this.deps.data.checkInput(reader, [revisionId]);
        this.verdicts.set(key, verdict);
        // A check that could not run is tried again next time.
        verdict.catch(() => this.verdicts.delete(key));
      }
      check = toSourceCheck(await verdict, source);
      // A Combined input reads the sources before this one first: the run's error names its own
      // line, which the check of the sources together gives (review L2, 2026-11-20).
      const before = this.sourcesBefore(role, source);
      if (check.state === 'ok' && check.readerError !== undefined && before.length > 0) {
        const ids = [...(await Promise.all(before.map((s) => this.revisionFor(s)))), revisionId];
        const joined = await this.deps.data.checkInput(reader, ids).catch(() => undefined);
        if (joined?.ok === true && joined.readerError !== undefined) check = { ...check, readerError: joined.readerError };
      }
    } catch (error) {
      check = { state: 'error', message: messageOf(error), ...(isUnavailable(error) ? { unavailable: true as const } : {}) };
    }
    if (this.generations.get(task) !== generation || this.state.get().program !== program || this.state.get().dataLoaders !== dataLoaders) return;
    const checkMs = this.deps.now() - started;
    this.updateSource(role, key, (s) =>
      s.excluded === source.excluded && s.base === source.base ? { ...s, check, checkMs } : s,
    );
    if (check.state !== 'error') void this.askAgain();
  }

  /**
   * The sources that a Combined run of the query reads before `source`: every source of the run's
   * input before it, whose bytes the run passes as they are, those without records too
   * (`DataService.runInputBytes`; review L-a, 2026-11-20). None otherwise, and none when a source
   * before ends the input with its reader's error, or `source` cannot follow them: a run then does
   * not read `source` after them, and its own check stands.
   */
  private sourcesBefore(role: InputRole, source: DraftSource): readonly DraftSource[] {
    const draft = this.role(role);
    if (role !== 'query' || draft.mode !== 'combined' || cannotFollow(draft.sources, source)) return [];
    const at = draft.sources.findIndex((s) => s.key === source.key);
    const before = draft.sources.slice(0, Math.max(0, at)).filter((s) => s.status === 'ready');
    return before.some((s) => endsInput(s.base)) ? [] : before;
  }

  /**
   * Checks again the query sources whose reader's error a run meets after other sources: the line
   * that its message names follows them (`check`, review L2, 2026-11-20).
   */
  private recheckReaderErrors(role: InputRole, except?: string): void {
    if (role !== 'query') return;
    for (const source of this.role(role).sources) {
      if (source.key !== except && source.status === 'ready' && source.base?.readerError !== undefined) {
        this.scheduleCheck(role, source.key, 0);
      }
    }
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
    this.debounce('validate', () => this.track(this.validateInTurn()));
  }

  /**
   * One validation at a time, the latest winning: one that comes while another runs waits for it, and
   * then validates only if no newer change came meanwhile. So a validation with the file of -conffile,
   * which the Data worker reads on an instance of its own, never runs beside another one that reads the
   * file again, and none starts for an argv that is already out of date.
   */
  /** Settles once no validation waits or runs; one that waits runs now, as idle runs it (askAgain). */
  private async validationSettled(): Promise<void> {
    for (;;) {
      if (this.timers.has('validate')) {
        clearTimeout(this.timers.get('validate'));
        this.timers.delete('validate');
        this.runTimer('validate');
      }
      const turn = this.validationTurn;
      await turn;
      if (!this.timers.has('validate') && turn === this.validationTurn) return;
    }
  }

  private validateInTurn(): Promise<void> {
    const generation = this.generations.get('validate');
    const turn = this.validationTurn.then(() => (this.generations.get('validate') === generation ? this.validate() : undefined));
    this.validationTurn = turn.catch(() => undefined);
    return turn;
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
      // No argv, so no start of one: the -conffile row goes with the words that named the file.
      this.set({ validation: { state: 'invalid', message: messageOf(error) }, start: undefined });
      return;
    }
    this.set({ validation: { state: 'checking' } });
    let validation: Validation;
    let start: StartInfo | undefined;
    let registry: ChosenRegistryFile | undefined;
    // The data loaders of the argv's run, where the engine found the argv valid.
    let dataLoaders: boolean | undefined;
    // Whether the engine answered (its start did), whatever the verdict.
    let answered = false;
    try {
      // NCBI's application start first (line RG E1w-3b): the name of the registry file that the argv
      // reads, under which the chosen file is checked. An argv that it cannot start (a program that
      // ABI v2 does not run) has the engine's message, as validate gives it.
      const started = await this.deps.engine.start(argv);
      answered = true;
      // A newer change came meanwhile: its validation follows, and this one reads no file.
      if (this.generations.get('validate') !== generation) return;
      if (!started.ok) {
        validation = { state: 'invalid', message: started.message };
      } else {
        start = started.start;
        registry = this.registryFile(program, start);
        const result = await this.deps.engine.validate(argv, registry);
        validation = result.ok ? { state: 'ok' } : { state: 'invalid', message: result.message };
        // Also from a refused argv whose registry the engine read (line WQA's D13 round 3, L1): the setting does not depend
        // on the option that it refuses.
        if (result.dataLoaders !== undefined) dataLoaders = result.dataLoaders;
      }
    } catch (error) {
      // The chosen file changed or went after it was chosen: the browser no longer reads it.
      const message =
        registry !== undefined && isFileReadError(error) ? unreadableConffile(registry.name, registry.file.name, error) : messageOf(error);
      // A part of LOSAT that could not be loaded is not a refusal of the options (line WQA's sweep, D1).
      validation = isUnavailable(error) ? { state: 'unavailable', message } : { state: 'invalid', message };
    }
    if (this.generations.get('validate') !== generation) return;
    this.validArgv = validation.state === 'ok' ? argv : undefined;
    this.set({ validation, start });
    if (dataLoaders !== undefined) this.setDataLoaders(dataLoaders);
    this.prepareEngine();
    if (answered) void this.askOptionsAgain();
  }

  /** The registry file of the argv's `-conffile`: the chosen file under the argv's name (`registryName`); undefined without either. */
  private registryFile(program: ProgramId, start: StartInfo): ChosenRegistryFile | undefined {
    const name = registryName(start);
    const file = this.state.get().conffiles[program];
    return name === undefined || file === undefined ? undefined : { name, file };
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
      inputBytes += Math.max(...groups.map((sources) => sources.reduce((sum, source) => sum + source.size, 0)));
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
    // Before the listeners, which may change the state again.
    this.holdFormSources(next as DraftState);
    this.state.set(next as DraftState);
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
      return `${accession} makes the fetched text larger than ${withUnit(FETCH_BYTE_LIMIT / 1_000_000, 'MB')}: download it from ${label} and open it with "Open FASTA files…".`;
    case 'unreachable':
      return `${label} could not be reached (${result.detail}).`;
    case 'refused':
      return result.status === 429
        ? `${label} asks to wait (too many requests): fetch again in a moment.`
        : `${label} answered HTTP ${result.status} for ${accession}.`;
  }
}

/** The registry file of a validation: the File chosen for `-conffile`, under the argv's name. */
interface ChosenRegistryFile extends RegistryFile {
  readonly file: File;
}

/** Whether the form holds inputs: a pasted text (accessions too) or a source of either role. */
export function holdsInputs(state: DraftState): boolean {
  return ROLES.some((role) => state[role].paste !== '' || state[role].sources.length > 0);
}

/** The included records of a ready source. */
export function includedOf(source: DraftSource): readonly DatasetRecord[] {
  return source.base === undefined ? [] : includedRecords({ ...source.base, excluded: source.excluded });
}

/**
 * The sources whose records a Combined run reads, in order: those with included records or the
 * reader's error, up to the first with the error, which ends the run's input
 * (`DataService.runInputBytes`, line WQA's D10, 2026-11-19).
 */
function combinedSources(sources: readonly DraftSource[]): readonly DraftSource[] {
  return readSources(sources.filter((source) => source.status === 'ready')).filter(
    (source) => includedOf(source).length > 0 || endsInput(source.base),
  );
}

/**
 * The sources that a run's input reads, in order: up to and with the first whose reading ends with the reader's error, which
 * ends the input, the engine reading no record after it (`DataService.runInputBytes`, line WQA's D10, 2026-11-19). The one
 * place in the draft that cuts a list of sources there: the readiness of a source after it (combinedSources), the name of the
 * input (nameOf) and the words that say what it is made of (fromOf, line WQA's D22) take their sources from it, so that the name
 * and the words name what the bytes hold. `sourcesBefore` and `queryEndBefore` ask the same question of one source, with
 * `endsInput`, the one test of it (domain/dataset.ts). `DataService.runInputBytes` is the data layer's own reading of the same
 * revisions: it needs the error's offset to keep the bytes after the records, and it is outside the draft (a layer, and held by
 * D17), so it tests the same field there.
 */
function readSources(sources: readonly DraftSource[]): readonly DraftSource[] {
  const stop = sources.findIndex((source) => endsInput(source.base));
  return stop < 0 ? sources : sources.slice(0, stop + 1);
}

/** Whether what a run reads of a source starts without a defline (review M2, 2026-11-20). */
function startsWithoutDefline(source: DraftSource): boolean {
  const first = includedOf(source)[0];
  return first !== undefined ? first.sequence_offset === 0 : source.base?.readerError?.defline === false;
}

/**
 * Whether `source` starts without a defline after a source whose records a Combined input of
 * `sources` reads: the record before it would read its first line, so no run takes it
 * (`SearchDraft`'s readiness; review M2, 2026-11-20), and its card says only that (review L-c).
 */
export function cannotFollow(sources: readonly DraftSource[], source: DraftSource): boolean {
  return combinedSources(sources).findIndex((s) => s.key === source.key) > 0 && startsWithoutDefline(source);
}

/** What a source that cannot follow another input starts with (`cannotFollow`), after its name. */
export function cannotFollowText(source: DraftSource): string {
  return includedOf(source).length > 0
    ? 'starts with a sequence without a defline (">"), which cannot follow another input.'
    : 'starts with a line that is not a defline (">"), which cannot follow another input. Correct that line or remove this input.';
}

/**
 * The source before `source` whose reader's error ends a Combined query input, after which the run
 * reads nothing (`DataService.runInputBytes`): `source` is not searched (review L2, 2026-11-20). A
 * source that cannot follow another input ends no run: no run starts (review L-c, 2026-11-20).
 */
export function queryEndBefore(sources: readonly DraftSource[], source: DraftSource): DraftSource | undefined {
  for (const before of sources) {
    if (before.key === source.key) return undefined;
    if (before.status === 'ready' && endsInput(before.base) && !cannotFollow(sources, before)) return before;
  }
  return undefined;
}

/**
 * Why what runs take of a source that is read holds no residue, if it holds none: `none`, the input has
 * no residue at all (an empty file, white space only, deflines alone); `excluded`, it has records with
 * residues, but they are left out and the records included have none. The engine reads such an input,
 * as NCBI's reader does, and runs search it as the command line does: it adds nothing to the search,
 * and NCBI BLAST+ 2.17.0 warns of it or ends the run with an error, by role, program and the other
 * inputs (an empty subject alone: `Empty CBlastQueryVector`; a lone subject record without residues:
 * `The average subject length is too short` in BLASTN and TBLASTX, a warning in the others; an empty
 * query alone: `Query is Empty!`, a warning; a lone query record without residues: `Sequence contains
 * no data`, an error). So its card shows no "Accepted" mark (line WQA's sweep, D4, 2026-11-19; D4
 * review). Undefined for a source with an included record that has residues, and for one whose
 * records are all left out: runs leave it out.
 */
export function noSequenceOf(source: DraftSource): 'none' | 'excluded' | undefined {
  if (source.status !== 'ready' || source.base === undefined) return undefined;
  const { records } = source.base;
  if (records.length === 0) return 'none';
  const included = includedOf(source);
  if (included.length === 0 || included.some((record) => record.length > 0)) return undefined;
  return records.some((record) => record.length > 0) ? 'excluded' : 'none';
}

/** Whether the sources are one file of the user's, with every record in: the bytes of the run's input are the file's own. */
function isFileAsGiven(sources: readonly DraftSource[]): boolean {
  // Of the sources that the run reads: a Combined input whose first file ends with the reader's error is that file.
  const read = readSources(sources);
  return read.length === 1 && read[0]!.origin === 'file' && read[0]!.excluded.length === 0;
}

/**
 * The name that an input has in the argv, where the command line writes the name of a file: the name of a file whose
 * bytes are what the run searches. The user's own file name, as the command line is typed, where the input is that file as
 * given: a file named `.gz`, `.bz2` or `.zst` is inflated for the engine (the card's `name` is the plain one), but the engine's
 * run takes the bytes of its inputs, not files, and writes the names of its argv where the command line writes the names of its
 * files (the `Input:` line of a subject set in outfmt 0 to 4 and 7, SAM's command line), so the argv names the file the user gave
 * (line WQA's D21, S3-3 and the review's F1, 2026-11-23). The user's file with records left out is not that file: its name would
 * make the command search other bytes than the run did, so the argv names the file of the included records, `included_` and the
 * name (the plain one: the bytes are not compressed), as `combined_query.fa` names several inputs together and `query.fa` pasted
 * text; the page offers these files with the run (line WQA's D22, 2026-11-23). `runNames` keeps the two names of a run apart.
 */
function nameOf(sources: readonly DraftSource[], role: InputRole): string {
  const read = readSources(sources);
  if (isFileAsGiven(read)) return read[0]!.file.name;
  if (read.length > 1) return COMBINED_NAMES[role];
  const [source] = read;
  return source === undefined || source.origin === 'paste' ? PASTED_NAMES[role] : `${INCLUDED_PREFIX}${source.name}`;
}

/** An input's name in the argv, and whether it names a file that LOSAT made (not the user's file as given). */
export interface NamedInput {
  readonly name: string;
  readonly made: boolean;
}

/**
 * The names of the two inputs of one run (line WQA's D22, review M1, 2026-11-23): a command that names two inputs by one name
 * cannot be run on the files that the page offers, which then hold other bytes under one name. A file that LOSAT made takes
 * its role into its name (`included_query_q.fa`) where the other input has the same name: both made (the same file in both
 * roles with other records left out), or the other the user's own file of that name (the role makes the name longer than the
 * other's). The user's own names stay as they are.
 */
export function runNames(query: NamedInput, subject: NamedInput): Record<InputRole, string> {
  if (query.name !== subject.name) return { query: query.name, subject: subject.name };
  // The same name: a file of LOSAT's takes its role, which makes it differ from the other's name; the user's own stays.
  return {
    query: query.made ? withRole(query.name, 'query') : query.name,
    subject: subject.made ? withRole(subject.name, 'subject') : subject.name,
  };
}

/** `included_q.fa` -> `included_query_q.fa`; `query.fa` -> `query_query.fa`. */
function withRole(name: string, role: InputRole): string {
  return name.startsWith(INCLUDED_PREFIX) ? `${INCLUDED_PREFIX}${role}_${name.slice(INCLUDED_PREFIX.length)}` : `${role}_${name}`;
}

/**
 * What an input made of `sources` is made of, for the run's snapshot; undefined where it is the user's file as given. Only the
 * sources that the run reads (readSources), and of the records left out the first few IDs and their number: a file with a
 * million records and ten kept has a million left out, which the page names three of.
 */
function fromOf(sources: readonly DraftSource[]): readonly InputPart[] | undefined {
  if (isFileAsGiven(sources)) return undefined;
  return readSources(sources).map((source) => ({
    // The user's file by its own name (a compressed file by the name that the user gave), pasted text as such.
    name: source.origin === 'paste' ? 'pasted text' : source.file.name,
    records: includedOf(source).length,
    leftOut: leftOutIdsOf(source),
    leftOutCount: source.excluded.length,
  }));
}

/**
 * The IDs of the first records that a source leaves out, as the page names them: a record without a title by its place in the
 * file (`record 3`), a record whose ID another record has too by its place (`r1 (#2)`, the number of the card's list).
 */
function leftOutIdsOf(source: DraftSource): readonly string[] {
  const table = source.base?.records ?? [];
  const named = source.excluded.slice(0, MOST_LEFT_OUT_IDS);
  // How many records have the ID of each, in one pass over the table.
  const counts = new Map<string, number>(named.map((index) => [table[index]?.id ?? '', 0]));
  for (const record of table) {
    const count = counts.get(record.id);
    if (count !== undefined) counts.set(record.id, count + 1);
  }
  return named.map((index) => {
    const id = table[index]?.id ?? '';
    if (id === '') return `record ${index + 1}`;
    return (counts.get(id) ?? 0) > 1 ? `${id} (#${index + 1})` : id;
  });
}

function describeSource(source: DraftSource): string {
  return source.origin === 'paste' ? 'pasted' : `file ${source.name}`;
}

/** The engine's verdict, with the record its message names mapped to the record table. */
function toSourceCheck(check: InputCheck, source: DraftSource): SourceCheck {
  if (check.ok) {
    const { readerError } = check;
    return { state: 'ok', records: check.records.length, ...(readerError === undefined ? {} : { readerError }) };
  }
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
