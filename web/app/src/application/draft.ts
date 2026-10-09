// The search draft (design §3.1, §3.2; plan §5.2-§5.4): the job being edited on the search
// screen. It is independent of the queue, so the next job can be edited while a search
// runs, and "Add to queue" freezes the draft's inputs and conditions into run snapshots
// (Coordinator.enqueueAll); the draft itself stays, so a subject can be kept while the
// query changes.
//
// Inputs: every pasted text and every file becomes a source of the data layer and is
// indexed when it is added (DatasetStore). The record table of a source is its
// `base` revision; the records that the user leaves out make other revisions, created on
// demand and reused for the same selection, so that a subject whose selection does not
// change is the same revision from search to search (the engine then keeps it, R1).
// Whether the engine reads a source is the engine's verdict (`checkInput`, the program's
// `register`), shown with its message. Several sources of a role are one combined search
// input, or separate searches (one run per source, queued as a group).
import { buildArgv, COMBINED_NAMES, PASTED_NAMES } from '../domain/argv';
import { includedRecords, type DatasetRecord, type DatasetRevision, type FastaParserKind } from '../domain/dataset';
import { formParameters, setField, type FieldValue, type FormValues } from '../domain/parameters';
import { PROGRAMS, programById, type InputRole, type ProgramId } from '../domain/programs';
import { REGION_FLAG, regionProblem, regionValue, type RegionText } from '../domain/region';
import type { DataGateway } from '../ports/data';
import type { EngineGateway, ProgramDescription } from '../ports/engine';
import type { InputCheck } from '../ports/input-check';
import type { EnqueueAllResult, SearchRequest } from './coordinator';
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

export interface DraftSource {
  /** Identifies the source within the draft. */
  readonly key: string;
  readonly origin: 'file' | 'paste';
  readonly name: string;
  readonly size: number;
  /** The File of the source (pasted text is a File too); a source is indexed again from it. */
  readonly file: File;
  readonly status: 'indexing' | 'ready' | 'failed';
  /** Why the source could not be indexed (the index scan's error). */
  readonly error?: string;
  /** The record table: the revision of the data layer that includes every record. */
  readonly base?: DatasetRevision;
  /** Indices of the records that runs leave out, ascending. */
  readonly excluded: readonly number[];
  /** The first lines of the source, as text, for its preview. */
  readonly head?: string;
  readonly check?: SourceCheck;
  /** Milliseconds that the index (scan and record hashes) took. */
  readonly indexMs?: number;
  /** Milliseconds that the latest engine check took (the run input and `register`). */
  readonly checkMs?: number;
}

export interface RoleDraft {
  /** The text of the paste box; its source is the role's paste source. */
  readonly paste: string;
  /** The paste source first, then the files in the order they were added. */
  readonly sources: readonly DraftSource[];
  readonly mode: InputMode;
  /**
   * The region of the role's only record (DW-9), and that record (`regionIdentity`). It
   * applies only while that record is still the role's only included record; another
   * record, or another count of records, leaves it aside.
   */
  readonly region?: { readonly text: RegionText; readonly record: string };
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
  readonly threads: number | 'auto';
  /** The engine's `describe` of the program, once it has arrived. */
  readonly description?: ProgramDescription;
  readonly descriptionError?: string;
  /** The engine's `validate` of the argv that the draft makes now. */
  readonly validation: Validation;
  readonly submitting: boolean;
  readonly message?: { readonly kind: 'info' | 'error'; readonly text: string };
}

export interface DraftDeps {
  readonly engine: Pick<EngineGateway, 'describe' | 'validate'>;
  readonly data: Pick<
    DataGateway,
    'addSource' | 'indexSource' | 'reviseDataset' | 'checkInput' | 'previewSource'
  >;
  readonly enqueueAll: (requests: readonly SearchRequest[]) => Promise<EnqueueAllResult>;
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
  /** The latest request of each asynchronous task; older results are dropped. */
  private readonly generations = new Map<string, number>();
  private readonly timers = new Map<string, ReturnType<typeof setTimeout>>();
  private readonly workers = new Map<string, () => void>();
  private readonly pending = new Set<Promise<unknown>>();
  private nextKey = 1;

  constructor(private readonly deps: DraftDeps) {
    const values = Object.fromEntries(PROGRAMS.map((program) => [program.id, Object.freeze({})])) as Record<
      ProgramId,
      FormValues
    >;
    this.state = new Store<DraftState>({
      program: 'blastn',
      query: emptyRole(),
      subject: emptyRole(),
      values,
      threads: 'auto',
      validation: { state: 'idle' },
      submitting: false,
    });
    this.loadDescription('blastn');
    this.scheduleValidation();
  }

  // --- program and form -------------------------------------------------------------------

  setProgram(program: ProgramId): void {
    if (this.state.get().program === program) return;
    const next = programById(program);
    const description = this.described.get(program);
    this.set({ program, description, descriptionError: undefined, message: undefined });
    this.loadDescription(program);
    for (const role of ROLES) {
      for (const source of this.role(role).sources) {
        if (next.unavailable !== undefined) {
          this.scheduleCheck(role, source.key);
          continue;
        }
        if (source.base !== undefined && source.base.parser !== indexParser(program)) {
          // A program with another FASTA reader needs another record table (BLASTX, plan TD-8).
          this.reindex(role, source.key);
        } else {
          this.scheduleCheck(role, source.key, 0);
        }
      }
    }
    this.scheduleValidation();
  }

  setField(flag: string, value: FieldValue): void {
    const { program, values } = this.state.get();
    this.set({ values: { ...values, [program]: setField(values[program], flag, value) }, message: undefined });
    this.scheduleValidation();
  }

  setThreads(threads: number | 'auto'): void {
    this.set({ threads, message: undefined });
  }

  // --- inputs -----------------------------------------------------------------------------

  /**
   * The paste box's text. Its source shows that it is being read at once, and is indexed
   * after a pause in typing.
   */
  setPaste(role: InputRole, text: string): void {
    this.updateRole(role, (draft) => ({ ...draft, paste: text }));
    this.set({ message: undefined });
    if (text === '') {
      this.removeSource(role, PASTE_KEY);
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
    this.debounce(`paste:${role}`, () => this.track(this.index(role, PASTE_KEY, file)));
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
    for (const file of files) {
      const key = `file-${this.nextKey++}`;
      this.updateRole(role, (draft) => ({
        ...draft,
        sources: [
          ...draft.sources,
          { key, origin: 'file', name: file.name, size: file.size, file, status: 'indexing', excluded: [] },
        ],
      }));
      this.track(this.index(role, key, file));
    }
    this.set({ message: undefined });
    this.scheduleValidation();
  }

  removeSource(role: InputRole, key: string): void {
    this.cancelTask(`${role}:${key}`);
    if (key === PASTE_KEY) this.cancelTask(`paste:${role}`);
    this.updateRole(role, (draft) => ({
      ...draft,
      paste: key === PASTE_KEY ? '' : draft.paste,
      sources: draft.sources.filter((source) => source.key !== key),
    }));
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
    this.updateRole(role, (draft) => {
      const next = { ...draft };
      if (region === undefined || record === undefined) delete next.region;
      else next.region = { text: region, record };
      return next;
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
    const region = this.role(role).region;
    return region !== undefined && region.record === this.regionIdentity(role) ? region.text : undefined;
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

  /** The name of the role's input in the argv of a combined search (plan §5.3). */
  inputName(role: InputRole): string {
    const sources = this.role(role).sources.filter((source) => source.status !== 'failed');
    if (sources.length === 1) return nameOf(sources[0]!, role);
    return sources.length === 0 ? PASTED_NAMES[role] : COMBINED_NAMES[role];
  }

  /** The runs that "Add to queue" would make now. */
  runCount(): number {
    return ROLES.map((role) => this.inputGroups(role).length).reduce((a, b) => a * Math.max(1, b), 1);
  }

  // --- submit -----------------------------------------------------------------------------

  /** Freezes the draft into queued runs (one, or a group of separate searches). */
  async submit(): Promise<EnqueueAllResult> {
    if (this.state.get().submitting) return { ok: false, message: 'The draft is being added already.' };
    this.set({ submitting: true, message: { kind: 'info', text: 'Preparing the inputs…' } });
    try {
      let result: EnqueueAllResult;
      try {
        result = await this.prepareAndEnqueue();
      } catch (error) {
        result = { ok: false, message: `The search could not be queued: ${messageOf(error)}` };
      }
      if (result.ok) {
        const count = result.runIds?.length ?? 0;
        this.set({ message: { kind: 'info', text: count > 1 ? `Added ${count} runs to the queue as a group.` : 'Added to the queue.' } });
      } else {
        this.set({ message: { kind: 'error', text: result.message } });
      }
      return result;
    } finally {
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

  private async prepareAndEnqueue(): Promise<EnqueueAllResult> {
    await this.idle();
    const { program, threads } = this.state.get();
    const descriptor = programById(program);
    if (descriptor.unavailable !== undefined) return { ok: false, message: descriptor.unavailable };
    for (const role of ROLES) {
      const problem = this.inputProblem(role);
      if (problem !== undefined) return { ok: false, message: problem };
    }
    // Everything from one state, before any wait.
    const groups = { query: this.inputGroups('query'), subject: this.inputGroups('subject') };
    const parameters = this.parameters();
    const resolved = {
      query: await Promise.all(groups.query.map((group) => this.resolveGroup(group))),
      subject: await Promise.all(groups.subject.map((group) => this.resolveGroup(group))),
    };
    const requests: SearchRequest[] = [];
    for (const query of resolved.query) {
      for (const subject of resolved.subject) {
        requests.push({ program, query: { dataset: query }, subject: { dataset: subject }, parameters, requestedThreads: threads });
      }
    }
    return this.deps.enqueueAll(requests);
  }

  /** Why the role cannot be searched as it is (the full message), or undefined. */
  private inputProblem(role: InputRole): string | undefined {
    return this.problemOf(role)?.message;
  }

  /**
   * What keeps the draft from being queued now, in short, for the line above the button;
   * undefined when nothing does. "Add to queue" reports the full message.
   */
  readiness(): string | undefined {
    if (programById(this.state.get().program).unavailable !== undefined) return undefined;
    for (const role of ROLES) {
      const problem = this.problemOf(role);
      if (problem !== undefined) return problem.short;
    }
    return undefined;
  }

  private problemOf(role: InputRole): { readonly message: string; readonly short: string } | undefined {
    const draft = this.role(role);
    const title = role === 'query' ? 'query' : 'subject';
    const label = role === 'query' ? 'Query' : 'Subject';
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
    return { name: group.name, revisionIds };
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

  private reindex(role: InputRole, key: string): void {
    const source = this.role(role).sources.find((s) => s.key === key);
    if (source === undefined) return;
    this.updateSource(role, key, (s) => {
      const next: DraftSource = { ...s, status: 'indexing', excluded: [] };
      delete (next as { base?: unknown }).base;
      delete (next as { check?: unknown }).check;
      return next;
    });
    this.track(this.index(role, key, source.file));
  }

  private async index(role: InputRole, key: string, file: File): Promise<void> {
    const task = `${role}:${key}`;
    const generation = this.nextGeneration(task);
    const current = () => this.generations.get(task) === generation;
    const started = this.deps.now();
    const parser = indexParser(this.state.get().program);
    try {
      const ref = await this.deps.data.addSource(file);
      const head = await this.deps.data.previewSource(ref.sourceId, HEAD_BYTES).then(headText, () => undefined);
      const base = await this.deps.data.indexSource(ref.sourceId, parser);
      if (!current()) return;
      this.updateSource(role, key, (source) => ({
        ...source,
        status: 'ready',
        base,
        excluded: [],
        check: { state: 'pending' },
        indexMs: this.deps.now() - started,
        ...(head === undefined ? {} : { head }),
      }));
      this.scheduleCheck(role, key, 0);
    } catch (error) {
      if (!current()) return;
      this.updateSource(role, key, (source) => ({ ...source, status: 'failed', error: messageOf(error) }));
    } finally {
      this.scheduleValidation();
    }
  }

  private scheduleCheck(role: InputRole, key: string, delay?: number): void {
    if (programById(this.state.get().program).unavailable !== undefined) {
      // Nothing is checked for a program that cannot be searched yet.
      this.cancelTask(`check:${role}:${key}`);
      this.updateSource(role, key, (source) => {
        const rest = { ...source };
        delete (rest as { check?: unknown }).check;
        return rest;
      });
      return;
    }
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
      revision = this.deps.data.reviseDataset(base.revisionId, source.excluded).then((r) => r.revisionId);
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
    if (programById(program).unavailable !== undefined) {
      this.set({ validation: { state: 'idle' } });
      return;
    }
    let argv: readonly string[];
    try {
      argv = buildArgv({
        program,
        queryName: this.inputName('query'),
        subjectName: this.inputName('subject'),
        parameters: this.parameters(),
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
    if (this.generations.get('validate') === generation) this.set({ validation });
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
  }
}

/**
 * The FASTA reader that indexes the sources of a program: its own, or for a program that
 * cannot be searched yet, that of BLASTN (the sources are indexed again when another
 * reader is needed).
 */
function indexParser(program: ProgramId): FastaParserKind {
  const descriptor = programById(program);
  return descriptor.unavailable === undefined ? descriptor.fastaParser : programById('blastn').fastaParser;
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
