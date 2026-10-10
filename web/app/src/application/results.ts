// The results screen's state (plan §5.7, design §11.1): which run, query, subject and HSP
// are selected, the view filters and sort orders (ViewState, never a new search), and what
// the screen shows for them. The selection is held by the HSP's identity (run, query,
// rank), never by a row number of a table, and a detail that arrives after the selection
// moved on is dropped. Values shown come from the stored outputs and the HSP records
// (docs/web/results_columns.md); the domain functions group, sort and filter them.
import type { RecordKey } from '../domain/dataset';
import { optionValue } from '../domain/argv';
import { translates, type Unit } from '../domain/coordinates';
import { frame, out0Range, out0SubjectRange, out6Range, subjectSpan } from '../domain/hsp-table';
import { splitOutfmt6Row, type Outfmt6Row } from '../domain/outfmt6';
import { programById, residueUnit, type ProgramId, type SequenceKind } from '../domain/programs';
import {
  buildResultIndex,
  filterQuery,
  hitLimits,
  NO_FILTERS,
  orientation,
  sortHsps,
  sortSubjects,
  type HitLimits,
  type HspSortKey,
  type Orientation,
  type ResultIndex,
  type SortSpec,
  type SubjectGroup,
  type SubjectSortKey,
  type ViewFilters,
} from '../domain/result-index';
import { verificationBadge, type Badge, type VerificationTable } from '../domain/verification';
import type { RunStore } from '../ports/data';
import type { ProgramDescription } from '../ports/engine';
import type { AppState, RunView } from './coordinator';
import { Store } from './store';

/** The identity of an HSP: its run, its query record and its rank among the query's HSPs. */
export interface HspId {
  readonly runId: string;
  readonly qIdx: number;
  readonly rank: number;
}

export const sameHsp = (a: HspId | undefined, b: HspId | undefined): boolean =>
  a !== undefined && b !== undefined && a.runId === b.runId && a.qIdx === b.qIdx && a.rank === b.rank;

export interface QueryEntry {
  readonly qIdx: number;
  readonly id: string;
  readonly length: number;
  readonly subjects: number;
  readonly hsps: number;
}

export interface SubjectEntry {
  readonly sIdx: number;
  /** The subject's ID in the run's record table. */
  readonly recordId: string;
  readonly length: number;
  /** Position of the subject in the engine's order for this query (1-based). */
  readonly order: number;
  /** The outfmt 6 row of the subject's first HSP. */
  readonly first: Outfmt6Row;
  /** HSPs of the subject in this query (before the view filters). */
  readonly hspCount: number;
  /** HSPs that pass the view filters, in the HSP sort order. */
  readonly rows: readonly number[];
  readonly inOutfmt0: boolean;
  /** The subject has as many HSPs as -max_hsps keeps. */
  readonly atHspLimit: boolean;
}

export interface HspEntry {
  readonly row: number;
  readonly id: HspId;
  readonly fields: Outfmt6Row;
  readonly orientation: Orientation;
  readonly queryFrame?: number;
  readonly subjectFrame?: number;
  readonly inOutfmt0: boolean;
}

/** An HSP of the selected subject as the Alignments show it: NCBI's "Range n: a to b". */
export interface RangeEntry {
  readonly id: HspId;
  readonly row: number;
  /** The HSP's position among its subject's HSPs in the engine's order (1-based), whatever the view filters hide. */
  readonly n: number;
  /** The smaller and the larger subject coordinate of the HSP record. */
  readonly from: number;
  readonly to: number;
  readonly inOutfmt0: boolean;
}

export interface Detail {
  readonly id: HspId;
  readonly state: 'loading' | 'ready' | 'failed';
  /** The outfmt 6 row, as written. */
  readonly row?: string;
  /** The subject heading of outfmt 0 (defline and Length=), as written. */
  readonly heading?: string;
  /** The HSP's section of outfmt 0 (score lines and alignment), as written. */
  readonly section?: string;
  readonly error?: string;
}

export interface Units {
  readonly query: Unit;
  readonly subject: Unit;
}

export interface LoadedRun {
  readonly run: RunView;
  readonly index: ResultIndex;
  readonly out6: Uint8Array;
  readonly description: ProgramDescription;
  readonly limits: HitLimits;
  readonly badge: Badge;
  readonly diagnostics: string;
  readonly kinds: { readonly query: SequenceKind; readonly subject: SequenceKind };
  readonly units: Units;
  /** The task that the search used: the argv's -task, or the engine's default; undefined for a program without tasks. */
  readonly task?: string;
}

export interface ResultsState {
  readonly runId?: string;
  readonly phase: 'none' | 'loading' | 'ready' | 'unavailable' | 'failed';
  /** Why the run has no results to show (not completed, cancelled, failed, or a read error). */
  readonly message?: string;
  readonly loaded?: LoadedRun;
  readonly filters: ViewFilters;
  readonly subjectSort: SortSpec<SubjectSortKey>;
  readonly hspSort: SortSpec<HspSortKey>;
  /** The queries listed, after the query filters. */
  readonly queries: readonly QueryEntry[];
  readonly qIdx?: number;
  /** Subjects of the query that pass the view filters, in the subject sort order. */
  readonly subjects: readonly SubjectEntry[];
  /** The selected query's subjects and HSPs before the view filters. */
  readonly queryTotals?: { readonly subjects: number; readonly hsps: number };
  readonly hidden: { readonly subjects: number; readonly hsps: number };
  /** The query has as many subjects as the hit list keeps. */
  readonly atSubjectLimit: boolean;
  /** Subjects of the query that outfmt 0 shows alignments for. */
  readonly outfmt0Subjects: number;
  readonly sIdx?: number;
  readonly hsps: readonly HspEntry[];
  /** The selected subject's HSPs that pass the view filters, in the engine's order (the Ranges of the Alignments). */
  readonly ranges: readonly RangeEntry[];
  readonly hsp?: HspId;
  readonly detail?: Detail;
  /** Subject headings of outfmt 0 read so far, by subject record index (the same in every query). */
  readonly headings: ReadonlyMap<number, string>;
}

export interface ResultsDeps {
  readonly data: Pick<RunStore, 'readHitTable' | 'readOutput' | 'readOutputRange' | 'readDiagnostics'>;
  readonly describe: (program: ProgramId) => Promise<ProgramDescription>;
  readonly runs: Store<AppState>;
  readonly verification: VerificationTable;
}

const INITIAL: ResultsState = Object.freeze<ResultsState>({
  phase: 'none',
  filters: NO_FILTERS,
  subjectSort: { key: 'order', descending: false },
  hspSort: { key: 'rank', descending: false },
  queries: [],
  subjects: [],
  hidden: { subjects: 0, hsps: 0 },
  atSubjectLimit: false,
  outfmt0Subjects: 0,
  hsps: [],
  ranges: [],
  headings: new Map(),
});

const decoder = new TextDecoder();

export class ResultsBrowser {
  readonly state = new Store<ResultsState>(INITIAL);
  /** Outfmt 6 rows already split, by table row, for the loaded run. */
  private rows = new Map<number, Outfmt6Row>();
  /**
   * Table rows of the loaded run by rank, per query, made when an HSP of the query is first looked
   * up: a map of every HSP built at each opening cost about 20 ms for 150,000 HSPs (W4b F1).
   */
  private rowsByRank = new Map<number, Int32Array>();
  private loadToken = 0;
  /** Subjects of the loaded run whose heading waits to be read. */
  private readonly headingQueue: number[] = [];
  /** The load token of the heading reader that runs, if one does. */
  private headingReader: number | undefined;
  /** Outfmt 0 sections of the loaded run read or being read, by `${qIdx}:${rank}`; another run drops them. */
  private sections = new Map<string, Promise<string | undefined>>();
  /** Outfmt 0 subject headings of the loaded run read or being read, by subject record index. */
  private headingReads = new Map<number, Promise<string | undefined>>();

  /** The table the verification badges come from (generated at build time). */
  get verification(): VerificationTable {
    return this.deps.verification;
  }

  constructor(private readonly deps: ResultsDeps) {
    // A run that the screen shows may finish, or be cancelled, while it is selected. A completed
    // run whose results could not be read ('failed') is not read again at each change of the runs.
    deps.runs.subscribe((app) => {
      const { runId, phase } = this.state.get();
      if (runId === undefined || phase !== 'unavailable') return;
      const run = app.runs.find((r) => r.snapshot.runId === runId);
      if (run?.status === 'completed') void this.open(runId);
      else if (run !== undefined) this.set({ message: unavailableMessage(run) });
    });
  }

  /**
   * Shows the results of a run: the run's HSP records and outfmt 6 text are read once. Another
   * run starts without view filters, as with the selection and the sort orders; the same run
   * (opened again, or completing while it is shown) keeps them.
   */
  async open(runId: string): Promise<void> {
    const token = ++this.loadToken;
    const run = this.deps.runs.get().runs.find((r) => r.snapshot.runId === runId);
    this.rows = new Map();
    this.rowsByRank = new Map();
    this.sections = new Map();
    this.headingReads = new Map();
    this.headingQueue.length = 0;
    if (run === undefined) return;
    const filters = this.state.get().runId === runId ? this.state.get().filters : NO_FILTERS;
    if (run.status !== 'completed') {
      this.state.set({ ...INITIAL, filters, runId, phase: 'unavailable', message: unavailableMessage(run) });
      return;
    }
    this.state.set({ ...INITIAL, filters, runId, phase: 'loading' });
    try {
      const [table, out6, description, diagnostics] = await Promise.all([
        this.deps.data.readHitTable(runId),
        this.deps.data.readOutput(runId, 6),
        this.deps.describe(run.snapshot.program),
        this.deps.data.readDiagnostics(runId),
      ]);
      if (token !== this.loadToken) return;
      const program = programById(run.snapshot.program);
      const kinds = { query: program.query, subject: program.subject };
      const loaded: LoadedRun = {
        run,
        index: buildResultIndex(table),
        out6,
        description,
        limits: hitLimits(run.snapshot.argv, description.parameters),
        badge: verificationBadge(
          {
            program: run.snapshot.program,
            argv: run.snapshot.argv,
            formats: description.formats,
            runtimePath: run.record.runtimePath,
            threads: run.record.threads,
            grammar: grammarOf(description),
          },
          this.deps.verification,
        ),
        diagnostics,
        kinds,
        units: { query: residueUnit(kinds.query), subject: residueUnit(kinds.subject) },
        ...taskOf(run.snapshot.argv, description),
      };
      const queries = this.queryEntries(loaded, this.state.get().filters);
      const first = queries.find((q) => q.hsps > 0) ?? queries[0];
      this.state.set({ ...this.state.get(), phase: 'ready', loaded, queries });
      if (first !== undefined) this.selectQuery(first.qIdx);
    } catch (error) {
      if (token !== this.loadToken) return;
      this.set({ phase: 'failed', message: `The results of run ${run.snapshot.number} could not be read: ${errorMessage(error)}` });
    }
  }

  selectQuery(qIdx: number): void {
    const state = this.state.get();
    if (state.loaded === undefined) return;
    this.set({ qIdx, sIdx: undefined, hsp: undefined, detail: undefined });
    this.refresh();
    const subject = this.state.get().subjects[0];
    if (subject !== undefined) this.selectSubject(subject.sIdx);
  }

  selectSubject(sIdx: number): void {
    this.set({ sIdx });
    this.refresh();
    const first = this.state.get().hsps[0];
    this.selectHsp(first?.id);
  }

  /**
   * Selects an HSP of the loaded run; its query and subject follow. Selecting the HSP whose detail
   * could not be read reads it again ("Try again"); one that is loading or ready is left alone.
   */
  selectHsp(id: HspId | undefined): void {
    const state = this.state.get();
    if (id === undefined || state.loaded === undefined || id.runId !== state.runId) {
      this.set({ hsp: undefined, detail: undefined });
      return;
    }
    const row = this.rowOf(state.loaded, id.qIdx, id.rank);
    if (row === undefined) return;
    const table = state.loaded.index.table;
    const sIdx = table.sIdx[row]!;
    if (state.qIdx !== id.qIdx || state.sIdx !== sIdx) {
      this.set({ qIdx: id.qIdx, sIdx });
      this.refresh();
    }
    const { hsp, detail } = this.state.get();
    if (sameHsp(hsp, id) && detail !== undefined && detail.state !== 'failed') return;
    this.set({ hsp: id });
    void this.loadDetail(state.loaded, id, row);
  }

  /**
   * Applies the view filters. A query filter that hides the selected query moves the selection
   * to the first query listed (as a subject filter does for subjects); none listed, none selected.
   */
  setFilters(filters: ViewFilters): void {
    const state = this.state.get();
    if (state.loaded === undefined) {
      this.set({ filters });
      return;
    }
    const queries = this.queryEntries(state.loaded, filters);
    this.set({ filters, queries });
    if (!queries.some((query) => query.qIdx === state.qIdx)) {
      const first = queries[0];
      if (first !== undefined) this.selectQuery(first.qIdx);
      else {
        this.set({ qIdx: undefined, sIdx: undefined, hsp: undefined, detail: undefined });
        this.refresh();
      }
      return;
    }
    this.refresh();
    this.keepSelection();
  }

  setSubjectSort(subjectSort: SortSpec<SubjectSortKey>): void {
    this.set({ subjectSort });
    this.refresh();
  }

  setHspSort(hspSort: SortSpec<HspSortKey>): void {
    this.set({ hspSort });
    this.refresh();
  }

  /** The outfmt 6 row of a table row of the loaded run, split into its fields. */
  row(row: number): Outfmt6Row {
    let fields = this.rows.get(row);
    if (fields === undefined) {
      const loaded = this.state.get().loaded!;
      const range = out6Range(loaded.index.table, row);
      if (range === undefined) throw new Error(`HSP row ${row} has no outfmt 6 row`);
      fields = splitOutfmt6Row(decoder.decode(loaded.out6.subarray(range[0], range[1])));
      this.rows.set(row, fields);
    }
    return fields;
  }

  /**
   * The outfmt 0 section of an HSP of the loaded run (its score lines and alignment), as written;
   * undefined for an HSP that outfmt 0 does not show, or of another run. Each section is read once
   * per run (the Alignments read the Ranges that come into view); a read that fails is tried again
   * at the next request.
   */
  readSection(id: HspId): Promise<string | undefined> {
    const state = this.state.get();
    const loaded = state.loaded;
    if (loaded === undefined || id.runId !== state.runId) return Promise.resolve(undefined);
    const key = `${id.qIdx}:${id.rank}`;
    let read = this.sections.get(key);
    if (read === undefined) {
      const sections = this.sections;
      const row = this.rowOf(loaded, id.qIdx, id.rank);
      const range = row === undefined ? undefined : out0Range(loaded.index.table, row);
      read =
        range === undefined
          ? Promise.resolve(undefined)
          : this.deps.data.readOutputRange(loaded.run.snapshot.runId, 0, range[0], range[1]).then((bytes) => decoder.decode(bytes));
      sections.set(key, read);
      read.catch(() => {
        if (sections.get(key) === read) sections.delete(key);
      });
    }
    return read;
  }

  /**
   * Reads the outfmt 0 headings of subjects of the selected query (the rows in view of the
   * subject list), once per subject. Subjects that outfmt 0 does not show have none.
   */
  requestHeadings(sIdxs: readonly number[]): void {
    const state = this.state.get();
    if (state.loaded === undefined) return;
    for (const sIdx of sIdxs) {
      if (!state.headings.has(sIdx) && !this.headingQueue.includes(sIdx)) this.headingQueue.push(sIdx);
    }
    // A reader of a run shown before stops at its next read, and leaves the queue to this one.
    if (this.headingReader !== this.loadToken) void this.readHeadings(this.loadToken);
  }

  private async readHeadings(token: number): Promise<void> {
    this.headingReader = token;
    try {
      while (this.headingQueue.length > 0 && token === this.loadToken) {
        const batch = this.headingQueue.splice(0, 40);
        const state = this.state.get();
        const loaded = state.loaded;
        const query = state.qIdx === undefined ? undefined : loaded?.index.queries.get(state.qIdx);
        if (loaded === undefined || query === undefined) break;
        const read = new Map<number, string>();
        for (const sIdx of batch) {
          if (token !== this.loadToken) break;
          const subject = query.subjects.find((s) => s.sIdx === sIdx);
          const heading = subject === undefined ? undefined : await this.readHeading(loaded, sIdx, subject.rows[0]!);
          if (heading !== undefined) read.set(sIdx, heading);
        }
        if (token !== this.loadToken) break;
        if (read.size > 0) this.set({ headings: new Map([...this.state.get().headings, ...read]) });
      }
    } catch {
      // A heading that cannot be read leaves the description empty; the detail reports read errors.
    } finally {
      if (this.headingReader === token) {
        this.headingReader = undefined;
        this.headingQueue.length = 0;
      }
    }
  }

  /** The table row of an HSP of the loaded run (its query's rows by rank, made once), or undefined. */
  private rowOf(loaded: LoadedRun, qIdx: number, rank: number): number | undefined {
    const table = loaded.index.table;
    let rows = this.rowsByRank.get(qIdx);
    if (rows === undefined) {
      const query = loaded.index.queries.get(qIdx);
      if (query === undefined) return undefined;
      rows = new Int32Array(query.hspCount).fill(-1);
      for (const subject of query.subjects) {
        for (const row of subject.rows) if (table.rank[row]! >= 0 && table.rank[row]! < rows.length) rows[table.rank[row]!] = row;
      }
      this.rowsByRank.set(qIdx, rows);
    }
    const row = rank >= 0 && rank < rows.length ? rows[rank]! : -1;
    return row >= 0 && table.qIdx[row] === qIdx && table.rank[row] === rank ? row : undefined;
  }

  /**
   * A subject's outfmt 0 heading (the same in every query), read once per run from the heading
   * before `row`'s section: the selected HSP's detail and the lists that show headings share it.
   * Undefined where outfmt 0 does not show the subject for that row's query.
   */
  private readHeading(loaded: LoadedRun, sIdx: number, row: number): Promise<string | undefined> {
    const range = out0SubjectRange(loaded.index.table, row);
    if (range === undefined) return Promise.resolve(undefined);
    const known = this.state.get().headings.get(sIdx);
    if (known !== undefined) return Promise.resolve(known);
    let read = this.headingReads.get(sIdx);
    if (read === undefined) {
      const reads = this.headingReads;
      read = this.deps.data.readOutputRange(loaded.run.snapshot.runId, 0, range[0], range[1]).then((bytes) => decoder.decode(bytes));
      reads.set(sIdx, read);
      // A read that fails is tried again at the next request.
      read.catch(() => {
        if (reads.get(sIdx) === read) reads.delete(sIdx);
      });
    }
    return read;
  }

  /** Recomputes the lists from the selection, the filters and the sort orders. */
  private refresh(): void {
    const state = this.state.get();
    const loaded = state.loaded;
    if (loaded === undefined || state.qIdx === undefined) {
      this.set({ subjects: [], hsps: [], ranges: [], hidden: { subjects: 0, hsps: 0 }, queryTotals: undefined, atSubjectLimit: false, outfmt0Subjects: 0 });
      return;
    }
    const { table } = loaded.index;
    const query = loaded.index.queries.get(state.qIdx);
    if (query === undefined) {
      this.set({
        subjects: [],
        hsps: [],
        ranges: [],
        hidden: { subjects: 0, hsps: 0 },
        queryTotals: { subjects: 0, hsps: 0 },
        atSubjectLimit: false,
        outfmt0Subjects: 0,
      });
      return;
    }
    const subjectRecords = loaded.run.snapshot.subject.records;
    const filtered = filterQuery(query, state.filters, table, (subject) => [
      subjectRecords[subject.sIdx]?.id ?? '',
      this.row(subject.rows[0]!).sseqid,
    ]);
    const order = new Map(query.subjects.map((subject, i) => [subject.sIdx, i + 1]));
    // The list shows each subject's first HSP and HSP count before the filters, and sorts by them.
    const unfiltered = new Map(query.subjects.map((subject) => [subject.sIdx, subject]));
    const shown = new Map(filtered.subjects.map((subject) => [subject.sIdx, subject]));
    const sorted = sortSubjects(
      filtered.subjects.map((subject) => unfiltered.get(subject.sIdx)!),
      state.subjectSort,
      table,
      (sIdx) => subjectRecords[sIdx]?.length ?? 0,
    );
    const subjects = sorted.map((all) => this.subjectEntry(loaded, shown.get(all.sIdx)!, all, order.get(all.sIdx)!, state.hspSort));
    const selected = subjects.find((s) => s.sIdx === state.sIdx);
    const hsps = selected === undefined ? [] : selected.rows.map((row) => this.hspEntry(loaded, row));
    const ranges = selected === undefined ? [] : this.rangeEntries(loaded, unfiltered.get(selected.sIdx)!, new Set(selected.rows));
    this.set({
      subjects,
      hsps,
      ranges,
      hidden: { subjects: filtered.hiddenSubjects, hsps: filtered.hiddenHsps },
      queryTotals: { subjects: query.subjects.length, hsps: query.hspCount },
      atSubjectLimit: loaded.limits.maxTargetSeqs !== undefined && query.subjects.length >= loaded.limits.maxTargetSeqs,
      outfmt0Subjects: query.subjects.filter((s) => out0SubjectRange(table, s.rows[0]!) !== undefined).length,
    });
  }

  /** After a filter change, keeps the selected HSP if it is still shown, else selects the first shown. */
  private keepSelection(): void {
    const state = this.state.get();
    if (state.hsps.some((hsp) => sameHsp(hsp.id, state.hsp))) return;
    const subject = state.subjects.find((s) => s.sIdx === state.sIdx) ?? state.subjects[0];
    if (subject === undefined) {
      this.set({ sIdx: undefined, hsps: [], ranges: [], hsp: undefined, detail: undefined });
      return;
    }
    this.selectSubject(subject.sIdx);
  }

  private subjectEntry(
    loaded: LoadedRun,
    shown: SubjectGroup,
    all: SubjectGroup,
    order: number,
    hspSort: SortSpec<HspSortKey>,
  ): SubjectEntry {
    const record: RecordKey | undefined = loaded.run.snapshot.subject.records[shown.sIdx];
    const first = all.rows[0]!;
    return {
      sIdx: shown.sIdx,
      recordId: record?.id ?? '',
      length: record?.length ?? 0,
      order,
      first: this.row(first),
      hspCount: all.rows.length,
      rows: sortHsps(shown.rows, hspSort, loaded.index.table),
      inOutfmt0: out0SubjectRange(loaded.index.table, first) !== undefined,
      atHspLimit: loaded.limits.maxHsps !== undefined && all.rows.length >= loaded.limits.maxHsps,
    };
  }

  /** The subject's HSPs in the engine's order, numbered among all of them, that the view filters show. */
  private rangeEntries(loaded: LoadedRun, subject: SubjectGroup, shown: ReadonlySet<number>): readonly RangeEntry[] {
    const { table } = loaded.index;
    const runId = loaded.run.snapshot.runId;
    return subject.rows.flatMap((row, i) => {
      if (!shown.has(row)) return [];
      const [from, to] = subjectSpan(table, row);
      return [{ id: { runId, qIdx: table.qIdx[row]!, rank: table.rank[row]! }, row, n: i + 1, from, to, inOutfmt0: out0Range(table, row) !== undefined }];
    });
  }

  private hspEntry(loaded: LoadedRun, row: number): HspEntry {
    const { table } = loaded.index;
    // Frames belong to the translated roles only (domain/coordinates.ts `translates`; TBLASTN
    // and TBLASTX): the engine's BLASTP records carry frame 1 for both sequences.
    const translated = (kind: SequenceKind) => translates(loaded.run.snapshot.program, kind);
    const queryFrame = translated(loaded.kinds.query) ? frame(table.queryFrame, row) : undefined;
    const subjectFrame = translated(loaded.kinds.subject) ? frame(table.subjectFrame, row) : undefined;
    return {
      row,
      id: { runId: loaded.run.snapshot.runId, qIdx: table.qIdx[row]!, rank: table.rank[row]! },
      fields: this.row(row),
      orientation: orientation(table, row, loaded.kinds),
      ...(queryFrame === undefined ? {} : { queryFrame }),
      ...(subjectFrame === undefined ? {} : { subjectFrame }),
      inOutfmt0: out0Range(table, row) !== undefined,
    };
  }

  private queryEntries(loaded: LoadedRun, filters: ViewFilters): readonly QueryEntry[] {
    const text = (filters.queryText ?? '').toLowerCase();
    const entries: QueryEntry[] = [];
    loaded.run.snapshot.query.records.forEach((record, qIdx) => {
      const group = loaded.index.queries.get(qIdx);
      if (filters.queriesWithHitsOnly && group === undefined) return;
      if (text !== '' && !record.id.toLowerCase().includes(text)) return;
      entries.push({ qIdx, id: record.id, length: record.length, subjects: group?.subjects.length ?? 0, hsps: group?.hspCount ?? 0 });
    });
    return entries;
  }

  /** Reads the HSP's row, subject heading and section; drops the answer if the selection moved on. */
  private async loadDetail(loaded: LoadedRun, id: HspId, row: number): Promise<void> {
    const table = loaded.index.table;
    const range6 = out6Range(table, row);
    const rowText = range6 === undefined ? undefined : decoder.decode(loaded.out6.subarray(range6[0], range6[1]));
    this.set({ detail: { id, state: 'loading', ...(rowText === undefined ? {} : { row: rowText }) } });
    try {
      const [heading, section] = await Promise.all([this.readHeading(loaded, table.sIdx[row]!, row), this.readSection(id)]);
      if (!sameHsp(this.state.get().hsp, id)) return;
      this.set({
        detail: {
          id,
          state: 'ready',
          ...(rowText === undefined ? {} : { row: rowText }),
          ...(heading === undefined ? {} : { heading }),
          ...(section === undefined ? {} : { section }),
        },
      });
    } catch (error) {
      if (!sameHsp(this.state.get().hsp, id)) return;
      this.set({ detail: { id, state: 'failed', error: errorMessage(error) } });
    }
  }

  private set(change: { [K in keyof ResultsState]?: ResultsState[K] | undefined }): void {
    const next = { ...this.state.get() } as Record<string, unknown>;
    for (const [key, value] of Object.entries(change)) {
      if (value === undefined) delete next[key];
      else next[key] = value;
    }
    this.state.set(next as unknown as ResultsState);
  }
}

/** The task of a run: the argv's -task, else the default that the engine describes (none for a program without tasks). */
function taskOf(argv: readonly string[], description: ProgramDescription): { task?: string } {
  const task = optionValue(argv, '-task') ?? description.parameters.find((p) => p.flag === '-task')?.defaultValue;
  return task === undefined ? {} : { task };
}

/** The option grammar of `describe` that the verification badge reads options with. */
function grammarOf(description: ProgramDescription) {
  const valued = new Set(description.parameters.filter((p) => p.takesValue).map((p) => p.flag));
  const task = description.parameters.find((p) => p.flag === '-task')?.defaultValue;
  return { takesValue: (flag: string) => valued.has(flag), ...(task === undefined ? {} : { defaultTask: task }) };
}

function unavailableMessage(run: RunView): string {
  const n = run.snapshot.number;
  switch (run.status) {
    case 'cancelled':
      return `Run ${n} was cancelled. A cancelled run keeps no results.`;
    case 'failed':
      return `Run ${n} failed: ${run.record.error ?? 'no reason was recorded'}. A failed run keeps no results.`;
    case 'completed':
      return '';
    default:
      return `Run ${n} has not finished yet. Its results appear when it completes.`;
  }
}

function errorMessage(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
