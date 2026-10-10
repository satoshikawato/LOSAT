// LOSAT Web's own files of the run that the results screen shows (S15 item 2; design §12.1,
// §12.3): CSV, JSON and the static HTML report (domain/hsp-export.ts, domain/hsp-report.ts), for
// the whole run, after the view filters, or for the marked subjects. They are kept apart from the
// compatibility outputs (outfmt 0/6/7 byte for byte, the whole run), which the view never changes.
//
// The scope is taken when the export starts: a filter or a mark changed while the file is written
// does not change it. Every file is written through `writeFile` in blocks, in order, and a failure
// at any point saves nothing (application/export-writer.ts). The values come from the stored
// outputs and the HSP records: the outfmt 6 text that the results screen already holds, the HSP
// records read from the Data worker in batches (never the whole run at once, never one read per
// HSP), and the outfmt 0 text read in windows of about 1 MiB that serve many sections each. The
// report is made in parts of about TEXT_CHARS characters and waits only for those parts and for
// the windows of outfmt 0: an `await` for each of its million small pieces made the report of
// 100,000 queries take 16-25 s in Firefox (fix round 2).
import { hspLabel } from '../domain/extraction';
import {
  allRows,
  csvHeader,
  csvLine,
  filteredRows,
  JSON_TAIL,
  jsonHead,
  jsonHsp,
  markedRows,
  SCOPE_LABELS,
  type ExportedHsp,
  type ExportRun,
  type ExportScope,
  type ExportScopeInfo,
} from '../domain/hsp-export';
import {
  REPORT_ALIGNMENTS_START,
  REPORT_QUERY_END,
  REPORT_TABLE_END,
  reportHead,
  reportHeading,
  reportNotInOutfmt0,
  reportQueryStart,
  reportSection,
  reportTableRow,
  reportTail,
} from '../domain/hsp-report';
import { frame, out0Range, out0SubjectRange, out6Range, type ByteRange } from '../domain/hsp-table';
import { splitOutfmt6Row, type Outfmt6Row } from '../domain/outfmt6';
import { programById } from '../domain/programs';
import type { RunStore } from '../ports/data';
import type { Downloader } from '../ports/download';
import { writeFile, type ExportWriter } from './export-writer';
import type { LoadedRun, ResultsState } from './results';
import { Store } from './store';

export type ExportFormat = 'csv' | 'json' | 'report';

/** File name ending and media type of each format. File names never hold an input's name. */
export const EXPORT_FILES: Readonly<Record<ExportFormat, { readonly ending: string; readonly mime: string; readonly label: string }>> = Object.freeze({
  csv: { ending: 'hsps.csv', mime: 'text/csv', label: 'CSV' },
  json: { ending: 'hsps.json', mime: 'application/json', label: 'JSON' },
  report: { ending: 'report.html', mime: 'text/html', label: 'report' },
});

export const exportFileName = (run: { readonly number: number; readonly program: string }, format: ExportFormat): string =>
  `losat-run${run.number}-${run.program}-${EXPORT_FILES[format].ending}`;

/** HSP records read from the Data worker per call. */
export const RECORD_BATCH = 1000;
/** A batch also ends when the HSPs' coordinates span this many residues, so that long aligned rows stay few per batch. */
export const RECORD_BATCH_RESIDUES = 4_000_000;
/** Rows of CSV handed to the Writer at once. */
const ROWS_PER_TEXT = 500;
/** Characters of a report collected before they are handed to the Writer. */
const TEXT_CHARS = 1 << 16;
/**
 * HSPs of a JSON batch formatted in one task. A batch of 1,000 with their aligned rows, formatted at
 * once, held Chromium's page for up to 300 ms where the garbage collector marked a large heap in
 * steps as the text was made (fix round 2); a pause after each group lets the page draw.
 */
const JSON_GROUP = 250;
/** Bytes of outfmt 0 read at once: one read serves the headings and sections that lie in it. */
export const OUTFMT0_WINDOW = 1 << 20;

export type ScopeCounts = Readonly<Record<ExportScope, number>>;

export interface ExportOptions {
  /** JSON: write the HSP records' aligned rows (default true). */
  readonly aligned?: boolean;
  /** Report: write the outfmt 0 headings and sections of the HSPs (default true). */
  readonly alignments?: boolean;
}

export interface ExportSummary {
  readonly format: ExportFormat;
  readonly scope: ExportScope;
  readonly fileName: string;
  readonly hsps: number;
  readonly bytes: number;
  /** A report saved without the outfmt 0 text of its alignments. */
  readonly withoutAlignments?: true;
}

export interface ExportState {
  /** The file being written. */
  readonly busy?: { readonly format: ExportFormat; readonly scope: ExportScope; readonly fileName: string };
  /** The last file saved. */
  readonly last?: ExportSummary;
  /** Why the last export saved nothing. */
  readonly error?: string;
}

export interface ResultExporterDeps {
  /** What the results screen shows: the loaded run, the view filters and the marks. */
  readonly results: Store<ResultsState>;
  readonly data: Pick<RunStore, 'readHspRecords' | 'readOutputRange'>;
  readonly downloader: Pick<Downloader, 'open'>;
  readonly now: () => number;
  /** Lets the page draw between two parts of a file's work (the composition gives the browser's next task). None: no pause. */
  readonly pause?: () => Promise<void>;
}

/** What one export reads: fixed when it starts. */
interface ExportJob {
  readonly loaded: LoadedRun;
  readonly rows: Int32Array;
  readonly scope: ExportScopeInfo;
  readonly exportedAt: number;
}

const decoder = new TextDecoder();

export class ResultExporter {
  readonly state = new Store<ExportState>({});
  private filtered?: { readonly key: readonly unknown[]; readonly rows: Int32Array };
  private marked?: { readonly key: readonly unknown[]; readonly rows: Int32Array };

  constructor(private readonly deps: ResultExporterDeps) {}

  /** The HSPs of each scope for the screen's state (none without a loaded run). */
  counts(results: ResultsState): ScopeCounts {
    const loaded = results.loaded;
    if (loaded === undefined) return { all: 0, filtered: 0, marked: 0 };
    return { all: loaded.index.table.count, filtered: this.rows(results, 'filtered').length, marked: this.rows(results, 'marked').length };
  }

  /** The table rows of a scope, in the engine's order. The filtered and marked rows are kept until their inputs change. */
  rows(results: ResultsState, scope: ExportScope): Int32Array {
    const loaded = results.loaded;
    if (loaded === undefined) return new Int32Array(0);
    const { table } = loaded.index;
    if (scope === 'all') return allRows(table);
    if (scope === 'filtered') {
      const key = [loaded, results.filters, results.queries];
      if (!sameKey(this.filtered?.key, key)) {
        const subjectRecords = loaded.run.snapshot.subject.records;
        const rows = filteredRows(
          loaded.index,
          results.queries.map((query) => query.qIdx),
          results.filters,
          (subject) => [subjectRecords[subject.sIdx]?.id ?? '', outfmt6Row(loaded, subject.rows[0]!)?.sseqid ?? ''],
        );
        this.filtered = { key, rows };
      }
      return this.filtered!.rows;
    }
    const key = [loaded, results.qIdx, results.subjects, results.marked];
    if (!sameKey(this.marked?.key, key)) this.marked = { key, rows: markedRows(table, results.subjects, results.marked) };
    return this.marked!.rows;
  }

  /**
   * Writes a file of the run that the results screen shows, for a scope of its HSPs, and saves it.
   * One export at a time: a call while one is written is ignored. Resolves with what was saved, or
   * undefined when nothing was (the state's error says why).
   */
  async export(format: ExportFormat, scope: ExportScope, options: ExportOptions = {}): Promise<ExportSummary | undefined> {
    if (this.state.get().busy !== undefined) return undefined;
    const results = this.deps.results.get();
    const loaded = results.loaded;
    if (results.phase !== 'ready' || loaded === undefined) {
      this.state.set({ error: 'No results are shown, so there is nothing to export.' });
      return undefined;
    }
    const rows = this.rows(results, scope);
    if (rows.length === 0) {
      this.state.set({ error: `${SCOPE_LABELS[scope]}: there are no HSPs to export.` });
      return undefined;
    }
    const job: ExportJob = { loaded, rows, scope: scopeInfo(loaded, results, scope, rows.length), exportedAt: this.deps.now() };
    const { snapshot } = loaded.run;
    const fileName = exportFileName(snapshot, format);
    const aligned = options.aligned ?? true;
    const alignments = options.alignments ?? true;
    this.state.set({ busy: { format, scope, fileName } });
    try {
      const bytes = await writeFile(this.deps.downloader, fileName, EXPORT_FILES[format].mime, (writer) => {
        switch (format) {
          case 'csv':
            return writeCsv(writer, job);
          case 'json':
            return this.writeJson(writer, job, aligned);
          case 'report':
            return this.writeReport(writer, job, alignments);
        }
      });
      const last: ExportSummary = {
        format,
        scope,
        fileName,
        hsps: rows.length,
        bytes,
        ...(format === 'report' && !alignments ? { withoutAlignments: true as const } : {}),
      };
      this.state.set({ last });
      return last;
    } catch (error) {
      this.state.set({
        error: `The ${EXPORT_FILES[format].label} of run ${snapshot.number} could not be written, so nothing was saved: ${errorMessage(error)}`,
      });
      return undefined;
    }
  }

  /** The JSON document: the run, the scope, then each HSP with its record, read in batches. */
  private async writeJson(writer: ExportWriter, job: ExportJob, aligned: boolean): Promise<void> {
    const { loaded, rows } = job;
    const { table } = loaded.index;
    const runId = loaded.run.snapshot.runId;
    await writer.text(jsonHead(exportRun(loaded), job.scope, { aligned, exportedAt: job.exportedAt }));
    let first = true;
    for (const batch of recordBatches(rows, table)) {
      const records = await this.deps.data.readHspRecords(runId, Array.from(batch, (row) => table.index[row]!));
      for (let start = 0; start < batch.length; start += JSON_GROUP) {
        if (start > 0) await this.deps.pause?.();
        const parts: string[] = [];
        for (let i = start; i < Math.min(batch.length, start + JSON_GROUP); i++) {
          const row = batch[i]!;
          const record = records[i];
          if (record === undefined || record.index !== table.index[row]) {
            throw new Error(`the Data worker returned another HSP record than HSP ${table.index[row]} of run ${loaded.run.snapshot.number}`);
          }
          parts.push(jsonHsp(exportedHsp(loaded, row), record, aligned, first));
          first = false;
        }
        await writer.texts(parts);
      }
    }
    await writer.text(JSON_TAIL);
  }

  /**
   * The report: the run, then each query of the scope with its HSP table and, unless they are left
   * out, the outfmt 0 headings and sections of the HSPs that outfmt 0 shows, then the run's warnings.
   */
  private async writeReport(writer: ExportWriter, job: ExportJob, alignments: boolean): Promise<void> {
    const { loaded, rows } = job;
    const { table } = loaded.index;
    const { snapshot } = loaded.run;
    const out0 = new RangeReader(
      (start, end) => this.deps.data.readOutputRange(snapshot.runId, 0, start, end),
      loaded.run.result?.byteLengths[0],
    );
    const head = { run: exportRun(loaded), formats: loaded.description.formats, scope: job.scope, exportedAt: job.exportedAt, alignments };
    let parts: string[] = [];
    let chars = 0;
    const add = (part: string) => {
      parts.push(part);
      chars += part.length;
    };
    const handOn = async () => {
      await writer.texts(parts);
      parts = [];
      chars = 0;
    };
    add(reportHead(head));
    for (const [qIdx, queryRows] of groupBy(rows, table.qIdx)) {
      const record = snapshot.query.records[qIdx];
      add(reportQueryStart({ position: qIdx, id: record?.id ?? '', length: record?.length ?? 0, unit: loaded.units.query, hsps: queryRows.length }));
      for (const row of queryRows) {
        add(reportTableRow(exportedHsp(loaded, row)));
        if (chars >= TEXT_CHARS) await handOn();
      }
      add(REPORT_TABLE_END);
      if (alignments) {
        const shown = queryRows.filter((row) => out0Range(table, row) !== undefined);
        if (shown.length > 0) add(REPORT_ALIGNMENTS_START);
        for (const [, subjectRows] of groupBy(shown, table.sIdx)) {
          const heading = subjectRows.map((row) => out0SubjectRange(table, row)).find((range) => range !== undefined);
          // A range in the window last read is taken without waiting.
          if (heading !== undefined) add(reportHeading(out0.cached(heading) ?? (await out0.text(heading))));
          for (const row of subjectRows) {
            const range = out0Range(table, row)!;
            add(reportSection(exportedHsp(loaded, row), out0.cached(range) ?? (await out0.text(range))));
            if (chars >= TEXT_CHARS) await handOn();
          }
        }
        add(reportNotInOutfmt0(queryRows.filter((row) => out0Range(table, row) === undefined).map((row) => hspLabel(qIdx, table.rank[row]!))));
      }
      add(REPORT_QUERY_END);
      if (chars >= TEXT_CHARS) await handOn();
    }
    add(reportTail(loaded.diagnostics));
    await handOn();
  }
}

/** The CSV: a header row, then one row per HSP. */
async function writeCsv(writer: ExportWriter, job: ExportJob): Promise<void> {
  await writer.text(csvHeader());
  const { loaded, rows } = job;
  for (let i = 0; i < rows.length; i += ROWS_PER_TEXT) {
    const parts: string[] = [];
    for (const row of rows.subarray(i, i + ROWS_PER_TEXT)) parts.push(csvLine(exportedHsp(loaded, row)));
    await writer.texts(parts);
  }
}

/** The identity of an HSP of the loaded run, its outfmt 6 fields as written, its frames and outfmt 0. */
function exportedHsp(loaded: LoadedRun, row: number): ExportedHsp {
  const { table } = loaded.index;
  const { snapshot } = loaded.run;
  const qIdx = table.qIdx[row]!;
  const sIdx = table.sIdx[row]!;
  return {
    run: snapshot.number,
    qIdx,
    sIdx,
    queryId: snapshot.query.records[qIdx]?.id ?? '',
    subjectId: snapshot.subject.records[sIdx]?.id ?? '',
    index: table.index[row]!,
    rank: table.rank[row]!,
    outfmt6: outfmt6Row(loaded, row),
    queryFrame: frame(table.queryFrame, row) ?? null,
    subjectFrame: frame(table.subjectFrame, row) ?? null,
    inOutfmt0: out0Range(table, row) !== undefined,
  };
}

/** The fields of an HSP's outfmt 6 row, from the outfmt 6 text that the results screen holds. */
function outfmt6Row(loaded: LoadedRun, row: number): Outfmt6Row | undefined {
  const range = out6Range(loaded.index.table, row);
  return range === undefined ? undefined : splitOutfmt6Row(decoder.decode(loaded.out6.subarray(range[0], range[1])));
}

/** What the JSON and the report say about the run (its snapshot, record and verification badge). */
function exportRun(loaded: LoadedRun): ExportRun {
  const { snapshot, record } = loaded.run;
  const input = (role: 'query' | 'subject') => ({
    name: snapshot[role].name,
    records: snapshot[role].records.length,
    // A run loaded from a session file has no input bytes here; its session gives their length.
    bytes: snapshot[role].bytes?.length ?? loaded.run.fromSession?.inputs[role].length ?? 0,
    sha256: snapshot[role].sha256,
  });
  const { badge } = loaded;
  return {
    number: snapshot.number,
    ...(snapshot.title === undefined ? {} : { title: snapshot.title }),
    program: snapshot.program,
    programLabel: programById(snapshot.program).label,
    ...(loaded.task === undefined ? {} : { task: loaded.task }),
    argv: snapshot.argv,
    requestedThreads: snapshot.requestedThreads,
    ...(record.engineBuild === undefined ? {} : { engineBuild: record.engineBuild }),
    ...(record.runtimePath === undefined ? {} : { runtimePath: record.runtimePath }),
    ...(record.threads === undefined ? {} : { threads: record.threads }),
    queuedAt: snapshot.queuedAt,
    ...(record.startedAt === undefined ? {} : { startedAt: record.startedAt }),
    ...(record.endedAt === undefined ? {} : { endedAt: record.endedAt }),
    ...(snapshot.group === undefined ? {} : { group: { position: snapshot.group.position, size: snapshot.group.size } }),
    query: input('query'),
    subject: input('subject'),
    verification: { level: badge.level, label: badge.label, details: badge.details, exceptions: badge.exceptions },
  };
}

/** The scope as the files describe it: the view filters, and the query and subjects marked. */
function scopeInfo(loaded: LoadedRun, results: ResultsState, scope: ExportScope, hsps: number): ExportScopeInfo {
  if (scope === 'all') return { scope, hsps };
  if (scope === 'filtered') return { scope, hsps, filters: results.filters };
  const { query, subject } = loaded.run.snapshot;
  const qIdx = results.qIdx!;
  return {
    scope,
    hsps,
    filters: results.filters,
    query: { position: qIdx, id: query.records[qIdx]?.id ?? '' },
    markedSubjects: results.subjects
      .filter((entry) => results.marked.has(entry.sIdx))
      .map((entry) => ({ position: entry.sIdx, id: subject.records[entry.sIdx]?.id ?? '' })),
  };
}

/**
 * Consecutive rows read with one `readHspRecords`: at most RECORD_BATCH, and fewer when their
 * coordinates span more than RECORD_BATCH_RESIDUES residues (the aligned rows are about as long).
 */
export function* recordBatches(
  rows: Int32Array,
  table: { readonly qStart: Float64Array; readonly qEnd: Float64Array; readonly sStart: Float64Array; readonly sEnd: Float64Array },
): Generator<Int32Array> {
  let start = 0;
  let residues = 0;
  for (let i = 0; i < rows.length; i++) {
    const row = rows[i]!;
    const span = Math.abs(table.qEnd[row]! - table.qStart[row]!) + Math.abs(table.sEnd[row]! - table.sStart[row]!) + 2;
    if (i > start && (i - start >= RECORD_BATCH || residues + span > RECORD_BATCH_RESIDUES)) {
      yield rows.subarray(start, i);
      start = i;
      residues = 0;
    }
    residues += span;
  }
  if (start < rows.length) yield rows.subarray(start);
}

/** Rows grouped by a column's value, the groups in the order of their first rows, each group's rows in order. */
function groupBy(rows: ArrayLike<number> & Iterable<number>, column: Int32Array): Map<number, number[]> {
  const groups = new Map<number, number[]>();
  for (const row of rows) {
    const key = column[row]!;
    let group = groups.get(key);
    if (group === undefined) groups.set(key, (group = []));
    group.push(row);
  }
  return groups;
}

/**
 * Reads byte ranges of a stored output that come mostly in increasing order (outfmt 0's headings
 * and sections in the engine's order): one read of OUTFMT0_WINDOW bytes, or of the range if it is
 * longer, serves every range that lies in it.
 */
export class RangeReader {
  private start = 0;
  private bytes: Uint8Array = new Uint8Array(0);

  constructor(
    private readonly read: (start: number, end: number) => Promise<Uint8Array>,
    /** The output's length, which no read passes; unknown, a read ends at its range. */
    private readonly length: number | undefined,
    private readonly window = OUTFMT0_WINDOW,
  ) {}

  /** The text of a range that the last read holds, without reading; undefined if it does not hold it. */
  cached(range: ByteRange): string | undefined {
    const [start, end] = range;
    if (start < this.start || end > this.start + this.bytes.length) return undefined;
    return decoder.decode(this.bytes.subarray(start - this.start, end - this.start));
  }

  async text(range: ByteRange): Promise<string> {
    const [start, end] = range;
    if (start < this.start || end > this.start + this.bytes.length) {
      const limit = Math.max(end, Math.min(this.length ?? end, start + this.window));
      this.bytes = await this.read(start, limit);
      this.start = start;
      if (this.bytes.length < end - start) throw new Error(`outfmt 0 bytes [${start}, ${end}) could not be read`);
    }
    return decoder.decode(this.bytes.subarray(start - this.start, end - this.start));
  }
}

const sameKey = (a: readonly unknown[] | undefined, b: readonly unknown[]): boolean =>
  a !== undefined && a.length === b.length && a.every((value, i) => value === b[i]);

function errorMessage(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
