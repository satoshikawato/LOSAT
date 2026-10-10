// The exporter of LOSAT Web's own files (S15 item 2) against a RunStore double and the results
// screen's state: the scopes and their counts, CSV, JSON and the report written in blocks through
// the Writer contract, the HSP records read in batches and outfmt 0 in windows, the scope fixed when
// the export starts, one export at a time, and nothing saved when a read fails.
import { describe, expect, it } from 'vitest';
import type { AppState, RunView } from '../../src/application/coordinator';
import { RECORD_BATCH, ResultExporter, exportFileName, recordBatches } from '../../src/application/result-export';
import { ResultsBrowser, type ResultsState } from '../../src/application/results';
import { Store } from '../../src/application/store';
import { CSV_COLUMNS } from '../../src/domain/hsp-export';
import { hspTable } from '../../src/domain/hsp-table';
import type { RunSnapshot } from '../../src/domain/run';
import type { HspRecord, ProgramDescription } from '../../src/ports/engine';
import { memoryDownloader, type SavedFile } from './support/memory-downloader';

const DESCRIPTION: ProgramDescription = { program: 'blastn', formats: [0, 6, 7], parameters: [] };
const encoder = new TextEncoder();
const decoder = new TextDecoder();

interface Spec {
  readonly q: number;
  readonly s: number;
  readonly bits: number;
  readonly e: number;
  /** Whether outfmt 0 shows the HSP (default true). */
  readonly shown?: boolean;
}

/**
 * A stored run of the specs, in their order (the engine's): outfmt 6 rows and outfmt 0 headings and
 * sections at the records' ranges, as the engine writes them.
 */
function makeRun(specs: readonly Spec[], ids: { readonly query: readonly string[]; readonly subject: readonly string[] }) {
  const out6: string[] = [];
  const out0: string[] = [];
  const end = { out6: 0, out0: 0 };
  const append = (stream: 'out6' | 'out0', text: string): readonly [number, number] => {
    (stream === 'out6' ? out6 : out0).push(text);
    const start = end[stream];
    end[stream] += encoder.encode(text).length;
    return [start, end[stream]];
  };
  append('out0', 'BLASTN header\n');
  const records: HspRecord[] = [];
  const ranks = new Map<number, number>();
  const headings = new Map<string, readonly [number, number]>();
  specs.forEach((spec, index) => {
    const qseqid = ids.query[spec.q] || `Query_${spec.q + 1}`;
    const sseqid = ids.subject[spec.s] || `Subject_${spec.s + 1}`;
    const row = append('out6', [qseqid, sseqid, '99.000', '50', '0', '0', 1, 50, 101, 150, String(spec.e), String(spec.bits)].join('\t') + '\n');
    const shown = spec.shown ?? true;
    const key = `${spec.q}:${spec.s}`;
    if (shown && !headings.has(key)) headings.set(key, append('out0', `> ${sseqid} subject\nLength=500\n\n`));
    const section = shown ? append('out0', ` Score = ${spec.bits} bits, Expect = ${spec.e} (HSP index ${index})\n\nQuery  1  ACGT  50\n\n`) : null;
    const rank = ranks.get(spec.q) ?? 0;
    ranks.set(spec.q, rank + 1);
    records.push({
      index,
      q_idx: spec.q,
      s_idx: spec.s,
      rank,
      raw_score: spec.bits * 2,
      bit_score: spec.bits,
      e_value: spec.e,
      q_start: 1,
      q_end: 50,
      s_start: 101,
      s_end: 150,
      query_frame: null,
      subject_frame: null,
      subject_length: 500,
      query_aligned: `QA${index}`,
      subject_aligned: `SA${index}`,
      out6: row,
      out0: section,
      out0_subject: shown ? headings.get(key)! : null,
    });
  });
  return { records, out6: encoder.encode(out6.join('')), out0: encoder.encode(out0.join('')) };
}

const IDS = { query: ['q1', '', 'q,3'], subject: ['s1', '<b>"s2"</b>', 's3,x'] };
const SPECS: readonly Spec[] = [
  { q: 0, s: 1, bits: 90, e: 1e-20 },
  { q: 0, s: 1, bits: 40, e: 1e-5 },
  { q: 0, s: 0, bits: 70, e: 1e-12 },
  { q: 0, s: 2, bits: 30, e: 0.5, shown: false },
  { q: 1, s: 0, bits: 60, e: 1e-9 },
  { q: 2, s: 2, bits: 20, e: 2 },
];

function snapshot(ids: typeof IDS): RunSnapshot {
  const input = (name: string, records: readonly string[]) => ({
    name,
    bytes: new Uint8Array(42),
    sha256: name.repeat(8),
    revisionIds: [],
    records: records.map((id) => ({ id, length: 500 })),
  });
  return {
    runId: 'r1',
    number: 7,
    program: 'blastn',
    title: 'Export test',
    argv: ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-evalue', '10'],
    query: input('q.fa', ids.query),
    subject: input('s.fa', ids.subject),
    requestedThreads: 'auto',
    queuedAt: Date.UTC(2026, 9, 10, 1, 0, 0),
  };
}

async function setup(options: { specs?: readonly Spec[]; ids?: typeof IDS } = {}) {
  const ids = options.ids ?? IDS;
  const run = makeRun(options.specs ?? SPECS, ids);
  const view: RunView = {
    snapshot: snapshot(ids),
    status: 'completed',
    record: { runtimePath: 'serial', threads: 1, engineBuild: 'test-build', startedAt: Date.UTC(2026, 9, 10, 1, 0, 1), endedAt: Date.UTC(2026, 9, 10, 1, 0, 2) },
    result: { runId: 'r1', byteLengths: { 0: run.out0.length, 6: run.out6.length, 7: 0 }, hitCount: run.records.length },
  };
  const runs = new Store<AppState>({ runs: [view] });
  const results = new ResultsBrowser({
    data: {
      readHitTable: async () => hspTable(run.records),
      readOutput: async (_id, format) => (format === 6 ? run.out6.slice() : run.out0.slice()),
      readOutputRange: async (_id, _format, start, end) => run.out0.slice(start, end),
      readDiagnostics: async () => 'Warning: <check> & "this"\n',
    },
    describe: async () => DESCRIPTION,
    runs,
    verification: { ncbi: '2.17.0', sources: [], programs: {} },
  });
  await results.open('r1');
  const saved: SavedFile[] = [];
  const recordReads: number[][] = [];
  const rangeReads: [number, number][] = [];
  let failRecordsAt: number | undefined;
  let failRanges = false;
  let hold: Promise<void> | undefined;
  const exporter = new ResultExporter({
    results: results.state,
    data: {
      readHspRecords: async (_id, indices) => {
        recordReads.push([...indices]);
        if (hold !== undefined) await hold;
        if (failRecordsAt !== undefined && recordReads.length >= failRecordsAt) throw new Error('the stored records are gone');
        return indices.map((index) => run.records[index]!);
      },
      readOutputRange: async (_id, format, start, end) => {
        expect(format).toBe(0);
        if (end > run.out0.length) throw new RangeError(`bytes [${start}, ${end}) are not in outfmt 0`);
        rangeReads.push([start, end]);
        if (failRanges) throw new Error('the stored output is gone');
        return run.out0.slice(start, end);
      },
    },
    downloader: memoryDownloader((file) => saved.push(file)),
    now: () => Date.UTC(2026, 9, 10, 2, 0, 0),
  });
  return {
    run,
    results,
    exporter,
    saved,
    recordReads,
    rangeReads,
    text: (i = saved.length - 1) => decoder.decode(saved[i]!.bytes),
    failRecords: (at: number) => (failRecordsAt = at),
    failRanges: () => (failRanges = true),
    holdRecords: () => {
      let release!: () => void;
      hold = new Promise((resolve) => (release = resolve));
      return () => {
        hold = undefined;
        release();
      };
    },
  };
}

/** The rows of a CSV text without its header (the fields have no line end of their own here). */
const csvRows = (text: string) => text.split('\r\n').slice(1, -1);

/** The records of a CSV text (RFC 4180), the header's first: a reader independent of the writer. */
function parseCsv(text: string): string[][] {
  const records: string[][] = [];
  let fields: string[] = [];
  let i = 0;
  while (i < text.length) {
    let field = '';
    if (text[i] === '"') {
      i++;
      for (;;) {
        const quote = text.indexOf('"', i);
        field += text.slice(i, quote);
        i = quote + 1;
        if (text[i] !== '"') break;
        field += '"';
        i++;
      }
    } else {
      while (i < text.length && text[i] !== ',' && text[i] !== '\r') field += text[i++];
    }
    fields.push(field);
    if (text[i] === ',') i++;
    else {
      expect(text.slice(i, i + 2)).toBe('\r\n');
      i += 2;
      records.push(fields);
      fields = [];
    }
  }
  return records;
}

/** The HSP labels of a CSV text's rows. */
const csvLabels = (text: string) => parseCsv(text).slice(1).map((fields) => fields[5]);

describe('the scopes and their counts', () => {
  it('counts the whole run, the HSPs after the view filters, and those of the marked subjects', async () => {
    const { results, exporter } = await setup();
    expect(exporter.counts(results.state.get())).toEqual({ all: 6, filtered: 6, marked: 0 });
    results.setFilters({ maxEValue: 1e-6 });
    expect(exporter.counts(results.state.get())).toEqual({ all: 6, filtered: 3, marked: 0 });
    results.setFilters({ maxEValue: 1e-6, queryText: 'q1' });
    expect(exporter.counts(results.state.get())).toEqual({ all: 6, filtered: 2, marked: 0 });
    // The marks are the selected query's listed subjects; the filters apply to their HSPs.
    results.markSubjects([1, 0], true);
    expect(exporter.counts(results.state.get())).toEqual({ all: 6, filtered: 2, marked: 2 });
    results.setFilters({});
    expect(exporter.counts(results.state.get()).marked).toBe(3);
  });

  it('counts nothing without a loaded run', () => {
    const exporter = new ResultExporter({
      results: new Store<ResultsState>({ phase: 'none' } as ResultsState),
      data: { readHspRecords: async () => [], readOutputRange: async () => new Uint8Array() },
      downloader: memoryDownloader(() => undefined),
      now: () => 0,
    });
    expect(exporter.counts({ phase: 'none' } as ResultsState)).toEqual({ all: 0, filtered: 0, marked: 0 });
  });
});

describe('CSV', () => {
  it('writes the whole run with the record IDs, the HSP labels and the outfmt 6 fields as written', async () => {
    const { exporter, saved, text, recordReads, rangeReads } = await setup();
    const summary = await exporter.export('csv', 'all');
    expect(saved).toHaveLength(1);
    expect(saved[0]).toMatchObject({ name: 'losat-run7-blastn-hsps.csv', mime: 'text/csv' });
    expect(summary).toEqual({ format: 'csv', scope: 'all', fileName: 'losat-run7-blastn-hsps.csv', hsps: 6, bytes: saved[0]!.bytes.length });
    const rows = csvRows(text());
    expect(rows).toEqual([
      '7,1,q1,2,"<b>""s2""</b>",1.1,0,0,q1,"<b>""s2""</b>",99.000,50,0,0,1,50,101,150,1e-20,90,,,true',
      '7,1,q1,2,"<b>""s2""</b>",1.2,1,1,q1,"<b>""s2""</b>",99.000,50,0,0,1,50,101,150,0.00001,40,,,true',
      '7,1,q1,1,s1,1.3,2,2,q1,s1,99.000,50,0,0,1,50,101,150,1e-12,70,,,true',
      '7,1,q1,3,"s3,x",1.4,3,3,q1,"s3,x",99.000,50,0,0,1,50,101,150,0.5,30,,,false',
      '7,2,,1,s1,2.1,4,0,Query_2,s1,99.000,50,0,0,1,50,101,150,1e-9,60,,,true',
      '7,3,"q,3",3,"s3,x",3.1,5,0,"q,3","s3,x",99.000,50,0,0,1,50,101,150,2,20,,,true',
    ]);
    expect(parseCsv(text())[0]).toEqual([...CSV_COLUMNS]);
    expect(parseCsv(text())[1]![4]).toBe('<b>"s2"</b>');
    // Everything comes from what the results screen holds: no read of the Data worker.
    expect(recordReads).toEqual([]);
    expect(rangeReads).toEqual([]);
    expect(exporter.state.get()).toEqual({ last: summary });
  });

  it('writes the filtered and the marked scopes in the engine order', async () => {
    const { results, exporter, text } = await setup();
    results.setFilters({ maxEValue: 1e-6 });
    await exporter.export('csv', 'filtered');
    expect(csvLabels(text())).toEqual(['1.1', '1.3', '2.1']);
    results.setFilters({});
    results.setSubjectSort({ key: 'bitScore', descending: false });
    results.markSubjects([2, 1], true);
    await exporter.export('csv', 'marked');
    expect(csvLabels(text())).toEqual(['1.1', '1.2', '1.4']);
    expect(parseCsv(text())[3]!.slice(3, 5)).toEqual(['3', 's3,x']);
  });

  it('writes a large scope in more than one block', async () => {
    const specs = Array.from({ length: 12_000 }, (_, i) => ({ q: Math.floor(i / 40), s: i % 3, bits: 50, e: 1e-5 }));
    const ids = { query: Array.from({ length: 300 }, (_, q) => `query_${q}`), subject: ['s1', 's2', 's3'] };
    const { exporter, saved, text } = await setup({ specs, ids });
    await exporter.export('csv', 'all');
    expect(saved[0]!.blocks).toBeGreaterThan(1);
    expect(csvRows(text())).toHaveLength(12_000);
  });
});

describe('JSON', () => {
  it('writes one valid document: the run, the scope, then each HSP with its record and aligned rows', async () => {
    const { results, exporter, text, saved } = await setup();
    results.markSubjects([1], true);
    await exporter.export('json', 'marked');
    expect(saved[0]).toMatchObject({ name: 'losat-run7-blastn-hsps.json', mime: 'application/json' });
    const doc = JSON.parse(text());
    expect(doc).toMatchObject({ format: 'LOSAT Web HSP export', schema: 1, exported_at: '2026-10-10T02:00:00.000Z' });
    expect(doc.note).toMatch(/not an NCBI BLAST output/);
    expect(doc.run).toMatchObject({
      number: 7,
      title: 'Export test',
      program: 'blastn',
      argv: ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-evalue', '10'],
      engine_build: 'test-build',
      runtime_path: 'serial',
      threads: 1,
      requested_threads: 'auto',
      started_at: '2026-10-10T01:00:01.000Z',
      query: { name: 'q.fa', records: 3, engine_input_bytes: 42, engine_input_sha256: 'q.fa'.repeat(8) },
      subject: { name: 's.fa', records: 3 },
    });
    expect(doc.scope).toEqual({
      name: 'marked',
      label: 'Marked subjects',
      hsps: 2,
      filters: {},
      query: { record: 1, id: 'q1' },
      marked_subjects: [{ record: 2, id: '<b>"s2"</b>' }],
      aligned_sequences: true,
    });
    expect(doc.hsps.map((hsp: { hsp: string }) => hsp.hsp)).toEqual(['1.1', '1.2']);
    expect(doc.hsps[1]).toMatchObject({
      run: 7,
      query_record: 1,
      query_id: 'q1',
      subject_record: 2,
      subject_id: '<b>"s2"</b>',
      index: 1,
      rank: 1,
      in_outfmt0: true,
      outfmt6: { qseqid: 'q1', sseqid: '<b>"s2"</b>', evalue: '0.00001', bitscore: '40' },
      record: { index: 1, q_idx: 0, s_idx: 1, rank: 1, raw_score: 80, bit_score: 40, e_value: 1e-5, subject_length: 500, query_aligned: 'QA1', subject_aligned: 'SA1' },
    });
  });

  it('leaves out the aligned rows when they are not chosen', async () => {
    const { exporter, text } = await setup();
    await exporter.export('json', 'all', { aligned: false });
    const doc = JSON.parse(text());
    expect(doc.scope).toMatchObject({ name: 'all', filters: null, aligned_sequences: false });
    expect(doc.hsps).toHaveLength(6);
    expect(doc.hsps[0].record.query_aligned).toBeUndefined();
  });

  it('reads the HSP records in batches, never one by one', async () => {
    const specs = Array.from({ length: 2_500 }, (_, i) => ({ q: Math.floor(i / 25), s: i % 3, bits: 50, e: 1e-5 }));
    const ids = { query: Array.from({ length: 100 }, (_, q) => `query_${q}`), subject: ['s1', 's2', 's3'] };
    const { exporter, text, recordReads } = await setup({ specs, ids });
    await exporter.export('json', 'all');
    expect(recordReads.map((batch) => batch.length)).toEqual([RECORD_BATCH, RECORD_BATCH, 500]);
    expect(recordReads.flat()).toEqual(Array.from({ length: 2_500 }, (_, i) => i));
    expect(JSON.parse(text()).hsps).toHaveLength(2_500);
  });

  it('saves nothing when a read fails, and says why', async () => {
    const specs = Array.from({ length: 2_500 }, (_, i) => ({ q: Math.floor(i / 25), s: i % 3, bits: 50, e: 1e-5 }));
    const ids = { query: Array.from({ length: 100 }, (_, q) => `query_${q}`), subject: ['s1', 's2', 's3'] };
    const { exporter, saved, failRecords } = await setup({ specs, ids });
    failRecords(2);
    expect(await exporter.export('json', 'all')).toBeUndefined();
    expect(saved).toEqual([]);
    expect(exporter.state.get().error).toBe('The JSON of run 7 could not be written, so nothing was saved: the stored records are gone');
    expect(exporter.state.get().busy).toBeUndefined();
  });
});

describe('the report', () => {
  it('writes the run, the scope, each query with its HSP table and outfmt 0 text, and the warnings, all escaped', async () => {
    const { results, exporter, text, saved, run } = await setup();
    results.setFilters({ maxEValue: 1 });
    await exporter.export('report', 'filtered');
    expect(saved[0]).toMatchObject({ name: 'losat-run7-blastn-report.html', mime: 'text/html' });
    const html = text();
    expect(html).toMatch(/^<!DOCTYPE html>/);
    expect(html).toContain("default-src 'none'; style-src 'unsafe-inline'");
    expect(html).toContain('not an NCBI BLAST report');
    expect(html).toContain('<dt>HSPs</dt><dd>After the view filters: 5 HSPs</dd>');
    expect(html).toContain('<dt>View filters</dt><dd>E value at most 1</dd>');
    expect(html).toContain('LOSAT blastn -query q.fa -subject s.fa -evalue 10 -outfmt 0');
    expect(html).toContain('Engine-supported, outside certified profile');
    expect(html).toContain('<h2>Query 1: q1</h2>');
    expect(html).toContain('<h2>Query 2</h2>');
    expect(html).not.toContain('Query 3');
    // The subject's heading and sections of outfmt 0, as written, escaped.
    const out0 = decoder.decode(run.out0);
    const heading = out0.slice(out0.indexOf('> <b>'), out0.indexOf('Length=500\n\n') + 'Length=500\n\n'.length);
    expect(html).toContain(`<pre class="heading">\n${heading.replaceAll('<', '&lt;').replaceAll('>', '&gt;').replaceAll('"', '&quot;')}</pre>`);
    expect(html).toContain('<td>&lt;b&gt;&quot;s2&quot;&lt;/b&gt;</td>');
    expect(html).not.toContain('<b>');
    expect(html.match(/<pre class="section">/g)).toHaveLength(4);
    expect(html).toContain('outfmt 0 does not show HSP 1.4.');
    expect(html).toContain('<pre>\nWarning: &lt;check&gt; &amp; &quot;this&quot;\n</pre>');
    expect(html.trimEnd().endsWith('</html>')).toBe(true);
  });

  it('reads outfmt 0 in windows that serve many sections, never past its end', async () => {
    const specs = Array.from({ length: 3_000 }, (_, i) => ({ q: Math.floor(i / 30), s: i % 3, bits: 50, e: 1e-5 }));
    const ids = { query: Array.from({ length: 100 }, (_, q) => `query_${q}`), subject: ['s1', 's2', 's3'] };
    const { exporter, text, rangeReads, run } = await setup({ specs, ids });
    await exporter.export('report', 'all');
    expect(text().match(/<pre class="section">/g)).toHaveLength(3_000);
    expect(rangeReads.length).toBeLessThanOrEqual(Math.ceil(run.out0.length / (1 << 20)) + 1);
    expect(Math.max(...rangeReads.map(([, end]) => end))).toBeLessThanOrEqual(run.out0.length);
  });

  it('saves nothing when outfmt 0 cannot be read', async () => {
    const { exporter, saved, failRanges } = await setup();
    failRanges();
    expect(await exporter.export('report', 'all')).toBeUndefined();
    expect(saved).toEqual([]);
    expect(exporter.state.get().error).toMatch(/^The report of run 7 could not be written, so nothing was saved: the stored output is gone/);
  });
});

describe('an export in progress', () => {
  it('keeps the scope it started with, and runs alone', async () => {
    const { results, exporter, text, holdRecords } = await setup();
    results.setFilters({ maxEValue: 1e-6 });
    const release = holdRecords();
    const writing = exporter.export('json', 'filtered');
    await Promise.resolve();
    expect(exporter.state.get().busy).toEqual({ format: 'json', scope: 'filtered', fileName: 'losat-run7-blastn-hsps.json' });
    expect(await exporter.export('csv', 'all')).toBeUndefined();
    results.setFilters({});
    release();
    await writing;
    const doc = JSON.parse(text());
    expect(doc.scope.filters).toEqual({ max_e_value: 1e-6 });
    expect(doc.hsps.map((hsp: { hsp: string }) => hsp.hsp)).toEqual(['1.1', '1.3', '2.1']);
    expect(exporter.state.get().last?.hsps).toBe(3);
  });

  it('refuses an empty scope with a message', async () => {
    const { exporter, saved } = await setup();
    expect(await exporter.export('csv', 'marked')).toBeUndefined();
    expect(saved).toEqual([]);
    expect(exporter.state.get().error).toBe('Marked subjects: there are no HSPs to export.');
  });
});

describe('file names and batches', () => {
  it('names files after the run and the program, never an input', () => {
    expect(exportFileName({ number: 12, program: 'tblastx' }, 'report')).toBe('losat-run12-tblastx-report.html');
  });

  it('ends a batch of records early when the HSPs span many residues', () => {
    const table = hspTable(
      Array.from({ length: 5 }, (_, i) => ({
        index: i,
        q_idx: 0,
        s_idx: 0,
        rank: i,
        raw_score: 1,
        bit_score: 1,
        e_value: 1,
        q_start: 1,
        q_end: 1_500_000,
        s_start: 1,
        s_end: 1_500_000,
        query_frame: null,
        subject_frame: null,
        out6: null,
        out0: null,
        out0_subject: null,
      })),
    );
    const batches = [...recordBatches(Int32Array.from([0, 1, 2, 3, 4]), table)].map((batch) => [...batch]);
    expect(batches).toEqual([[0], [1], [2], [3], [4]]);
  });
});
