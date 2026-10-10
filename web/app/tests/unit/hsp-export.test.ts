// LOSAT Web's own files of a run's HSPs (S15 item 2): the scopes' rows, CSV quoting and columns,
// the JSON pieces joined into one valid document, and the report's escaping, checked with hostile
// IDs and names against an allow-list of the report's own tags and attributes.
import { describe, expect, it } from 'vitest';
import {
  allRows,
  csvField,
  csvHeader,
  csvLine,
  csvRecord,
  CSV_COLUMNS,
  engineOrder,
  filteredRows,
  filterWords,
  JSON_FORMAT,
  JSON_TAIL,
  jsonHead,
  jsonHsp,
  markedRows,
  type ExportedHsp,
  type ExportedRecord,
  type ExportRun,
  type ExportScopeInfo,
} from '../../src/domain/hsp-export';
import {
  escapeHtml,
  REPORT_ALIGNMENTS_START,
  REPORT_CSP,
  REPORT_QUERY_END,
  REPORT_TABLE_END,
  reportHead,
  reportHeading,
  reportNotInOutfmt0,
  reportQueryStart,
  reportSection,
  reportTableRow,
  reportTail,
} from '../../src/domain/hsp-report';
import { hspTable, type HspSummary } from '../../src/domain/hsp-table';
import { splitOutfmt6Row } from '../../src/domain/outfmt6';
import { buildResultIndex } from '../../src/domain/result-index';

/** IDs and names that a page must never read as markup or a URL. */
const HOSTILE = [
  '<script>alert(1)</script>',
  '"><img src=x onerror=alert(2)>',
  "'><svg/onload=alert(3)>",
  'javascript:alert(4)',
  '</style><style>body{background:url(https://example.invalid/x)}</style>',
  '</pre></table><a href="https://example.invalid/">x</a>',
  'a&amp;b',
];

function summary(index: number, q: number, s: number, rank: number, e = 1e-5, bits = 50): HspSummary {
  return {
    index,
    q_idx: q,
    s_idx: s,
    rank,
    raw_score: bits * 2,
    bit_score: bits,
    e_value: e,
    q_start: 1,
    q_end: 10,
    s_start: 1,
    s_end: 10,
    query_frame: null,
    subject_frame: null,
    out6: [0, 1],
    out0: null,
    out0_subject: null,
  };
}

function hsp(fields: Partial<ExportedHsp> = {}): ExportedHsp {
  return {
    run: 3,
    qIdx: 0,
    sIdx: 1,
    queryId: 'q1',
    subjectId: 's2',
    index: 0,
    rank: 0,
    outfmt6: splitOutfmt6Row('q1\ts2\t98.684\t76\t1\t0\t1\t76\t1\t76\t2.37e-38\t138\n'),
    queryFrame: null,
    subjectFrame: null,
    inOutfmt0: true,
    ...fields,
  };
}

const RECORD: ExportedRecord = {
  index: 0,
  q_idx: 0,
  s_idx: 1,
  rank: 0,
  raw_score: 276,
  bit_score: 138.0291,
  e_value: 2.3672e-38,
  q_start: 1,
  q_end: 76,
  s_start: 1,
  s_end: 76,
  query_frame: null,
  subject_frame: null,
  subject_length: 500,
  query_aligned: 'ACGT-A',
  subject_aligned: 'ACGTTA',
};

function run(fields: Partial<ExportRun> = {}): ExportRun {
  return {
    number: 3,
    title: 'Run title',
    program: 'blastn',
    programLabel: 'BLASTN',
    task: 'megablast',
    argv: ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-evalue', '1e-5'],
    requestedThreads: 'auto',
    engineBuild: 'build-1',
    runtimePath: 'threaded',
    threads: 4,
    queuedAt: Date.UTC(2026, 9, 10, 1, 2, 3),
    startedAt: Date.UTC(2026, 9, 10, 1, 2, 4),
    endedAt: Date.UTC(2026, 9, 10, 1, 2, 5),
    query: { name: 'q.fa', records: 2, bytes: 120, sha256: 'a'.repeat(64) },
    subject: { name: 's.fa', records: 3, bytes: 300, sha256: 'b'.repeat(64) },
    verification: { level: 'outside', label: 'Engine-supported, outside certified profile', details: ['Why.'], exceptions: [] },
    ...fields,
  };
}

describe('the scopes of an export', () => {
  // Rows in another order than the HSP index; query 0 has subjects 1 and 0, query 1 subject 2.
  const records = [summary(2, 0, 0, 2, 1e-3, 30), summary(0, 0, 1, 0, 1e-20, 90), summary(3, 1, 2, 0, 1e-10, 60), summary(1, 0, 1, 1, 1e-5, 40)];
  const table = hspTable(records);
  const index = buildResultIndex(table);
  const indices = (rows: Int32Array) => Array.from(rows, (row) => table.index[row]);

  it('puts rows in the engine order', () => {
    expect(indices(allRows(table))).toEqual([0, 1, 2, 3]);
    expect(indices(engineOrder(table, [2, 0, 1]))).toEqual([0, 2, 3]);
  });

  it('keeps the HSPs that the view filters keep, in the listed queries', () => {
    const ids = (subject: { sIdx: number }) => [`s${subject.sIdx}`];
    expect(indices(filteredRows(index, [0, 1], {}, ids))).toEqual([0, 1, 2, 3]);
    expect(indices(filteredRows(index, [0, 1], { maxEValue: 1e-4 }, ids))).toEqual([0, 1, 3]);
    expect(indices(filteredRows(index, [0, 1], { minBitScore: 50 }, ids))).toEqual([0, 3]);
    expect(indices(filteredRows(index, [0, 1], { subjectText: 'S1' }, ids))).toEqual([0, 1]);
    // A query that the query filters do not list keeps none of its HSPs; one without HSPs adds none.
    expect(indices(filteredRows(index, [1, 5], {}, ids))).toEqual([3]);
  });

  it('takes the marked subjects of the list, in the engine order whatever the list sort', () => {
    const rowOf = (i: number) => [...table.index].indexOf(i);
    const listed = [
      { sIdx: 0, rows: [rowOf(2)] },
      { sIdx: 1, rows: [rowOf(1), rowOf(0)] },
    ];
    expect(indices(markedRows(table, listed, new Set([1])))).toEqual([0, 1]);
    expect(indices(markedRows(table, listed, new Set([0, 1])))).toEqual([0, 1, 2]);
    expect(indices(markedRows(table, listed, new Set([7])))).toEqual([]);
  });

  it('names the filters in words', () => {
    expect(filterWords({})).toEqual([]);
    expect(filterWords({ maxEValue: 1e-5, minBitScore: 40, subjectText: 'abc', queriesWithHitsOnly: true, queryText: '' })).toEqual([
      'E value at most 0.00001',
      'Bit score at least 40',
      'Subject ID contains "abc"',
      'Queries with hits only',
    ]);
  });
});

describe('CSV', () => {
  it('quotes a field with a comma, a quote, CR or LF, and doubles its quotes (RFC 4180)', () => {
    expect(csvField('plain')).toBe('plain');
    expect(csvField('a,b')).toBe('"a,b"');
    expect(csvField('say "hi"')).toBe('"say ""hi"""');
    expect(csvField('line\nbreak')).toBe('"line\nbreak"');
    expect(csvField('cr\rhere')).toBe('"cr\rhere"');
    expect(csvField('')).toBe('');
    expect(csvRecord(['a', 'b,c', ''])).toBe('a,"b,c",\r\n');
  });

  it('writes a header row and one row per HSP with the outfmt 6 fields as written', () => {
    expect(csvHeader()).toBe(
      'run,query_record,query_id,subject_record,subject_id,hsp,index,rank,qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,query_frame,subject_frame,in_outfmt0\r\n',
    );
    expect(csvLine(hsp({ qIdx: 1, rank: 4, index: 9, queryFrame: 2, subjectFrame: -1, inOutfmt0: false }))).toBe(
      '3,2,q1,2,s2,2.5,9,4,q1,s2,98.684,76,1,0,1,76,1,76,2.37e-38,138,2,-1,false\r\n',
    );
    // An HSP that outfmt 6 does not show has empty fields there.
    expect(csvLine(hsp({ outfmt6: undefined })).split(',').slice(8, 20)).toEqual(Array(12).fill(''));
  });

  it('writes IDs as they are, quoted where they must be, never altered', () => {
    const line = csvLine(hsp({ queryId: '=HYPERLINK("x")', subjectId: 'a,"b"\nc' }));
    expect(line.startsWith('3,1,"=HYPERLINK(""x"")",2,"a,""b""\nc",1.1,')).toBe(true);
    expect(CSV_COLUMNS).toHaveLength(23);
  });
});

describe('JSON', () => {
  const scope: ExportScopeInfo = {
    scope: 'marked',
    hsps: 2,
    filters: { maxEValue: 1e-5, subjectText: 'x' },
    query: { position: 0, id: 'q1' },
    markedSubjects: [{ position: 1, id: HOSTILE[0]! }],
  };

  it('joins its pieces into one valid document that says what it is', () => {
    const text =
      jsonHead(run({ title: HOSTILE[1]! }), scope, { aligned: true, exportedAt: Date.UTC(2026, 9, 10) }) +
      jsonHsp(hsp({ queryId: HOSTILE[2]! }), RECORD, true, true) +
      jsonHsp(hsp({ index: 1, rank: 1, outfmt6: undefined }), { ...RECORD, index: 1, rank: 1, query_aligned: null, subject_aligned: null }, true, false) +
      JSON_TAIL;
    const doc = JSON.parse(text);
    expect(doc.format).toBe(JSON_FORMAT);
    expect(doc.schema).toBe(1);
    expect(doc.note).toMatch(/LOSAT Web/);
    expect(doc.note).toMatch(/not an NCBI BLAST output/);
    expect(doc.exported_at).toBe('2026-10-10T00:00:00.000Z');
    expect(doc.run).toMatchObject({
      number: 3,
      title: HOSTILE[1],
      program: 'blastn',
      task: 'megablast',
      argv: run().argv,
      requested_threads: 'auto',
      engine_build: 'build-1',
      threads: 4,
      queued_at: '2026-10-10T01:02:03.000Z',
      ended_at: '2026-10-10T01:02:05.000Z',
      query: { name: 'q.fa', records: 2, engine_input_bytes: 120, engine_input_sha256: 'a'.repeat(64) },
      verification: { label: 'Engine-supported, outside certified profile', exceptions: [] },
    });
    expect(doc.scope).toEqual({
      name: 'marked',
      label: 'Marked subjects',
      hsps: 2,
      filters: { max_e_value: 1e-5, subject_id_contains: 'x' },
      query: { record: 1, id: 'q1' },
      marked_subjects: [{ record: 2, id: HOSTILE[0] }],
      aligned_sequences: true,
    });
    expect(doc.hsps).toHaveLength(2);
    expect(doc.hsps[0]).toEqual({
      run: 3,
      query_record: 1,
      query_id: HOSTILE[2],
      subject_record: 2,
      subject_id: 's2',
      hsp: '1.1',
      index: 0,
      rank: 0,
      in_outfmt0: true,
      outfmt6: splitOutfmt6Row('q1\ts2\t98.684\t76\t1\t0\t1\t76\t1\t76\t2.37e-38\t138'),
      record: { ...RECORD },
    });
    expect(doc.hsps[1].outfmt6).toBeNull();
    expect(doc.hsps[1].record.query_aligned).toBeNull();
    // One HSP per line after the head.
    expect(text.split('\n')).toHaveLength(5);
  });

  it('leaves the aligned rows out when they are not chosen, and writes an empty scope as an empty array', () => {
    const doc = JSON.parse(jsonHead(run(), { scope: 'all', hsps: 1 }, { aligned: false, exportedAt: 0 }) + jsonHsp(hsp(), RECORD, false, true) + JSON_TAIL);
    expect(doc.scope).toMatchObject({ name: 'all', filters: null, aligned_sequences: false });
    expect(Object.keys(doc.hsps[0].record)).not.toContain('query_aligned');
    expect(doc.hsps[0].record.e_value).toBe(2.3672e-38);
    expect(JSON.parse(jsonHead(run(), { scope: 'all', hsps: 0 }, { aligned: true, exportedAt: 0 }) + JSON_TAIL).hsps).toEqual([]);
  });
});

// --- the report ------------------------------------------------------------------------------------

/** The report's own tags (an allow-list): nothing that loads, runs, links or submits. */
const ALLOWED_TAGS = new Set([
  'html', 'head', 'meta', 'title', 'style', 'body', 'header', 'footer', 'section', 'h1', 'h2', 'h3', 'p', 'strong', 'span', 'br',
  'dl', 'dt', 'dd', 'ul', 'li', 'code', 'pre', 'div', 'table', 'thead', 'tbody', 'tr', 'th', 'td',
]);
const ALLOWED_ATTRIBUTES = new Set(['lang', 'charset', 'http-equiv', 'content', 'name', 'class']);

interface Tag {
  readonly name: string;
  readonly attributes: readonly (readonly [string, string])[];
}

/** Every tag of an HTML text, with its attributes (the report quotes every attribute value). */
function tagsOf(html: string): Tag[] {
  const tags: Tag[] = [];
  for (const match of html.matchAll(/<(\/?)([A-Za-z!][^\s/>]*)([^>]*)>/g)) {
    const attributes = [...match[3]!.matchAll(/([^\s=/]+)(?:\s*=\s*"([^"]*)")?/g)].map((a) => [a[1]!.toLowerCase(), a[2] ?? ''] as const);
    tags.push({ name: match[2]!.toLowerCase(), attributes });
  }
  return tags;
}

function report(ids: { query: string; subject: string; name: string; title: string; argv: string; out: string; warnings: string }): string {
  const hostileRun = run({
    title: ids.title,
    argv: ['blastn', '-query', ids.name, '-subject', 's.fa', '-entrez_query', ids.argv],
    query: { name: ids.name, records: 1, bytes: 1, sha256: ids.name },
    engineBuild: ids.out,
    verification: { level: 'outside', label: ids.title, details: [ids.argv], exceptions: [ids.out] },
  });
  const row = hsp({
    queryId: ids.query,
    subjectId: ids.subject,
    outfmt6: splitOutfmt6Row(`${ids.query}\t${ids.subject}\t98.684\t76\t1\t0\t1\t76\t1\t76\t2.37e-38\t138`),
  });
  return [
    reportHead({
      run: hostileRun,
      formats: [0, 6, 7],
      scope: { scope: 'marked', hsps: 1, filters: { subjectText: ids.subject, queryText: ids.query }, query: { position: 0, id: ids.query }, markedSubjects: [{ position: 1, id: ids.subject }] },
      exportedAt: 0,
    }),
    reportQueryStart({ position: 0, id: ids.query, length: 76, unit: 'nt', hsps: 1 }),
    reportTableRow(row),
    REPORT_TABLE_END,
    REPORT_ALIGNMENTS_START,
    reportHeading(`> ${ids.subject} ${ids.out}\nLength=500\n`),
    reportSection(row, ` Score = 138 bits\n${ids.out}\n`),
    reportNotInOutfmt0(['1.2']),
    REPORT_QUERY_END,
    reportTail(`Warning: ${ids.warnings}\n`),
  ].join('');
}

describe('the report', () => {
  it('escapes & < > " and \'', () => {
    expect(escapeHtml(`<a href="x" title='y'>&amp;</a>`)).toBe('&lt;a href=&quot;x&quot; title=&#39;y&#39;&gt;&amp;amp;&lt;/a&gt;');
    expect(escapeHtml('plain text')).toBe('plain text');
  });

  it('is a static page that loads nothing and says what it is', () => {
    const html = report({ query: 'q1', subject: 's2', name: 'q.fa', title: 'Title', argv: '-x', out: 'out', warnings: 'none' });
    expect(html.startsWith('<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta http-equiv="Content-Security-Policy"')).toBe(true);
    expect(html).toContain(`content="${REPORT_CSP}"`);
    expect(REPORT_CSP).toMatch(/^default-src 'none'; style-src 'unsafe-inline'/);
    expect(html).toContain('<h1>LOSAT Web report</h1>');
    expect(html).toMatch(/An application format of LOSAT Web/);
    expect(html).toMatch(/not an NCBI BLAST report/);
    expect(html).toContain('LOSAT blastn -query q.fa -subject s.fa -entrez_query -x -outfmt 6');
    expect(html).not.toMatch(/url\(|@import|https?:|<script|<img|<link|<a\s|<form|<iframe/i);
  });

  it('writes every text of the run as text, whatever it holds', () => {
    const benign = report({ query: 'q1', subject: 's2', name: 'q.fa', title: 'Title', argv: '-x', out: 'out', warnings: 'none' });
    for (const text of HOSTILE) {
      const html = report({ query: text, subject: text, name: text, title: text, argv: text, out: text, warnings: text });
      const tags = tagsOf(html);
      for (const tag of tags) {
        expect(ALLOWED_TAGS.has(tag.name) || tag.name === '!doctype', `tag <${tag.name}> of ${text}`).toBe(true);
        for (const [name, value] of tag.attributes) {
          if (tag.name === '!doctype') continue;
          expect(ALLOWED_ATTRIBUTES.has(name), `attribute ${name} of ${text}`).toBe(true);
          expect(name === 'src' || name === 'href' || name === 'style' || name.startsWith('on')).toBe(false);
          expect(value).not.toMatch(/url\(|javascript:|https?:/i);
        }
      }
      // The hostile texts add no tag and no attribute: the page's markup is that of benign IDs.
      expect(tags.map((tag) => `${tag.name}${JSON.stringify(tag.attributes)}`)).toEqual(
        tagsOf(benign).map((tag) => `${tag.name}${JSON.stringify(tag.attributes)}`),
      );
      // The text itself is kept: an escaped copy for each place it was written.
      expect(html).toContain(escapeHtml(text));
      expect(html.match(/<style>/g)).toHaveLength(1);
      expect(html.match(/<\/style>/g)).toHaveLength(1);
    }
  });

  it('keeps a text that starts with a line end as written inside <pre> (the parser drops one after the tag)', () => {
    expect(reportSection(hsp(), '\n Score = 1\n')).toBe('<p class="hsp-label">HSP 1.1</p>\n<pre class="section">\n\n Score = 1\n</pre>\n');
    expect(reportHeading('> s1\nLength=5\n')).toBe('<pre class="heading">\n&gt; s1\nLength=5\n</pre>\n');
  });

  it('lists the HSPs that outfmt 0 does not show, and says when there are no warnings', () => {
    expect(reportNotInOutfmt0([])).toBe('');
    expect(reportNotInOutfmt0(['1.2', '1.3'])).toContain('outfmt 0 does not show HSPs 1.2, 1.3.');
    expect(reportTail('')).toContain('The engine wrote no warnings.');
  });
});
