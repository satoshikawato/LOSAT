// The HSP records of a run and the rows and sections of its outputs (S13, plan §4.4-§4.5,
// TD-3): for every outfmt 0 fixture of the programs that LOSAT Web runs
// (LOSAT/tests/outfmt0_manifest.tsv; BLASTX joins in SX), the serial reactor runs the
// search in Node, and
// - the stream 0, 6 or 7 bytes are NCBI's frozen bytes of the fixture (a fixture of an
//   outfmt 6 or 7 field list only adds its search; for the approved
//   -db_gencode exception, NCBI's database search outside the lines that name the database,
//   as docs/evidence/losat_web_e2a/check_losat.py compares them),
// - the HSP records' `out6` ranges tile the outfmt 6 text in the order of the records, and
//   each row's coordinates are the record's,
// - outfmt 0 has as many sections and subject headings as the records point to,
// - each `out0` range is the HSP's section (its score lines and its alignment, whose first
//   and last Query and Sbjct coordinates are the record's), in the order of the records,
// - each `out0_subject` range is the subject's heading, the same for the HSPs of one subject
//   in one query, and before their sections,
// - an HSP has no section exactly when its subject has no heading, and such subjects come
//   after the subjects that outfmt 0 shows (-num_alignments),
// - the results screen's index (domain/result-index.ts) groups the same HSPs.
// It needs LOSAT_WEB_REACTORS; the gate runs it (docs/evidence/losat_web_w4/run_gate.sh).
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { findReactors } from '../../build/reactors';
import { hspTable } from '../../src/domain/hsp-table';
import { splitOutfmt6Row } from '../../src/domain/outfmt6';
import { buildResultIndex } from '../../src/domain/result-index';
import { ROLE_QUERY, ROLE_SUBJECT, type ReactorAbi } from '../../src/infra/reactor/abi';
import { instantiateSerial } from '../../src/infra/reactor/instance';
import type { HspRecord } from '../../src/ports/engine';

const ENGINE = resolve(import.meta.dirname, '../../../../LOSAT');
const PROGRAMS = new Set(['blastn', 'blastp', 'tblastn', 'tblastx']);

interface Fixture {
  readonly id: string;
  readonly program: string;
  readonly words: readonly string[];
  readonly query: string;
  readonly subject: string;
  /**
   * The stream that the fixture froze (0, 6 or 7) and NCBI's SHA-256 of it; for the approved
   * -db_gencode exception, the file of NCBI's database search (`<fixture_id>.db.out`).
   * Undefined for an outfmt 6 or 7 field list, which the reactor's streams do not have.
   */
  readonly frozen: { readonly stream: 0 | 6 | 7; readonly sha256: string; readonly database?: string } | undefined;
}

/** The fixtures of the manifest, one search for each (several fixtures may freeze the same search). */
function fixtures(): Map<string, Fixture[]> {
  const lines = readFileSync(join(ENGINE, 'tests/outfmt0_manifest.tsv'), 'utf8').split('\n').filter((l) => l !== '' && !l.startsWith('#'));
  const header = lines[0]!.split('\t');
  const searches = new Map<string, Fixture[]>();
  for (const line of lines.slice(1)) {
    const row = Object.fromEntries(line.split('\t').map((value, i) => [header[i]!, value])) as Record<string, string>;
    const program = row['program']!;
    const contract = row['contract'] ?? '';
    if (!PROGRAMS.has(program) || (contract !== '' && contract !== 'approved_db_gencode_deviation')) continue;
    if ((row['extra_args'] ?? '').includes("'") || (row['extra_args'] ?? '').includes('"')) continue;
    const words = [...(row['task'] ? ['-task', row['task']] : []), ...(row['extra_args'] ?? '').split(/\s+/).filter((w) => w !== '')];
    if (words.includes('-outfmt') || words.includes('-out')) continue;
    const stream = ({ '': 0, '0': 0, '6': 6, '7': 7 } as const)[row['outfmt'] ?? ''];
    const fixture: Fixture = {
      id: row['fixture_id']!,
      program,
      words,
      query: row['query']!,
      subject: row['subject']!,
      // The approved -db_gencode exception compares with NCBI's -db oracle (AGENTS.md).
      frozen:
        stream === undefined
          ? undefined
          : contract === ''
          ? { stream, sha256: row['stdout_sha256']! }
          : { stream, sha256: row['db_stdout_sha256']!, database: `tests/fixtures/outfmt0/${row['fixture_id']!}.db.out` },
    };
    const key = JSON.stringify([program, fixture.query, fixture.subject, words]);
    searches.set(key, [...(searches.get(key) ?? []), fixture]);
  }
  return searches;
}

interface Outputs {
  readonly streams: Map<number, Uint8Array>;
  readonly records: HspRecord[];
}

function search(abi: ReactorAbi, fixture: Fixture): Outputs {
  const parts = new Map<number, Uint8Array[]>();
  const query = abi.register(fixture.program, ROLE_QUERY, readFileSync(join(ENGINE, fixture.query)));
  const subject = abi.register(fixture.program, ROLE_SUBJECT, readFileSync(join(ENGINE, fixture.subject)));
  try {
    const argv = [fixture.program, '-query', fixture.query, '-subject', fixture.subject, ...fixture.words, '-num_threads', '1'];
    abi.run(argv, query.handle, subject.handle, (stream, bytes) => parts.set(stream, [...(parts.get(stream) ?? []), bytes.slice()]));
  } finally {
    abi.release(query.handle);
    abi.release(subject.handle);
  }
  const streams = new Map<number, Uint8Array>();
  for (const [stream, chunks] of parts) streams.set(stream, Buffer.concat(chunks));
  const text = new TextDecoder().decode(streams.get(1) ?? new Uint8Array());
  const records = text.split('\n').filter((l) => l.trim() !== '').map((l) => JSON.parse(l) as HspRecord);
  return { streams, records };
}

const sha256 = (bytes: Uint8Array) => createHash('sha256').update(bytes).digest('hex');

/**
 * check_losat.py's normalize_database_lines: a report without the lines that name the
 * subjects (the database title of the outfmt 0 prolog up to the counts line and of the
 * epilog through `Posted date:`, outfmt 7's `# Database:` line) and with `> ` in every
 * subject heading (NCBI writes a database subject as `>title`).
 */
function withoutDatabaseLines(report: string): string {
  const lines = report.split('\n');
  const kept: string[] = [];
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i]!;
    if (line.startsWith('Database: ')) {
      while (i < lines.length && !/^ +\S+ sequences; /.test(lines[i]!)) i++;
      kept.push('<database>', lines[i] ?? '');
    } else if (line.startsWith('  Database: ')) {
      while (i < lines.length && !lines[i]!.startsWith('    Posted date:')) i++;
      kept.push('<database>');
    } else if (line.startsWith('# Database: ')) {
      kept.push('<database>');
    } else if (line.startsWith('>') && !line.startsWith('> ')) {
      kept.push(`> ${line.slice(1)}`);
    } else {
      kept.push(line);
    }
  }
  return kept.join('\n');
}

/** Whether a search's stream is the fixture's frozen NCBI bytes (see the file comment); what it compared. */
function expectFrozen(fixture: Fixture, outputs: Outputs): 'bytes' | 'database' | 'none' {
  if (fixture.frozen === undefined) return 'none';
  const { stream, sha256: frozen, database } = fixture.frozen;
  const bytes = outputs.streams.get(stream) ?? new Uint8Array();
  const label = `${fixture.id}: NCBI's frozen outfmt ${stream}`;
  if (database === undefined) {
    expect(sha256(bytes), label).toBe(frozen);
    return 'bytes';
  }
  const oracle = readFileSync(join(ENGINE, database));
  // run_oracle.py hashes the file with the time makeblastdb ran replaced (without_posted_date).
  const posted = Buffer.from(oracle.toString('latin1').replace(/^( {4}Posted date: {2}).*$/gm, '$1(the time makeblastdb ran)'), 'latin1');
  expect(sha256(posted), `${fixture.id}: ${database}`).toBe(frozen);
  const decoder = new TextDecoder();
  expect(withoutDatabaseLines(decoder.decode(bytes)), `${label} (database search, -db_gencode exception)`).toBe(
    withoutDatabaseLines(decoder.decode(oracle)),
  );
  return 'database';
}

/** The first and last coordinates of the alignment lines that start with `label`. */
function alignmentEnds(section: string, label: 'Query' | 'Sbjct'): [number, number] {
  const lines = section.split('\n').filter((line) => line.startsWith(`${label} `));
  const parse = (line: string) => line.trim().split(/\s+/);
  const first = parse(lines[0]!);
  const last = parse(lines.at(-1)!);
  return [Number(first[1]), Number(last.at(-1))];
}

/** Checks one search and returns what it saw, for the summary. */
function check(fixture: Fixture, outputs: Outputs): { hsps: number; sections: number; unshown: number } {
  const decoder = new TextDecoder();
  const out6 = outputs.streams.get(6) ?? new Uint8Array();
  const out0 = outputs.streams.get(0) ?? new Uint8Array();
  const records = [...outputs.records].sort((a, b) => a.index - b.index);
  const where = (record: HspRecord) => `${fixture.id} HSP ${record.index}`;

  // outfmt 6: one row per record, in the order of the records, covering the text exactly.
  const rows = decoder.decode(out6).split('\n').filter((l) => l !== '');
  expect(records.length, `${fixture.id}: records and outfmt 6 rows`).toBe(rows.length);
  let at = 0;
  for (const record of records) {
    expect(record.out6, where(record)).not.toBeNull();
    const [start, end] = record.out6!;
    expect(start, `${where(record)} out6 start`).toBe(at);
    const fields = splitOutfmt6Row(decoder.decode(out6.subarray(start, end)));
    expect([fields.qstart, fields.qend, fields.sstart, fields.send], `${where(record)} coordinates`).toEqual(
      [record.q_start, record.q_end, record.s_start, record.s_end].map(String),
    );
    at = end;
  }
  expect(at, `${fixture.id}: the rows end with the text`).toBe(out6.length);

  // outfmt 0: sections in the order of the records; headings before their subjects' sections.
  let previousEnd = 0;
  const headings = new Map<string, readonly [number, number] | null>();
  let sections = 0;
  let unshown = 0;
  for (const record of records) {
    const pair = `${record.q_idx}:${record.s_idx}`;
    expect(record.out0 === null, `${where(record)} section and heading`).toBe(record.out0_subject === null);
    if (headings.has(pair)) expect(record.out0_subject, `${where(record)} heading`).toEqual(headings.get(pair));
    else headings.set(pair, record.out0_subject);
    if (record.out0 === null) {
      unshown++;
      continue;
    }
    sections++;
    const [start, end] = record.out0;
    expect(start, `${where(record)} section order`).toBeGreaterThanOrEqual(previousEnd);
    const [hStart, hEnd] = record.out0_subject!;
    expect(hEnd, `${where(record)} heading before section`).toBeLessThanOrEqual(start);
    const heading = decoder.decode(out0.subarray(hStart, hEnd));
    expect(heading.startsWith('>'), `${where(record)} heading text`).toBe(true);
    expect(heading, `${where(record)} heading Length=`).toMatch(/\nLength=\d+\n/);
    const section = decoder.decode(out0.subarray(start, end));
    expect(section.startsWith(' Score ='), `${where(record)} section text`).toBe(true);
    expect(alignmentEnds(section, 'Query'), `${where(record)} Query lines`).toEqual([record.q_start, record.q_end]);
    expect(alignmentEnds(section, 'Sbjct'), `${where(record)} Sbjct lines`).toEqual([record.s_start, record.s_end]);
    previousEnd = end;
  }
  // The text has as many sections and subject headings as the records point to: no section or
  // heading of outfmt 0 lacks its records.
  const text = decoder.decode(out0);
  expect(text.match(/^ Score =/gm)?.length ?? 0, `${fixture.id}: outfmt 0 sections and records with one`).toBe(sections);
  const pointed = new Set([...headings.values()].flatMap((range) => (range === null ? [] : [range[0]])));
  expect(text.match(/^>/gm)?.length ?? 0, `${fixture.id}: outfmt 0 subject headings and the records' headings`).toBe(pointed.size);

  // Subjects without alignments in outfmt 0 come after those with them, in each query.
  const table = hspTable(outputs.records);
  const index = buildResultIndex(table);
  let grouped = 0;
  for (const query of index.queries.values()) {
    const shown = query.subjects.map((subject) => table.out0SubjectStart[subject.rows[0]!]! >= 0);
    const firstUnshown = shown.indexOf(false);
    if (firstUnshown >= 0) expect(shown.slice(firstUnshown).every((s) => !s), `${fixture.id} q${query.qIdx}: unshown subjects last`).toBe(true);
    for (const subject of query.subjects) {
      grouped += subject.rows.length;
      const ranks = subject.rows.map((row) => table.rank[row]!);
      expect(ranks, `${fixture.id} q${query.qIdx} s${subject.sIdx}: rank order`).toEqual([...ranks].sort((a, b) => a - b));
    }
  }
  expect(grouped, `${fixture.id}: every HSP is in the index`).toBe(records.length);
  return { hsps: records.length, sections, unshown };
}

const reactors = findReactors();

describe.skipIf(reactors === undefined)('HSP records and the rows and sections of outfmt 6 and 0 (every outfmt 0 fixture)', () => {
  const searches = reactors === undefined ? new Map<string, Fixture[]>() : fixtures();
  const module = reactors === undefined ? undefined : new WebAssembly.Module(reactors.serial.bytes as BufferSource);
  const summary = { searches: 0, hsps: 0, sections: 0, unshown: 0, programs: new Set<string>(), frozen: { bytes: 0, database: 0, none: 0 } };

  it(
    'agree for every search',
    async () => {
      let abi: ReactorAbi | undefined;
      let runs = 0;
      for (const group of searches.values()) {
        // A fresh instance now and then keeps the instance's memory small.
        if (abi === undefined || runs++ % 20 === 0) abi = (await instantiateSerial(module!)).abi;
        const fixture = group[0]!;
        const outputs = search(abi, fixture);
        for (const frozen of group) summary.frozen[expectFrozen(frozen, outputs)]++;
        const seen = check(fixture, outputs);
        summary.searches++;
        summary.hsps += seen.hsps;
        summary.sections += seen.sections;
        summary.unshown += seen.unshown;
        summary.programs.add(fixture.program);
      }
      console.log(
        `HSP correspondence: ${summary.searches} searches (${[...summary.programs].sort().join(', ')}), ${summary.hsps} HSPs, ` +
          `${summary.sections} outfmt 0 sections, ${summary.unshown} HSPs that outfmt 0 does not show; fixtures equal to NCBI: ` +
          `${summary.frozen.bytes} byte for byte, ${summary.frozen.database} as the database search (-db_gencode exception), ` +
          `${summary.frozen.none} field lists not compared`,
      );
      expect([...summary.programs].sort()).toEqual(['blastn', 'blastp', 'tblastn', 'tblastx']);
      expect(summary.unshown, 'the fixtures include HSPs that outfmt 0 does not show').toBeGreaterThan(0);
      expect(summary.frozen.bytes, 'fixtures compared byte for byte').toBeGreaterThan(0);
      expect(summary.frozen.database, 'fixtures of the -db_gencode exception compared').toBeGreaterThan(0);
    },
    1_800_000,
  );
});
