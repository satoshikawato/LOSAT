// The HSP records of a run and the rows and sections of its outputs (S13, plan §4.4-§4.5,
// TD-3): for every outfmt 0 fixture of the programs that LOSAT Web runs
// (LOSAT/tests/outfmt0_manifest.tsv; BLASTX joins in SX), the serial reactor runs the
// search in Node, and
// - the stream 0, 6 or 7 bytes are NCBI's frozen bytes of the fixture in LOSAT's words (LOSAT
//   states its own version, LOSAT/src/wording.rs; the frozen hash judges LOSAT's output read back
//   in NCBI's words, `losat_wording.py hash-form`) (a fixture of an outfmt 6 or 7 field list only
//   adds its search; for the approved -db_gencode exception, NCBI's database search in LOSAT's
//   words outside the lines that name the database, as docs/evidence/losat_web_e2a/check_losat.py
//   compares them),
// - the HSP records' `out6` ranges tile the outfmt 6 text in the order of the records, and
//   each row's coordinates are the record's,
// - outfmt 0 has as many sections and subject headings as the records point to, except in a run
//   of two searches (-num_descriptions / -num_alignments that give outfmt 0 another hit list size
//   than outfmt 6 and 7, abi_v2.md §7): the records are the tabular list, and the sections and
//   headings that no record points to belong to pairs whose records are absent or differ,
// - each `out0` range is the HSP's section (its score lines and its alignment, whose first
//   and last Query and Sbjct coordinates are the record's), in the order of the records,
// - each `out0_subject` range is the subject's heading, the same for the HSPs of one subject
//   in one query, and before their sections,
// - an HSP has no section exactly when its subject has no heading, and such subjects come
//   after the subjects that outfmt 0 shows (-num_alignments; in a run of two searches a pair is
//   linked only where the outfmt 0 search has exactly its HSPs, so there is no such order),
// - the results screen's index (domain/result-index.ts) groups the same HSPs.
// It needs LOSAT_WEB_REACTORS; the gate runs it (docs/evidence/losat_web_w4/run_gate.sh).
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { instantiateSerial } from 'losat';
import type { HspRecord } from 'losat/ports';
import { findReactors } from '../../build/reactors';
import { hspTable } from '../../src/domain/hsp-table';
import { htmlShownText } from '../../src/domain/outfmt0-html';
import { splitOutfmt6Row } from '../../src/domain/outfmt6';
import { buildResultIndex } from '../../src/domain/result-index';
import { ROLE_QUERY, ROLE_SUBJECT, type ReactorAbi } from '../../../packages/losat/src/infra/reactor/abi';

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
    // LOSAT Web reads FASTA inputs only: a PSSM query (-in_pssm, the manifest's `pssm:` and `pssm-gz:`
    // cells) is the CLI's.
    if (/^pssm(-gz)?:/.test(row['query']!)) continue;
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
  /** The ids of the registered records (the `id` of the register responses). */
  readonly queryIds: readonly string[];
  readonly subjectIds: readonly string[];
  /** The IDs that the outputs show for the subjects, where the engine read the deflines' ids (-parse_deflines). */
  readonly shownSubjectIds?: readonly string[];
}

/**
 * The hit list size that the engine derives for the output format `format` from the options of a run
 * (NCBI `blast_args.cpp:2894-2973`, the adapter's `hit_list_size`): written out here, not read from the
 * app's `hitLimits`, which the results screen uses, so that a wrong `separate` there cannot select the
 * weaker rule below. outfmt 0 (any format up to 4): -max_target_seqs, or the larger of -num_descriptions
 * (500) and -num_alignments (250) where either is given; the tabular formats: -max_target_seqs, else
 * -num_alignments, else 500.
 */
function hitListSize(words: readonly string[], format: 0 | 6): number {
  const value = (flag: string): number | undefined => {
    const at = words.lastIndexOf(flag);
    return at >= 0 ? Number(words[at + 1]) : undefined;
  };
  const [maxTargetSeqs, descriptions, alignments] = [value('-max_target_seqs'), value('-num_descriptions'), value('-num_alignments')];
  if (format === 0) {
    return descriptions !== undefined || alignments !== undefined ? Math.max(descriptions ?? 500, alignments ?? 250) : (maxTargetSeqs ?? 500);
  }
  return maxTargetSeqs ?? alignments ?? 500;
}

/** Whether the run searched twice: outfmt 0 has a hit list size of its own (the sizes of the formats differ). */
const searchedTwice = (fixture: Fixture): boolean => hitListSize(fixture.words, 0) !== hitListSize(fixture.words, 6);

function search(abi: ReactorAbi, fixture: Fixture): Outputs {
  const parts = new Map<number, Uint8Array[]>();
  const argv = [fixture.program, '-query', fixture.query, '-subject', fixture.subject, ...fixture.words, '-num_threads', '1'];
  const [queryBytes, subjectBytes] = [readFileSync(join(ENGINE, fixture.query)), readFileSync(join(ENGINE, fixture.subject))];
  const started = abi.start(argv);
  // A run with -parse_deflines reads its inputs together for its argv (line WQA's sweep, D7, 2026-11-19), as the app's runs do.
  const together = started.ok && started.start.registerInputs ? abi.registerInputs(argv, queryBytes, subjectBytes) : undefined;
  const query = together?.query ?? abi.register(fixture.program, ROLE_QUERY, true, queryBytes);
  const subject = together?.subject ?? abi.register(fixture.program, ROLE_SUBJECT, true, subjectBytes);
  try {
    abi.run(argv, query.handle, subject.handle, (stream, bytes) => parts.set(stream, [...(parts.get(stream) ?? []), bytes.slice()]));
  } finally {
    abi.release(query.handle);
    abi.release(subject.handle);
  }
  const streams = new Map<number, Uint8Array>();
  for (const [stream, chunks] of parts) streams.set(stream, Buffer.concat(chunks));
  const text = new TextDecoder().decode(streams.get(1) ?? new Uint8Array());
  const records = text.split('\n').filter((l) => l.trim() !== '').map((l) => JSON.parse(l) as HspRecord);
  return {
    streams,
    records,
    queryIds: query.records.map((r) => r.id),
    subjectIds: subject.records.map((r) => r.id),
    ...(together === undefined ? {} : { shownSubjectIds: together.shownIds.subject }),
  };
}

const sha256 = (bytes: Uint8Array) => createHash('sha256').update(bytes).digest('hex');

/** The one mapping between NCBI's words and LOSAT's (LOSAT/tests/losat_wording.py, the port of
 * LOSAT/src/wording.rs): `map-output` gives LOSAT's bytes of NCBI's output, `hash-form` LOSAT's
 * output read back in NCBI's words (`ncbi-words: ...` when it still writes NCBI's words). */
function wording(command: 'map-output' | 'hash-form', bytes: Uint8Array): Uint8Array | string {
  const run = spawnSync('python3', ['-I', join(ENGINE, 'tests/losat_wording.py'), command, '-'], {
    input: bytes,
    maxBuffer: 1 << 30,
  });
  if (run.error) throw run.error;
  return run.status === 0 ? new Uint8Array(run.stdout) : `ncbi-words: ${run.stderr.toString()}`;
}

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
    const read = wording('hash-form', bytes);
    expect(typeof read === 'string' ? read : sha256(read), label).toBe(frozen);
    return 'bytes';
  }
  const oracle = readFileSync(join(ENGINE, database));
  // run_oracle.py hashes the file with the time makeblastdb ran replaced (without_posted_date).
  const posted = Buffer.from(oracle.toString('latin1').replace(/^( {4}Posted date: {2}).*$/gm, '$1(the time makeblastdb ran)'), 'latin1');
  expect(sha256(posted), `${fixture.id}: ${database}`).toBe(frozen);
  const decoder = new TextDecoder();
  const inLosatWords = wording('map-output', oracle);
  if (typeof inLosatWords === 'string') throw new Error(inLosatWords);
  expect(withoutDatabaseLines(decoder.decode(bytes)), `${label} (database search, -db_gencode exception)`).toBe(
    withoutDatabaseLines(decoder.decode(inLosatWords)),
  );
  return 'database';
}

/** The first and last coordinates of the alignment lines that start with `label`. */
/**
 * The first and last coordinates of the label's lines of one section. With -line_length 1 NCBI writes
 * the first Query line without its start (`Query       T`, the frozen bytes of the t1.*.ll1 rows); the
 * start is then undefined and only the end can be compared. The section of an HSP of length 0
 * (BLASTP `-ungapped -window_size 0`, the ungapped.blastp.len0.* rows) has its score lines and no
 * alignment lines: undefined.
 */
function alignmentEnds(section: string, label: 'Query' | 'Sbjct'): [number | undefined, number] | undefined {
  const lines = section.split('\n').filter((line) => line.startsWith(`${label} `));
  if (lines.length === 0) return undefined;
  const parse = (line: string) => line.trim().split(/\s+/);
  const first = parse(lines[0]!);
  const last = parse(lines.at(-1)!);
  const start = /^\d+$/.test(first[1] ?? '') ? Number(first[1]) : undefined;
  return [start, Number(last.at(-1))];
}

/** The record's coordinates as the section shows them: the start only where the section writes it. */
function shownEnds(ends: [number | undefined, number], start: number, end: number): [number | undefined, number] {
  return [ends[0] === undefined ? undefined : start, end];
}

/** The raw score, the first and last coordinates of one HSP section of outfmt 0, as a comparable text. */
function sectionKey(section: string): string {
  const raw = /^ Score = +[\d.]+ bits \((\d+)\)/.exec(section)?.[1];
  return JSON.stringify([Number(raw), ...(alignmentEnds(section, 'Query') ?? []), ...(alignmentEnds(section, 'Sbjct') ?? [])]);
}

/**
 * The subject headings of outfmt 0 in order with their offsets, the query they belong to (the index
 * among the registered queries, matched by id in order) and the sections of their HSPs. `raw` is the
 * bytes as written (the offsets are the records' coordinates; under `-html` a raw `>` at a line's start
 * is always a heading, NCBI writes a title's `>` as `&gt;`); `shown` gives the text of a range (the ids,
 * the sections' lines, and the query lines: `Query= ` as shown, which under `-html` NCBI writes
 * `<b>Query=</b> `).
 */
function headingsOf(
  raw: string,
  queryIds: readonly string[],
  shown: (text: string) => string,
): Array<{ qIdx: number; offset: number; id: string; hsps: string[] }> {
  const found: Array<{ qIdx: number; offset: number; id: string; hsps: string[] }> = [];
  let qIdx = -1;
  let offset = 0;
  let current: { qIdx: number; offset: number; id: string; hsps: string[] } | undefined;
  let section: string[] | undefined;
  const close = () => {
    if (section !== undefined && current !== undefined) current.hsps.push(sectionKey(shown(section.join('\n'))));
    section = undefined;
  };
  for (const line of raw.split('\n')) {
    const text = shown(line);
    if (text.startsWith('Query= ')) {
      close();
      current = undefined;
      const id = text.slice(7).trim().split(/\s+/)[0]!;
      do qIdx++;
      while (qIdx < queryIds.length && queryIds[qIdx] !== id);
    } else if (line.startsWith('>') && qIdx >= 0) {
      close();
      current = { qIdx, offset, id: shown(line).slice(1).trim().split(/\s+/)[0]!, hsps: [] };
      found.push(current);
    } else if (text.startsWith(' Score =') && current !== undefined) {
      close();
      section = [line];
    } else if (section !== undefined) {
      section.push(line);
    }
    offset += Buffer.byteLength(line) + 1;
  }
  close();
  return found;
}

/** Checks one search and returns what it saw, for the summary. */
function check(fixture: Fixture, outputs: Outputs): { hsps: number; sections: number; unshown: number; unpointed: number } {
  const decoder = new TextDecoder();
  const out6 = outputs.streams.get(6) ?? new Uint8Array();
  const out0 = outputs.streams.get(0) ?? new Uint8Array();
  const records = [...outputs.records].sort((a, b) => a.index - b.index);
  const where = (record: HspRecord) => `${fixture.id} HSP ${record.index}`;
  // A run with -html wrote NCBI's HTML: the ranges are read as the app shows them (`htmlShownText`).
  const html = fixture.words.includes('-html');
  const shownText = (text: string): string => (html ? htmlShownText(text) : text);
  const outfmt0Text = (bytes: Uint8Array): string => shownText(decoder.decode(bytes));
  /** The name of each subject's heading, by subject (-parse_deflines). */
  const headingNames = new Map<number, string>();

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
  // -sorthits and -sorthsps sort the pairwise output only (NCBI: "Not applicable for outfmt > 4";
  // blast_format.cpp:1523-1530), so with them the sections follow that order, not the records'.
  const sortsPairwise = fixture.words.includes('-sorthits') || fixture.words.includes('-sorthsps');
  // BLASTN's matrix-only scoring (`-task rmblastn`, reward 0 and penalty 0) opens its report with
  // RMBLASTN and shows the raw score alone (` Score = 27204`; blast_format.cpp:168-177,
  // showalign.cpp:3586-3591,3605); the other searches show `<bits> bits (<raw>),`.
  const rawScoreOnly = decoder.decode(out0.subarray(0, 9)) === 'RMBLASTN ';
  let previousEnd = 0;
  const headings = new Map<string, readonly [number, number] | null>();
  let sections = 0;
  let unshown = 0;
  // A subject whose title NCBI cannot decode (pd.undec.*, under -parse_deflines) has an empty heading, and
  // NCBI's line stands in place of each of its alignments (c++/src/objtools/align_format/showalign.cpp:1976-1981; LOSAT's
  // report/pairwise.rs write_skipped_alignment; line WQA's sweep, D7, 2026-11-19).
  const skippedLine = /^Sequence with id \S+ no longer exists in database\.\.\.alignment skipped\n$/;
  const skippedHeadings = new Set<number>();
  let skipped = 0;
  for (const record of records) {
    const pair = `${record.q_idx}:${record.s_idx}`;
    expect(record.out0 === null, `${where(record)} section and heading`).toBe(record.out0_subject === null);
    if (headings.has(pair)) expect(record.out0_subject, `${where(record)} heading`).toEqual(headings.get(pair));
    else headings.set(pair, record.out0_subject);
    if (record.out0 === null) {
      unshown++;
      continue;
    }
    const [start, end] = record.out0;
    if (!sortsPairwise) expect(start, `${where(record)} section order`).toBeGreaterThanOrEqual(previousEnd);
    const [hStart, hEnd] = record.out0_subject!;
    expect(hEnd, `${where(record)} heading before section`).toBeLessThanOrEqual(start);
    const section = outfmt0Text(out0.subarray(start, end));
    if (skippedLine.test(section)) {
      expect(hEnd, `${where(record)} the empty heading of a skipped alignment`).toBe(hStart);
      skippedHeadings.add(hStart);
      skipped++;
      previousEnd = end;
      continue;
    }
    sections++;
    const heading = outfmt0Text(out0.subarray(hStart, hEnd));
    expect(heading.startsWith('>'), `${where(record)} heading text`).toBe(true);
    expect(heading, `${where(record)} heading Length=`).toMatch(/\nLength=\d+\n/);
    // The heading is the record's subject (not another subject with the same coordinates); a record without an ID has
    // no name in the heading, or NCBI's `unnamed`.
    const named = heading.slice(1).split('\n')[0]!.trim().split(/\s+/)[0]!;
    const subjectId = outputs.subjectIds[record.s_idx]!;
    if (outputs.shownSubjectIds === undefined) {
      expect(subjectId === '' ? ['', 'unnamed'] : [subjectId], `${where(record)} heading names the record's subject`).toContain(named);
    } else {
      // Under -parse_deflines NCBI's heading names a subject by an ID of its choice (`NP_000006.1` of
      // `gi|6|ref|NP_000006.1|`, `2DEF_B` of `pdb|2DEF|B`): each subject's headings have one name, which
      // no subject of another shown ID has (checked below).
      expect(named, `${where(record)} heading names the record's subject`).toBe(headingNames.get(record.s_idx) ?? named);
      headingNames.set(record.s_idx, named);
    }
    expect(section.startsWith(' Score ='), `${where(record)} section text`).toBe(true);
    if (rawScoreOnly) expect(section, `${where(record)} section raw score`).toMatch(new RegExp(`^ Score = ${record.raw_score}\n`));
    else expect(section, `${where(record)} section raw score`).toContain(`bits (${record.raw_score}),`);
    const queryEnds = alignmentEnds(section, 'Query');
    const subjectEnds = alignmentEnds(section, 'Sbjct');
    if (queryEnds === undefined || subjectEnds === undefined) {
      // An HSP of length 0: no alignment lines, and the record's ranges end one before they start.
      expect([queryEnds, subjectEnds], `${where(record)} a section without alignment lines`).toEqual([undefined, undefined]);
      expect(section, `${where(record)} the section of an HSP of length 0`).toMatch(/^ Identities = 0\/0 \(100%\),/m);
      expect([record.q_end - record.q_start, record.s_end - record.s_start], `${where(record)} length 0`).toEqual([-1, -1]);
    } else {
      expect(queryEnds, `${where(record)} Query lines`).toEqual(shownEnds(queryEnds, record.q_start, record.q_end));
      expect(subjectEnds, `${where(record)} Sbjct lines`).toEqual(shownEnds(subjectEnds, record.s_start, record.s_end));
    }
    previousEnd = end;
  }
  const subjectsNamed = new Map<string, Set<string>>();
  for (const [sIdx, name] of headingNames) subjectsNamed.set(name, (subjectsNamed.get(name) ?? new Set()).add(outputs.shownSubjectIds![sIdx]!));
  for (const [name, shown] of subjectsNamed) expect([...shown], `${fixture.id}: the heading ${name} names one subject`).toHaveLength(1);
  // The text has as many sections and subject headings as the records point to: no section or
  // heading of outfmt 0 lacks its records. In a run of two searches outfmt 0 has a hit list of its own:
  // a section or heading that no record points to is a pair whose records are absent or differ.
  // The headings are counted and located on the raw bytes (the records' offsets); the shown text serves the contents.
  const raw = decoder.decode(out0);
  const text = outfmt0Text(out0);
  const pointed = new Set([...headings.values()].flatMap((range) => (range === null || skippedHeadings.has(range[0]) ? [] : [range[0]])));
  const twice = searchedTwice(fixture);
  let unpointed = 0;
  if (!twice) {
    expect(text.match(/^ Score =/gm)?.length ?? 0, `${fixture.id}: outfmt 0 sections and records with one`).toBe(sections);
    expect(text.match(/^Sequence with id \S+ no longer exists in database\.\.\.alignment skipped$/gm)?.length ?? 0, `${fixture.id}: skipped alignments`).toBe(skipped);
    expect(raw.match(/^>/gm)?.length ?? 0, `${fixture.id}: outfmt 0 subject headings and the records' headings`).toBe(pointed.size);
  } else {
    const recordsOf = new Map<string, string[]>();
    for (const record of records) {
      const key = `${record.q_idx}:${record.s_idx}`;
      recordsOf.set(key, [...(recordsOf.get(key) ?? []), JSON.stringify([record.raw_score, record.q_start, record.q_end, record.s_start, record.s_end])]);
    }
    let seenPointed = 0;
    for (const heading of headingsOf(raw, outputs.queryIds, shownText)) {
      if (pointed.has(heading.offset)) {
        seenPointed++;
        continue;
      }
      unpointed += heading.hsps.length;
      // The ids of the registered subjects that the heading can name; the pair differs for each of them.
      for (const [sIdx, id] of outputs.subjectIds.entries()) {
        if (id !== heading.id) continue;
        const ours = recordsOf.get(`${heading.qIdx}:${sIdx}`) ?? [];
        expect(ours, `${fixture.id}: q${heading.qIdx} ${heading.id}: outfmt 0's HSPs are not the records' (the pair would be linked)`).not.toEqual(heading.hsps);
      }
    }
    expect(seenPointed, `${fixture.id}: every heading that the records point to is in outfmt 0`).toBe(pointed.size);
    expect(text.match(/^ Score =/gm)?.length ?? 0, `${fixture.id}: outfmt 0 sections = linked + unlinked`).toBe(sections + unpointed);
  }

  // Subjects without alignments in outfmt 0 come after those with them, in each query.
  const table = hspTable(outputs.records);
  const index = buildResultIndex(table);
  let grouped = 0;
  for (const query of index.queries.values()) {
    const shown = query.subjects.map((subject) => table.out0SubjectStart[subject.rows[0]!]! >= 0);
    const firstUnshown = shown.indexOf(false);
    if (firstUnshown >= 0 && !twice) expect(shown.slice(firstUnshown).every((s) => !s), `${fixture.id} q${query.qIdx}: unshown subjects last`).toBe(true);
    for (const subject of query.subjects) {
      grouped += subject.rows.length;
      const ranks = subject.rows.map((row) => table.rank[row]!);
      expect(ranks, `${fixture.id} q${query.qIdx} s${subject.sIdx}: rank order`).toEqual([...ranks].sort((a, b) => a - b));
    }
  }
  expect(grouped, `${fixture.id}: every HSP is in the index`).toBe(records.length);
  return { hsps: records.length, sections, unshown, unpointed };
}

const reactors = findReactors();

describe.skipIf(reactors === undefined)('HSP records and the rows and sections of outfmt 6 and 0 (every outfmt 0 fixture)', () => {
  const searches = reactors === undefined ? new Map<string, Fixture[]>() : fixtures();
  const module = reactors === undefined ? undefined : new WebAssembly.Module(reactors.serial.bytes as BufferSource);
  const summary = {
    searches: 0,
    twice: 0,
    hsps: 0,
    sections: 0,
    unshown: 0,
    unpointed: 0,
    parsed: 0,
    refused: 0,
    programs: new Set<string>(),
    frozen: { bytes: 0, database: 0, none: 0 },
  };

  it(
    'agree for every search',
    async () => {
      let abi: ReactorAbi | undefined;
      let runs = 0;
      for (const group of searches.values()) {
        // A fresh instance now and then keeps the instance's memory small.
        if (abi === undefined || runs++ % 20 === 0) abi = (await instantiateSerial(module!)).abi;
        const fixture = group[0]!;
        let outputs: Outputs;
        try {
          outputs = search(abi, fixture);
        } catch (error) {
          // The search of these fixtures is outfmt 0 (the fixture's other format only adds its search):
          // LOSAT refuses outfmt 0 of a record that is both a query and a subject and has masks, as on the
          // command line (the pd.self rows of formats other than 0).
          const shared = /is both a query and a subject under -parse_deflines and has masks/.test(String(error));
          if (!shared || !fixture.words.includes('-parse_deflines')) throw error;
          summary.refused++;
          continue;
        }
        if (outputs.shownSubjectIds !== undefined) summary.parsed++;
        for (const frozen of group) summary.frozen[expectFrozen(frozen, outputs)]++;
        const seen = check(fixture, outputs);
        summary.searches++;
        summary.hsps += seen.hsps;
        summary.sections += seen.sections;
        summary.unshown += seen.unshown;
        summary.unpointed += seen.unpointed;
        if (searchedTwice(fixture)) summary.twice++;
        summary.programs.add(fixture.program);
      }
      console.log(
        `HSP correspondence: ${summary.searches} searches (${[...summary.programs].sort().join(', ')}), ${summary.hsps} HSPs, ` +
          `${summary.sections} outfmt 0 sections, ${summary.unshown} HSPs that outfmt 0 does not show; ${summary.twice} searches of two hit list sizes, ${summary.unpointed} outfmt 0 sections no record points to; fixtures equal to NCBI: ` +
          `${summary.frozen.bytes} byte for byte, ${summary.frozen.database} as the database search (-db_gencode exception), ` +
          `${summary.frozen.none} field lists not compared; ${summary.parsed} searches with -parse_deflines, ${summary.refused} refused as on the command line`,
      );
      expect([...summary.programs].sort()).toEqual(['blastn', 'blastp', 'tblastn', 'tblastx']);
      expect(summary.unshown, 'the fixtures include HSPs that outfmt 0 does not show').toBeGreaterThan(0);
      expect(summary.twice, 'the fixtures include searches of two hit list sizes').toBeGreaterThan(0);
      expect(summary.unpointed, 'the split fixture has an outfmt 0 section that no record points to').toBeGreaterThan(0);
      expect(summary.frozen.bytes, 'fixtures compared byte for byte').toBeGreaterThan(0);
      expect(summary.frozen.database, 'fixtures of the -db_gencode exception compared').toBeGreaterThan(0);
      expect(summary.parsed, 'searches with -parse_deflines').toBeGreaterThan(0);
    },
    1_800_000,
  );
});
