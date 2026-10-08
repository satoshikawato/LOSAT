// The verification table of the results screen (plan §6.1), generated at build time from
// docs/web/verification_cells.tsv (the browser runtime) and the certification records of the
// repository (the option sets compared byte for byte with NCBI BLAST+ 2.17.0). Nobody writes
// the table by hand: a record that cannot be read as options, or that is not a byte
// comparison (a rejection, an environment variable, an approved NCBI crash), adds nothing.
// The plugin serves the table to the application as the module `virtual:losat-verification`.
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import type { Plugin } from 'vite';
import { optionKey, type BrowserVerification, type OptionGrammar, type ProgramVerification, type VerificationTable } from '../src/domain/verification.ts';
import type { OutputFormat } from '../src/domain/output-format.ts';
import type { ProgramId } from '../src/domain/programs.ts';

const VIRTUAL_ID = 'virtual:losat-verification';
const RESOLVED_ID = `\0${VIRTUAL_ID}`;
const NCBI_VERSION = '2.17.0';
const PROGRAMS = ['blastn', 'blastp', 'tblastn', 'tblastx'] as const satisfies readonly ProgramId[];
const DESCRIBE = fileURLToPath(new URL('../src/infra/fake/describe.json', import.meta.url));

type Program = (typeof PROGRAMS)[number];
type Source = VerificationTable['sources'][number];

/** The options of a program in the engine's `describe`. */
interface Described {
  readonly takesValue: ReadonlyMap<string, boolean>;
  readonly defaultTask: string | undefined;
}

/** One comparison of a certification source: the program, its option words and the format compared. */
interface Comparison {
  readonly program: string;
  readonly words: readonly string[];
  readonly format: string;
}

/** The table, from the repository root (the directory with docs/ and LOSAT/). */
export function generateVerificationTable(repository: string): VerificationTable {
  const describe = JSON.parse(readFileSync(DESCRIBE, 'utf8')) as Dictionary<{
    parameters: { flag: string; takes_value: boolean; default?: string }[];
  }>;
  const described = new Map<Program, Described>(
    PROGRAMS.map((program) => {
      const parameters = describe[program]?.parameters ?? [];
      return [
        program,
        {
          takesValue: new Map(parameters.map((p) => [p.flag, p.takes_value])),
          defaultTask: parameters.find((p) => p.flag === '-task')?.default,
        },
      ];
    }),
  );
  const sets = new Map<Program, Map<string, Set<OutputFormat>>>(PROGRAMS.map((program) => [program, new Map()]));
  const sources: Source[] = [];

  for (const [path, read] of SOURCES) {
    const bytes = readFileSync(join(repository, path));
    let accepted = 0;
    for (const record of read(bytes.toString('utf8'))) {
      if (!isProgram(record.program)) continue;
      const format = Number(record.format);
      if (format !== 0 && format !== 6 && format !== 7) continue;
      const key = keyOf(record.program, record.words, described);
      if (key === undefined) continue;
      const formats = sets.get(record.program)!;
      formats.set(key, (formats.get(key) ?? new Set()).add(format));
      accepted++;
    }
    sources.push({ path, sha256: createHash('sha256').update(bytes).digest('hex'), optionSets: accepted });
  }

  const browsers = browserVerification(readTable(readFileSync(join(repository, 'docs/web/verification_cells.tsv'), 'utf8')));
  const programs: Partial<Dictionary<ProgramVerification>> = {};
  for (const program of PROGRAMS) {
    const optionSets: Dictionary<readonly OutputFormat[]> = {};
    for (const key of [...sets.get(program)!.keys()].sort()) optionSets[key] = [...sets.get(program)!.get(key)!].sort((a, b) => a - b);
    const browser = browsers.get(program);
    programs[program] = browser === undefined ? { optionSets } : { optionSets, browser };
  }
  return { ncbi: NCBI_VERSION, sources, programs };
}

/** The key of an option set, or undefined when a flag is not one of the program's described options. */
function keyOf(program: Program, words: readonly string[], described: ReadonlyMap<Program, Described>): string | undefined {
  const { takesValue, defaultTask } = described.get(program)!;
  let unknown = false;
  const grammar: OptionGrammar = {
    takesValue(flag) {
      const value = takesValue.get(flag);
      if (value === undefined) unknown = true;
      return value ?? false;
    },
    ...(defaultTask === undefined ? {} : { defaultTask }),
  };
  const key = optionKey([program, ...words], grammar);
  return unknown ? undefined : key;
}

function isProgram(program: string): program is Program {
  return (PROGRAMS as readonly string[]).includes(program);
}

// Each source file and how its records are read, in the order of the sources of the table.
const SOURCES: readonly (readonly [string, (text: string) => Iterable<Comparison>])[] = [
  ['LOSAT/tests/outfmt0_manifest.tsv', outfmt0Records],
  ['LOSAT/tests/fixtures/blastn_regression/manifest.tsv', (text) => regressionRecords(text, 'blastn')],
  ['LOSAT/tests/fixtures/tblastx_regression/manifest.tsv', (text) => regressionRecords(text, 'tblastx')],
  ['LOSAT/tests/fixtures/range_regression/manifest.tsv', (text) => regressionRecords(text, undefined)],
  ['LOSAT/tests/blastp_v010_parity_manifest.tsv', blastpGateRecords],
  ['LOSAT/tests/tblastx_v010_parity_manifest.tsv', tblastxGateRecords],
  ['docs/evidence/tlosan_stage_g/matrix_162.jsonl', tblastnMatrixRecords],
  ...(['blastp', 'tblastn', 'tblastx'] as const).map(
    (program) => [`docs/evidence/losat_web_e2e/run-20261004T163746Z/sweeps/after-${program}.tsv`, (text: string) => sweepRecords(text, program)] as const,
  ),
];

/** outfmt0_manifest.tsv: NCBI's frozen reports, and the approved -db_gencode exception against NCBI's -db oracle. */
function* outfmt0Records(text: string): Iterable<Comparison> {
  for (const row of readTable(text)) {
    if (row.contract !== '' && row.contract !== 'approved_db_gencode_deviation') continue;
    const extra = row.extra_args ?? '';
    if (hasQuote(extra)) continue;
    const task = row.task ?? '';
    yield { program: row.program ?? '', words: [...(task === '' ? [] : ['-task', task]), ...split(extra)], format: row.outfmt === '' ? '0' : (row.outfmt ?? '') };
  }
}

/**
 * The frozen regression searches. The program is the one the fixture script runs: fixed for
 * the BLASTN and TBLASTX manifests, the first word of the argv for the range manifest.
 */
function* regressionRecords(text: string, fixed: Program | undefined): Iterable<Comparison> {
  for (const row of readTable(text)) {
    const argv = row.argv ?? '';
    if (hasQuote(argv) || row.env !== '' || row.losat_extra !== '' || row.exit !== '0') continue;
    const words = split(argv);
    const program = fixed ?? words.shift() ?? '';
    const at = words.indexOf('-outfmt');
    const format = words[at + 1];
    if (at < 0 || format === undefined) continue;
    yield { program, words, format };
  }
}

/** Gate A, BLASTP: the certified frozen bytes of outfmt 6. */
function* blastpGateRecords(text: string): Iterable<Comparison> {
  for (const row of readTable(text)) {
    const words: string[] = [];
    if (row.max_hsps_per_subject !== '') words.push('-max_hsps', row.max_hsps_per_subject ?? '');
    if (row.max_target_seqs !== '') words.push('-max_target_seqs', row.max_target_seqs ?? '');
    yield { program: 'blastp', words, format: row.outfmt ?? '' };
  }
}

/** Gate A, TBLASTX: the certified frozen bytes of outfmt 6, the genetic codes written when explicit. */
function* tblastxGateRecords(text: string): Iterable<Comparison> {
  for (const row of readTable(text)) {
    if (row.contract !== 'parity') continue;
    const words = row.gencode_args === 'explicit' ? ['-query_gencode', row.query_gencode ?? '', '-db_gencode', row.db_gencode ?? ''] : [];
    yield { program: 'tblastx', words, format: row.outfmt ?? '' };
  }
}

/** Stage G, TBLASTN: the matrix of 27 genetic codes, equal to the frozen bytes. */
function* tblastnMatrixRecords(text: string): Iterable<Comparison> {
  for (const line of text.split('\n')) {
    if (line.trim() === '') continue;
    const row = JSON.parse(line) as { equal?: boolean; exit?: number; losat_command?: string[] };
    const [, program, ...words] = row.losat_command ?? [];
    if (row.equal !== true || row.exit !== 0 || program === undefined) continue;
    const at = words.indexOf('-outfmt');
    const format = words[at + 1];
    if (at < 0 || format === undefined) continue;
    yield { program, words, format };
  }
}

/** The sweeps of E2e: the options whose outputs were the same as NCBI's, or the approved TBLASTN genetic-code exception. */
function* sweepRecords(text: string, program: Program): Iterable<Comparison> {
  for (const row of readTable(text)) {
    if (row.result !== 'same' && row.result !== 'exception-gencode') continue;
    const options = row.options ?? '';
    if (hasQuote(options)) continue;
    yield { program, words: split(options), format: row.outfmt ?? '' };
  }
}

/** The cells of the browser runtime that were checked, merged per program. */
function browserVerification(rows: readonly Row[]): Map<Program, BrowserVerification> {
  const merged = new Map<Program, { formats: Set<OutputFormat>; paths: Set<'serial' | 'threaded'>; threads: Set<number>; browsers: Set<string> }>();
  for (const row of rows) {
    const runtime = row.runtime_path ?? '';
    const outfmt = row.outfmt ?? '';
    if (!runtime.startsWith('browser') || row.status !== 'checked' || outfmt.includes('custom')) continue;
    const text = `${row.profile ?? ''} ${runtime}`;
    for (const program of (row.program ?? '').split('/').map((p) => p.trim())) {
      if (!isProgram(program)) continue;
      const cell = merged.get(program) ?? { formats: new Set(), paths: new Set(), threads: new Set(), browsers: new Set() };
      for (const [format] of outfmt.matchAll(/\b[067]\b/g)) cell.formats.add(Number(format) as OutputFormat);
      if (runtime.includes('serial')) cell.paths.add('serial');
      if (runtime.includes('threaded')) cell.paths.add('threaded');
      for (const [, count] of text.matchAll(/threads ([0-9]+(?:\/[0-9]+)*)/g)) {
        for (const n of count!.split('/')) cell.threads.add(Number(n));
      }
      for (const [name] of text.matchAll(/Chromium|Firefox|WebKit/g)) cell.browsers.add(name);
      merged.set(program, cell);
    }
  }
  return new Map(
    [...merged].map(([program, cell]) => [
      program,
      {
        formats: [...cell.formats].sort((a, b) => a - b),
        paths: [...cell.paths].sort(),
        threads: [...cell.threads].sort((a, b) => a - b),
        browsers: [...cell.browsers],
      },
    ]),
  );
}

type Dictionary<T> = { [key: string]: T };
type Row = Dictionary<string>;

/** A TSV with '#' comment lines and a header row, as rows by column name. */
function readTable(text: string): Row[] {
  const lines = text.split(/\r?\n/).filter((line) => line !== '' && !line.startsWith('#'));
  const header = (lines.shift() ?? '').split('\t');
  return lines.map((line) => {
    const cells = line.split('\t');
    return Object.fromEntries(header.map((name, i) => [name, cells[i] ?? '']));
  });
}

function split(text: string): string[] {
  return text.split(/\s+/).filter((word) => word !== '');
}

/** Words with a quote cannot be split on white space; those records are left out. */
function hasQuote(text: string): boolean {
  return /['"]/.test(text);
}

/** The plugin: the generated table as the module `virtual:losat-verification`. */
export function losatVerification(): Plugin {
  const repository = fileURLToPath(new URL('../../..', import.meta.url));
  return {
    name: 'losat-verification',
    resolveId(id) {
      return id === VIRTUAL_ID ? RESOLVED_ID : undefined;
    },
    load(id) {
      if (id !== RESOLVED_ID) return undefined;
      return `export const VERIFICATION_TABLE = Object.freeze(${JSON.stringify(generateVerificationTable(repository))});\n`;
    },
  };
}
