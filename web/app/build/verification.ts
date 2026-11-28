// The verification table of the results screen (plan §6.1), generated at build time from
// docs/web/verification_cells.tsv (the browser runtime) and the certification records of the
// repository (the option sets compared byte for byte with NCBI BLAST+ 2.17.0). Nobody writes
// the table by hand: a record that cannot be read as options, or that is not a byte
// comparison (a rejection, an environment variable, an approved NCBI crash), adds nothing.
// The plugin serves the table to the application as the module `virtual:losat-verification`.
import { createHash } from 'node:crypto';
import { existsSync, readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import type { OutputFormat, ProgramId } from 'losat/domain';
import type { Plugin } from 'vite';
import { parseToolDescription } from '../../packages/losat/src/domain/tool-description.ts';
import {
  optionKey,
  toolOptionGrammar,
  toolOptionKey,
  type BrowserVerification,
  type OptionGrammar,
  type ProgramVerification,
  type ToolOptionGrammar,
  type VerificationTable,
} from '../src/domain/verification.ts';

const VIRTUAL_ID = 'virtual:losat-verification';
const RESOLVED_ID = `\0${VIRTUAL_ID}`;
const NCBI_VERSION = '2.17.0';
const PROGRAMS = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'] as const satisfies readonly ProgramId[];
const DESCRIBE = fileURLToPath(new URL('../src/infra/fake/describe.json', import.meta.url));
const WINDOWMASKER_XMLHELP = fileURLToPath(new URL('../src/infra/fake/windowmasker-xmlhelp.json', import.meta.url));
const WINDOWMASKER_CASES_PATH = 'LOSAT/tests/windowmasker/cases.tsv';
const BLAST_FORMATTER_XMLHELP = fileURLToPath(new URL('../src/infra/fake/blast_formatter-xmlhelp.json', import.meta.url));
const BLAST_FORMATTER_CASES_PATH = 'LOSAT/tests/blast_formatter/cases.tsv';
const MAKEBLASTDB_XMLHELP = fileURLToPath(new URL('../src/infra/fake/makeblastdb-xmlhelp.json', import.meta.url));
const TOOL_SKIP_FLAGS = new Set(['-h', '-help', '-xmlhelp', '-version', '-version-full', '-version-full-xml', '-version-full-json', '-help-full', '-dryrun']);

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

  const wm = casesToolVerification(repository, 'windowmasker', WINDOWMASKER_XMLHELP, WINDOWMASKER_CASES_PATH, 'LOSAT/tests/windowmasker/ncbi');
  sources.push(wm.source);
  const bf = casesToolVerification(repository, 'blast_formatter', BLAST_FORMATTER_XMLHELP, BLAST_FORMATTER_CASES_PATH, 'LOSAT/tests/blast_formatter/ncbi');
  sources.push(bf.source);
  const mk = makeblastdbVerification(repository);
  sources.push(...mk.sources);

  const browsers = browserVerification(readTable(readFileSync(join(repository, 'docs/web/verification_cells.tsv'), 'utf8')));
  const programs: Partial<Dictionary<ProgramVerification>> = {};
  for (const program of PROGRAMS) {
    const optionSets: Dictionary<readonly OutputFormat[]> = {};
    for (const key of [...sets.get(program)!.keys()].sort()) optionSets[key] = [...sets.get(program)!.get(key)!].sort((a, b) => a - b);
    const browser = browsers.get(program);
    programs[program] = browser === undefined ? { optionSets } : { optionSets, browser };
  }
  return {
    ncbi: NCBI_VERSION,
    sources,
    programs,
    tools: {
      windowmasker: { optionSets: wm.optionSets },
      blast_formatter: { optionSets: bf.optionSets },
      makeblastdb: { optionSets: mk.optionSets },
    },
  };
}

/**
 * Generates a CLI tool's option sets from `casesPath` (`cases.tsv`) and its NCBI-frozen
 * `ncbi/<name>.rc` files (matching `masker_fixtures.rs`: last row per name wins,
 * shell-word quoting, skipping non-zero exits, environment variables, stdin, help and `-dryrun`).
 */
function casesToolVerification(
  repository: string,
  tool: 'windowmasker' | 'blast_formatter',
  xmlhelpPath: string,
  casesPath: string,
  ncbiRelDir: string,
): { readonly source: Source; readonly optionSets: readonly string[] } {
  const casesBytes = readFileSync(join(repository, casesPath));
  const cases = new Map<string, string>();
  for (const line of casesBytes.toString('utf8').split(/\r?\n/)) {
    if (line === '' || line.startsWith('#')) continue;
    const tab = line.indexOf('\t');
    if (tab < 0) continue;
    cases.set(line.slice(0, tab), line.slice(tab + 1));
  }

  const ncbiDir = join(repository, ncbiRelDir);
  const commands: (readonly string[])[] = [];
  for (const [name, args] of cases) {
    if (name.startsWith('refused_')) continue;
    const rcPath = join(ncbiDir, `${name}.rc`);
    if (!existsSync(rcPath) || readFileSync(rcPath, 'utf8').trim() !== '0') continue;
    const parsed = parseShellWords(args);
    if (parsed !== undefined && parsed.stdin === undefined) commands.push(parsed.words);
  }
  const sets = toolOptionSets(tool, xmlhelpPath, commands);
  return {
    source: { path: casesPath, sha256: createHash('sha256').update(casesBytes).digest('hex'), optionSets: sets.accepted },
    optionSets: sets.optionSets,
  };
}

/**
 * makeblastdb's option sets (line WW, WW4b2 R1), from every ledger of makeblastdb runs that NCBI BLAST+ 2.17.0 froze and
 * that a LOSAT test compares with LOSAT's run: the exit status, the standard output and error, and every file the run
 * writes, byte for byte (dates, times and run directories masked). Each ledger's row shape becomes the words after
 * `makeblastdb` in its reader of `MAKEBLASTDB_LEDGERS`; a row counts when NCBI ended it with exit status 0, without
 * standard input or an environment variable. The ledgers whose test leaves a written file uncompared add nothing:
 * `tests/seqdb/mask_writer` and `makeblastdb_app_cases.rs` (the metadata file), `tests/seqdb/mask_parse` and
 * `mask_ids` (the masks read back, not the files), `gi_mask_read` (the copies read back).
 */
function makeblastdbVerification(repository: string): { readonly sources: readonly Source[]; readonly optionSets: readonly string[] } {
  const sources: Source[] = [];
  const keys = new Set<string>();
  for (const [path, read] of MAKEBLASTDB_LEDGERS) {
    const bytes = readFileSync(join(repository, path));
    const sets = toolOptionSets('makeblastdb', MAKEBLASTDB_XMLHELP, read(bytes.toString('utf8'), dirname(join(repository, path))));
    for (const key of sets.optionSets) keys.add(key);
    sources.push({ path, sha256: createHash('sha256').update(bytes).digest('hex'), optionSets: sets.accepted });
  }
  return { sources, optionSets: [...keys].sort() };
}

/** A ledger of makeblastdb runs: its rows as the words after `makeblastdb`, read from its text and its directory. */
type MakeblastdbLedger = (text: string, directory: string) => Iterable<readonly string[]>;

// Each ledger of makeblastdb runs and how its rows are read, in the order of the sources of the table.
const MAKEBLASTDB_LEDGERS: readonly (readonly [string, MakeblastdbLedger])[] = [
  ['LOSAT/tests/dropin_corpus_manifest.tsv', dropinMakeblastdbRuns],
  ['LOSAT/tests/fixtures/asn1/makeblastdb/ncbi_2_17_0_manifest.tsv', asn1MakeblastdbRuns],
  ['LOSAT/tests/seqdb/gi_mask/cases.tsv', giMaskMakeblastdbRuns],
  ['LOSAT/tests/non_utf8_names/cases.tsv', nonUtf8MakeblastdbRuns],
  ['LOSAT/src/writedb/tests/makeblastdb_v5_cases.rs', v5MakeblastdbRuns],
  ['LOSAT/src/writedb/tests/makeblastdb_parse_cases.rs', parseMakeblastdbRuns],
  ['LOSAT/src/writedb/tests/makeblastdb_fix_cases.rs', fixMakeblastdbRuns],
];

/**
 * The drop-in ledger (`dropin_corpus.py check`): the `makeblastdb ...` rows checked as `bytes` or `bytes-dir` (stdout
 * and every file the run writes) where NCBI wrote nothing on stderr (otherwise `check` only wants some text there),
 * with no standard input.
 */
function* dropinMakeblastdbRuns(text: string): Iterable<readonly string[]> {
  for (const row of readTable(text)) {
    const command = row['command'] ?? '';
    if (!command.startsWith('makeblastdb ') || (row['check'] !== 'bytes' && row['check'] !== 'bytes-dir')) continue;
    if (row['rc'] !== '0' || row['stderr_empty'] !== 'yes' || (row['stdin'] ?? '') !== '') continue;
    const parsed = parseShellWords(command.slice('makeblastdb '.length));
    if (parsed !== undefined && parsed.stdin === undefined) yield parsed.words;
  }
}

/**
 * The ASN.1 inputs of makeblastdb (`makeblastdb_asn1_tests.rs`): `-in <in_spec> -input_type <input_type> -dbtype
 * <dbtype> -title T -out <dir>/db` and the words of `extra_args`; `env_spec` names an environment variable and
 * `stdin_file` the standard input.
 */
function* asn1MakeblastdbRuns(text: string): Iterable<readonly string[]> {
  for (const row of readTable(text)) {
    if (row['rc'] !== '0' || row['env_spec'] !== '' || row['stdin_file'] !== '') continue;
    const fixed = ['-in', row['in_spec'] ?? '', '-input_type', row['input_type'] ?? '', '-dbtype', row['dbtype'] ?? ''];
    yield [...fixed, '-title', 'T', '-out', '@D/db', ...split(row['extra_args'] ?? '')];
  }
}

/**
 * `-gi_mask` (`tests/makeblastdb_gi_mask.rs`): name, `losat`, setup, then the words, one per field; NCBI's exit status
 * in `expected/<name>.rc`. Only the rows of `losat` = `same`. A single-volume database built with `-gi_mask` is
 * NCBI's run with the change of PD-LOSAT-NCBI-DEFECTS approved exception 25 (NCBI ended it with exit status 2,
 * `make_gi_mask_fixtures.py`), so a `-gi_mask` row counts only when its files hold numbered volumes.
 */
function* giMaskMakeblastdbRuns(text: string, directory: string): Iterable<readonly string[]> {
  for (const line of text.split('\n')) {
    const [name, losat, , ...words] = line.split('\t');
    if (name === undefined || name === '' || losat !== 'same') continue;
    if (readFileSync(join(directory, 'expected', `${name}.rc`), 'utf8').trim() !== '0') continue;
    const files = readFileSync(join(directory, 'expected', `${name}.files`), 'utf8');
    if (words.includes('-gi_mask') && !/\.\d{2,}\.[np][a-z]{2}\t/.test(files)) continue;
    yield words;
  }
}

/**
 * Names that are not UTF-8 (`tests/non_utf8_names.rs`): case, kind, env, setup, then the command line split at single
 * spaces with `\xHH`, `\n`, `\s` and `\\` in its words; NCBI's exit status in `expected/<case>.rc`.
 */
function* nonUtf8MakeblastdbRuns(text: string, directory: string): Iterable<readonly string[]> {
  for (const line of text.split('\n')) {
    if (line === '' || line.startsWith('#')) continue;
    const [name, kind, env, , argv] = line.split('\t');
    if (kind !== 'ncbi' || env !== '-' || argv === undefined || !argv.startsWith('makeblastdb ')) continue;
    if (readFileSync(join(directory, 'expected', `${name ?? ''}.rc`), 'utf8').trim() !== '0') continue;
    yield argv.split(' ').slice(1).map(unescapeManifestWord);
  }
}

/** A word of `non_utf8_names/cases.tsv`: `\xHH` as the code point HH, `\n`, `\s` (a space) and `\\`. */
function unescapeManifestWord(word: string): string {
  return word.replace(/\\(x[0-9a-fA-F]{2}|n|s|\\)/g, (_, escape: string) =>
    escape === 'n' ? '\n' : escape === 's' ? ' ' : escape === '\\' ? '\\' : String.fromCharCode(Number.parseInt(escape.slice(1), 16)),
  );
}

/**
 * The version 5 runs (`run_v5_case` of `LOSAT/src/writedb/tests.rs`): each case runs its words in one directory `{D}`,
 * then the files directly in `{D}` are compared. A run counts when its `-out` is a file of `{D}` that no later run of
 * the case writes again (the runs that build a source database in a folder of `{D}` are not compared).
 */
function* v5MakeblastdbRuns(text: string): Iterable<readonly string[]> {
  for (const block of rustBlocks(text, 'V5Case')) {
    const runs = rustRuns(block);
    for (const [at, run] of runs.entries()) {
      const out = valueOf(run.words, '-out');
      if (run.status !== 0 || out === undefined || !/^\{D\}\/[^/]+$/.test(out)) continue;
      if (runs.slice(at + 1).some((later) => valueOf(later.words, '-out') === out)) continue;
      yield run.words;
    }
  }
}

/** The `-parse_seqids` runs (`run_parse_case`): `-in {D}/in.fa` and the case's words. */
function* parseMakeblastdbRuns(text: string): Iterable<readonly string[]> {
  for (const block of rustBlocks(text, 'ParseCase')) {
    for (const run of rustRuns(block)) if (run.status === 0) yield ['-in', '{D}/in.fa', ...run.words];
  }
}

/** The review-fix runs (`run_fix_case`): every word after `makeblastdb`, and a standard input. */
function* fixMakeblastdbRuns(text: string): Iterable<readonly string[]> {
  for (const block of rustBlocks(text, 'FixCase')) {
    if (!/\n\s*stdin: b"",\n/.test(block)) continue;
    for (const run of rustRuns(block)) if (run.status === 0) yield run.words;
  }
}

/** The cases of a generated Rust table (`    <Struct> {` at the start of a line), each as its text. */
function rustBlocks(text: string, struct: string): string[] {
  return text.split(`\n    ${struct} {\n`).slice(1);
}

/** The runs of a case of a generated Rust table: each `words: &[...]` (string literals) and the `status` after it. */
function rustRuns(block: string): { readonly words: readonly string[]; readonly status: number }[] {
  const runs: { words: string[]; status: number }[] = [];
  for (const match of block.matchAll(/words: &\[((?:[^\]"]|"(?:[^"\\]|\\.)*")*)\],\s*status: (-?\d+),/g)) {
    runs.push({ words: [...match[1]!.matchAll(/"((?:[^"\\]|\\.)*)"/g)].map(([, literal]) => rustString(literal!)), status: Number(match[2]) });
  }
  return runs;
}

/** A Rust string literal's text: `\"`, `\\`, `\n`, `\r`, `\t`, `\0`, `\xHH` and `\u{H...}`. */
function rustString(literal: string): string {
  return literal.replace(/\\(x[0-9a-fA-F]{2}|u\{[0-9a-fA-F]+\}|.)/g, (_, escape: string) => {
    if (escape.startsWith('x')) return String.fromCharCode(Number.parseInt(escape.slice(1), 16));
    if (escape.startsWith('u{')) return String.fromCodePoint(Number.parseInt(escape.slice(2, -1), 16));
    return ({ n: '\n', r: '\r', t: '\t', '0': '\0' } as Dictionary<string>)[escape] ?? escape;
  });
}

/** The value of the last `flag` of words, or undefined. */
function valueOf(words: readonly string[], flag: string): string | undefined {
  const at = words.lastIndexOf(flag);
  return at < 0 ? undefined : words[at + 1];
}

/**
 * The option sets of a tool's command lines (the words after the tool's name): skipping environment variables, help,
 * `-dryrun` and flags that the tool's `-xmlhelp` does not list.
 */
function toolOptionSets(
  tool: string,
  xmlhelpPath: string,
  commands: Iterable<readonly string[]>,
): { readonly accepted: number; readonly optionSets: readonly string[] } {
  const xmlDoc = JSON.parse(readFileSync(xmlhelpPath, 'utf8')) as { xml: string };
  const description = parseToolDescription(xmlDoc.xml);
  const baseGrammar = toolOptionGrammar(description);
  const knownFlags = new Set(['-logfile', '-conffile', ...description.parameters.map((param) => param.flag)]);
  const keys = new Set<string>();
  let accepted = 0;
  for (const words of commands) {
    if (words.length === 0 || hasEnvAssignment(words[0]!)) continue;
    if (words.some((word) => TOOL_SKIP_FLAGS.has(word.split('=', 1)[0]!))) continue;

    let unknown = false;
    const grammar: ToolOptionGrammar = {
      takesValue: (flag) => {
        if (!knownFlags.has(flag)) unknown = true;
        return baseGrammar.takesValue(flag);
      },
      ...(baseGrammar.isBoolean === undefined ? {} : { isBoolean: baseGrammar.isBoolean }),
      ...(baseGrammar.defaultValue === undefined ? {} : { defaultValue: baseGrammar.defaultValue }),
    };
    const key = toolOptionKey([tool, ...words], grammar);
    if (unknown || key === undefined) continue;
    keys.add(key);
    accepted++;
  }
  return { accepted, optionSets: [...keys].sort() };
}

/** Splits `cases.tsv` shell words with single/double quotes and `< file` stdin (matching `masker_fixtures.rs`). */
function parseShellWords(line: string): { readonly words: readonly string[]; readonly stdin?: string } | undefined {
  const words: string[] = [];
  let word: string | undefined;
  let quote: '"' | "'" | undefined;
  for (const ch of line) {
    if (quote === undefined && (ch === '"' || ch === "'")) {
      quote = ch;
      word ??= '';
    } else if (quote !== undefined && ch === quote) {
      quote = undefined;
    } else if (quote === undefined && (ch === ' ' || ch === '\t')) {
      if (word !== undefined) {
        words.push(word);
        word = undefined;
      }
    } else {
      word = (word ?? '') + ch;
    }
  }
  if (quote !== undefined) return undefined;
  if (word !== undefined) words.push(word);
  const stdinAt = words.indexOf('<');
  if (stdinAt >= 0) {
    const stdin = words[stdinAt + 1];
    if (stdin === undefined) return undefined;
    words.splice(stdinAt, 2);
    return { words, stdin };
  }
  return { words };
}

function hasEnvAssignment(word: string): boolean {
  const eq = word.indexOf('=');
  if (eq <= 0) return false;
  return /^[A-Z_]+$/.test(word.slice(0, eq));
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
