// Single command-line builders for the tools (windowmasker: WW3b DESIGN §2.5, brief §4.4;
// blast_formatter: RG-F3p7; makeblastdb: line WW WW4b2): the same argv array is validated by `validateTool`, run by `runTool`,
// displayed beside "Run LOSAT" as `losat windowmasker ...` and shown in Run details. Defaults are
// omitted; Make counts and Convert counts write `-out wmasker.<sformat>` until the user edits it.
import type { ToolDescription, ToolParameterDescription } from 'losat/domain';
import { splitOptions } from './argv';
import { writtenValue, type FieldValue, type FormValues } from './parameters';
import {
  allStageParameters,
  BLAST_FORMATTER_BLOCK_FLAGS,
  blastFormatterOptionParameters,
  defaultCountsOutName,
  makeblastdbOptionParameters,
  makeblastdbParameters,
  stageOptionParameters,
  toToolParameterField,
  WINDOWMASKER_BLOCK_FLAGS,
  type ToolId,
  type ToolStage,
  type WindowMaskerStage,
} from './tools';

export interface WindowMaskerArgvInput {
  readonly stage: WindowMaskerStage;
  /** Value for `-in` (FASTA file name, combined_input.fa, counts file for Convert, or blastdb name). */
  readonly inputName?: string;
  /** Value for `-ustat` (Mask stage only). */
  readonly ustatName?: string;
  /**
   * Value for `-out`: undefined uses the stage default (`wmasker.<sformat>` for `mk_counts`/`convert`,
   * `-` for `mask`). Typing `'-'` sends output to standard output (the CLI default, omitted).
   */
  readonly outName?: string;
  /** Option values keyed by `-flag` (block selects `-infmt`, `-sformat`, `-outfmt`, plus WindowMasker options). */
  readonly values: FormValues;
  /** Parsed `-xmlhelp` description (provides defaults and usage order). */
  readonly description?: ToolDescription;
  /** Extra words parsed from "Additional options" (`splitOptions`). */
  readonly extra?: readonly string[];
}

function toolWrittenValue(param: ToolParameterDescription, value: FieldValue | undefined): string | boolean | undefined {
  const field = toToolParameterField(param);
  const known = param.defaultValue === undefined ? undefined : { defaultValue: param.defaultValue };
  return writtenValue(field, value, known);
}

/**
 * Appends the words of the options set in the form: in `ordered`'s order with the defaults left out, or,
 * before `-xmlhelp` arrives (`ordered` empty), every value set except those of `skip`, so the argv stays faithful.
 */
function pushOptionWords(
  argv: string[],
  ordered: readonly ToolParameterDescription[],
  values: FormValues,
  skip: ReadonlySet<string>,
): void {
  if (ordered.length > 0) {
    for (const param of ordered) {
      const written = toolWrittenValue(param, values[param.flag]);
      if (written === true) argv.push(param.flag);
      else if (typeof written === 'string') argv.push(param.flag, written);
    }
    return;
  }
  for (const [flag, val] of Object.entries(values)) {
    if (skip.has(flag)) continue;
    if (val === true) argv.push(flag);
    else if (typeof val === 'string' && val.trim() !== '') argv.push(flag, val.trim());
  }
}

/** Appends the words of Additional options; throws for `''` and a word with NUL, which the search form refuses too (DW-42, `domain/argv.ts`). */
function pushExtraWords(argv: string[], extra: readonly string[] | undefined): void {
  for (const word of extra ?? []) {
    const unread = word === '' ? "''" : word.includes('\0') ? 'a NUL character' : undefined;
    if (unread !== undefined) throw new Error(`Additional options: ${unread} is not accepted`);
    argv.push(word);
  }
}

/** The count of the options in `params` that differ from their defaults, plus one for Additional options. */
function changedOptionsCount(params: readonly ToolParameterDescription[], values: FormValues, extra: string): number {
  let count = 0;
  for (const param of params) {
    if (toolWrittenValue(param, values[param.flag]) !== undefined) count++;
  }
  if (extra.trim() !== '') count++;
  return count;
}

/**
 * Builds `['windowmasker', stageFlag, ...]` for validation, execution and display.
 * Throws for a word of `extra` that Additional options refuse (`pushExtraWords`).
 */
export function buildWindowMaskerArgv(input: WindowMaskerArgvInput): readonly string[] {
  const argv: string[] = ['windowmasker'];
  const { stage, values, description } = input;

  if (stage === 'mk_counts') {
    argv.push('-mk_counts');
  } else if (stage === 'convert') {
    argv.push('-convert');
  } else {
    const ustat = input.ustatName?.trim() ?? '';
    if (ustat !== '') argv.push('-ustat', ustat);
  }

  const inputName = input.inputName?.trim() ?? '';
  if (inputName !== '' && (stage !== 'mask' || inputName !== '-')) {
    argv.push('-in', inputName);
  }

  const sformat = typeof values['-sformat'] === 'string' ? values['-sformat'].trim() : '';
  const effectiveOut =
    input.outName !== undefined
      ? input.outName.trim()
      : stage === 'mk_counts' || stage === 'convert'
        ? defaultCountsOutName(sformat)
        : '';
  if (effectiveOut !== '' && effectiveOut !== '-') {
    argv.push('-out', effectiveOut);
  }

  // Remaining block selects (-infmt, -sformat, -outfmt) and WindowMasker options in usage order.
  const ordered = allStageParameters(description, stage).filter((param) => !WINDOWMASKER_NAMED_FLAGS.has(param.flag));
  pushOptionWords(argv, ordered, values, WINDOWMASKER_NAMED_FLAGS);
  pushExtraWords(argv, input.extra);
  return Object.freeze(argv);
}

/** windowmasker's options whose words the builder writes from the blocks (the file names). */
const WINDOWMASKER_NAMED_FLAGS: ReadonlySet<string> = new Set(['-in', '-out', '-ustat']);

export interface BlastFormatterArgvInput {
  /** Archive file name in the run directory (`-archive`), omitted when empty. */
  readonly archiveName?: string;
  /** Output format string (`-outfmt`), omitted when empty or NCBI's default `'0'`. */
  readonly outfmt?: string;
  /** Output file name (`-out`), omitted when empty or `'-'` (standard output). */
  readonly outName?: string;
  /** Formatting option values keyed by `-flag` (`-max_target_seqs`, `-parse_deflines`, etc.). */
  readonly values: FormValues;
  /** Parsed `-xmlhelp` description (provides defaults and usage order). */
  readonly description?: ToolDescription;
  /** Extra words parsed from "Additional options" (`splitOptions`). */
  readonly extra?: readonly string[];
}

/**
 * Builds `['blast_formatter', ...]` for validation, execution and display.
 * Throws for a word of `extra` that Additional options refuse (`pushExtraWords`).
 */
export function buildBlastFormatterArgv(input: BlastFormatterArgvInput): readonly string[] {
  const argv: string[] = ['blast_formatter'];
  const archiveName = input.archiveName?.trim() ?? '';
  if (archiveName !== '') argv.push('-archive', archiveName);

  const outfmt = input.outfmt?.trim() ?? '';
  if (outfmt !== '' && outfmt !== '0') argv.push('-outfmt', outfmt);

  const outName = input.outName?.trim() ?? '';
  if (outName !== '' && outName !== '-') argv.push('-out', outName);

  pushOptionWords(argv, blastFormatterOptionParameters(input.description), input.values, BLAST_FORMATTER_BLOCK_FLAGS);
  pushExtraWords(argv, input.extra);
  return Object.freeze(argv);
}

/** Count of changed options inside the collapsible "Formatting options" section (plus Additional options). */
export function changedBlastFormatterOptionsCount(
  values: FormValues,
  extra: string,
  description: ToolDescription | undefined,
): number {
  return changedOptionsCount(blastFormatterOptionParameters(description), values, extra);
}

export interface MakeblastdbArgvInput {
  /**
   * `-in`'s value as the command line takes it, omitted when empty: the files of the run as one list
   * (`seqdbNameList`), or a database's name or list as the user typed it.
   */
  readonly input: string;
  /** Option values keyed by `-flag` (`-input_type`, `-dbtype`, `-out`, and the options of the groups). */
  readonly values: FormValues;
  /** Parsed `-xmlhelp` description (provides defaults and usage order). */
  readonly description?: ToolDescription;
  /** Extra words parsed from "Additional options" (`splitOptions`). */
  readonly extra?: readonly string[];
}

/**
 * Builds `['makeblastdb', '-in', ..., ...]` for validation, execution and display: `-in` first, then the
 * options set in the form in NCBI's `-help` usage order, then Additional options.
 * Throws for a word of `extra` that Additional options refuse (`pushExtraWords`).
 */
export function buildMakeblastdbArgv(input: MakeblastdbArgvInput): readonly string[] {
  const argv: string[] = ['makeblastdb'];
  if (input.input.trim() !== '') argv.push('-in', input.input);
  const ordered = makeblastdbParameters(input.description).filter((param) => param.flag !== '-in');
  pushOptionWords(argv, ordered, input.values, MAKEBLASTDB_NAMED_FLAGS);
  pushExtraWords(argv, input.extra);
  return Object.freeze(argv);
}

/** makeblastdb's option whose words the builder writes from the Input block (the input's names). */
const MAKEBLASTDB_NAMED_FLAGS: ReadonlySet<string> = new Set(['-in']);

/** Count of changed options inside the collapsible "makeblastdb options" section (plus Additional options). */
export function changedMakeblastdbOptionsCount(values: FormValues, extra: string, description: ToolDescription | undefined): number {
  return changedOptionsCount(makeblastdbOptionParameters(description), values, extra);
}

/** The names of a comma-separated list of files (`-mask_data a.asnb,b.asnb`), without empty entries. */
export function fileListNames(value: unknown): readonly string[] {
  if (typeof value !== 'string') return [];
  return value
    .split(',')
    .map((name) => name.trim())
    .filter((name) => name !== '');
}

/** The command words as shown beside Run LOSAT and in Run details (`losat windowmasker ...`). */
export function toolCommandWords(argv: readonly string[]): readonly string[] {
  return ['losat', ...argv];
}

/** Count of changed options inside the collapsible "WindowMasker options" section (plus Additional options). */
export function changedToolOptionsCount(
  stage: WindowMaskerStage,
  values: FormValues,
  extra: string,
  description: ToolDescription | undefined,
): number {
  const params = stageOptionParameters(description, stage).filter((param) => !WINDOWMASKER_BLOCK_FLAGS.has(param.flag));
  return changedOptionsCount(params, values, extra);
}

/** Normalizes simple relative prefixes (`./foo` -> `foo`) so `-out ./x -logfile x` matches `x`. */
function normalizeRelName(name: string): string {
  let out = name.trim();
  while (out.startsWith('./')) out = out.slice(2);
  return out;
}

/** Reads the last value of `flag` in `argv` (supporting both `-flag val` and `-flag=val`). */
export function toolOptionValue(argv: readonly string[], flag: string): string | undefined {
  let found: string | undefined;
  for (let i = 1; i < argv.length; i++) {
    const word = argv[i]!;
    if (word === flag && i + 1 < argv.length) {
      found = argv[++i];
    } else if (word.startsWith(`${flag}=`)) {
      found = word.slice(flag.length + 1);
    }
  }
  return found;
}

function namedOutPath(argv: readonly string[]): string | undefined {
  const out = toolOptionValue(argv, '-out');
  return out !== undefined && out !== '-' ? normalizeRelName(out) : undefined;
}

/** Identifies which option(s) (`-out`, `-logfile`) named an output file produced by a tool run. */
export function toolFileOptionLabel(argv: readonly string[], fileName: string): string | undefined {
  const target = normalizeRelName(fileName);
  const log = toolOptionValue(argv, '-logfile');
  const matchesOut = namedOutPath(argv) === target;
  const matchesLog = log !== undefined && normalizeRelName(log) === target;
  if (matchesOut && matchesLog) return '-out, -logfile';
  if (matchesOut) return '-out';
  if (matchesLog) return '-logfile';
  return undefined;
}

/** Whether `fileName` is the counts output file of WindowMasker `-mk_counts` or `-convert` (eligible for "Use in BLASTN" / "Mask with these counts"). */
export function isCountsOutputFile(stage: ToolStage, argv: readonly string[], fileName: string, tool: ToolId = 'windowmasker'): boolean {
  if (tool !== 'windowmasker' || (stage !== 'mk_counts' && stage !== 'convert')) return false;
  return namedOutPath(argv) === normalizeRelName(fileName);
}

/**
 * The `-outfmt` values of the masking tools that write mask data for makeblastdb's `-mask_data`, with an extension for a
 * file of each. makeblastdb reads ASN.1 binary and text only: a file in another encoding (`maskinfo_xml`) is refused with
 * "Unknown encoding for mask data in '<file>'" (exit 1; NCBI makeblastdb.cpp:160-181, LOSAT's
 * `LOSAT/src/writedb/mask_info/reader.rs`), so no hand-over is offered for it.
 */
const MASK_DATA_FORMATS: Readonly<Record<string, string>> = Object.freeze({
  maskinfo_asn1_bin: 'asnb',
  maskinfo_asn1_text: 'asnt',
});

/** The tool runs that mask with a `-outfmt`: WindowMasker's Mask (`-ustat`); the masking tools of line WW WW4c join here. */
function masks(tool: ToolId, stage: ToolStage): boolean {
  return tool === 'windowmasker' && stage === 'mask';
}

/**
 * The name under which an output of a tool run goes to makeblastdb's `-mask_data` ("Use as -mask_data in makeblastdb"),
 * or undefined for an output that is not mask data: the `-out` file of a masking run whose `-outfmt` writes mask data
 * that makeblastdb reads (`maskinfo_asn1_bin`, `maskinfo_asn1_text`; not `maskinfo_xml`, see `MASK_DATA_FORMATS`), or
 * its standard output when it has no `-out`, named for the run (`windowmasker-run3.asnb`).
 */
export function maskDataName(
  run: { readonly tool: ToolId; readonly stage: ToolStage; readonly argv: readonly string[]; readonly number: number },
  entry: { readonly key: string; readonly name: string },
): string | undefined {
  if (!masks(run.tool, run.stage)) return undefined;
  const extension = MASK_DATA_FORMATS[toolOptionValue(run.argv, '-outfmt') ?? ''];
  if (extension === undefined) return undefined;
  const out = namedOutPath(run.argv);
  if (out === undefined) return entry.key === 'stdout' ? `${run.tool}-run${run.number}.${extension}` : undefined;
  return entry.key !== 'stdout' && entry.key !== 'stderr' && normalizeRelName(entry.name) === out ? entry.name : undefined;
}

/**
 * Extracts the concise error line from a `validateTool` refusal message.
 * NCBI `ncbi_app::args` errors prefix the error line with a multi-line `USAGE ... ========`
 * block; stripping everything up to the `========` rule leaves the verbatim CLI error sentence.
 */
export function toolValidationSummary(message: string): string {
  const sep = /^={8,}$/m.exec(message);
  const tail = sep === null ? message : message.slice(sep.index + sep[0].length);
  const first = tail
    .split(/\r?\n/)
    .map((line) => line.trim())
    .find((line) => line !== '');
  return first ?? message.trim();
}

/**
 * Maps a `validateTool` refusal message to the form fields and/or "Additional options" it names.
 * Handles both NCBI `ncbi_app::args` phrasing (`Argument "unit"`, `Incompatible with argument: \`dust'`,
 * `Unknown argument: "foo"`, `defined already: window`) and `-flag` tokens.
 */
export function refusedToolFields(
  message: string,
  extra: string,
): { readonly extra: boolean; readonly flags: ReadonlySet<string> } {
  const summary = toolValidationSummary(message);
  const named = new Set<string>();
  for (const match of summary.matchAll(/(?:^|[^\w-])(-[A-Za-z][A-Za-z0-9_]*)/g)) named.add(match[1]!);
  for (const match of summary.matchAll(/\bArgument\s+"([A-Za-z][A-Za-z0-9_]*)"/g)) named.add(`-${match[1]!}`);
  for (const match of summary.matchAll(/\bIncompatible with argument:\s+`([A-Za-z][A-Za-z0-9_]*)'/g)) named.add(`-${match[1]!}`);
  for (const match of summary.matchAll(/\bdefined already:\s+([A-Za-z][A-Za-z0-9_]*)\b/g)) named.add(`-${match[1]!}`);
  for (const match of summary.matchAll(/\bUnknown argument:\s+"([^"]+)"/g)) {
    const raw = match[1]!;
    named.add(raw.startsWith('-') ? raw : `-${raw}`);
  }

  let words: readonly string[];
  try {
    words = splitOptions(extra);
  } catch {
    words = [];
  }
  const typed = words.map((word) => word.split('=', 1)[0]!);
  const extraMatchesWord = (word: string): boolean => {
    if (word === '') return false;
    const bare = word.startsWith('-') ? word.slice(1) : word;
    return (
      named.has(word) ||
      named.has(`-${bare}`) ||
      summary.includes(`"${word}"`) ||
      summary.includes(`"${bare}"`) ||
      summary.includes(`'${word}'`) ||
      summary.includes(`\`${bare}'`)
    );
  };
  const extraRefused =
    message.startsWith('Additional options') ||
    typed.some(extraMatchesWord) ||
    words.some(extraMatchesWord) ||
    (summary.includes('Registry file') && typed.includes('-conffile'));
  return {
    extra: extraRefused,
    flags: new Set([...named].filter((flag) => !typed.includes(flag))),
  };
}
