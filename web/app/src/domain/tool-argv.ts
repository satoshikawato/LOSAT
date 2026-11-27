// Single command-line builders for the tools (windowmasker: WW3b DESIGN §2.5, brief §4.4;
// blast_formatter: RG-F3p7): the same argv array is validated by `validateTool`, run by `runTool`,
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
