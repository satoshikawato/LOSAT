// Single command-line builder for windowmasker (WW3b DESIGN §2.5, brief §4.4):
// the same argv array is validated by `validateTool`, run by `runTool`, displayed beside
// "Run LOSAT" as `losat windowmasker ...` and shown in Run details. Defaults are omitted;
// Make counts and Convert counts write `-out wmasker.<sformat>` until the user edits it.
import type { ToolDescription } from 'losat/domain';
import { splitOptions } from './argv';
import { writtenValue, type FormValues } from './parameters';
import {
  allStageParameters,
  defaultCountsOutName,
  stageOptionParameters,
  toToolParameterField,
  WINDOWMASKER_BLOCK_FLAGS,
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

/**
 * Builds `['windowmasker', stageFlag, ...]` for validation, execution and display.
 * Throws if `extra` contains a word that ABI's NUL-joined argv cannot transmit (`''` or NUL).
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
  const ordered = allStageParameters(description, stage).filter(
    (param) => param.flag !== '-in' && param.flag !== '-out' && param.flag !== '-ustat',
  );
  if (ordered.length > 0) {
    for (const param of ordered) {
      const field = toToolParameterField(param);
      const known = param.defaultValue === undefined ? undefined : { defaultValue: param.defaultValue };
      const written = writtenValue(field, values[param.flag], known);
      if (written === true) argv.push(param.flag);
      else if (typeof written === 'string') argv.push(param.flag, written);
    }
  } else {
    // Before `-xmlhelp` arrives, still write any explicitly set values so argv stays faithful.
    for (const [flag, val] of Object.entries(values)) {
      if (flag === '-in' || flag === '-out' || flag === '-ustat') continue;
      if (val === true) argv.push(flag);
      else if (typeof val === 'string' && val.trim() !== '') argv.push(flag, val.trim());
    }
  }

  for (const word of input.extra ?? []) {
    const unread = word === '' ? "''" : word.includes('\0') ? 'a NUL character' : undefined;
    if (unread !== undefined) throw new Error(`Additional options: ${unread} is not accepted`);
    argv.push(word);
  }

  return Object.freeze(argv);
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
  let count = 0;
  for (const param of stageOptionParameters(description, stage)) {
    if (WINDOWMASKER_BLOCK_FLAGS.has(param.flag)) continue;
    const field = toToolParameterField(param);
    const known = param.defaultValue === undefined ? undefined : { defaultValue: param.defaultValue };
    if (writtenValue(field, values[param.flag], known) !== undefined) count++;
  }
  if (extra.trim() !== '') count++;
  return count;
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

/** Identifies which option(s) (`-out`, `-logfile`) named an output file produced by a tool run. */
export function toolFileOptionLabel(argv: readonly string[], fileName: string): string | undefined {
  const target = normalizeRelName(fileName);
  const out = toolOptionValue(argv, '-out');
  const log = toolOptionValue(argv, '-logfile');
  const matchesOut = out !== undefined && out !== '-' && normalizeRelName(out) === target;
  const matchesLog = log !== undefined && normalizeRelName(log) === target;
  if (matchesOut && matchesLog) return '-out, -logfile';
  if (matchesOut) return '-out';
  if (matchesLog) return '-logfile';
  return undefined;
}

/** Whether `fileName` is the counts output file of `-mk_counts` or `-convert` (eligible for "Use in BLASTN" / "Mask with these counts"). */
export function isCountsOutputFile(stage: WindowMaskerStage, argv: readonly string[], fileName: string): boolean {
  if (stage !== 'mk_counts' && stage !== 'convert') return false;
  const out = toolOptionValue(argv, '-out');
  return out !== undefined && out !== '-' && normalizeRelName(out) === normalizeRelName(fileName);
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
