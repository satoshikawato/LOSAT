import type { OutputFormat, ProgramId, ToolDescription } from 'losat/domain';
import type { ToolId } from './tools';
// The verification badge of a run (plan §6.1, design §10.3, PD-LOSAT-WEB-APP-BOUNDARY
// "Compatibility contract"). A run matches a certified profile only when its program and
// options were compared byte for byte with NCBI BLAST+ for every output format (the records
// compare the formats that every run writes, so a run that also writes the output of its own
// -outfmt with fields is outside them, W10), and its
// runtime path in the browser was checked against the native engine with its thread count.
// Auto's count depends on the user's processors (serial below 20,000 bytes of input, then half
// of the processors, at most 8 below 1 MB and 32 from 1 MB), so V-BR runs every count from 1 to
// 32 for each program (WA-RS-03).
// The table that says what was compared is generated from docs/web/verification_cells.tsv
// and the certification records at build time (build/verification.ts); nobody writes it
// by hand. The badge speaks about the engine, the options and the runtime path, never
// about the user's inputs: no run is compared with NCBI when it runs.

/** What the records show for one program. */
export interface ProgramVerification {
  /**
   * Option sets (keys of `optionKey`) whose outputs were compared byte for byte with NCBI
   * BLAST+, with the formats compared (a record that checked an approved exception with
   * NCBI's oracle counts, with its note).
   */
  readonly optionSets: Readonly<Record<string, readonly OutputFormat[]>>;
  /** The browser runtime checked against the native engine (V-BR), if it was. */
  readonly browser?: BrowserVerification;
}

export interface BrowserVerification {
  readonly formats: readonly OutputFormat[];
  readonly paths: readonly ('serial' | 'threaded')[];
  readonly threads: readonly number[];
  readonly browsers: readonly string[];
}

/** What the records show for one CLI tool (`cases.tsv` frozen against NCBI BLAST+). */
export interface ToolVerification {
  /** Sorted canonical keys (`toolOptionKey`) that have an NCBI-frozen exit-0 row. */
  readonly optionSets: readonly string[];
}

export interface VerificationTable {
  /** NCBI BLAST+ version of the comparisons. */
  readonly ncbi: string;
  /** The records the table was generated from, with their SHA-256. */
  readonly sources: readonly { readonly path: string; readonly sha256: string; readonly optionSets: number }[];
  readonly programs: Readonly<Partial<Record<ProgramId, ProgramVerification>>>;
  /** CLI tools besides BLAST searches (windowmasker, blast_formatter, dustmasker, segmasker, convert2blastmask). */
  readonly tools?: Readonly<Partial<Record<ToolId, ToolVerification>>>;
}

/**
 * Flags that do not belong to a profile: the inputs, the outputs and the threads, the file of
 * `-logfile`, which names a file of the run's messages as `-out` names its output, and the search
 * strategy file of `-export_search_strategy`, a file of the run's requests: the search's outputs do
 * not change with either (line RG E1w-3b; line RG F1w). `-conffile` and `-dryrun` change what the
 * run reads or does, and stay options (no record compares them).
 */
const NOT_OPTIONS = new Set(['-query', '-subject', '-out', '-outfmt', '-num_threads', '-logfile', '-export_search_strategy']);
/** Flags whose value is a coordinate range of the inputs (S11): a profile has the flag, not the range. */
const RANGE_OPTIONS = new Set(['-query_loc', '-subject_loc']);
export const RANGE_VALUE = '<range>';
/** Placeholder for the values of TOOL_FILE_OPTIONS in tool option keys. */
export const FILE_VALUE = '<file>';
/** Tool flags that name only input/output/log files and do not change the option set. */
const TOOL_NOT_OPTIONS = new Set(['-in', '-archive', '-out', '-logfile']);
/**
 * Tool flags whose values are names: files in the run directory (`-ustat`, `-ids`, `-exclude_ids`, makeblastdb's
 * `-mask_data` and `-taxid_map`) and makeblastdb's database `-title`. The flag stays in the key, the name does not.
 */
const TOOL_FILE_OPTIONS = new Set(['-ustat', '-ids', '-exclude_ids', '-mask_data', '-taxid_map', '-title']);

export interface OptionGrammar {
  /** Whether the flag takes a value (the engine's `describe`). */
  takesValue(flag: string): boolean;
  /** The program's default task, which an argv may leave out or write. */
  readonly defaultTask?: string;
}

export interface ToolOptionGrammar extends OptionGrammar {
  /** Whether the flag takes an NCBI Boolean value (`T`/`true`/`yes` vs `F`/`false`/`no`). */
  isBoolean?(flag: string): boolean;
  /** Default value of the flag in `-xmlhelp`, omitted from the key when written explicitly. */
  defaultValue?(flag: string): string | undefined;
}

/** Builds a `ToolOptionGrammar` from a tool's parsed `-xmlhelp` description. */
export function toolOptionGrammar(description: ToolDescription | undefined): ToolOptionGrammar {
  const valued = new Set(['-logfile', '-conffile']);
  const booleans = new Set<string>();
  const defaults = new Map<string, string>();
  for (const param of description?.parameters ?? []) {
    if (param.takesValue) valued.add(param.flag);
    if (param.kind === 'Boolean') booleans.add(param.flag);
    if (param.defaultValue !== undefined && param.defaultValue !== '') defaults.set(param.flag, param.defaultValue);
  }
  return {
    takesValue: (flag) => valued.has(flag),
    isBoolean: (flag) => booleans.has(flag),
    defaultValue: (flag) => defaults.get(flag),
  };
}

function normalizeToolBool(value: string): string {
  const lower = value.trim().toLowerCase();
  if (lower === 't' || lower === 'true' || lower === 'y' || lower === 'yes' || lower === '1') return 'true';
  if (lower === 'f' || lower === 'false' || lower === 'n' || lower === 'no' || lower === '0') return 'false';
  return value;
}

/**
 * The canonical key of a tool argv (with or without the tool name first): its stage flag and
 * options other than `-in`, `-out`, `-logfile` and file names, with boolean spellings and explicit
 * defaults normalized. Returns undefined when the words cannot be read as options.
 */
export function toolOptionKey(words: readonly string[], grammar: ToolOptionGrammar): string | undefined {
  const pairs: [string, string][] = [];
  let i = words[0] !== undefined && !words[0].startsWith('-') ? 1 : 0;
  while (i < words.length) {
    const word = words[i]!;
    if (!word.startsWith('-')) return undefined;
    const eq = word.indexOf('=');
    let flag: string;
    let value: string;
    if (eq >= 0) {
      flag = word.slice(0, eq);
      value = word.slice(eq + 1);
      i++;
    } else {
      flag = word;
      const takesValue = TOOL_NOT_OPTIONS.has(flag) || grammar.takesValue(flag);
      if (takesValue && i + 1 >= words.length) return undefined;
      value = takesValue ? words[i + 1]! : '';
      i += takesValue ? 2 : 1;
    }
    if (TOOL_NOT_OPTIONS.has(flag)) continue;
    if (TOOL_FILE_OPTIONS.has(flag)) {
      pairs.push([flag, FILE_VALUE]);
      continue;
    }
    const isBool = grammar.isBoolean?.(flag) === true;
    const normalized = isBool ? normalizeToolBool(value) : value;
    const def = grammar.defaultValue?.(flag);
    if (def !== undefined && def !== '') {
      const normalizedDef = isBool ? normalizeToolBool(def) : def;
      if (normalized === normalizedDef) continue;
    }
    pairs.push([flag, normalized]);
  }
  // No option besides the inputs and outputs is the tool's default option set (`dustmasker -in a.fa`), a set of its own.
  pairs.sort(([fa, va], [fb, vb]) => (fa < fb ? -1 : fa > fb ? 1 : va < vb ? -1 : va > vb ? 1 : 0));
  return JSON.stringify(pairs);
}

/**
 * The key of the option set of an argv (with or without the program name first): its
 * options other than the inputs, outputs and threads, in flag order, with the default
 * task written and ranges replaced by RANGE_VALUE. Returns undefined when the words cannot
 * be read as options (an unknown flag or a missing value), so they never match.
 */
export function optionKey(words: readonly string[], grammar: OptionGrammar): string | undefined {
  const pairs: [string, string][] = [];
  let i = words[0] !== undefined && !words[0].startsWith('-') ? 1 : 0;
  while (i < words.length) {
    const flag = words[i]!;
    if (!flag.startsWith('-')) return undefined;
    // `-outfmt=6 qseqid` is the command line's form of `-outfmt '6 qseqid'`.
    if (NOT_OPTIONS.has(flag.split('=', 1)[0]!) && flag.includes('=')) {
      i++;
      continue;
    }
    const takesValue = NOT_OPTIONS.has(flag) || grammar.takesValue(flag);
    if (takesValue && i + 1 >= words.length) return undefined;
    const value = takesValue ? words[i + 1]! : '';
    i += takesValue ? 2 : 1;
    if (NOT_OPTIONS.has(flag)) continue;
    pairs.push([flag, RANGE_OPTIONS.has(flag) ? RANGE_VALUE : value]);
  }
  if (grammar.defaultTask !== undefined && !pairs.some(([flag]) => flag === '-task')) pairs.push(['-task', grammar.defaultTask]);
  pairs.sort(([fa, va], [fb, vb]) => (fa < fb ? -1 : fa > fb ? 1 : va < vb ? -1 : va > vb ? 1 : 0));
  return JSON.stringify(pairs);
}

export interface BadgeInput {
  readonly program: ProgramId;
  readonly argv: readonly string[];
  /** The formats the run wrote. */
  readonly formats: readonly OutputFormat[];
  /** The -outfmt value of the output of the run's own -outfmt (`customOutfmt`), which no record compares; undefined without one. */
  readonly custom?: string | undefined;
  readonly runtimePath: 'threaded' | 'serial' | 'fake' | undefined;
  readonly threads: number | undefined;
  readonly grammar: OptionGrammar;
}

export interface Badge {
  readonly level: 'certified' | 'outside' | 'development';
  readonly label: string;
  /** Why the run is outside the certified profiles, or what the certification covers. */
  readonly details: readonly string[];
  /** Approved exceptions that apply to the run's options (AGENTS.md). */
  readonly exceptions: readonly string[];
}

// The words of the badge (the maintainer, 2026-11-02: no "certified", which reads as an official
// certification; the badge says whether the run is supported, DW-36; 2026-11-04: no "verification"
// either, and the two marks are Supported and Not yet supported). The levels keep their names.
export const SUPPORTED_LABEL = 'Supported';
export const NOT_SUPPORTED_LABEL = 'Not yet supported';

export function verificationBadge(input: BadgeInput, table: VerificationTable): Badge {
  const exceptions = approvedExceptions(input.program, input.argv);
  if (input.runtimePath === 'fake') {
    return {
      level: 'development',
      label: 'Development build',
      details: ['The fake engine wrote these outputs. They are not search results.'],
      exceptions: [],
    };
  }
  const verification = table.programs[input.program];
  const key = optionKey(input.argv, input.grammar);
  const compared = key === undefined ? undefined : verification?.optionSets[key];
  const outside: string[] = [];
  const custom = input.custom === undefined ? [] : [input.custom];
  const missing = [...input.formats.filter((format) => !(compared ?? []).includes(format)), ...custom];
  if (compared === undefined) {
    outside.push(`These options of ${input.program} were not compared with NCBI BLAST+ ${table.ncbi} on LOSAT's test fixtures.`);
  } else if (missing.length > 0) {
    outside.push(`These options were compared with NCBI BLAST+ ${table.ncbi} for outfmt ${list(compared)}, not for outfmt ${list(missing)}.`);
  }
  const browser = verification?.browser;
  if (browser === undefined) {
    outside.push(`The browser runtime of ${input.program} was not checked against the native engine.`);
  } else {
    const formats = [...input.formats.filter((format) => !browser.formats.includes(format)), ...custom];
    if (formats.length > 0) outside.push(`The browser runtime was not checked for outfmt ${list(formats)}.`);
    if (input.runtimePath === undefined || !browser.paths.includes(input.runtimePath)) {
      outside.push(`The ${input.runtimePath ?? 'unknown'} runtime path was not checked in the browser.`);
    }
    if (input.threads === undefined || !browser.threads.includes(input.threads)) {
      outside.push(
        `This run used ${input.threads ?? 'an unknown number of'} ${input.threads === 1 ? 'thread' : 'threads'}; ` +
          `the browser runtime was checked with ${counts(browser.threads)} threads.`,
      );
    }
  }
  if (outside.length > 0) return { level: 'outside', label: NOT_SUPPORTED_LABEL, details: outside, exceptions };
  // The counts of the run's path: the serial module runs one thread, the threaded module two or more.
  const pathThreads = browser!.threads.filter((n) => (input.runtimePath === 'serial' ? n === 1 : n > 1));
  return {
    level: 'certified',
    label: SUPPORTED_LABEL,
    details: [
      `The program and these options were compared byte for byte with NCBI BLAST+ ${table.ncbi} for outfmt ${list(input.formats)} on LOSAT's test fixtures, apart from the line that names the program's version (LOSAT names its own).`,
      `The ${input.runtimePath} runtime in the browser was checked against the native engine with ${counts(pathThreads)} ` +
        `${pathThreads.length === 1 && pathThreads[0] === 1 ? 'thread' : 'threads'} in ${list(browser!.browsers)}.`,
      'This is a statement about the engine build and the options, not a comparison of your inputs with NCBI.',
    ],
    exceptions,
  };
}

export interface ToolBadgeInput {
  readonly tool: ToolId;
  readonly argv: readonly string[];
  readonly runtimePath: 'serial' | 'fake' | undefined;
  readonly grammar: ToolOptionGrammar;
}

export function toolVerificationBadge(input: ToolBadgeInput, table: VerificationTable): Badge {
  if (input.runtimePath === 'fake') {
    return {
      level: 'development',
      label: 'Development build',
      details: ['The fake engine wrote these outputs. They are not tool results.'],
      exceptions: [],
    };
  }
  const key = toolOptionKey(input.argv, input.grammar);
  const optionSets = table.tools?.[input.tool]?.optionSets;
  if (key !== undefined && optionSets !== undefined && optionSets.includes(key)) {
    return {
      level: 'certified',
      label: SUPPORTED_LABEL,
      details: [
        `This stage and these options of ${input.tool} were compared byte for byte with NCBI BLAST+ ${table.ncbi} on LOSAT's test fixtures.`,
        'This is a statement about the engine build and the options, not a comparison of your inputs with NCBI.',
      ],
      exceptions: [],
    };
  }
  return {
    level: 'outside',
    label: NOT_SUPPORTED_LABEL,
    details: [
      `This stage and these options of ${input.tool} were not compared with NCBI BLAST+ ${table.ncbi} on LOSAT's test fixtures.`,
    ],
    exceptions: [],
  };
}

/**
 * The approved exceptions of AGENTS.md that change a run's results: a local subject searched
 * with a non-default subject genetic code (TBLASTN `PD-TLOSAN-LOCAL-GENCODE-32`, TBLASTX).
 */
export function approvedExceptions(program: ProgramId, argv: readonly string[]): readonly string[] {
  // The approved -db_gencode exception applies only to local -subject searches, not -db database searches.
  if (argv.some((word) => word === '-db' || word.startsWith('-db='))) return [];
  const at = argv.lastIndexOf('-db_gencode');
  const code = at >= 0 ? argv[at + 1]?.trim() : undefined;
  if (code === undefined || code === '1' || (program !== 'tblastn' && program !== 'tblastx')) return [];
  const decision = program === 'tblastn' ? ' (PD-TLOSAN-LOCAL-GENCODE-32)' : '';
  return [
    `Approved exception${decision}: LOSAT translates the local subjects with -db_gencode ${code}. NCBI BLAST+ ` +
      'treats -db_gencode differently for local -subject sequences, so its results can differ for this genetic code.',
  ];
}

/** Thread counts in words, a run of three or more consecutive counts as a range: "1, 2 and 4", "2 to 32". */
function counts(numbers: readonly number[]): string {
  const words: string[] = [];
  for (let i = 0; i < numbers.length; ) {
    let j = i;
    while (j + 1 < numbers.length && numbers[j + 1] === numbers[j]! + 1) j++;
    if (j - i >= 2) words.push(`${numbers[i]} to ${numbers[j]}`);
    else for (let k = i; k <= j; k++) words.push(String(numbers[k]));
    i = j + 1;
  }
  return list(words);
}

function list(items: readonly (string | number)[]): string {
  if (items.length <= 1) return items.join('');
  return `${items.slice(0, -1).join(', ')} and ${items.at(-1)}`;
}
