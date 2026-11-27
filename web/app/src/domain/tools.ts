// Tool descriptors and helpers for CLI tools besides BLAST searches (WW3b DESIGN §1-§2, RG-F3p7):
// windowmasker's three stages (-mk_counts, -ustat, -convert), blast_formatter (-archive),
// their option filters driven by `-xmlhelp`, mapping into ParameterField / ParameterDescription,
// and output file naming. Nothing here lists option defaults, value lists or ranges: all
// of those come from `ToolEngine.describeTool(...)`.
import type { PrintedKind, ToolDescription, ToolParameterDescription } from 'losat/domain';
import { exportPrefix, FORMAT_NAMES } from './output-format';
import type { ParameterField } from './programs';

/** Structural mirror of `ParameterDescription` so the domain layer stays independent of `ports`. */
export interface ToolEngineOption {
  readonly flag: string;
  readonly help: string;
  readonly takesValue: boolean;
  readonly defaultValue?: string;
  readonly choices?: readonly string[];
}

export type ToolId = 'windowmasker' | 'blast_formatter';

export type WindowMaskerStage = 'mk_counts' | 'mask' | 'convert';

export type ToolStage = WindowMaskerStage | 'format';

/** Human-readable tool label for Results selector and Run details. */
export function toolLabel(tool: ToolId): string {
  return tool === 'windowmasker' ? 'WindowMasker' : 'Format Archive';
}

/** Stage label for Results selector and Run details. */
export function toolStageLabel(tool: ToolId, stage: ToolStage): string {
  if (tool === 'blast_formatter') return 'blast_formatter';
  return windowMaskerStageById(stage as WindowMaskerStage).label;
}

export interface ToolStageDescriptor {
  readonly id: WindowMaskerStage;
  readonly label: string;
  /** The CLI flag identifying this stage: `-mk_counts`, `-ustat`, or `-convert`. */
  readonly flag: string;
  readonly fallbackSummary: string;
}

export const WINDOWMASKER_STAGES: readonly ToolStageDescriptor[] = Object.freeze([
  { id: 'mk_counts', label: 'Make counts', flag: '-mk_counts', fallbackSummary: 'generate frequency counts for a database' },
  { id: 'mask', label: 'Mask', flag: '-ustat', fallbackSummary: 'Window based sequence masker' },
  { id: 'convert', label: 'Convert counts', flag: '-convert', fallbackSummary: 'convert counts between different formats' },
]);

/** The stage descriptor for `stage`. */
export function windowMaskerStageById(stage: WindowMaskerStage): ToolStageDescriptor {
  const found = WINDOWMASKER_STAGES.find((entry) => entry.id === stage);
  if (found === undefined) throw new Error(`unknown windowmasker stage: ${stage}`);
  return found;
}

/** NCBI's sentence for a stage, read from `-xmlhelp` once available. */
export function windowMaskerStageSummary(stage: WindowMaskerStage, description: ToolDescription | undefined): string {
  const descriptor = windowMaskerStageById(stage);
  if (description === undefined) return descriptor.fallbackSummary;
  if (stage === 'mask') return description.description || descriptor.fallbackSummary;
  return description.parameters.find((param) => param.flag === descriptor.flag)?.help || descriptor.fallbackSummary;
}

/** Options placed in the Input, Counts and Output blocks rather than the "WindowMasker options" section. */
export const WINDOWMASKER_BLOCK_FLAGS: ReadonlySet<string> = new Set(['-in', '-infmt', '-ustat', '-out', '-sformat', '-outfmt']);

/** Stage-selector flags (one is always chosen by the stage switcher). */
export const WINDOWMASKER_STAGE_FLAGS: ReadonlySet<string> = new Set(['-mk_counts', '-ustat', '-convert']);

/** String-typed options in windowmasker.xml whose values name files in the run directory. */
export const WINDOWMASKER_FILE_STRING_FLAGS: ReadonlySet<string> = new Set(['-ustat', '-ids', '-exclude_ids']);

/** Help and version flags, entered via Additional options. */
export const WINDOWMASKER_HELP_FLAGS: ReadonlySet<string> = new Set(['-h', '-help', '-xmlhelp']);
export const WINDOWMASKER_VERSION_FLAGS: ReadonlySet<string> = new Set(['-version', '-version-full', '-version-full-xml', '-version-full-json']);

/** Options placed in the Archive and Output blocks of Format Archive rather than "Formatting options". */
export const BLAST_FORMATTER_BLOCK_FLAGS: ReadonlySet<string> = new Set(['-archive', '-rid', '-outfmt', '-out']);

/** NCBI `-help` usage order for blast_formatter's options (excluding `-rid`, which LOSAT omits from `-xmlhelp` and refuses). */
export const BLAST_FORMATTER_USAGE_ORDER: readonly string[] = Object.freeze([
  '-archive',
  '-outfmt',
  '-show_gis',
  '-num_descriptions',
  '-num_alignments',
  '-line_length',
  '-html',
  '-sorthits',
  '-sorthsps',
  '-max_target_seqs',
  '-out',
  '-parse_deflines',
]);

/**
 * The `-outfmt` choices of the Format Archive tab: the formats that LOSAT's `blast_formatter` writes
 * (`docs/manual/src/programs/blast_formatter.md`); "Other / custom format…" takes any other value.
 * Each PR that makes `blast_formatter` write another format adds it here, and the e2e test runs every
 * choice on the engine (`tests/e2e/blast-formatter-tool.spec.ts`).
 */
export const BLAST_FORMATTER_OUTFMT_CHOICES: readonly { readonly value: string; readonly label: string }[] = Object.freeze(
  ([8, 9, 12] as const).map((num) => ({
    value: String(num),
    label: `${num} — ${FORMAT_NAMES[num]}`,
  })),
);

/** Default `-outfmt` selected in the Format Archive tab (the first of the choices). */
export const BLAST_FORMATTER_DEFAULT_OUTFMT = BLAST_FORMATTER_OUTFMT_CHOICES[0]!.value;

/** All group-defined non-hidden options of `blast_formatter` in NCBI `-help` usage order. */
export function allBlastFormatterParameters(description: ToolDescription | undefined): readonly ToolParameterDescription[] {
  if (description === undefined) return [];
  return sortByUsage(
    description.parameters.filter((param) => param.group !== undefined && !param.hidden),
    BLAST_FORMATTER_USAGE_ORDER,
  );
}

/** Options of `blast_formatter` shown under the collapsible "Formatting options" section. */
export function blastFormatterOptionParameters(description: ToolDescription | undefined): readonly ToolParameterDescription[] {
  return allBlastFormatterParameters(description).filter((param) => !BLAST_FORMATTER_BLOCK_FLAGS.has(param.flag));
}

/** NCBI `-help` usage order for windowmasker's options; any unknown future option sorts after these in XML order. */
export const WINDOWMASKER_USAGE_ORDER: readonly string[] = Object.freeze([
  '-ustat',
  '-in',
  '-out',
  '-checkdup',
  '-fa_list',
  '-mem',
  '-meta',
  '-unit',
  '-genome_size',
  '-window',
  '-t_extend',
  '-t_thres',
  '-set_t_high',
  '-set_t_low',
  '-parse_seqids',
  '-outfmt',
  '-t_high',
  '-t_high_pct',
  '-t_thres_pct',
  '-t_extend_pct',
  '-t_low',
  '-t_low_pct',
  '-infmt',
  '-exclude_ids',
  '-ids',
  '-text_match',
  '-sformat',
  '-smem',
  '-dust',
  '-dust_level',
  '-mk_counts',
  '-convert',
]);

/** Sorts parameters in a tool's NCBI `-help` usage order (`order`), keeping relative XML order for any option not in it. */
function sortByUsage(parameters: readonly ToolParameterDescription[], order: readonly string[]): readonly ToolParameterDescription[] {
  const rank = (flag: string) => {
    const index = order.indexOf(flag);
    return index < 0 ? order.length : index;
  };
  return parameters
    .map((param, index) => ({ param, index, rank: rank(param.flag) }))
    .sort((a, b) => (a.rank !== b.rank ? a.rank - b.rank : a.index - b.index))
    .map((entry) => entry.param);
}

/** Flags excluded by `stage` according to `-xmlhelp`'s `<dependencies>`. */
export function stageExcludedFlags(description: ToolDescription, stage: WindowMaskerStage): ReadonlySet<string> {
  const stageFlag = windowMaskerStageById(stage).flag;
  const excluded = new Set<string>();
  for (const [first, second] of description.excludes) {
    if (first === stageFlag) excluded.add(second);
    else if (second === stageFlag) excluded.add(first);
  }
  return excluded;
}

/** Sorts windowmasker's parameters in NCBI `-help` usage order, keeping relative XML order for any new option. */
export function sortToolParameters(parameters: readonly ToolParameterDescription[]): readonly ToolParameterDescription[] {
  return sortByUsage(parameters, WINDOWMASKER_USAGE_ORDER);
}

type MaybeToolDescription = ToolDescription | undefined;

function filterStageParameters(
  description: MaybeToolDescription,
  stage: WindowMaskerStage,
  includeBlock: boolean,
): readonly ToolParameterDescription[] {
  if (description === undefined) return [];
  const excluded = stageExcludedFlags(description, stage);
  const allowed = description.parameters.filter(
    (param) =>
      param.group !== undefined &&
      !param.hidden &&
      !WINDOWMASKER_STAGE_FLAGS.has(param.flag) &&
      (includeBlock || !WINDOWMASKER_BLOCK_FLAGS.has(param.flag)) &&
      !excluded.has(param.flag),
  );
  return sortToolParameters(allowed);
}

/**
 * Every option allowed in `stage` that belongs under "WindowMasker options": all group-1 options
 * from `-xmlhelp` that are not stage flags, not block options, and not excluded by the stage flag,
 * sorted in NCBI `-help` usage order.
 */
export function stageOptionParameters(description: MaybeToolDescription, stage: WindowMaskerStage): readonly ToolParameterDescription[] {
  return filterStageParameters(description, stage, false);
}

/** Every non-stage group-1 option allowed in `stage` (including block options like `-infmt`, `-sformat`, `-outfmt`), in usage order. */
export function allStageParameters(description: MaybeToolDescription, stage: WindowMaskerStage): readonly ToolParameterDescription[] {
  return filterStageParameters(description, stage, true);
}

/** Looks up a block option (`-infmt`, `-sformat`, `-outfmt`, etc.) if allowed in `stage`. */
export function stageBlockParameter(description: MaybeToolDescription, stage: WindowMaskerStage, flag: string): ToolParameterDescription | undefined {
  if (description === undefined) return undefined;
  if (stageExcludedFlags(description, stage).has(flag)) return undefined;
  return description.parameters.find((param) => param.flag === flag);
}

/**
 * Maps an `-xmlhelp` parameter to a `ParameterField` for the shared ParameterField component:
 * label is the option name (`label: ''`, rendered with OptionName), with NCBI's sentence as `hint`:
 * the first line of its description (DW-55, DW-57). NCBI's `-help` writes some descriptions over
 * several lines (blast_formatter's `-sorthits`: "Sorting option for hits:", then one line a value);
 * the "?" holds every line (`toolParameterHelpText`), as the Search tab's "?" holds the CLI's help.
 */
export function toToolParameterField(param: ToolParameterDescription): ParameterField {
  const hint = (param.help.split('\n')[0] ?? '').trim();
  if (param.kind === 'Flag') {
    return { flag: param.flag, label: '', kind: 'flag', hint };
  }
  if (param.kind === 'Boolean') {
    return { flag: param.flag, label: '', kind: 'boolean', select: true, hint };
  }
  if (param.choices !== undefined && param.choices.length > 0) {
    return { flag: param.flag, label: '', kind: 'choice', choices: param.choices, hint };
  }
  if (WINDOWMASKER_FILE_STRING_FLAGS.has(param.flag) || param.kind === 'File_In') {
    return { flag: param.flag, label: '', kind: 'file-text', hint };
  }
  if (param.kind === 'Integer' || param.kind === 'Int8') {
    return { flag: param.flag, label: '', kind: 'text', inputMode: 'numeric', hint };
  }
  if (param.kind === 'Real') {
    return { flag: param.flag, label: '', kind: 'text', inputMode: 'decimal', hint };
  }
  return { flag: param.flag, label: '', kind: 'text', hint };
}

/** Formats the "?" HelpTip text for a tool parameter, including numeric range / constraint details. */
export function toolParameterHelpText(param: ToolParameterDescription): string {
  const parts: string[] = [];
  if (param.help) parts.push(param.help);
  if (param.constraintDescription) {
    parts.push(`Allowed: ${param.constraintDescription}`);
  } else if (param.range !== undefined) {
    parts.push(`Allowed: ${param.range.min}..${param.range.max}`);
  }
  if (param.defaultValue !== undefined && param.defaultValue !== '') {
    parts.push(`Default: ${param.defaultValue}`);
  }
  return parts.join('\n');
}

/** Maps an `-xmlhelp` parameter into the option description shape that `FormModel` reads. */
export function toToolEngineOption(param: ToolParameterDescription): ToolEngineOption {
  return {
    flag: param.flag,
    help: toolParameterHelpText(param),
    takesValue: param.takesValue,
    ...(param.defaultValue === undefined ? {} : { defaultValue: param.defaultValue }),
    ...(param.choices === undefined ? {} : { choices: param.choices }),
  };
}

/** Flags accepted by NCBI C++ Toolkit's `CNcbiApplication` that are not always listed in `windowmasker.xml`. */
const TOOL_COMMON_VALUED_FLAGS: ReadonlySet<string> = new Set(['-logfile', '-conffile']);

/** Grammar predicate for `CommandText` on a tool argv. */
export function toolTakesValue(description: ToolDescription | undefined): (flag: string) => boolean {
  const valued = new Set(TOOL_COMMON_VALUED_FLAGS);
  for (const param of description?.parameters ?? []) {
    if (param.takesValue) valued.add(param.flag);
  }
  return (flag: string) => valued.has(flag);
}

/** Whether the argv asks the tool to print help or version instead of running a stage. */
export function toolPrintsKind(argv: readonly string[]): PrintedKind | null {
  for (const word of argv.slice(1)) {
    const flag = word.split('=', 1)[0]!;
    if (WINDOWMASKER_HELP_FLAGS.has(flag)) return 'help';
    if (WINDOWMASKER_VERSION_FLAGS.has(flag)) return 'version';
  }
  return null;
}

/** Download file name for a tool run output or file (`runFileName` pattern). */
export function toolOutputFileName(
  run: { readonly number: number; readonly tool: ToolId; readonly stage: ToolStage; readonly queuedAt: number },
  name: string,
): string {
  const safe = name.replace(/\/+$/, '').replaceAll('/', '_');
  return `${exportPrefix({ number: run.number, program: run.tool, task: run.stage, queuedAt: run.queuedAt })}.${safe}`;
}

/** Default counts output file name for Make counts and Convert counts (`wmasker.<sformat>`). */
export function defaultCountsOutName(sformat: string | undefined): string {
  const format = sformat?.trim();
  return `wmasker.${format ? format : 'ascii'}`;
}

/** Checks whether bytes decode as UTF-8 with no NUL byte (DESIGN §3). */
export function decodeUtf8Text(bytes: Uint8Array): string | undefined {
  for (let i = 0; i < bytes.length; i++) {
    if (bytes[i] === 0) return undefined;
  }
  try {
    return new TextDecoder('utf-8', { fatal: true }).decode(bytes);
  } catch {
    return undefined;
  }
}

/**
 * Decodes a prefix slice of a larger file as UTF-8 without failing on a multi-byte sequence split at
 * the slice boundary: trims up to 3 trailing incomplete UTF-8 bytes when `truncated` is true.
 */
export function decodeUtf8Preview(bytes: Uint8Array, truncated: boolean): string | undefined {
  let end = bytes.length;
  if (truncated && end > 0) {
    // Walk back over trailing UTF-8 continuation bytes (10xxxxxx) to inspect the lead byte.
    let i = end - 1;
    while (i >= Math.max(0, end - 4) && (bytes[i]! & 0xc0) === 0x80) i--;
    if (i >= 0) {
      const lead = bytes[i]!;
      const expected =
        (lead & 0x80) === 0
          ? 1
          : (lead & 0xe0) === 0xc0
            ? 2
            : (lead & 0xf0) === 0xe0
              ? 3
              : (lead & 0xf8) === 0xf0
                ? 4
                : 0;
      if (expected > 0 && end - i < expected) end = i;
    }
  }
  return decodeUtf8Text(bytes.subarray(0, end));
}
