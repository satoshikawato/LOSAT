// Tool descriptors and helpers for CLI tools besides BLAST searches (WW3b DESIGN §1-§2, RG-F3p7, line WW
// WW4b2 and WW4c): windowmasker's three stages (-mk_counts, -ustat, -convert), blast_formatter (-archive),
// makeblastdb (its options in the groups of NCBI's `-help`), the Masking tab's dustmasker, segmasker and
// convert2blastmask, their option filters driven by `-xmlhelp`, mapping into ParameterField /
// ParameterDescription, and output file naming. Nothing here lists option defaults, value lists or
// ranges: all of those come from `ToolEngine.describeTool(...)`.
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

/** The tools of the Masking tab (WW4c), in the order of its tool switcher. */
export type MaskingTool = 'dustmasker' | 'segmasker' | 'convert2blastmask';

export type ToolId = 'windowmasker' | 'blast_formatter' | 'makeblastdb' | MaskingTool;

export type WindowMaskerStage = 'mk_counts' | 'mask' | 'convert';

/**
 * A run's stage: windowmasker's three, `format` for blast_formatter, `build` for makeblastdb; a run of the
 * Masking tab is a `mask` run.
 */
export type ToolStage = WindowMaskerStage | 'format' | 'build';

export interface MaskingToolDescriptor {
  readonly id: MaskingTool;
  readonly label: string;
  /** The sentence shown until `-xmlhelp` has arrived (its `<description>`, as the tool writes it). */
  readonly fallbackSummary: string;
}

/** The Masking tab's tool switcher, DustMasker first (selected by default). */
export const MASKING_TOOLS: readonly MaskingToolDescriptor[] = Object.freeze([
  { id: 'dustmasker', label: 'DustMasker', fallbackSummary: 'Low complexity region masker based on Symmetric DUST algorithm' },
  { id: 'segmasker', label: 'SegMasker', fallbackSummary: 'Low complexity region masker based on the SEG algorithm' },
  {
    id: 'convert2blastmask',
    label: 'convert2blastmask',
    fallbackSummary: 'Convert masking information in lower-case masked FASTA input to file formats suitable for makeblastdb',
  },
]);

/** NCBI's sentence for a tool of the Masking tab, read from `-xmlhelp` once available. */
export function maskingToolSummary(tool: MaskingTool, description: ToolDescription | undefined): string {
  return description?.description || MASKING_TOOLS.find((entry) => entry.id === tool)!.fallbackSummary;
}

/** Whether `tool` is one of the Masking tab's tools. */
export function isMaskingTool(tool: ToolId): tool is MaskingTool {
  return MASKING_TOOLS.some((entry) => entry.id === tool);
}

/** Human-readable tool label for Results selector and Run details. */
export function toolLabel(tool: ToolId): string {
  if (tool === 'windowmasker') return 'WindowMasker';
  if (tool === 'blast_formatter') return 'Format Archive';
  if (tool === 'makeblastdb') return 'makeblastdb';
  return MASKING_TOOLS.find((entry) => entry.id === tool)!.label;
}

/** Stage label for Results selector and Run details: windowmasker's stage, else the command-line tool's name. */
export function toolStageLabel(tool: ToolId, stage: ToolStage): string {
  if (tool !== 'windowmasker') return tool;
  return windowMaskerStageById(stage as WindowMaskerStage).label;
}

/**
 * What Run details and the run's line in Results call the tool of a run: WindowMasker with its stage
 * ("WindowMasker · Mask"), Format Archive, whose tab is named otherwise than its tool, with the tool's
 * name ("Format Archive (blast_formatter)"), and the other tools by the name of their tab or tool
 * ("makeblastdb", "DustMasker").
 */
export function toolRunTitle(tool: ToolId, stage: ToolStage): string {
  if (tool === 'windowmasker') return `${toolLabel(tool)} · ${toolStageLabel(tool, stage)}`;
  if (tool === 'blast_formatter') return `${toolLabel(tool)} (${toolStageLabel(tool, stage)})`;
  return toolLabel(tool);
}

/**
 * Whether the tool writes its files into the run's directory (makeblastdb's database, ABI v2 version 17), rather than
 * the files that its options name as it writes them (`-out`, `-logfile`): Results lists them, with a download of all.
 */
export function toolWritesDirectory(tool: ToolId): boolean {
  return tool === 'makeblastdb';
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

/**
 * String-typed options whose values name files in the run directory: windowmasker's `-ustat`, `-ids`,
 * `-exclude_ids`, and makeblastdb's `-mask_data` (a comma-separated list of them, `TOOL_FILE_LIST_FLAGS`).
 */
export const TOOL_FILE_STRING_FLAGS: ReadonlySet<string> = new Set(['-ustat', '-ids', '-exclude_ids', '-mask_data']);

/** Options whose value is a comma-separated list of files of the run directory (makeblastdb's `-mask_data`). */
export const TOOL_FILE_LIST_FLAGS: ReadonlySet<string> = new Set(['-mask_data']);

/**
 * The help flags of every tool of LOSAT Web (`-xmlhelp` prints the XML, the others the usage). `-help-full` is one of
 * them for every tool, as NCBI BLAST+ 2.17.0 and LOSAT print each tool's help for it with exit status 0, makeblastdb's
 * listing its hidden `-randomize`, which `-help` does not list (review of WW4c, F4, probed 2026-11-28). Additional
 * options take them and the version flags. The one list of them: `toolPrintsKind` reads it, and every reader of a tool
 * run's help (the drafts, FakeEngine) asks `toolPrintsKind`.
 */
export const WINDOWMASKER_HELP_FLAGS: ReadonlySet<string> = new Set(['-h', '-help', '-help-full', '-xmlhelp']);

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

/** Options placed in the Input and Output blocks of the Masking tab rather than "<Tool> options". */
export const MASKING_BLOCK_FLAGS: ReadonlySet<string> = new Set(['-in', '-infmt', '-out', '-outfmt']);

/** NCBI `-help` usage order of each Masking tool's options (the tool's `USAGE` line); an unknown future option sorts after these in XML order. */
export const MASKING_USAGE_ORDER: Readonly<Record<MaskingTool, readonly string[]>> = Object.freeze({
  dustmasker: Object.freeze(['-in', '-out', '-window', '-level', '-linker', '-infmt', '-outfmt', '-parse_seqids', '-hard_masking']),
  segmasker: Object.freeze(['-in', '-out', '-infmt', '-parse_seqids', '-outfmt', '-window', '-locut', '-hicut']),
  convert2blastmask: Object.freeze(['-in', '-out', '-outfmt', '-parse_seqids', '-masking_algorithm', '-masking_options']),
});

/** Every option of a Masking tool from `-xmlhelp` but help and version (which Additional options take), in NCBI `-help` usage order. */
export function allMaskingParameters(description: ToolDescription | undefined, tool: MaskingTool): readonly ToolParameterDescription[] {
  if (description === undefined) return [];
  return sortByUsage(
    description.parameters.filter(
      (param) => !param.hidden && !WINDOWMASKER_HELP_FLAGS.has(param.flag) && !WINDOWMASKER_VERSION_FLAGS.has(param.flag),
    ),
    MASKING_USAGE_ORDER[tool],
  );
}

/** The options that the command line requires (convert2blastmask's `-masking_algorithm` and `-masking_options`): a block of their own. */
export function maskingRequiredParameters(description: ToolDescription | undefined, tool: MaskingTool): readonly ToolParameterDescription[] {
  return allMaskingParameters(description, tool).filter((param) => param.mandatory && !MASKING_BLOCK_FLAGS.has(param.flag));
}

/** The options of a Masking tool shown under its collapsible "<Tool> options" section. */
export function maskingOptionParameters(description: ToolDescription | undefined, tool: MaskingTool): readonly ToolParameterDescription[] {
  return allMaskingParameters(description, tool).filter((param) => !param.mandatory && !MASKING_BLOCK_FLAGS.has(param.flag));
}

/** A block option of a Masking tool (`-infmt`, `-outfmt`), if its `-xmlhelp` lists it. */
export function maskingBlockParameter(description: ToolDescription | undefined, flag: string): ToolParameterDescription | undefined {
  return description?.parameters.find((param) => param.flag === flag);
}

/** NCBI `-help` usage order for makeblastdb's options; any unknown future option sorts after these in XML order. */
export const MAKEBLASTDB_USAGE_ORDER: readonly string[] = Object.freeze([
  '-in',
  '-input_type',
  '-dbtype',
  '-title',
  '-parse_seqids',
  '-hash_index',
  '-mask_data',
  '-mask_id',
  '-mask_desc',
  '-gi_mask',
  '-gi_mask_name',
  '-out',
  '-blastdb_version',
  '-max_file_sz',
  '-metadata_output_prefix',
  '-logfile',
  '-taxid',
  '-taxid_map',
  '-oid_masks',
]);

/**
 * The headings of makeblastdb's option groups, as NCBI's `-help` prints them (` *** Input options`): `-xmlhelp` gives
 * each option its group's number only. A group that a later version adds is shown under "Other options".
 */
export const MAKEBLASTDB_GROUPS: Readonly<Record<string, string>> = Object.freeze({
  '1': 'Input options',
  '2': 'Configuration options',
  '3': 'Sequence masking options',
  '4': 'Output options',
  '5': 'Taxonomy options',
});

/** makeblastdb's options placed in the Input and Output blocks of its tab rather than in its option groups. */
export const MAKEBLASTDB_BLOCK_FLAGS: ReadonlySet<string> = new Set(['-in', '-input_type', '-dbtype', '-out']);

/**
 * Every option of makeblastdb that its `-help` lists, in its usage order: those of a group, and the required
 * `-dbtype` (NCBI's "REQUIRED ARGUMENTS", which have no group); not the help and version flags, and not hidden
 * options (`-randomize`, which `-help-full` lists; Additional options take it).
 */
export function makeblastdbParameters(description: ToolDescription | undefined): readonly ToolParameterDescription[] {
  if (description === undefined) return [];
  return sortByUsage(
    description.parameters.filter((param) => (param.group !== undefined || param.mandatory) && !param.hidden),
    MAKEBLASTDB_USAGE_ORDER,
  );
}

/** The sections of makeblastdb's options outside the Input and Output blocks: one for each group of `-help`, in its order. */
export function makeblastdbOptionSections(
  description: ToolDescription | undefined,
): readonly { readonly title: string; readonly parameters: readonly ToolParameterDescription[] }[] {
  const sections = new Map<string, ToolParameterDescription[]>();
  for (const param of makeblastdbParameters(description)) {
    if (MAKEBLASTDB_BLOCK_FLAGS.has(param.flag)) continue;
    const title = (param.group === undefined ? undefined : MAKEBLASTDB_GROUPS[param.group]) ?? 'Other options';
    const list = sections.get(title);
    if (list === undefined) sections.set(title, [param]);
    else list.push(param);
  }
  return [...sections].map(([title, parameters]) => ({ title, parameters }));
}

/** makeblastdb's options shown in its option groups (every option of `makeblastdbOptionSections`). */
export function makeblastdbOptionParameters(description: ToolDescription | undefined): readonly ToolParameterDescription[] {
  return makeblastdbOptionSections(description).flatMap((section) => section.parameters);
}

/** The input types of makeblastdb whose `-in` names files that the run directory holds as picked (not a database). */
export type MakeblastdbInputType = 'fasta' | 'asn1_bin' | 'asn1_txt' | 'blastdb';

/** The input type that `-input_type` gives (`fasta`, its default, for an empty or unknown value). */
export function makeblastdbInputType(value: unknown): MakeblastdbInputType {
  const text = typeof value === 'string' ? value.trim() : '';
  return text === 'asn1_bin' || text === 'asn1_txt' || text === 'blastdb' ? text : 'fasta';
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
/**
 * The tool options whose text the tool writes or stores as given, with no list of values: makeblastdb's database
 * `-title` and the masking algorithms' `-mask_id` and `-mask_desc`, convert2blastmask's `-masking_algorithm` (compared
 * as given, case aside) and `-masking_options` (its `algo-options`), windowmasker's `-meta` (the counts' note). The
 * tool forms write them as typed, white space and all (`tool-argv.ts`; review of WW4c, F3); the search forms and the
 * tools' names, numbers and keywords drop the white space around a value.
 */
export const TOOL_WRITTEN_TEXT_FLAGS: ReadonlySet<string> = new Set([
  '-title',
  '-mask_id',
  '-mask_desc',
  '-masking_algorithm',
  '-masking_options',
  '-meta',
]);

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
  if (TOOL_FILE_STRING_FLAGS.has(param.flag) || param.kind === 'File_In') {
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

/**
 * Whether the argv asks the tool to print its version or its help instead of running a stage, as the command line
 * does: a version flag anywhere wins (NCBI BLAST+ 2.17.0 and LOSAT print the version for `-help -version`), else a
 * help flag (`toolHelpFlag`).
 */
export function toolPrintsKind(argv: readonly string[]): PrintedKind | null {
  const flags = argv.slice(1).map((word) => word.split('=', 1)[0]!);
  if (flags.some((flag) => WINDOWMASKER_VERSION_FLAGS.has(flag))) return 'version';
  return toolHelpFlag(argv) === undefined ? null : 'help';
}

/** The first help flag of a tool's argv (`-xmlhelp` prints the XML, the others the usage); undefined when there is none. */
export function toolHelpFlag(argv: readonly string[]): string | undefined {
  return argv
    .slice(1)
    .map((word) => word.split('=', 1)[0]!)
    .find((flag) => WINDOWMASKER_HELP_FLAGS.has(flag));
}

/** Download file name for a tool run output or file (`runFileName` pattern). */
export function toolOutputFileName(
  run: { readonly number: number; readonly tool: ToolId; readonly stage: ToolStage; readonly queuedAt: number },
  name: string,
): string {
  const safe = name.replace(/\/+$/, '').replaceAll('/', '_');
  return `${exportPrefix({ number: run.number, program: run.tool, task: run.stage, queuedAt: run.queuedAt })}.${safe}`;
}

/**
 * The name under which an output of a tool run downloads: a file that the tool wrote into its run directory keeps its
 * own name (makeblastdb's database files name each other: an alias file its volumes), else the run's start of names and
 * the output's name (`stdout.txt` and `stderr.txt` for the streams).
 */
export function toolOutputDownloadName(
  run: { readonly number: number; readonly tool: ToolId; readonly stage: ToolStage; readonly queuedAt: number },
  entry: { readonly key: string; readonly name: string },
): string {
  if (entry.key === 'stdout') return toolOutputFileName(run, 'stdout.txt');
  if (entry.key === 'stderr') return toolOutputFileName(run, 'stderr.txt');
  if (toolWritesDirectory(run.tool)) return entry.name.slice(entry.name.lastIndexOf('/') + 1);
  return toolOutputFileName(run, entry.name);
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
