// The Masking form draft (WW4c): one draft for the Masking tab's three tools, dustmasker, segmasker
// and convert2blastmask, on the tool drafts' shared core (`ToolDraftCore`: `-xmlhelp`, files of the
// run, check, submit) and with WindowMasker's FASTA input (`FastaToolDraft`: paste, files inflated into
// session OPFS, Combined or Separate runs, the `-in` name of a BLAST database). The input is shared by
// the tools; each tool keeps its own option values, Additional options and `-out` name (empty: no
// `-out` word, the standard output, as on the command line). Each tool's `-xmlhelp` is read when the
// tool is first shown.
import type { DirectoryFile, PrintedKind } from 'losat/domain';
import { splitOptions } from '../domain/argv';
import type { FieldValue, FormValues } from '../domain/parameters';
import { buildMaskingArgv, namesOption, runFileNames } from '../domain/tool-argv';
import { MASKING_TOOLS, maskingOptionParameters, toolPrintsKind, type MaskingTool, type ToolId } from '../domain/tools';
import type { ToolRequest } from './coordinator';
import { EMPTY_FASTA_INPUT, FastaToolDraft, type FastaInputState, type FastaToolDraftDeps } from './fasta-tool-draft';
import type { DatabaseOption } from './tool-draft-core';

export interface MaskingDraftState extends FastaInputState {
  /** The tool chosen in the tool switcher. */
  readonly tool: MaskingTool;
  readonly values: Readonly<Record<MaskingTool, FormValues>>;
  readonly extra: Readonly<Record<MaskingTool, string>>;
  /** Each tool's `-out` name: empty (or `-`) writes no `-out` word (the standard output). */
  readonly outNames: Readonly<Record<MaskingTool, string>>;
}

export type MaskingDraftDeps = FastaToolDraftDeps;

/** The stage of every run of the Masking tab (`ToolStage`): each tool has one. */
const STAGE = 'mask';

function perTool<T>(value: T): Readonly<Record<MaskingTool, T>> {
  return Object.freeze(Object.fromEntries(MASKING_TOOLS.map((entry) => [entry.id, value])) as Record<MaskingTool, T>);
}

export class MaskingDraft extends FastaToolDraft<MaskingDraftState> {
  constructor(deps: MaskingDraftDeps) {
    super('dustmasker', deps, {
      tool: 'dustmasker',
      values: perTool<FormValues>(Object.freeze({})),
      extra: perTool(''),
      outNames: perTool(''),
      ...EMPTY_FASTA_INPUT,
      runFiles: [],
      validation: { state: 'idle' },
      submitting: false,
    });
  }

  // --- the tool and its options ---------------------------------------------------------------

  /** Shows another tool: its `-xmlhelp` (read once), its options, Additional options and `-out` name. */
  setTool(tool: MaskingTool): void {
    if (this.state.get().tool === tool) return;
    this.set({ tool, message: undefined, runFilesError: undefined });
    this.showDescription();
    this.scheduleValidation();
  }

  protected optionValue(state: MaskingDraftState, flag: string): FieldValue | undefined {
    return state.values[state.tool][flag];
  }

  setField(flag: string, value: FieldValue): void {
    const { tool, values } = this.state.get();
    const next: Record<string, FieldValue> = { ...values[tool] };
    if (value === undefined || value === '' || value === false) delete next[flag];
    else next[flag] = value;
    this.set({ values: { ...values, [tool]: Object.freeze(next) }, message: undefined });
    this.scheduleValidation();
  }

  setOutName(text: string): void {
    const { tool, outNames } = this.state.get();
    this.set({ outNames: { ...outNames, [tool]: text }, message: undefined });
    this.scheduleValidation();
  }

  setExtra(text: string): void {
    const { tool, extra } = this.state.get();
    this.set({ extra: { ...extra, [tool]: text }, message: undefined });
    this.scheduleValidation();
  }

  /** Resets the options of the "<Tool> options" section and Additional options of the tool shown. */
  resetOptions(): void {
    const { tool, values, extra, description } = this.state.get();
    const next: Record<string, FieldValue> = { ...values[tool] };
    for (const param of maskingOptionParameters(description, tool)) delete next[param.flag];
    this.set({
      values: { ...values, [tool]: Object.freeze(next) },
      extra: { ...extra, [tool]: '' },
      message: undefined,
    });
    this.scheduleValidation();
  }

  // --- the command line -----------------------------------------------------------------------

  /** Whether the tool reads `-in` from the FASTA box: always for convert2blastmask, else unless `-infmt` names a BLAST database. */
  usesFastaBox(state: MaskingDraftState = this.state.get()): boolean {
    const infmt = state.values[state.tool]['-infmt'];
    const text = typeof infmt === 'string' ? infmt.trim() : '';
    return text === '' || text === 'fasta';
  }

  /** The argv shown beside Run LOSAT; a quote left open in Additional options leaves their words out (the check reports it). */
  previewArgv(state: MaskingDraftState = this.state.get()): readonly string[] {
    let extra: readonly string[];
    try {
      extra = this.parseExtra(state);
    } catch {
      extra = [];
    }
    const inputName = this.inputNameOf(state, extra);
    try {
      return this.argvOf(state, extra, inputName);
    } catch {
      return [state.tool];
    }
  }

  holdsInputs(state: MaskingDraftState = this.state.get()): boolean {
    return this.holdsFastaInputs(state);
  }

  // --- what the core asks -------------------------------------------------------------------

  protected override shownTool(state: MaskingDraftState): ToolId {
    return state.tool;
  }

  protected parseExtra(state: MaskingDraftState): readonly string[] {
    return splitOptions(state.extra[state.tool]);
  }

  /**
   * The argv of the check: an input not given yet is named `input.fa` (none where the words name `-in` themselves). A
   * database of the files of the run that `-in` cannot name (`databaseRefusal`) is the check's refusal, unless the words
   * print help or the version.
   */
  protected validationArgv(state: MaskingDraftState, extra: readonly string[]): readonly string[] {
    const prints = this.printsKindOf(state, extra) !== null;
    const refusal = prints ? undefined : this.databaseRefusal(state);
    if (refusal !== undefined) throw new Error(refusal);
    const placeholder = prints || namesOption(extra, '-in') ? undefined : 'input.fa';
    return this.argvOf(state, extra, this.inputNameOf(state, extra) ?? placeholder);
  }

  /**
   * The tool's `-in` with its `-infmt` (`isSeqdbList`: dustmasker's and segmasker's with `-infmt blastdb`), and the kind
   * of the database that it reads: protein for segmasker, nucleotide for dustmasker.
   */
  protected override inDatabase(state: MaskingDraftState): DatabaseOption {
    const infmt = state.values[state.tool]['-infmt'];
    return {
      option: { command: state.tool, flag: '-in', inputFormat: typeof infmt === 'string' ? infmt : 'fasta' },
      kind: state.tool === 'segmasker' ? 'protein' : 'nucleotide',
    };
  }

  /** The check's directory holds the input that the run reads, under its name (none for help and the version). */
  protected override validationInputs(state: MaskingDraftState): readonly DirectoryFile[] {
    return this.printsKindOf(state, this.parseExtra(state)) === null ? this.fastaCheckInputs(state) : [];
  }

  /** A help or version word of Additional options (`toolPrintsKind`): the words that the form does not build. */
  protected printsKindOf(state: MaskingDraftState, extra: readonly string[]): PrintedKind | null {
    return toolPrintsKind([state.tool, ...extra]);
  }

  protected async buildRequests(state: MaskingDraftState, extra: readonly string[]): Promise<readonly ToolRequest[]> {
    const prints = this.printsKindOf(state, extra);
    if (prints !== null) {
      // Help and version read no input; the files of the run are there, as in the directory where the command line runs.
      return [{ tool: state.tool, stage: STAGE, argv: this.argvOf(state, extra, undefined), directory: this.runDirectoryOf(state, []), prints }];
    }
    const inputs = await this.runInputs(state);
    return inputs.map(
      (input): ToolRequest => ({
        tool: state.tool,
        stage: STAGE,
        argv: this.argvOf(state, extra, input.name),
        directory: this.runDirectoryOf(state, input.file === undefined ? [] : [input.file]),
        ...(input.sessionTokens.length === 0 ? {} : { sessionTokens: Object.freeze([...input.sessionTokens]) }),
      }),
    );
  }

  // --- internals ----------------------------------------------------------------------------

  /** The `-in` name of the argv shown: none for a run that prints help or the version. */
  private inputNameOf(state: MaskingDraftState, extra: readonly string[]): string | undefined {
    return this.printsKindOf(state, extra) === null ? this.fastaInputName(state) : undefined;
  }

  private argvOf(state: MaskingDraftState, extra: readonly string[], inputName: string | undefined): readonly string[] {
    const { tool, description } = state;
    return buildMaskingArgv({
      tool,
      ...(inputName === undefined ? {} : { inputName }),
      outName: state.outNames[tool],
      values: state.values[tool],
      ...(description === undefined ? {} : { description }),
      extra,
      fileNames: runFileNames(state.runFiles),
    });
  }
}
