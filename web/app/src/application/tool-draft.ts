// The WindowMasker form draft (WW3b DESIGN §2, §4): independent of queued and completed runs,
// driven by `ToolEngine.describeTool('windowmasker')` (-xmlhelp) and validated by `validateTool`
// through the tool drafts' shared core (`ToolDraftCore`: files of the run, check, submit). Its
// FASTA input (paste, files inflated into session OPFS, Combined or Separate runs, the `-in` name of
// a BLAST database) is the one that the Masking tab reads too (`FastaToolDraft`).
import type { DirectoryFile } from 'losat/domain';
import { splitOptions } from '../domain/argv';
import { normalizeRunPath } from '../domain/blast-db';
import { countsName } from '../domain/counts';
import type { FieldValue, FormValues } from '../domain/parameters';
import { buildWindowMaskerArgv, handedOverName, namesOption } from '../domain/tool-argv';
import {
  defaultCountsOutName,
  stageOptionParameters,
  toolPrintsKind,
  type WindowMaskerStage,
} from '../domain/tools';
import type { DataGateway } from '../ports/data';
import type { ToolDirectoryFile, ToolRequest } from './coordinator';
import { EMPTY_FASTA_INPUT, FastaToolDraft, type FastaInputState, type FastaToolDraftDeps } from './fasta-tool-draft';

import type { DatabaseOption } from './tool-draft-core';
export type { ToolExtraFile } from './tool-draft-core';
export type { ToolInputFile } from './fasta-tool-draft';

/** Counts file selected for Mask (`-ustat`) or Convert counts (`-in`). */
export type ToolCountsChoice =
  | { readonly kind: 'none' }
  | {
      readonly kind: 'session';
      readonly runId: string;
      readonly runNumber: number;
      readonly key: string;
      readonly name: string;
      readonly size: number;
    }
  | { readonly kind: 'saved'; readonly taxid: number; readonly name: string; readonly size: number }
  | { readonly kind: 'file'; readonly file: File; readonly name: string; readonly size: number };

export interface ToolDraftState extends FastaInputState {
  readonly stage: WindowMaskerStage;
  readonly values: Readonly<Record<WindowMaskerStage, FormValues>>;
  readonly extra: Readonly<Record<WindowMaskerStage, string>>;
  readonly outNames: Readonly<Record<WindowMaskerStage, string>>;
  readonly outEdited: Readonly<Record<WindowMaskerStage, boolean>>;
  readonly ustatSource: ToolCountsChoice;
  readonly convertSource: ToolCountsChoice;
}

export interface ToolDraftDeps extends FastaToolDraftDeps {
  readonly data: Pick<DataGateway, 'inflateToSessionFile' | 'dropSessionFile' | 'getToolOutputFile'>;
  /** Reads a counts file saved in this browser for a taxonomy ID (OPFS-backed File). */
  readonly readSavedCounts?: (taxid: number) => Promise<File | undefined>;
}

const STAGES: readonly WindowMaskerStage[] = ['mk_counts', 'mask', 'convert'];

export class ToolDraft extends FastaToolDraft<ToolDraftState> {
  constructor(private readonly deps: ToolDraftDeps) {
    const emptyValues = Object.fromEntries(STAGES.map((stage) => [stage, Object.freeze({})])) as Record<WindowMaskerStage, FormValues>;
    const emptyExtra = Object.fromEntries(STAGES.map((stage) => [stage, ''])) as Record<WindowMaskerStage, string>;
    super('windowmasker', deps, {
      stage: 'mk_counts',
      values: emptyValues,
      extra: emptyExtra,
      outNames: { mk_counts: defaultCountsOutName(undefined), mask: '', convert: defaultCountsOutName(undefined) },
      outEdited: { mk_counts: false, mask: false, convert: false },
      ...EMPTY_FASTA_INPUT,
      ustatSource: { kind: 'none' },
      convertSource: { kind: 'none' },
      runFiles: [],
      validation: { state: 'idle' },
      submitting: false,
    });
  }

  // --- stage and options ------------------------------------------------------------------

  setStage(stage: WindowMaskerStage): void {
    if (this.state.get().stage === stage) return;
    this.set({ stage, message: undefined, runFilesError: undefined });
    this.scheduleValidation();
  }

  private patchStageValues(
    values: ToolDraftState['values'],
    stage: WindowMaskerStage,
    mutate: (next: Record<string, FieldValue>) => void,
  ): ToolDraftState['values'] {
    const next: Record<string, FieldValue> = { ...values[stage] };
    mutate(next);
    return { ...values, [stage]: Object.freeze(next) };
  }

  setField(flag: string, value: FieldValue): void {
    const state = this.state.get();
    const { stage, values, outNames, outEdited } = state;
    const nextValues = this.patchStageValues(values, stage, (next) => {
      if (value === undefined || value === '' || value === false) delete next[flag];
      else next[flag] = value;
    });

    let nextOutNames = outNames;
    if (flag === '-sformat' && (stage === 'mk_counts' || stage === 'convert') && !outEdited[stage]) {
      nextOutNames = {
        ...outNames,
        [stage]: defaultCountsOutName(typeof value === 'string' ? value : undefined),
      };
    }
    this.set({
      values: nextValues,
      outNames: nextOutNames,
      message: undefined,
    });
    this.scheduleValidation();
  }

  setOutName(text: string): void {
    const { stage, outNames, outEdited } = this.state.get();
    this.set({
      outNames: { ...outNames, [stage]: text },
      outEdited: { ...outEdited, [stage]: true },
      message: undefined,
    });
    this.scheduleValidation();
  }

  setExtra(text: string): void {
    const { stage, extra } = this.state.get();
    this.set({ extra: { ...extra, [stage]: text }, message: undefined });
    this.scheduleValidation();
  }

  /** Resets options inside the collapsible "WindowMasker options" section and Additional options for the active stage. */
  resetOptions(): void {
    const state = this.state.get();
    const { stage, values, extra, description, runFiles } = state;
    const optionFlags = new Set(stageOptionParameters(description, stage).map((param) => param.flag));
    const clearedValues = this.patchStageValues(values, stage, (next) => {
      for (const flag of optionFlags) delete next[flag];
    });
    const nextRunFiles = runFiles.filter((entry) => entry.flag === undefined || !optionFlags.has(entry.flag));
    this.set({
      values: clearedValues,
      extra: { ...extra, [stage]: '' },
      runFiles: Object.freeze(nextRunFiles),
      runFilesError: undefined,
      message: undefined,
    });
    this.scheduleValidation();
  }

  // --- counts and run directory files ------------------------------------------------------

  setUstatSession(runId: string, runNumber: number, key: string, name: string, size: number): void {
    this.set({
      ustatSource: sessionCountsChoice(runId, runNumber, key, name, size),
      message: undefined,
    });
    this.scheduleValidation();
  }

  setUstatSaved(taxid: number, size: number): void {
    this.set({
      ustatSource: savedCountsChoice(taxid, size),
      message: undefined,
    });
    this.scheduleValidation();
  }

  setUstatFile(file: File | undefined): void {
    this.set({
      ustatSource: fileCountsChoice(file),
      message: undefined,
    });
    this.scheduleValidation();
  }

  clearUstat(): void {
    this.set({ ustatSource: { kind: 'none' }, message: undefined });
    this.scheduleValidation();
  }

  setConvertSession(runId: string, runNumber: number, key: string, name: string, size: number): void {
    this.set({
      convertSource: sessionCountsChoice(runId, runNumber, key, name, size),
      message: undefined,
    });
    this.scheduleValidation();
  }

  setConvertSaved(taxid: number, size: number): void {
    this.set({
      convertSource: savedCountsChoice(taxid, size),
      message: undefined,
    });
    this.scheduleValidation();
  }

  setConvertFile(file: File | undefined): void {
    this.set({
      convertSource: fileCountsChoice(file),
      message: undefined,
    });
    this.scheduleValidation();
  }

  clearConvert(): void {
    this.set({ convertSource: { kind: 'none' }, message: undefined });
    this.scheduleValidation();
  }

  /** Switches to Mask and preselects a counts file produced by an earlier run in this session. */
  useCountsForMask(runId: string, runNumber: number, key: string, name: string, size: number): void {
    this.ensureLoaded();
    const ustatSource = sessionCountsChoice(runId, runNumber, key, name, size);
    this.set({ stage: 'mask', ustatSource, message: undefined });
    this.scheduleValidation();
  }

  // --- command line & validation ----------------------------------------------------------

  /** Whether the active stage currently reads `-in` as FASTA files/paste rather than a named blastdb/seqids/fa_list. */
  usesFastaBox(state: ToolDraftState = this.state.get()): boolean {
    if (state.stage === 'convert') return false;
    const values = state.values[state.stage];
    const infmt = typeof values['-infmt'] === 'string' ? values['-infmt'].trim() : 'fasta';
    const faList = values['-fa_list'];
    const isFaList = faList === true || (typeof faList === 'string' && /^(?:true|t|1)$/i.test(faList.trim()));
    return (infmt === '' || infmt === 'fasta') && !isFaList;
  }

  /** Preview argv built by the single command-line builder (`buildWindowMaskerArgv`). */
  previewArgv(state: ToolDraftState = this.state.get()): readonly string[] {
    let extraWords: readonly string[] = [];
    try {
      extraWords = splitOptions(state.extra[state.stage]);
    } catch {
      // Unclosed quote in Additional options: omit broken extra words from preview; validation reports it.
      extraWords = [];
    }
    const inputName = this.previewInputName(state);
    const ustatName = state.stage === 'mask' && state.ustatSource.kind !== 'none' ? state.ustatSource.name : undefined;
    try {
      return this.buildArgv(state, extraWords, { inputName, ustatName });
    } catch {
      return ['windowmasker'];
    }
  }

  /** Whether the WindowMasker tab holds user inputs or files that a reload would lose. */
  holdsInputs(state: ToolDraftState = this.state.get()): boolean {
    return this.holdsFastaInputs(state) || state.ustatSource.kind !== 'none' || state.convertSource.kind !== 'none';
  }

  // --- what the core asks -----------------------------------------------------------------

  protected parseExtra(state: ToolDraftState): readonly string[] {
    return splitOptions(state.extra[state.stage]);
  }

  protected optionValue(state: ToolDraftState, flag: string): FieldValue | undefined {
    return state.values[state.stage][flag];
  }

  /**
   * The argv of the check: the inputs not chosen yet are named `input.fa` and `counts.ascii` (`-in` none where the words
   * name it themselves). A database of the files of the run that `-in` cannot name (`databaseRefusal`) is the check's
   * refusal, unless the words print help or the version.
   */
  protected validationArgv(state: ToolDraftState, extra: readonly string[]): readonly string[] {
    const placeholder = namesOption(extra, '-in') ? undefined : state.stage === 'convert' ? 'counts.ascii' : 'input.fa';
    const inputName = this.previewInputName(state) ?? placeholder;
    const ustatName = state.stage === 'mask' ? (state.ustatSource.kind === 'none' ? 'counts.ascii' : state.ustatSource.name) : undefined;
    const argv = this.buildArgv(state, extra, { inputName, ustatName });
    const refusal = this.databaseRefusal(state);
    if (refusal !== undefined && toolPrintsKind(argv) === null) throw new Error(refusal);
    return argv;
  }

  /**
   * WindowMasker's `-in` with the stage's `-infmt` (`isSeqdbList`), which reads a nucleotide database; Convert's `-in`
   * is its counts file.
   */
  protected override inDatabase(state: ToolDraftState): DatabaseOption | undefined {
    if (state.stage === 'convert') return undefined;
    const infmt = state.values[state.stage]['-infmt'];
    return {
      option: { command: 'windowmasker', flag: '-in', inputFormat: typeof infmt === 'string' ? infmt : 'fasta' },
      kind: 'nucleotide',
    };
  }

  /**
   * The check's directory holds the inputs that the run reads, under their names (review of WW4c, F2): the FASTA
   * input and, for Mask, the `-ustat` counts; for Convert the counts. A counts file that the form reads only when the
   * run starts (saved or made in the session) stands in by its name, as the check does not read the counts.
   */
  protected override validationInputs(state: ToolDraftState): readonly DirectoryFile[] {
    if (toolPrintsKind(this.buildArgv(state, this.parseExtra(state))) !== null) return [];
    const counts = (choice: ToolCountsChoice): readonly DirectoryFile[] =>
      choice.kind === 'none' ? [] : [{ name: choice.name, file: choice.kind === 'file' ? choice.file : new Blob([]) }];
    if (state.stage === 'convert') return counts(state.convertSource);
    return [...this.fastaCheckInputs(state), ...(state.stage === 'mask' ? counts(state.ustatSource) : [])];
  }

  protected async buildRequests(state: ToolDraftState, extra: readonly string[]): Promise<readonly ToolRequest[]> {
    const prints = toolPrintsKind(this.buildArgv(state, extra));
    return prints !== null ? [this.buildPrintRequest(state, extra, prints)] : this.buildStageRequests(state, extra);
  }

  // --- internals --------------------------------------------------------------------------

  private buildArgv(
    state: ToolDraftState,
    extra: readonly string[],
    names: { readonly inputName?: string | undefined; readonly ustatName?: string | undefined } = {},
  ): readonly string[] {
    const { stage, outNames, values, description } = state;
    return buildWindowMaskerArgv({
      stage,
      ...(names.inputName === undefined ? {} : { inputName: names.inputName }),
      ...(names.ustatName === undefined ? {} : { ustatName: names.ustatName }),
      outName: outNames[stage],
      values: values[stage],
      ...(description === undefined ? {} : { description }),
      extra,
    });
  }

  private buildPrintRequest(state: ToolDraftState, extraWords: readonly string[], prints: 'help' | 'version'): ToolRequest {
    const argv = this.buildArgv(state, extraWords);
    return { tool: 'windowmasker', stage: state.stage, argv, directory: [], prints };
  }

  private req(stage: WindowMaskerStage, argv: readonly string[], directory: readonly ToolDirectoryFile[]): ToolRequest[] {
    return [{ tool: 'windowmasker', stage, argv, directory }];
  }

  private async buildStageRequests(state: ToolDraftState, extraWords: readonly string[]): Promise<ToolRequest[]> {
    const { stage } = state;

    if (stage === 'convert') {
      const counts = await this.resolveCountsChoice(state.convertSource, 'Choose a counts file to convert.');
      const dir = this.runDirectoryOf(state, [counts]);
      const argv = this.buildArgv(state, extraWords, { inputName: counts.name });
      return this.req(stage, argv, dir);
    }

    const ustatFile =
      stage === 'mask' ? await this.resolveCountsChoice(state.ustatSource, 'Choose a counts file for -ustat.') : undefined;
    const inputs = await this.runInputs(state);
    return inputs.map((input): ToolRequest => {
      const directory = this.runDirectoryOf(state, [
        ...(input.file === undefined ? [] : [input.file]),
        ...(ustatFile === undefined ? [] : [ustatFile]),
      ]);
      const argv = this.buildArgv(state, extraWords, { inputName: input.name, ustatName: ustatFile?.name });
      return {
        tool: 'windowmasker',
        stage,
        argv,
        directory,
        ...(input.sessionTokens.length === 0 ? {} : { sessionTokens: Object.freeze([...input.sessionTokens]) }),
      };
    });
  }

  private async resolveCountsChoice(choice: ToolCountsChoice, emptyMessage: string): Promise<ToolDirectoryFile> {
    if (choice.kind === 'none') throw new Error(emptyMessage);
    if (choice.kind === 'file') {
      return { name: choice.name, file: choice.file, origin: { kind: 'picked' } };
    }
    if (choice.kind === 'saved') {
      const file = await this.deps.readSavedCounts?.(choice.taxid);
      if (file === undefined) {
        throw new Error(`The saved counts for taxonomy ID ${choice.taxid} are no longer available.`);
      }
      return { name: choice.name, file, origin: { kind: 'saved', taxid: choice.taxid } };
    }
    const file = await this.deps.data.getToolOutputFile(choice.runId, choice.key, choice.name);
    return { name: choice.name, file, origin: { kind: 'made', runNumber: choice.runNumber } };
  }

  private previewInputName(state: ToolDraftState): string | undefined {
    if (state.stage === 'convert') {
      return state.convertSource.kind === 'none' ? undefined : state.convertSource.name;
    }
    return this.fastaInputName(state);
  }
}

/** Counts made by a run of this session, under their name in this run's directory (`handedOverName`: `../c.counts` is `c.counts`). */
function sessionCountsChoice(runId: string, runNumber: number, key: string, name: string, size: number): ToolCountsChoice {
  return { kind: 'session', runId, runNumber, key, name: handedOverName(name), size };
}

function savedCountsChoice(taxid: number, size: number): ToolCountsChoice {
  return { kind: 'saved', taxid, name: countsName(taxid), size };
}

function fileCountsChoice(file: File | undefined): ToolCountsChoice {
  return file === undefined
    ? { kind: 'none' }
    : { kind: 'file', file, name: normalizeRunPath(file.name) ?? file.name, size: file.size };
}
