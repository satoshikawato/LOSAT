// The WindowMasker form draft (WW3b DESIGN §2, §4): independent of queued and completed runs,
// driven by `ToolEngine.describeTool('windowmasker')` (-xmlhelp) and validated by `validateTool`.
// Large inputs stay as disk/OPFS-backed Files: compressed inputs inflate into session OPFS via
// `DataGateway.inflateToSessionFile`, and combined FASTA files join Blob references by checking
// only the final byte of each non-final file so page JS heap never holds the input sequences.
import type { DirectoryFile, ToolDescription } from 'losat/domain';
import { isUnavailable, type ToolEngine, type ToolValidation } from 'losat/ports';
import { splitOptions } from '../domain/argv';
import { countsName } from '../domain/counts';
import type { FieldValue, FormValues } from '../domain/parameters';
import { buildWindowMaskerArgv, toolCommandWords, toolValidationSummary } from '../domain/tool-argv';
import {
  defaultCountsOutName,
  detectBlastDbName,
  stageOptionParameters,
  toolPrintsKind,
  type WindowMaskerStage,
} from '../domain/tools';
import type { DataGateway } from '../ports/data';
import { failureText } from './availability';
import type { EnqueueAllResult, ToolDirectoryFile, ToolRequest } from './coordinator';
import type { InputMode, Validation } from './draft';
import { Store } from './store';

/** One FASTA file added to WindowMasker's Input block (inflated into OPFS when compressed). */
export interface ToolInputFile {
  readonly key: string;
  readonly name: string;
  readonly size: number;
  readonly file: File;
  readonly status: 'inflating' | 'ready' | 'failed';
  readonly inflatedFrom?: string;
  readonly sessionToken?: string;
  readonly error?: string;
  readonly unavailable?: true;
}

/** One extra file added to "Files of the run" (for `-ids`, `-exclude_ids`, `-conffile`, `-fa_list`, `-infmt blastdb`). */
export interface ToolExtraFile {
  readonly key: string;
  readonly name: string;
  readonly size: number;
  readonly file: File;
  /** Option flag that picked this file via `Choose file…` (`-ids`, `-exclude_ids`), if any. */
  readonly flag?: string;
}

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

export interface ToolDraftState {
  readonly stage: WindowMaskerStage;
  readonly values: Readonly<Record<WindowMaskerStage, FormValues>>;
  readonly extra: Readonly<Record<WindowMaskerStage, string>>;
  readonly outNames: Readonly<Record<WindowMaskerStage, string>>;
  readonly outEdited: Readonly<Record<WindowMaskerStage, boolean>>;
  readonly fastaPaste: string;
  readonly fastaFiles: readonly ToolInputFile[];
  readonly fastaMode: InputMode;
  readonly blastDbIn: string;
  readonly blastDbInEdited: boolean;
  readonly ustatSource: ToolCountsChoice;
  readonly convertSource: ToolCountsChoice;
  readonly runFiles: readonly ToolExtraFile[];
  readonly runFilesError?: string;
  readonly description?: ToolDescription;
  readonly descriptionError?: string;
  readonly descriptionUnavailable?: true;
  readonly validation: Validation;
  readonly submitting: boolean;
  readonly message?: { readonly kind: 'info' | 'error'; readonly text: string };
}

export interface ToolDraftDeps {
  readonly engine: Pick<ToolEngine, 'describeTool' | 'validateTool'>;
  readonly data: Pick<DataGateway, 'inflateToSessionFile' | 'dropSessionFile' | 'getToolOutputFile'>;
  readonly enqueueToolAll: (requests: readonly ToolRequest[]) => Promise<EnqueueAllResult>;
  /** Reads a counts file saved in this browser for a taxonomy ID (OPFS-backed File). */
  readonly readSavedCounts?: (taxid: number) => Promise<File | undefined>;
  /** True while the Engine worker is running another job (so `validateTool` would wait behind it). */
  readonly isBusy?: () => boolean;
  readonly debounceMs?: number;
}

const STAGES: readonly WindowMaskerStage[] = ['mk_counts', 'mask', 'convert'];
const PASTE_NAME = 'input.fa';
const COMBINED_NAME = 'combined_input.fa';
const LF = 0x0a;

export class ToolDraft {
  readonly state: Store<ToolDraftState>;
  private descriptionRequest: Promise<ToolDescription> | undefined;
  private validationTimer: ReturnType<typeof setTimeout> | undefined;
  private validationGeneration = 0;
  private validationTurn: Promise<void> = Promise.resolve();
  private validatedKey: string | undefined;
  private nextKey = 1;
  private readonly asking = new Set<Promise<unknown>>();

  constructor(private readonly deps: ToolDraftDeps) {
    const emptyValues = Object.fromEntries(STAGES.map((stage) => [stage, Object.freeze({})])) as Record<WindowMaskerStage, FormValues>;
    const emptyExtra = Object.fromEntries(STAGES.map((stage) => [stage, ''])) as Record<WindowMaskerStage, string>;
    this.state = new Store<ToolDraftState>({
      stage: 'mk_counts',
      values: emptyValues,
      extra: emptyExtra,
      outNames: { mk_counts: defaultCountsOutName(undefined), mask: '', convert: defaultCountsOutName(undefined) },
      outEdited: { mk_counts: false, mask: false, convert: false },
      fastaPaste: '',
      fastaFiles: [],
      fastaMode: 'combined',
      blastDbIn: '',
      blastDbInEdited: false,
      ustatSource: { kind: 'none' },
      convertSource: { kind: 'none' },
      runFiles: [],
      validation: { state: 'idle' },
      submitting: false,
    });
  }

  /** Loads `-xmlhelp` once when the WindowMasker tab is first opened. */
  ensureLoaded(): void {
    const state = this.state.get();
    if (state.description === undefined && this.descriptionRequest === undefined) {
      void this.loadDescription();
    }
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

  /** Attaches a file chosen for a file-typed option (`-ids`, `-exclude_ids`) and puts it in Files of the run. */
  setOptionFile(flag: string, file: File | undefined): void {
    const state = this.state.get();
    const kept = state.runFiles.filter((entry) => entry.flag !== flag);
    if (file === undefined) {
      this.set({ runFiles: Object.freeze(kept), runFilesError: undefined });
      this.setField(flag, '');
      return;
    }
    const name = normalizeEntryName(file.name);
    if (kept.some((entry) => entry.name === name)) {
      this.set({ runFilesError: `A file named "${name}" is already in the run directory.` });
      return;
    }
    const entry: ToolExtraFile = { key: `rf-${this.nextKey++}`, name, size: file.size, file, flag };
    this.set({ runFiles: Object.freeze([...kept, entry]), runFilesError: undefined });
    this.setField(flag, name);
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

  // --- FASTA, blastdb, counts and run directory files -------------------------------------

  setFastaPaste(text: string): void {
    this.set({ fastaPaste: text, message: undefined });
    this.scheduleValidation();
  }

  setFastaMode(mode: InputMode): void {
    this.set({ fastaMode: mode, message: undefined });
    this.scheduleValidation();
  }

  addFastaFiles(files: readonly File[]): void {
    if (files.length === 0) return;
    const existing = new Set(this.state.get().fastaFiles.map((entry) => entry.name));
    const added: ToolInputFile[] = [];
    for (const file of files) {
      const entryName = normalizeEntryName(file.name);
      if (existing.has(entryName)) {
        this.reportError(`A file named "${entryName}" is already added.`);
        continue;
      }
      existing.add(entryName);
      const key = `fa-${this.nextKey++}`;
      added.push({ key, name: entryName, size: file.size, file, status: 'inflating' });
    }
    if (added.length === 0) return;
    this.set({
      fastaFiles: Object.freeze([...this.state.get().fastaFiles, ...added]),
      message: undefined,
    });
    for (const entry of added) void this.inflateFastaEntry(entry.key, entry.file);
  }

  removeFastaFile(key: string): void {
    const target = this.state.get().fastaFiles.find((entry) => entry.key === key);
    if (target === undefined) return;
    this.releaseSessionToken(target.sessionToken);
    this.set({
      fastaFiles: Object.freeze(this.state.get().fastaFiles.filter((entry) => entry.key !== key)),
      message: undefined,
    });
    this.scheduleValidation();
  }

  setBlastDbIn(name: string): void {
    this.set({ blastDbIn: name, blastDbInEdited: true, message: undefined });
    this.scheduleValidation();
  }

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

  /** Adds files (or a folder picked with `webkitdirectory`) to "Files of the run"; refuses duplicate names. */
  addRunFiles(files: readonly File[], relativeNames?: readonly string[]): void {
    if (files.length === 0) return;
    const state = this.state.get();
    const existing = new Set(state.runFiles.map((entry) => entry.name));
    const added: ToolExtraFile[] = [];
    let duplicate: string | undefined;
    for (let i = 0; i < files.length; i++) {
      const file = files[i]!;
      const rawName = relativeNames?.[i] || (file as File & { webkitRelativePath?: string }).webkitRelativePath || file.name;
      const name = normalizeEntryName(rawName);
      if (name === '') continue;
      if (existing.has(name)) {
        duplicate ??= name;
        continue;
      }
      existing.add(name);
      added.push({ key: `rf-${this.nextKey++}`, name, size: file.size, file });
    }
    const nextRunFiles = Object.freeze([...state.runFiles, ...added]);
    let nextBlastDbIn = state.blastDbIn;
    if (!state.blastDbInEdited || state.blastDbIn.trim() === '') {
      const detected = detectBlastDbName(nextRunFiles.map((entry) => entry.name));
      if (detected !== undefined) nextBlastDbIn = detected;
    }
    this.set({
      runFiles: nextRunFiles,
      blastDbIn: nextBlastDbIn,
      runFilesError: duplicate === undefined ? undefined : `A file named "${duplicate}" is already in the run directory.`,
      message: undefined,
    });
    this.scheduleValidation();
  }

  removeRunFile(key: string): void {
    const state = this.state.get();
    const target = state.runFiles.find((entry) => entry.key === key);
    if (target === undefined) return;
    const nextRunFiles = Object.freeze(state.runFiles.filter((entry) => entry.key !== key));
    this.set({ runFiles: nextRunFiles, runFilesError: undefined, message: undefined });
    if (target.flag !== undefined && state.values[state.stage][target.flag] === target.name) {
      this.setField(target.flag, '');
    } else {
      this.scheduleValidation();
    }
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

  previewCommandWords(state: ToolDraftState = this.state.get()): readonly string[] {
    return toolCommandWords(this.previewArgv(state));
  }

  /** Asks for `-xmlhelp` again and re-runs validation when LOSAT became unavailable (Availability.recheck). */
  async askAgain(): Promise<void> {
    const state = this.state.get();
    const tasks: Array<Promise<unknown>> = [...this.asking];
    if (state.descriptionUnavailable === true) {
      const task = this.loadDescription(true);
      this.trackAsking(task);
      tasks.push(task);
    }
    if (state.validation.state === 'unavailable') {
      const task = this.runValidationNow();
      this.trackAsking(task);
      tasks.push(task);
    }
    await Promise.allSettled(tasks);
  }

  /** Submits the current WindowMasker job(s) to the shared Coordinator queue. */
  async submit(): Promise<EnqueueAllResult> {
    const state = this.state.get();
    if (state.submitting) return { ok: false, message: 'A submission is already in progress.' };
    this.set({ submitting: true, message: undefined });
    try {
      let extraWords: readonly string[];
      try {
        extraWords = this.parseExtra(state);
      } catch (error) {
        const message = errorMessage(error);
        this.reportError(message);
        return { ok: false, message };
      }

      const preview = this.buildArgv(state, extraWords);
      const prints = toolPrintsKind(preview);

      const busy = this.deps.isBusy?.() === true;
      if (!busy) {
        if (this.validationTimer !== undefined) {
          clearTimeout(this.validationTimer);
          this.validationTimer = undefined;
          await this.runValidationNow();
        } else {
          await this.validationTurn;
        }
      }
      const currentKey = this.validationKey(this.state.get());
      const validated = this.state.get().validation;
      if (validated.state === 'invalid' && (!busy || this.validatedKey === currentKey)) {
        const text = toolValidationSummary(validated.message);
        this.reportError(text);
        return { ok: false, message: text };
      }

      const requests = prints !== null
        ? [this.buildPrintRequest(state, extraWords, prints)]
        : await this.buildStageRequests(state, extraWords);

      const result = await this.deps.enqueueToolAll(requests);
      if (!result.ok) {
        const text = toolValidationSummary(result.message);
        this.set({
          validation: { state: 'invalid', message: result.message },
          message: { kind: 'error', text },
        });
        return { ok: false, message: text };
      }
      return result;
    } catch (error) {
      const message = failureText(error);
      this.reportError(message);
      return { ok: false, message };
    } finally {
      this.set({ submitting: false });
    }
  }

  // --- internals --------------------------------------------------------------------------

  private parseExtra(state: ToolDraftState): readonly string[] {
    return splitOptions(state.extra[state.stage]);
  }

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

  private reportError(text: string): void {
    this.set({ message: { kind: 'error', text } });
  }

  private markInvalid(error: unknown): void {
    this.set({ validation: { state: 'invalid', message: errorMessage(error) } });
  }

  private releaseSessionToken(token: string | undefined): void {
    if (token !== undefined) void this.deps.data.dropSessionFile(token).catch(() => undefined);
  }

  private buildPrintRequest(state: ToolDraftState, extraWords: readonly string[], prints: 'help' | 'version'): ToolRequest {
    const argv = this.buildArgv(state, extraWords);
    return { tool: 'windowmasker', stage: state.stage, argv, directory: [], prints };
  }

  private req(stage: WindowMaskerStage, argv: readonly string[], directory: readonly ToolDirectoryFile[]): ToolRequest[] {
    return [{ tool: 'windowmasker', stage, argv, directory }];
  }

  private async buildStageRequests(state: ToolDraftState, extraWords: readonly string[]): Promise<ToolRequest[]> {
    const { stage, runFiles } = state;
    const baseExtraDir: ToolDirectoryFile[] = runFiles.map((entry) => ({
      name: entry.name,
      file: entry.file,
      origin: { kind: 'picked' },
    }));

    if (stage === 'convert') {
      const counts = await this.resolveCountsChoice(state.convertSource, 'Choose a counts file to convert.');
      const dir = assembleDirectory([counts, ...baseExtraDir]);
      const argv = this.buildArgv(state, extraWords, { inputName: counts.name });
      return this.req(stage, argv, dir);
    }

    const ustatFile =
      stage === 'mask' ? await this.resolveCountsChoice(state.ustatSource, 'Choose a counts file for -ustat.') : undefined;

    if (!this.usesFastaBox(state)) {
      const inputName = state.blastDbIn.trim();
      if (inputName === '') throw new Error('Enter an input name or open files.');
      const dir = assembleDirectory([...(ustatFile === undefined ? [] : [ustatFile]), ...baseExtraDir]);
      const argv = this.buildArgv(state, extraWords, { inputName, ustatName: ustatFile?.name });
      return this.req(stage, argv, dir);
    }

    if (state.fastaFiles.some((entry) => entry.status === 'inflating')) {
      throw new Error('Wait until the opened files finish reading.');
    }
    const failed = state.fastaFiles.find((entry) => entry.status === 'failed');
    if (failed !== undefined) {
      throw new Error(`${failed.name}: ${failed.error ?? 'could not be read'}`);
    }

    const parts: ToolDirectoryFile[] = [];
    if (state.fastaPaste.trim() !== '') {
      parts.push({
        name: PASTE_NAME,
        file: new File([state.fastaPaste], PASTE_NAME, { type: 'text/plain' }),
        origin: { kind: 'pasted' },
      });
    }
    for (const entry of state.fastaFiles) {
      parts.push({
        name: entry.name,
        file: entry.file,
        origin: entry.inflatedFrom === undefined ? { kind: 'picked' } : { kind: 'inflated', from: entry.inflatedFrom },
      });
    }
    if (parts.length === 0) {
      throw new Error('Enter sequences or open FASTA files.');
    }

    const makeRequest = (inputEntry: ToolDirectoryFile): ToolRequest => {
      const directory = assembleDirectory([inputEntry, ...(ustatFile === undefined ? [] : [ustatFile]), ...baseExtraDir]);
      const argv = this.buildArgv(state, extraWords, { inputName: inputEntry.name, ustatName: ustatFile?.name });
      return { tool: 'windowmasker', stage, argv, directory };
    };

    if (parts.length === 1) return [makeRequest(parts[0]!)];
    if (state.fastaMode === 'separate') return parts.map(makeRequest);
    const combined = await combineFastaFiles(parts);
    return [
      makeRequest({
        name: COMBINED_NAME,
        file: combined,
        origin: { kind: 'combined', parts: Object.freeze(parts.map((part) => part.name)) },
      }),
    ];
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
    if (!this.usesFastaBox(state)) {
      const name = state.blastDbIn.trim();
      return name === '' ? undefined : name;
    }
    const count = (state.fastaPaste.trim() !== '' ? 1 : 0) + state.fastaFiles.length;
    if (count === 0) return undefined;
    if (count === 1) return state.fastaPaste.trim() !== '' ? PASTE_NAME : state.fastaFiles[0]!.name;
    if (state.fastaMode === 'combined') return COMBINED_NAME;
    return state.fastaPaste.trim() !== '' ? PASTE_NAME : state.fastaFiles[0]!.name;
  }

  private async inflateFastaEntry(key: string, file: File): Promise<void> {
    try {
      const ref = await this.deps.data.inflateToSessionFile(file);
      const current = this.state.get().fastaFiles;
      const stillPresent = current.some((entry) => entry.key === key);
      if (!stillPresent) {
        this.releaseSessionToken(ref.token);
        return;
      }
      const conflict = current.some((entry) => entry.key !== key && entry.name === ref.name);
      if (conflict) {
        this.releaseSessionToken(ref.token);
        this.updateFastaEntry(key, { status: 'failed', error: `A file named "${ref.name}" is already added.` });
        return;
      }
      this.updateFastaEntry(key, {
        name: ref.name,
        size: ref.size,
        file: ref.file,
        status: 'ready',
        ...(ref.inflatedFrom === undefined ? {} : { inflatedFrom: ref.inflatedFrom }),
        ...(ref.token === undefined ? {} : { sessionToken: ref.token }),
      });
      this.scheduleValidation();
    } catch (error) {
      this.updateFastaEntry(key, {
        status: 'failed',
        error: failureText(error),
        ...(isUnavailable(error) ? { unavailable: true as const } : {}),
      });
    }
  }

  private updateFastaEntry(key: string, patch: Partial<ToolInputFile>): void {
    this.set({
      fastaFiles: Object.freeze(
        this.state.get().fastaFiles.map((entry) => (entry.key === key ? { ...entry, ...patch } : entry)),
      ),
    });
  }

  private loadDescription(force = false): Promise<ToolDescription | undefined> {
    if (!force && this.descriptionRequest !== undefined) return this.descriptionRequest;
    const req = this.deps.engine.describeTool('windowmasker');
    this.descriptionRequest = req;
    return req.then(
      (description) => {
        this.set({ description, descriptionError: undefined, descriptionUnavailable: undefined });
        this.scheduleValidation();
        return description;
      },
      (error: unknown) => {
        if (this.descriptionRequest === req) this.descriptionRequest = undefined;
        const unavailable = isUnavailable(error);
        this.set({
          descriptionError: errorMessage(error),
          ...(unavailable ? { descriptionUnavailable: true as const } : {}),
        });
        return undefined;
      },
    );
  }

  private scheduleValidation(): void {
    if (this.validationTimer !== undefined) clearTimeout(this.validationTimer);
    const delay = this.deps.debounceMs ?? 150;
    this.validationTimer = setTimeout(() => {
      this.validationTimer = undefined;
      void this.runValidationNow();
    }, delay);
  }

  /** Runs `validateTool` immediately (serialized so older answers are dropped). */
  runValidationNow(): Promise<void> {
    const generation = ++this.validationGeneration;
    const stale = () => generation !== this.validationGeneration;
    const run = async () => {
      const state = this.state.get();
      if (stale()) return;
      let extraWords: readonly string[];
      try {
        extraWords = this.parseExtra(state);
      } catch (error) {
        this.markInvalid(error);
        return;
      }
      let argv: readonly string[];
      try {
        const inputName = this.previewInputName(state) ?? (state.stage === 'convert' ? 'counts.ascii' : 'input.fa');
        const ustatName = state.stage === 'mask' ? (state.ustatSource.kind === 'none' ? 'counts.ascii' : state.ustatSource.name) : undefined;
        argv = this.buildArgv(state, extraWords, { inputName, ustatName });
      } catch (error) {
        this.markInvalid(error);
        return;
      }
      this.set({ validation: { state: 'checking' } });
      const directory: DirectoryFile[] = state.runFiles.map((entry) => ({ name: entry.name, file: entry.file }));
      let verdict: ToolValidation;
      try {
        verdict = await this.deps.engine.validateTool(argv, { directory });
      } catch (error) {
        if (stale()) return;
        const message = errorMessage(error);
        this.set({
          validation: isUnavailable(error)
            ? { state: 'unavailable', message }
            : { state: 'invalid', message },
        });
        return;
      }
      if (stale()) return;
      this.validatedKey = this.validationKey(state);
      this.set({
        validation: verdict.ok ? { state: 'ok' } : { state: 'invalid', message: verdict.message },
      });
    };
    const next = this.validationTurn.then(run, run);
    this.validationTurn = next.then(() => undefined, () => undefined);
    return next;
  }

  private validationKey(state: ToolDraftState): string {
    return JSON.stringify([
      state.stage,
      state.values[state.stage],
      state.extra[state.stage],
      state.outNames[state.stage],
      state.runFiles.map((entry) => entry.name),
    ]);
  }

  private trackAsking(promise: Promise<unknown>): void {
    this.asking.add(promise);
    void promise.finally(() => this.asking.delete(promise));
  }

  private set(change: { readonly [K in keyof ToolDraftState]?: ToolDraftState[K] | undefined }): void {
    const next = { ...this.state.get(), ...change };
    for (const key of Object.keys(change) as Array<keyof ToolDraftState>) {
      if (change[key] === undefined) delete next[key];
    }
    this.state.set(next as ToolDraftState);
  }
}

/** Whether the WindowMasker tab holds user inputs or files that a reload would lose. */
export function holdsToolInputs(state: ToolDraftState): boolean {
  return (
    state.fastaPaste.trim() !== '' ||
    state.fastaFiles.length > 0 ||
    state.runFiles.length > 0 ||
    state.ustatSource.kind !== 'none' ||
    state.convertSource.kind !== 'none'
  );
}

function errorMessage(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}

function sessionCountsChoice(runId: string, runNumber: number, key: string, name: string, size: number): ToolCountsChoice {
  return { kind: 'session', runId, runNumber, key, name, size };
}

function savedCountsChoice(taxid: number, size: number): ToolCountsChoice {
  return { kind: 'saved', taxid, name: countsName(taxid), size };
}

function fileCountsChoice(file: File | undefined): ToolCountsChoice {
  return file === undefined ? { kind: 'none' } : { kind: 'file', file, name: normalizeEntryName(file.name), size: file.size };
}

/** Normalizes a picked or folder-relative entry name for the run directory. */
export function normalizeEntryName(raw: string): string {
  let name = raw.replaceAll('\\', '/').trim();
  while (name.startsWith('./')) name = name.slice(2);
  while (name.startsWith('/')) name = name.slice(1);
  return name;
}

/** Ensures no two files in a tool run directory share a name. */
function assembleDirectory(entries: readonly ToolDirectoryFile[]): readonly ToolDirectoryFile[] {
  const seen = new Set<string>();
  for (const entry of entries) {
    if (seen.has(entry.name)) {
      throw new Error(`Two files in the run directory have the same name: "${entry.name}".`);
    }
    seen.add(entry.name);
  }
  return Object.freeze(entries);
}

/**
 * Joins FASTA files into one `combined_input.fa` Blob without copying their bytes into page memory:
 * reads only the final 1 byte of each non-final file to insert `\n` when missing.
 */
async function combineFastaFiles(parts: readonly ToolDirectoryFile[]): Promise<File> {
  const blobs: BlobPart[] = [];
  for (let i = 0; i < parts.length; i++) {
    const part = parts[i]!;
    blobs.push(part.file);
    if (i + 1 < parts.length && part.file.size > 0) {
      const tail = new Uint8Array(await part.file.slice(part.file.size - 1, part.file.size).arrayBuffer());
      if (tail[0] !== LF) blobs.push('\n');
    }
  }
  return new File(blobs, COMBINED_NAME, { type: 'text/plain' });
}
