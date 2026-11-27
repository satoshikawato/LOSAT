// The makeblastdb form draft (line WW, WW4b2; design DESIGN_ww4b.md §3): independent of queued and
// completed runs, on the tool drafts' shared core (`ToolDraftCore`: `-xmlhelp`, files of the run, files
// of options, check, submit). What is its own: the input that `-in` names, by `-input_type` (FASTA
// pasted or opened, ASN.1 files, or a BLAST database of the files of the run), every option in `values`,
// and the hand-over of a masking run's mask data to `-mask_data` (`addMaskData`).
//
// The run's directory holds the input files as the user opened them: makeblastdb reads `.gz`, `.bz2`
// and `.zst` itself (its `-in` names the compressed file, as on the command line), and several files
// are one `-in` list, as the command line takes them (`SeqDB_SplitQuoted`), not one file made of them.
import type { DirectoryFile } from 'losat/domain';
import { splitOptions } from '../domain/argv';
import { duplicateRunFileMessage, normalizeRunPath, primaryBlastDbName } from '../domain/blast-db';
import type { FieldValue, FormValues } from '../domain/parameters';
import type { ToolFileOrigin } from '../domain/run';
import { buildMakeblastdbArgv, fileListNames, makeblastdbInputList } from '../domain/tool-argv';
import { makeblastdbInputType, makeblastdbOptionParameters, toolPrintsKind, type MakeblastdbInputType } from '../domain/tools';
import type { ToolDirectoryFile, ToolRequest } from './coordinator';
import { ToolDraftCore, type StatePatch, type ToolDraftCoreDeps, type ToolExtraFile, type ToolFormState } from './tool-draft-core';

export interface MakeblastdbDraftState extends ToolFormState {
  /** Option values keyed by `-flag` (`-input_type`, `-dbtype`, `-out`, and the options of the groups). */
  readonly values: FormValues;
  /** Words typed in "Additional options". */
  readonly extra: string;
  /** FASTA pasted in the Input block (`input.fa` in the run directory). */
  readonly fastaPaste: string;
  /** The files that `-in` names for FASTA and ASN.1 input, in the run directory under their names. */
  readonly inputFiles: readonly ToolExtraFile[];
  /** `-in` for `-input_type blastdb`: the database's name (its files are files of the run). */
  readonly dbIn: string;
  /** Whether the user typed `dbIn`; until then the files of the run name it. */
  readonly dbInEdited: boolean;
}

export type MakeblastdbDraftDeps = ToolDraftCoreDeps;

/** The name of pasted FASTA in the run directory. */
export const MAKEBLASTDB_PASTE_NAME = 'input.fa';

export class MakeblastdbDraft extends ToolDraftCore<MakeblastdbDraftState> {
  constructor(deps: MakeblastdbDraftDeps) {
    super('makeblastdb', deps, {
      values: Object.freeze({}),
      extra: '',
      fastaPaste: '',
      inputFiles: Object.freeze([]),
      dbIn: '',
      dbInEdited: false,
      runFiles: Object.freeze([]),
      validation: { state: 'idle' },
      submitting: false,
    });
  }

  // --- options --------------------------------------------------------------------------------

  setField(flag: string, value: FieldValue): void {
    const next: Record<string, FieldValue> = { ...this.state.get().values };
    if (value === undefined || value === '' || value === false) delete next[flag];
    else next[flag] = value;
    this.set({ values: Object.freeze(next), message: undefined });
    this.scheduleValidation();
  }

  setExtra(text: string): void {
    this.set({ extra: text, message: undefined });
    this.scheduleValidation();
  }

  /** Resets the options of "makeblastdb options" and Additional options; the files chosen for them go too. */
  resetOptions(): void {
    const { values, description, runFiles } = this.state.get();
    const flags = new Set(makeblastdbOptionParameters(description).map((param) => param.flag));
    const next: Record<string, FieldValue> = { ...values };
    for (const flag of flags) delete next[flag];
    this.set({
      values: Object.freeze(next),
      extra: '',
      runFiles: Object.freeze(runFiles.filter((entry) => entry.flag === undefined || !flags.has(entry.flag))),
      runFilesError: undefined,
      message: undefined,
    });
    this.scheduleValidation();
  }

  /**
   * "Use as -mask_data in makeblastdb" on a masking run's mask data (WindowMasker's Mask with `-outfmt maskinfo_*`; the
   * masking tools of line WW WW4c the same way): the file joins `-mask_data`'s list, under `name`, as a file of the run.
   * The mask data of an earlier run handed over under the same name is replaced in its place in the list (as
   * WindowMasker's "Mask with these counts" replaces the counts chosen). A name that the list cannot carry (a comma
   * separates its files) or that another file of the run holds is refused in the tab's message, where the user is.
   */
  addMaskData(file: File, name: string, origin: ToolFileOrigin): void {
    this.ensureLoaded();
    const state = this.state.get();
    const path = normalizeRunPath(name);
    const earlier = state.runFiles.find((entry) => entry.name === path && entry.flag === '-mask_data' && entry.origin?.kind === 'made');
    if (path !== undefined && earlier !== undefined) {
      const replaced: ToolExtraFile = { ...earlier, key: `rf-${this.nextKey++}`, size: file.size, file, origin };
      const names = fileListNames(state.values['-mask_data']);
      this.set({
        runFiles: Object.freeze(state.runFiles.map((entry) => (entry === earlier ? replaced : entry))),
        runFilesError: undefined,
        message: undefined,
      });
      this.setField('-mask_data', (names.includes(path) ? names : [...names, path]).join(','));
      return;
    }
    this.setOptionFile('-mask_data', file, { name, origin });
    const error = this.state.get().runFilesError;
    if (error !== undefined) this.set({ runFilesError: undefined, message: { kind: 'error', text: error } });
  }

  // --- the input ------------------------------------------------------------------------------

  setFastaPaste(text: string): void {
    this.set({ fastaPaste: text, message: undefined });
    this.scheduleValidation();
  }

  /** Adds files that `-in` names (FASTA or ASN.1); refuses a name already there. */
  addInputFiles(files: readonly File[]): void {
    if (files.length === 0) return;
    const state = this.state.get();
    const taken = new Set([...state.inputFiles.map((entry) => entry.name), ...(state.fastaPaste.trim() === '' ? [] : [MAKEBLASTDB_PASTE_NAME])]);
    const added: ToolExtraFile[] = [];
    let error: string | undefined;
    for (const file of files) {
      const name = normalizeRunPath(file.name);
      if (name === undefined) {
        error ??= `Invalid file name: "${file.name}".`;
        continue;
      }
      if (taken.has(name)) {
        error ??= duplicateRunFileMessage(name);
        continue;
      }
      taken.add(name);
      added.push({ key: `in-${this.nextKey++}`, name, size: file.size, file });
    }
    this.set({
      inputFiles: Object.freeze([...state.inputFiles, ...added]),
      message: error === undefined ? undefined : { kind: 'error', text: error },
    });
    this.scheduleValidation();
  }

  removeInputFile(key: string): void {
    const state = this.state.get();
    if (!state.inputFiles.some((entry) => entry.key === key)) return;
    this.set({ inputFiles: Object.freeze(state.inputFiles.filter((entry) => entry.key !== key)), message: undefined });
    this.scheduleValidation();
  }

  /**
   * `-in` for `-input_type blastdb`, as typed: the command line's value (`SeqDB_SplitQuoted` reads a list of names, a
   * name with a space between double quotes), written as it is.
   */
  setDbIn(name: string): void {
    this.set({ dbIn: name, dbInEdited: true, message: undefined });
    this.scheduleValidation();
  }

  /** The input type that `-input_type` chooses. */
  inputType(state: MakeblastdbDraftState = this.state.get()): MakeblastdbInputType {
    return makeblastdbInputType(state.values['-input_type']);
  }

  // --- the command line -----------------------------------------------------------------------

  /** The argv shown beside Run LOSAT; a quote left open in Additional options leaves their words out (the check reports it). */
  previewArgv(state: MakeblastdbDraftState = this.state.get()): readonly string[] {
    let extra: readonly string[];
    try {
      extra = this.parseExtra(state);
    } catch {
      extra = [];
    }
    try {
      return this.argvOf(state, extra);
    } catch {
      return this.argvOf(state, []);
    }
  }

  holdsInputs(state: MakeblastdbDraftState = this.state.get()): boolean {
    return state.fastaPaste.trim() !== '' || state.inputFiles.length > 0 || state.runFiles.length > 0;
  }

  // --- what the core asks -------------------------------------------------------------------

  protected parseExtra(state: MakeblastdbDraftState): readonly string[] {
    return splitOptions(state.extra);
  }

  protected optionValue(state: MakeblastdbDraftState, flag: string): FieldValue | undefined {
    return state.values[flag];
  }

  /** The `-in` name that files added to the run directory name, until the user edits it (a BLAST database). */
  protected override runFilesAdded(state: MakeblastdbDraftState, runFiles: readonly ToolExtraFile[]): StatePatch<MakeblastdbDraftState> {
    if (state.dbInEdited && state.dbIn.trim() !== '') return {};
    const detected = primaryBlastDbName(runFiles.map((entry) => entry.name));
    return detected === undefined ? {} : { dbIn: detected, dbInEdited: false };
  }

  /**
   * The argv of the check; none before an input is given, unless the words name `-in` themselves or print help or the
   * version (makeblastdb would then read its standard input, which a run on the Web has empty).
   */
  protected validationArgv(state: MakeblastdbDraftState, extra: readonly string[]): readonly string[] | undefined {
    const argv = this.argvOf(state, extra);
    if (this.inputValue(state) !== '' || toolPrintsKind(argv) !== null) return argv;
    return namesInput(argv) ? argv : undefined;
  }

  protected override validationInputs(state: MakeblastdbDraftState): readonly DirectoryFile[] {
    return this.inputEntries(state).map((entry) => ({ name: entry.name, file: entry.file }));
  }

  protected async buildRequests(state: MakeblastdbDraftState, extra: readonly string[]): Promise<readonly ToolRequest[]> {
    const argv = this.argvOf(state, extra);
    const prints = toolPrintsKind(argv);
    if (prints !== null) {
      return [{ tool: 'makeblastdb', stage: 'build', argv, directory: this.runDirectoryOf(state, []), prints }];
    }
    if (this.inputValue(state) === '' && !namesInput(argv)) throw new Error(MISSING_INPUT[this.inputType(state)]);
    return [{ tool: 'makeblastdb', stage: 'build', argv, directory: this.runDirectoryOf(state, this.inputEntries(state)) }];
  }

  // --- internals ----------------------------------------------------------------------------

  private argvOf(state: MakeblastdbDraftState, extra: readonly string[]): readonly string[] {
    return buildMakeblastdbArgv({
      input: this.inputValue(state),
      values: state.values,
      ...(state.description === undefined ? {} : { description: state.description }),
      extra,
    });
  }

  /**
   * `-in`'s value for the input type chosen: a database's name as the user typed it (a list of databases, quotes and
   * all), or the name that the files of the run give, quoted when it holds a space; the files of FASTA and ASN.1 input
   * as one list.
   */
  private inputValue(state: MakeblastdbDraftState): string {
    if (this.inputType(state) === 'blastdb') {
      const name = state.dbIn.trim();
      return state.dbInEdited ? name : makeblastdbInputList([name]);
    }
    return makeblastdbInputList(this.inputEntries(state).map((entry) => entry.name));
  }

  /** The files that `-in` names for FASTA and ASN.1 input: the pasted FASTA (`input.fa`), then the files opened. */
  private inputEntries(state: MakeblastdbDraftState): readonly ToolDirectoryFile[] {
    const type = this.inputType(state);
    if (type === 'blastdb') return [];
    const entries: ToolDirectoryFile[] = [];
    if (type === 'fasta' && state.fastaPaste.trim() !== '') {
      entries.push({ name: MAKEBLASTDB_PASTE_NAME, file: new Blob([state.fastaPaste], { type: 'text/plain' }), origin: { kind: 'pasted' } });
    }
    for (const entry of state.inputFiles) entries.push({ name: entry.name, file: entry.file, origin: { kind: 'picked' } });
    return entries;
  }
}

/** Whether the words name `-in` themselves. */
function namesInput(argv: readonly string[]): boolean {
  return argv.slice(1).some((word) => /^-in(?:=|$)/.test(word));
}

/** What Run LOSAT says when the input type chosen has no input yet. */
const MISSING_INPUT: Readonly<Record<MakeblastdbInputType, string>> = {
  fasta: 'Enter sequences or open FASTA files.',
  asn1_bin: 'Open the files of -in.',
  asn1_txt: 'Open the files of -in.',
  blastdb: 'Enter the database name of -in, or open its files.',
};
