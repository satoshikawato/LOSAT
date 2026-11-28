// The FASTA input of the tool tabs that read sequences (WindowMasker's `ToolDraft`, Masking's
// `MaskingDraft`; WW4c, one code path for both): pasted FASTA (`input.fa`), FASTA files opened or
// dropped, compressed ones (`.gz`, `.bz2`, `.zst`) inflated into session OPFS via
// `DataGateway.inflateToSessionFile` and held while the form or a run reads them (`SessionFiles`),
// several inputs as one run (`combined_input.fa`, joined as Blob references, checking only the last
// byte of each file so page JS heap never holds the sequences) or one run each, and, for an input
// that is not FASTA (`-infmt blastdb`), the `-in` name with the files of the run that it names: typed, or
// worked out from the files of the run, written as SeqDB reads it and refused when it cannot be (WW4d).
import type { DirectoryFile, PrintedKind } from 'losat/domain';
import { isUnavailable } from 'losat/ports';
import { normalizeRunPath } from '../domain/blast-db';
import type { FormValues } from '../domain/parameters';
import type { ToolId } from '../domain/tools';
import type { DataGateway } from '../ports/data';
import { failureText } from './availability';
import type { ToolDirectoryFile } from './coordinator';
import type { InputMode } from './draft';
import { SessionFiles } from './session-files';
import { isSeqdbList } from '../domain/list-flags';
import { namesOption } from '../domain/tool-argv';
import {
  foundDatabase,
  ToolDraftCore,
  type DatabaseOption,
  type StatePatch,
  type ToolDraftCoreDeps,
  type ToolFormState,
} from './tool-draft-core';

/** One FASTA file added to a tool tab's Input block (inflated into OPFS when compressed). */
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

/** The state of a tool tab's FASTA input (and of its `-in` name for an input that is not FASTA). */
export interface FastaInputState extends ToolFormState {
  readonly fastaPaste: string;
  readonly fastaFiles: readonly ToolInputFile[];
  readonly fastaMode: InputMode;
  /** The Input name typed (`-in` of an input that is not FASTA), as the command line takes it. */
  readonly blastDbIn: string;
  /**
   * Whether the user typed the Input name; until then, and while what is typed is blank, the files of the run name it
   * (`inName`). The field shows what was typed, a blank one with the name found as its placeholder.
   */
  readonly blastDbInEdited: boolean;
}

/** The FASTA input of a new form: empty, several inputs combined. */
export const EMPTY_FASTA_INPUT: Pick<FastaInputState, 'fastaPaste' | 'fastaFiles' | 'fastaMode' | 'blastDbIn' | 'blastDbInEdited'> =
  Object.freeze({ fastaPaste: '', fastaFiles: Object.freeze([]), fastaMode: 'combined', blastDbIn: '', blastDbInEdited: false });

export interface FastaToolDraftDeps extends ToolDraftCoreDeps {
  readonly data: Pick<DataGateway, 'inflateToSessionFile' | 'dropSessionFile'>;
  /** Shared holder tracker for inflated session files; defaults to a draft-local instance in unit tests. */
  readonly sessionFiles?: SessionFiles;
}

/**
 * The input of one run: the name that `-in` gives (none where words of Additional options name `-in`: `namesOption`), the
 * file of the run directory under it (none for a named BLAST database), and the session files that it reads.
 */
export interface FastaRunInput {
  readonly name?: string;
  readonly file?: ToolDirectoryFile;
  readonly sessionTokens: readonly string[];
}

/** What the Input block of a tool tab asks of its draft (`ToolFastaInputs.vue`). */
export interface FastaInputForm {
  usesFastaBox(state?: FastaInputState): boolean;
  setFastaPaste(text: string): void;
  setFastaMode(mode: InputMode): void;
  addFastaFiles(files: readonly File[]): void;
  removeFastaFile(key: string): void;
  setBlastDbIn(name: string): void;
  inName(state?: FastaInputState): string;
  databaseRefusal(state?: FastaInputState): string | undefined;
  addRunFiles(files: readonly File[], relativeNames?: readonly string[]): void;
  dropInputFiles(files: readonly File[]): void;
}

/**
 * The input of the check of a form (`FastaToolDraft.checkInput`): the name that `-in` gives (none where the words name
 * `-in`), the files of the check's directory that hold it, and whether the check reads it as FASTA in place of the form's
 * BLAST database (a stand-in for a database not named yet).
 */
export interface CheckInput {
  readonly name: string | undefined;
  readonly files: readonly DirectoryFile[];
  readonly asFasta: boolean;
}

const PASTE_NAME = 'input.fa';
/** The name of the input that stands in for one that the form does not have yet in the check (`checkInput`): a paste's. */
const STAND_IN_NAME = PASTE_NAME;
const COMBINED_NAME = 'combined_input.fa';
const LF = 0x0a;

export abstract class FastaToolDraft<S extends FastaInputState> extends ToolDraftCore<S> implements FastaInputForm {
  protected readonly sessionFiles: SessionFiles;
  private readonly fastaData: FastaToolDraftDeps['data'];
  /** The file of the pasted FASTA (`input.fa`), made once for each text. */
  private pasted: { readonly text: string; readonly file: File } | undefined;

  protected constructor(tool: ToolId, deps: FastaToolDraftDeps, initial: S) {
    const sessionFiles =
      deps.sessionFiles ?? new SessionFiles((token) => void deps.data.dropSessionFile(token).catch(() => undefined));
    super(tool, { ...deps, sessionFiles }, initial);
    this.sessionFiles = sessionFiles;
    this.fastaData = deps.data;
  }

  /** Whether the form reads `-in` from the FASTA box (paste and files) rather than as a name (a BLAST database). */
  abstract usesFastaBox(state?: S): boolean;

  /**
   * What a press of `state` with the words `extra` of Additional options prints in place of a run (help or the version),
   * or null when it runs the tool: the one answer that the press (`buildRequests`, which then starts one print request)
   * and the count of runs (`runCount`) read.
   */
  protected abstract printsKindOf(state: S, extra: readonly string[]): PrintedKind | null;

  /**
   * The `-in` option of the input that is not FASTA, with the `-infmt` that decides whether it is a SeqDB name list
   * (`isSeqdbList`: with `-infmt blastdb`), and the molecule kind of the database that the tool reads; undefined when
   * `-in` names no input of this block (WindowMasker's Convert).
   */
  protected abstract inDatabase(state: S): DatabaseOption | undefined;

  // --- the FASTA box and the -in name ---------------------------------------------------------

  setFastaPaste(text: string): void {
    this.setInput({ fastaPaste: text, message: undefined });
    this.scheduleValidation();
  }

  setFastaMode(mode: InputMode): void {
    this.setInput({ fastaMode: mode, message: undefined });
    this.scheduleValidation();
  }

  addFastaFiles(files: readonly File[]): void {
    if (files.length === 0) return;
    const existing = new Set(this.state.get().fastaFiles.map((entry) => entry.name));
    const added: ToolInputFile[] = [];
    for (const file of files) {
      const entryName = normalizeRunPath(file.name);
      if (entryName === undefined) {
        this.reportError(`Invalid file name: "${file.name}".`);
        continue;
      }
      if (existing.has(entryName)) {
        this.reportError(`A file named "${entryName}" is already added.`);
        continue;
      }
      existing.add(entryName);
      const key = `fa-${this.nextKey++}`;
      added.push({ key, name: entryName, size: file.size, file, status: 'inflating' });
    }
    if (added.length === 0) return;
    this.setInput({
      fastaFiles: Object.freeze([...this.state.get().fastaFiles, ...added]),
      message: undefined,
    });
    for (const entry of added) void this.inflateFastaEntry(entry.key, entry.file);
  }

  removeFastaFile(key: string): void {
    const target = this.state.get().fastaFiles.find((entry) => entry.key === key);
    if (target === undefined) return;
    this.sessionFiles.release(this.holder('form', key));
    this.setInput({
      fastaFiles: Object.freeze(this.state.get().fastaFiles.filter((entry) => entry.key !== key)),
      message: undefined,
    });
    this.scheduleValidation();
  }

  setBlastDbIn(name: string): void {
    this.setInput({ blastDbIn: name, blastDbInEdited: true, message: undefined });
    this.scheduleValidation();
  }

  /**
   * The Input name (`-in` when the form does not read the FASTA box): the name typed, or, while none is (a blank one is
   * none), the BLAST database of the tool's kind that the files of the run hold when `-in` is a SeqDB name list
   * (`foundDatabase`: `-infmt blastdb`, `"my dbs/wm"`; line WW, WW4d), as makeblastdb's `dbInName` is; none when the list
   * cannot carry its name (`databaseRefusal`) or words of Additional options name `-in`. With `-infmt seqids` or
   * `-fa_list`, `-in` is one file's name, which the user types.
   */
  inName(state: S = this.state.get()): string {
    if (typedName(state)) return state.blastDbIn;
    return this.foundIn(state)?.written ?? '';
  }

  /**
   * Why the BLAST database that the files of the run hold cannot be `-in` (a name with a double quote), while it would
   * name `-in`: `-infmt blastdb`, its files in the run, no name typed, no `-in` in Additional options. The Input block,
   * the check and "Run LOSAT" say it for as long as that lasts (line WW, WW4d R1 and R2).
   */
  databaseRefusal(state: S = this.state.get()): string | undefined {
    return typedName(state) ? undefined : this.foundIn(state)?.refusal;
  }

  /** Files dropped on the Input block: FASTA inputs while the form reads the FASTA box, else files of the run. */
  dropInputFiles(files: readonly File[]): void {
    if (this.usesFastaBox()) this.addFastaFiles(files);
    else this.addRunFiles(files);
  }

  /**
   * How many runs Run LOSAT starts, for every tab with a FASTA box (WindowMasker, Masking): one, or one for each FASTA
   * input when Separate is chosen, and one for a help or version word, as the press starts one print request
   * (`printsKindOf`; line WQA's D33, S4-8, R1 M3 and R2 H1). Words that cannot be read start no run; the count is then
   * the inputs' (the check says why).
   */
  runCount(state: S = this.state.get()): number {
    if (!this.usesFastaBox(state)) return 1;
    const count = (state.fastaPaste.trim() !== '' ? 1 : 0) + state.fastaFiles.length;
    if (state.fastaMode !== 'separate' || count <= 1) return 1;
    let extra: readonly string[];
    try {
      extra = this.parseExtra(state);
    } catch {
      return count;
    }
    return this.printsKindOf(state, extra) === null ? count : 1;
  }

  // --- for the tools --------------------------------------------------------------------------

  /** The `-in` name that the FASTA box or the `-in` name gives the argv shown, if there is an input. */
  protected fastaInputName(state: S): string | undefined {
    if (!this.usesFastaBox(state)) {
      const name = this.inName(state).trim();
      return name === '' ? undefined : name;
    }
    const count = (state.fastaPaste.trim() !== '' ? 1 : 0) + state.fastaFiles.length;
    if (count === 0) return undefined;
    if (count === 1) return state.fastaPaste.trim() !== '' ? PASTE_NAME : state.fastaFiles[0]!.name;
    if (state.fastaMode === 'combined') return COMBINED_NAME;
    return state.fastaPaste.trim() !== '' ? PASTE_NAME : state.fastaFiles[0]!.name;
  }

  /**
   * The FASTA input of the check's directory, under the name that the run gives it (review of WW4c, F2): the
   * pasted FASTA or the opened file of the first run, or for several inputs combined `combined_input.fa` (the
   * files one after another: the check opens the input and reads its start, as the run does before `-out`, which
   * the line ends that the run puts between files do not change); none for an input that is not FASTA (its files
   * are files of the run). So the check refuses an `-out` that names the input, as the run does.
   */
  protected fastaCheckInputs(state: S): readonly DirectoryFile[] {
    if (!this.usesFastaBox(state)) return [];
    const parts: DirectoryFile[] = [];
    if (state.fastaPaste.trim() !== '') parts.push({ name: PASTE_NAME, file: this.pasteFile(state.fastaPaste) });
    for (const entry of state.fastaFiles) parts.push({ name: entry.name, file: entry.file });
    if (parts.length <= 1 || state.fastaMode === 'separate') return parts.slice(0, 1);
    return [{ name: COMBINED_NAME, file: new Blob(parts.map((part) => part.file), { type: 'text/plain' }) }];
  }

  /**
   * The input of the check (`validationArgv`, `validationInputs`): the form's (`fastaInputName`, `fastaCheckInputs`), or,
   * while it has none and the words do not name `-in` themselves, a stand-in that the run's first steps open and read as
   * an empty input — the file `input.fa` of the check's directory, empty (`standIns`), read as FASTA: `asFasta` asks the
   * argv of the check to leave the form's `-infmt blastdb` out, as no BLAST database is empty. The check, which opens
   * `-in` as the run does (line NP WW4f), so reports the form's other errors and not the input that it lacks, which Run
   * LOSAT asks for ("Enter sequences or open FASTA files.", "Enter an input name or open files."); an empty FASTA box with
   * `-infmt blastdb` in Additional options is checked as an empty paste is (review R0 of WW4f, F1).
   */
  protected checkInput(state: S, extra: readonly string[]): CheckInput {
    const name = this.fastaInputName(state);
    if (name !== undefined) return { name, files: this.fastaCheckInputs(state), asFasta: false };
    if (namesOption(extra, '-in')) return { name: undefined, files: [], asFasta: false };
    const database = this.inDatabase(state);
    return {
      name: STAND_IN_NAME,
      files: this.standIns(state, [STAND_IN_NAME]),
      asFasta: database !== undefined && isSeqdbList(database.option),
    };
  }

  /**
   * The empty files of the check's directory that stand in for the inputs `names` that the check's argv names and the
   * form does not have yet (`checkInput`, WindowMasker's counts), none for a name that a file of the run holds or a folder
   * of them has (the run's file is then the input, as for the run).
   */
  protected standIns(state: S, names: readonly string[]): readonly DirectoryFile[] {
    const claimed = (name: string) => state.runFiles.some((entry) => entry.name === name || entry.name.startsWith(`${name}/`));
    return names.filter((name) => !claimed(name)).map((name) => ({ name, file: new Blob([]) }));
  }

  /** Whether the FASTA input or the files of the run hold anything that a reload would lose. */
  protected holdsFastaInputs(state: S): boolean {
    return state.fastaPaste.trim() !== '' || state.fastaFiles.length > 0 || state.runFiles.length > 0;
  }

  /**
   * The inputs of the runs of a submit: the `-in` name (a BLAST database, whose files are files of
   * the run), or the FASTA inputs: one, one for each (Separate), or `combined_input.fa` (Combined).
   * Throws the message to show when the inputs are not ready (a database that `-in` cannot name: `databaseRefusal`).
   */
  protected async runInputs(state: S): Promise<readonly FastaRunInput[]> {
    // Words of Additional options that name -in give the input themselves: the form asks for none (makeblastdb's rule).
    const wordsGiveIn = this.wordsNameIn(state);
    if (!this.usesFastaBox(state)) {
      const refusal = this.databaseRefusal(state);
      if (refusal !== undefined) throw new Error(refusal);
      const name = this.inName(state).trim();
      if (name === '') {
        if (wordsGiveIn) return [{ sessionTokens: [] }];
        throw new Error('Enter an input name or open files.');
      }
      return [{ name, sessionTokens: [] }];
    }

    if (state.fastaFiles.some((entry) => entry.status === 'inflating')) {
      throw new Error('Wait until the opened files finish reading.');
    }
    const failed = state.fastaFiles.find((entry) => entry.status === 'failed');
    if (failed !== undefined) {
      throw new Error(`${failed.name}: ${failed.error ?? 'could not be read'}`);
    }

    const parts: Array<ToolDirectoryFile & { readonly sessionToken?: string }> = [];
    if (state.fastaPaste.trim() !== '') {
      parts.push({ name: PASTE_NAME, file: this.pasteFile(state.fastaPaste), origin: { kind: 'pasted' } });
    }
    for (const entry of state.fastaFiles) {
      parts.push({
        name: entry.name,
        file: entry.file,
        origin: entry.inflatedFrom === undefined ? { kind: 'picked' } : { kind: 'inflated', from: entry.inflatedFrom },
        ...(entry.sessionToken === undefined ? {} : { sessionToken: entry.sessionToken }),
      });
    }
    if (parts.length === 0) {
      if (wordsGiveIn) return [{ sessionTokens: [] }];
      throw new Error('Enter sequences or open FASTA files.');
    }

    const inputOf = (part: ToolDirectoryFile & { readonly sessionToken?: string }): FastaRunInput => ({
      name: part.name,
      file: { name: part.name, file: part.file, origin: part.origin },
      sessionTokens: part.sessionToken === undefined ? [] : [part.sessionToken],
    });
    if (parts.length === 1 || state.fastaMode === 'separate') return parts.map(inputOf);
    const combined = await combineFastaFiles(parts);
    return [
      {
        name: COMBINED_NAME,
        file: {
          name: COMBINED_NAME,
          file: combined,
          origin: { kind: 'combined', parts: Object.freeze(parts.map((part) => part.name)) },
        },
        sessionTokens: parts.flatMap((part) => (part.sessionToken === undefined ? [] : [part.sessionToken])),
      },
    ];
  }

  protected override heldTokens(state: S): readonly string[] {
    return state.fastaFiles.flatMap((entry) => (entry.sessionToken === undefined ? [] : [entry.sessionToken]));
  }

  protected override inputsUnavailable(state: S): string | undefined {
    for (const file of state.fastaFiles) {
      if (file.status === 'failed' && file.unavailable === true && file.error !== undefined) return file.error;
    }
    return undefined;
  }

  // --- internals --------------------------------------------------------------------------------

  /** Whether the words of Additional options name `-in` themselves (`namesOption`): they give the input. */
  protected wordsNameIn(state: S): boolean {
    return namesOption(this.wordsOf(state), '-in');
  }

  /** The BLAST database of the files of the run, for `-in` (`foundDatabase`; none when `-in` names no input here). */
  private foundIn(state: S): ReturnType<typeof foundDatabase> {
    const input = this.inDatabase(state);
    return input === undefined ? undefined : foundDatabase(input, state.runFiles, this.wordsOf(state));
  }

  /** The file of the pasted FASTA `text` (`input.fa`), the same for the check and the runs of one text. */
  private pasteFile(text: string): File {
    if (this.pasted?.text !== text) this.pasted = { text, file: new File([text], PASTE_NAME, { type: 'text/plain' }) };
    return this.pasted.file;
  }

  /** Changes the state of the FASTA input (a key set to `undefined` is removed). */
  private setInput(change: StatePatch<FastaInputState>): void {
    this.set(change as StatePatch<S>);
  }

  private releaseSessionToken(token: string | undefined): void {
    if (token !== undefined) void this.fastaData.dropSessionFile(token).catch(() => undefined);
  }

  private async inflateFastaEntry(key: string, file: File): Promise<void> {
    try {
      const ref = await this.fastaData.inflateToSessionFile(file);
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
      if (ref.token !== undefined) this.sessionFiles.set(this.holder('form', key), [ref.token]);
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
    this.setInput({
      fastaFiles: Object.freeze(
        this.state.get().fastaFiles.map((entry) => (entry.key === key ? { ...entry, ...patch } : entry)),
      ),
    });
  }
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

/**
 * `values` without `-infmt`, for the check's stand-in of a BLAST database not named yet (`CheckInput.asFasta`): FASTA,
 * which is `-infmt`'s default for each tool that reads a database.
 */
export function withoutInfmt(values: FormValues): FormValues {
  return Object.fromEntries(Object.entries(values).filter(([flag]) => flag !== '-infmt'));
}

/** Whether the Input name was typed and is not blank (a blank one is none: the files of the run name `-in` then). */
function typedName(state: FastaInputState): boolean {
  return state.blastDbInEdited && state.blastDbIn.trim() !== '';
}
