// The FASTA input of the tool tabs that read sequences (WindowMasker's `ToolDraft`, Masking's
// `MaskingDraft`; WW4c, one code path for both): pasted FASTA (`input.fa`), FASTA files opened or
// dropped, compressed ones (`.gz`, `.bz2`, `.zst`) inflated into session OPFS via
// `DataGateway.inflateToSessionFile` and held while the form or a run reads them (`SessionFiles`),
// several inputs as one run (`combined_input.fa`, joined as Blob references, checking only the last
// byte of each file so page JS heap never holds the sequences) or one run each, and, for an input
// that is not FASTA (`-infmt blastdb`), the `-in` name with the files of the run that it names.
import { isUnavailable } from 'losat/ports';
import { normalizeRunPath, primaryBlastDbName } from '../domain/blast-db';
import type { ToolId } from '../domain/tools';
import type { DataGateway } from '../ports/data';
import { failureText } from './availability';
import type { ToolDirectoryFile } from './coordinator';
import type { InputMode } from './draft';
import { SessionFiles } from './session-files';
import {
  ToolDraftCore,
  type StatePatch,
  type ToolDraftCoreDeps,
  type ToolExtraFile,
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
  readonly blastDbIn: string;
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

/** The input of one run: the name that `-in` gives, the file of the run directory under it (none for a named BLAST database), and the session files that it reads. */
export interface FastaRunInput {
  readonly name: string;
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
  addRunFiles(files: readonly File[], relativeNames?: readonly string[]): void;
  dropInputFiles(files: readonly File[]): void;
}

const PASTE_NAME = 'input.fa';
const COMBINED_NAME = 'combined_input.fa';
const LF = 0x0a;

export abstract class FastaToolDraft<S extends FastaInputState> extends ToolDraftCore<S> implements FastaInputForm {
  protected readonly sessionFiles: SessionFiles;
  private readonly fastaData: FastaToolDraftDeps['data'];

  protected constructor(tool: ToolId, deps: FastaToolDraftDeps, initial: S) {
    const sessionFiles =
      deps.sessionFiles ?? new SessionFiles((token) => void deps.data.dropSessionFile(token).catch(() => undefined));
    super(tool, { ...deps, sessionFiles }, initial);
    this.sessionFiles = sessionFiles;
    this.fastaData = deps.data;
  }

  /** Whether the form reads `-in` from the FASTA box (paste and files) rather than as a name (a BLAST database). */
  abstract usesFastaBox(state?: S): boolean;

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

  /** Files dropped on the Input block: FASTA inputs while the form reads the FASTA box, else files of the run. */
  dropInputFiles(files: readonly File[]): void {
    if (this.usesFastaBox()) this.addFastaFiles(files);
    else this.addRunFiles(files);
  }

  /** How many runs Run LOSAT starts: one, or one for each FASTA input when Separate is chosen. */
  runCount(state: S = this.state.get()): number {
    if (!this.usesFastaBox(state)) return 1;
    const count = (state.fastaPaste.trim() !== '' ? 1 : 0) + state.fastaFiles.length;
    return state.fastaMode === 'separate' && count > 1 ? count : 1;
  }

  // --- for the tools --------------------------------------------------------------------------

  /** The `-in` name that the FASTA box or the `-in` name gives the argv shown, if there is an input. */
  protected fastaInputName(state: S): string | undefined {
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

  /** Whether the FASTA input or the files of the run hold anything that a reload would lose. */
  protected holdsFastaInputs(state: S): boolean {
    return state.fastaPaste.trim() !== '' || state.fastaFiles.length > 0 || state.runFiles.length > 0;
  }

  /**
   * The inputs of the runs of a submit: the `-in` name (a BLAST database, whose files are files of
   * the run), or the FASTA inputs: one, one for each (Separate), or `combined_input.fa` (Combined).
   * Throws the message to show when the inputs are not ready.
   */
  protected async runInputs(state: S): Promise<readonly FastaRunInput[]> {
    if (!this.usesFastaBox(state)) {
      const name = state.blastDbIn.trim();
      if (name === '') throw new Error('Enter an input name or open files.');
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
        ...(entry.sessionToken === undefined ? {} : { sessionToken: entry.sessionToken }),
      });
    }
    if (parts.length === 0) {
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

  /** The `-in` name that files added to the run directory name, until the user edits it (a BLAST database). */
  protected override runFilesAdded(state: S, runFiles: readonly ToolExtraFile[]): StatePatch<S> {
    if (state.blastDbInEdited && state.blastDbIn.trim() !== '') return {};
    const detected = primaryBlastDbName(runFiles.map((entry) => entry.name));
    return detected === undefined ? {} : ({ blastDbIn: detected } as StatePatch<S>);
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
