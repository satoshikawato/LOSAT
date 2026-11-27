// The Format Archive (`blast_formatter`) form draft (RG-F3p7): independent of queued and
// completed runs, mirroring ToolDraft's architecture. Re-formats `-outfmt 11` BLAST archives
// from a local file, pasted ASN.1 text, or a completed `-outfmt 11` search output of this
// session, with optional run directory files (e.g. a local BLAST database folder for `-db`
// archives, or `-conffile`).
import { type DirectoryFile, type ToolDescription } from 'losat/domain';
import { isUnavailable, type ToolEngine } from 'losat/ports';
import { splitOptions } from '../domain/argv';
import type { FieldValue, FormValues } from '../domain/parameters';
import { buildBlastFormatterArgv, toolCommandWords, toolValidationSummary } from '../domain/tool-argv';
import {
  BLAST_FORMATTER_DEFAULT_OUTFMT,
  blastFormatterOptionParameters,
  toolPrintsKind,
} from '../domain/tools';
import { failureText } from './availability';
import type { EnqueueAllResult, SessionArchiveEntry, ToolDirectoryFile, ToolRequest } from './coordinator';
import { Store } from './store';

export interface BlastFormatterRunFile {
  readonly key: string;
  readonly name: string;
  readonly size: number;
  readonly file: File;
}

export type BlastFormatterValidationState =
  | { readonly state: 'idle' }
  | { readonly state: 'validating' }
  | { readonly state: 'valid' }
  | { readonly state: 'invalid'; readonly message: string }
  | { readonly state: 'unavailable'; readonly message: string };

export type ArchiveSourceKind = 'file' | 'paste' | 'session';

export interface ArchivePickedFile {
  readonly file: File;
  readonly name: string;
  readonly size: number;
}

export interface BlastFormatterDraftState {
  readonly archiveSource: ArchiveSourceKind;
  readonly archiveFile?: ArchivePickedFile;
  readonly archivePaste: string;
  readonly sessionArchive?: SessionArchiveEntry;
  /** Output format string (`-outfmt`), e.g. `'8'`, `'9'`, `'12'`, `'6 qseqid sseqid'`, or `'0'`. */
  readonly outfmt: string;
  /** Output file name (`-out`), empty or `'-'` for standard output. */
  readonly outName: string;
  /** Formatting option values keyed by `-flag` (`-max_target_seqs`, `-parse_deflines`, etc.). */
  readonly values: FormValues;
  /** Words typed in "Additional options". */
  readonly extra: string;
  /** Additional files placed in the run directory (local BLAST database files for `-db` archives, `-conffile`, etc.). */
  readonly runFiles: readonly BlastFormatterRunFile[];
  readonly runFilesError?: string;
  readonly description?: ToolDescription;
  readonly descriptionError?: string;
  readonly descriptionUnavailable?: boolean;
  readonly validation: BlastFormatterValidationState;
  readonly submitting: boolean;
  readonly message?: { readonly kind: 'error' | 'info'; readonly text: string };
}

type StatePatch = { [K in keyof BlastFormatterDraftState]?: BlastFormatterDraftState[K] | undefined };

export interface BlastFormatterDraftDeps {
  readonly engine: ToolEngine;
  readonly enqueueToolAll: (requests: readonly ToolRequest[]) => Promise<EnqueueAllResult>;
  readonly getSessionArchiveFile: (runId: string, fileName: string) => Promise<File>;
  readonly isBusy?: () => boolean;
  readonly debounceMs?: number;
}

/** True when the Format Archive tab holds user inputs that reloading the page would lose. */
export function holdsBlastFormatterInputs(state: BlastFormatterDraftState): boolean {
  return (
    state.archiveFile !== undefined ||
    state.archivePaste.trim() !== '' ||
    state.runFiles.length > 0
  );
}

function normalizeRunPath(raw: string): string | undefined {
  const cleaned = raw.replaceAll('\\', '/').replace(/^\/+/, '').replace(/\/+$/, '');
  if (cleaned === '' || cleaned.includes('\0')) return undefined;
  const parts = cleaned.split('/').filter((segment) => segment !== '' && segment !== '.');
  if (parts.length === 0 || parts.some((segment) => segment === '..')) return undefined;
  return parts.join('/');
}

function fileRunPath(file: File): string | undefined {
  const relative = (file as File & { readonly webkitRelativePath?: string }).webkitRelativePath;
  return normalizeRunPath(relative && relative.trim() !== '' ? relative : file.name);
}

const PASTED_ARCHIVE_NAME = 'archive.asn';

export class BlastFormatterDraft {
  readonly state: Store<BlastFormatterDraftState>;
  private readonly debounceMs: number;
  private describePending: Promise<void> | undefined;
  private validateTimer: ReturnType<typeof setTimeout> | undefined;
  private validateToken = 0;
  private nextKey = 1;

  constructor(private readonly deps: BlastFormatterDraftDeps) {
    this.debounceMs = deps.debounceMs ?? 200;
    this.state = new Store<BlastFormatterDraftState>({
      archiveSource: 'file',
      archivePaste: '',
      outfmt: BLAST_FORMATTER_DEFAULT_OUTFMT,
      outName: '',
      values: {},
      extra: '',
      runFiles: Object.freeze([]),
      validation: { state: 'idle' },
      submitting: false,
    });
  }

  /** Loads `blast_formatter -xmlhelp` lazily when the Format Archive tab is first opened. */
  ensureLoaded(): void {
    const { description, descriptionError } = this.state.get();
    if (description !== undefined || descriptionError !== undefined || this.describePending !== undefined) return;
    void this.loadDescription();
  }

  /** Re-asks the engine for anything that failed with an UnavailableError (`Availability.recheck`). */
  async askAgain(): Promise<void> {
    const state = this.state.get();
    if (state.descriptionUnavailable === true) await this.loadDescription();
    if (this.state.get().validation.state === 'unavailable') await this.runValidation();
  }

  private loadDescription(): Promise<void> {
    const pending = (async () => {
      try {
        const description = await this.deps.engine.describeTool('blast_formatter');
        this.set({
          description,
          descriptionError: undefined,
          descriptionUnavailable: undefined,
        });
        this.scheduleValidation();
      } catch (error) {
        this.set({
          descriptionError: error instanceof Error ? error.message : String(error),
          descriptionUnavailable: isUnavailable(error),
        });
      }
    })();
    this.describePending = pending;
    void pending.finally(() => {
      if (this.describePending === pending) this.describePending = undefined;
    });
    return pending;
  }

  // --- Archive source & output settings ------------------------------------------------------

  setArchiveSource(source: ArchiveSourceKind): void {
    if (this.state.get().archiveSource === source) return;
    this.set({ archiveSource: source, message: undefined });
    this.scheduleValidation();
  }

  setArchiveFile(file: File | undefined): void {
    if (file === undefined) {
      this.set({ archiveFile: undefined, message: undefined });
      this.scheduleValidation();
      return;
    }
    const name = normalizeRunPath(file.name);
    if (name === undefined) {
      this.set({ message: { kind: 'error', text: `Invalid file name: "${file.name}".` } });
      return;
    }
    this.set({
      archiveSource: 'file',
      archiveFile: { file, name, size: file.size },
      message: undefined,
    });
    this.scheduleValidation();
  }

  setArchivePaste(text: string): void {
    this.set({ archiveSource: 'paste', archivePaste: text, message: undefined });
    this.scheduleValidation();
  }

  selectSessionArchive(entry: SessionArchiveEntry | undefined): void {
    this.set({ sessionArchive: entry, message: undefined });
    this.scheduleValidation();
  }

  /** Called from the `"Format archive…"` shortcut on a completed `-outfmt 11` search output. */
  useSessionArchive(entry: SessionArchiveEntry): void {
    this.set({
      archiveSource: 'session',
      sessionArchive: entry,
      message: undefined,
    });
    this.ensureLoaded();
    this.scheduleValidation();
  }

  setOutfmt(value: string): void {
    this.set({ outfmt: value, message: undefined });
    this.scheduleValidation();
  }

  setOutName(value: string): void {
    this.set({ outName: value, message: undefined });
    this.scheduleValidation();
  }

  setField(flag: string, value: FieldValue): void {
    const next: Record<string, FieldValue> = { ...this.state.get().values };
    if (value === '' || value === false) delete next[flag];
    else next[flag] = value;
    this.set({ values: next, message: undefined });
    this.scheduleValidation();
  }

  setExtra(text: string): void {
    this.set({ extra: text, message: undefined });
    this.scheduleValidation();
  }

  resetOptions(): void {
    const { values, description } = this.state.get();
    const optionFlags = new Set(blastFormatterOptionParameters(description).map((param) => param.flag));
    const next: Record<string, FieldValue> = { ...values };
    for (const flag of optionFlags) delete next[flag];
    this.set({
      values: next,
      extra: '',
      message: undefined,
    });
    this.scheduleValidation();
  }

  // --- Run directory files (for local `-db` archives & `-conffile`) --------------------------

  addRunFiles(files: readonly File[]): void {
    if (files.length === 0) return;
    const existing = new Map(this.state.get().runFiles.map((entry) => [entry.name, entry]));
    let error: string | undefined;
    for (const file of files) {
      const name = fileRunPath(file);
      if (name === undefined) {
        error = `Invalid file name: "${file.name}".`;
        continue;
      }
      const key = `bf-rf-${this.nextKey++}`;
      existing.set(name, { key, name, size: file.size, file });
    }
    this.set({
      runFiles: Object.freeze([...existing.values()]),
      runFilesError: error,
      message: undefined,
    });
    this.scheduleValidation();
  }

  removeRunFile(key: string): void {
    this.set({
      runFiles: Object.freeze(this.state.get().runFiles.filter((entry) => entry.key !== key)),
      runFilesError: undefined,
      message: undefined,
    });
    this.scheduleValidation();
  }

  clearRunFiles(): void {
    this.set({
      runFiles: Object.freeze([]),
      runFilesError: undefined,
      message: undefined,
    });
    this.scheduleValidation();
  }

  // --- Argv preview & validation -------------------------------------------------------------

  archiveNameOf(state: BlastFormatterDraftState = this.state.get()): string | undefined {
    if (state.archiveSource === 'file') return state.archiveFile?.name;
    if (state.archiveSource === 'paste') return state.archivePaste.trim() !== '' ? PASTED_ARCHIVE_NAME : undefined;
    return state.sessionArchive?.name;
  }

  private extraWordsOf(state: BlastFormatterDraftState): { readonly ok: true; readonly words: readonly string[] } | { readonly ok: false; readonly message: string } {
    if (state.extra.trim() === '') return { ok: true, words: [] };
    try {
      return { ok: true, words: splitOptions(state.extra) };
    } catch (error) {
      return { ok: false, message: `Additional options: ${error instanceof Error ? error.message : String(error)}` };
    }
  }

  private argvInputOf(state: BlastFormatterDraftState, archiveName: string | undefined, extra?: readonly string[]) {
    return {
      ...(archiveName === undefined ? {} : { archiveName }),
      outfmt: state.outfmt,
      outName: state.outName,
      values: state.values,
      ...(state.description === undefined ? {} : { description: state.description }),
      ...(extra === undefined ? {} : { extra }),
    };
  }

  private buildArgv(state: BlastFormatterDraftState): { readonly ok: true; readonly argv: readonly string[] } | { readonly ok: false; readonly message: string } {
    const extra = this.extraWordsOf(state);
    if (!extra.ok) return extra;
    const prints = toolPrintsKind(['blast_formatter', ...extra.words]);
    const archiveName = prints !== null ? undefined : this.archiveNameOf(state);
    try {
      return { ok: true, argv: buildBlastFormatterArgv(this.argvInputOf(state, archiveName, extra.words)) };
    } catch (error) {
      return { ok: false, message: error instanceof Error ? error.message : String(error) };
    }
  }

  previewCommandWords(state: BlastFormatterDraftState = this.state.get()): readonly string[] {
    const built = this.buildArgv(state);
    if (built.ok) return toolCommandWords(built.argv);
    return toolCommandWords(buildBlastFormatterArgv(this.argvInputOf(state, this.archiveNameOf(state))));
  }

  private nextValidateToken(): number {
    this.validateToken += 1;
    return this.validateToken;
  }

  private scheduleValidation(): void {
    if (this.validateTimer !== undefined) clearTimeout(this.validateTimer);
    const token = this.nextValidateToken();
    this.validateTimer = setTimeout(() => {
      this.validateTimer = undefined;
      if (token === this.validateToken) void this.runValidation();
    }, this.debounceMs);
  }

  private async runValidation(): Promise<void> {
    const state = this.state.get();
    const token = this.nextValidateToken();
    const checked = this.buildArgv(state);
    if (!checked.ok) {
      this.set({ validation: { state: 'invalid', message: checked.message } });
      return;
    }
    const prints = toolPrintsKind(checked.argv);
    // Avoid flashing "Argument rid has no value" before the user has provided an archive.
    if (prints === null && this.archiveNameOf(state) === undefined && !checked.argv.some((w) => w === '-archive' || w.startsWith('-archive=') || w === '-rid' || w.startsWith('-rid='))) {
      this.set({ validation: { state: 'idle' } });
      return;
    }
    if (this.deps.isBusy?.() === true) return;
    this.set({ validation: { state: 'validating' } });
    const directory: DirectoryFile[] = [];
    for (const entry of state.runFiles) directory.push({ name: entry.name, file: entry.file });
    if (state.archiveSource === 'file' && state.archiveFile !== undefined) {
      directory.push({ name: state.archiveFile.name, file: state.archiveFile.file });
    } else if (state.archiveSource === 'paste' && state.archivePaste.trim() !== '') {
      directory.push({ name: PASTED_ARCHIVE_NAME, file: new Blob([state.archivePaste], { type: 'text/plain' }) });
    } else if (state.archiveSource === 'session' && state.sessionArchive !== undefined) {
      directory.push({ name: state.sessionArchive.name, file: new Blob([]) });
    }
    let nextValidation: BlastFormatterValidationState;
    try {
      const verdict = await this.deps.engine.validateTool(checked.argv, { directory });
      nextValidation = verdict.ok ? { state: 'valid' } : { state: 'invalid', message: verdict.message };
    } catch (error) {
      const message = error instanceof Error ? error.message : String(error);
      nextValidation = isUnavailable(error) ? { state: 'unavailable', message } : { state: 'invalid', message };
    }
    if (token !== this.validateToken) return;
    this.set({ validation: nextValidation });
  }

  // --- Submit --------------------------------------------------------------------------------

  private failSubmit(text: string, validation?: BlastFormatterValidationState): EnqueueAllResult {
    this.set({
      ...(validation === undefined ? {} : { validation }),
      message: { kind: 'error', text },
    });
    return { ok: false, reason: 'argv', message: text };
  }

  async submit(): Promise<EnqueueAllResult> {
    const state = this.state.get();
    if (state.submitting) return { ok: false, reason: 'argv', message: 'Submission already in progress.' };
    this.set({ submitting: true, message: undefined });
    try {
      const result = this.buildArgv(state);
      if (!result.ok) {
        return this.failSubmit(result.message, { state: 'invalid', message: result.message });
      }
      const printed = toolPrintsKind(result.argv);
      const directory: ToolDirectoryFile[] = [];
      const seen = new Set<string>();

      if (printed === null) {
        if (state.archiveSource === 'file') {
          if (state.archiveFile === undefined) {
            return this.failSubmit('Choose a BLAST archive file (-outfmt 11), or paste one.');
          }
          directory.push({ name: state.archiveFile.name, file: state.archiveFile.file, origin: { kind: 'picked' } });
          seen.add(state.archiveFile.name);
        } else if (state.archiveSource === 'paste') {
          if (state.archivePaste.trim() === '') {
            return this.failSubmit('Paste a BLAST archive (-outfmt 11), or choose a file.');
          }
          directory.push({
            name: PASTED_ARCHIVE_NAME,
            file: new Blob([state.archivePaste], { type: 'text/plain' }),
            origin: { kind: 'pasted' },
          });
          seen.add(PASTED_ARCHIVE_NAME);
        } else {
          if (state.sessionArchive === undefined) {
            return this.failSubmit('Choose a completed -outfmt 11 search output from this session, or choose a file.');
          }
          const archiveFile = await this.deps.getSessionArchiveFile(state.sessionArchive.runId, state.sessionArchive.name);
          directory.push({
            name: state.sessionArchive.name,
            file: archiveFile,
            origin: { kind: 'archive', runNumber: state.sessionArchive.runNumber },
          });
          seen.add(state.sessionArchive.name);
        }
      }

      for (const entry of state.runFiles) {
        if (seen.has(entry.name)) {
          return this.failSubmit(`Run directory file "${entry.name}" conflicts with the archive file name.`);
        }
        seen.add(entry.name);
        directory.push({ name: entry.name, file: entry.file, origin: { kind: 'picked' } });
      }

      const verdict = await this.deps.engine.validateTool(result.argv, {
        directory: directory.map((entry) => ({ name: entry.name, file: entry.file })),
      });
      if (!verdict.ok) {
        this.set({
          validation: { state: 'invalid', message: verdict.message },
          message: { kind: 'error', text: toolValidationSummary(verdict.message) },
        });
        return { ok: false, reason: 'argv', message: verdict.message };
      }

      const request: ToolRequest = {
        tool: 'blast_formatter',
        stage: 'format',
        argv: result.argv,
        directory,
        ...(printed === null ? {} : { prints: printed }),
      };
      const queued = await this.deps.enqueueToolAll([request]);
      if (!queued.ok) {
        this.set({ message: { kind: 'error', text: queued.message } });
      } else {
        this.set({ validation: { state: 'valid' }, message: undefined });
      }
      return queued;
    } catch (error) {
      return this.failSubmit(failureText(error));
    } finally {
      this.set({ submitting: false });
    }
  }

  private set(patch: StatePatch): void {
    const next: Record<string, unknown> = { ...this.state.get() };
    for (const [key, value] of Object.entries(patch)) {
      if (value === undefined) delete next[key];
      else next[key] = value;
    }
    this.state.set(next as unknown as BlastFormatterDraftState);
  }
}
