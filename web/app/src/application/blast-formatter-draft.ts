// The Format Archive (`blast_formatter`) form draft (RG-F3p7): independent of queued and completed
// runs, on the tool drafts' shared core (`ToolDraftCore`: `-xmlhelp`, files of the run, check, submit).
// What is its own: the archive of `-archive` (a file, pasted ASN.1 text, or the `-outfmt 11` output
// of a search of this session) and `-outfmt` / `-out`.
import type { DirectoryFile } from 'losat/domain';
import { splitOptions } from '../domain/argv';
import { normalizeRunPath } from '../domain/blast-db';
import type { FieldValue, FormValues } from '../domain/parameters';
import { buildBlastFormatterArgv } from '../domain/tool-argv';
import { BLAST_FORMATTER_DEFAULT_OUTFMT, blastFormatterOptionParameters, toolPrintsKind } from '../domain/tools';
import type { SessionArchiveEntry, ToolDirectoryFile, ToolRequest } from './coordinator';
import { ToolDraftCore, type ToolDraftCoreDeps, type ToolFormState } from './tool-draft-core';

export type ArchiveSourceKind = 'file' | 'paste' | 'session';

export interface ArchivePickedFile {
  readonly file: File;
  readonly name: string;
  readonly size: number;
}

export interface BlastFormatterDraftState extends ToolFormState {
  readonly archiveSource: ArchiveSourceKind;
  readonly archiveFile?: ArchivePickedFile;
  readonly archivePaste: string;
  readonly sessionArchive?: SessionArchiveEntry;
  /** Output format string (`-outfmt`), e.g. `'8'`, `'9'`, `'12'`, `'6 qseqid sseqid'`. */
  readonly outfmt: string;
  /** Output file name (`-out`), empty or `'-'` for standard output. */
  readonly outName: string;
  /** Formatting option values keyed by `-flag` (`-max_target_seqs`, `-parse_deflines`, etc.). */
  readonly values: FormValues;
  /** Words typed in "Additional options". */
  readonly extra: string;
}

export interface BlastFormatterDraftDeps extends ToolDraftCoreDeps {
  /** Reads the `-outfmt 11` output of a search of this session as the archive file. */
  readonly getSessionArchiveFile: (runId: string, fileName: string) => Promise<File>;
}

const PASTED_ARCHIVE_NAME = 'archive.asn';

/** The archive that `-archive` names, as the run directory holds it; the bytes of a session's archive are read at the run only. */
type ArchiveInput = { readonly name: string; readonly origin: ToolDirectoryFile['origin'] } & (
  | { readonly file: Blob }
  | { readonly runId: string }
);

export class BlastFormatterDraft extends ToolDraftCore<BlastFormatterDraftState> {
  constructor(private readonly deps: BlastFormatterDraftDeps) {
    super('blast_formatter', deps, {
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

  // --- the archive and the output -------------------------------------------------------------

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
      this.reportError(`Invalid file name: "${file.name}".`);
      return;
    }
    this.set({ archiveSource: 'file', archiveFile: { file, name, size: file.size }, message: undefined });
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

  /** "Format archive…" on the `-outfmt 11` output of a completed search: that output becomes the archive. */
  useSessionArchive(entry: SessionArchiveEntry): void {
    this.set({ archiveSource: 'session', sessionArchive: entry, message: undefined });
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
    if (value === undefined || value === '' || value === false) delete next[flag];
    else next[flag] = value;
    this.set({ values: Object.freeze(next), message: undefined });
    this.scheduleValidation();
  }

  setExtra(text: string): void {
    this.set({ extra: text, message: undefined });
    this.scheduleValidation();
  }

  /** Resets the options of "Formatting options" and Additional options. */
  resetOptions(): void {
    const { values, description } = this.state.get();
    const next: Record<string, FieldValue> = { ...values };
    for (const param of blastFormatterOptionParameters(description)) delete next[param.flag];
    this.set({ values: Object.freeze(next), extra: '', message: undefined });
    this.scheduleValidation();
  }

  // --- the command line -----------------------------------------------------------------------

  /** The argv shown beside Run LOSAT; a quote left open in Additional options leaves their words out (the check reports it). */
  previewArgv(state: BlastFormatterDraftState = this.state.get()): readonly string[] {
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

  holdsInputs(state: BlastFormatterDraftState = this.state.get()): boolean {
    return state.archiveFile !== undefined || state.archivePaste.trim() !== '' || state.runFiles.length > 0;
  }

  // --- what the core asks -------------------------------------------------------------------

  protected parseExtra(state: BlastFormatterDraftState): readonly string[] {
    return splitOptions(state.extra);
  }

  /**
   * The argv of the check; none before an archive is chosen, unless the words name `-archive` or
   * `-rid` themselves or print help or the version (NCBI's error then would be `-rid`'s missing value).
   */
  protected validationArgv(state: BlastFormatterDraftState, extra: readonly string[]): readonly string[] | undefined {
    const argv = this.argvOf(state, extra);
    if (this.archiveInput(state) !== undefined || toolPrintsKind(argv) !== null) return argv;
    return argv.some((word) => /^-(?:archive|rid)(?:=|$)/.test(word)) ? argv : undefined;
  }

  protected override validationInputs(state: BlastFormatterDraftState): readonly DirectoryFile[] {
    const archive = this.archiveInput(state);
    // The check does not open -archive: a session's archive stands in by its name.
    return archive === undefined ? [] : [{ name: archive.name, file: 'file' in archive ? archive.file : new Blob([]) }];
  }

  protected async buildRequests(state: BlastFormatterDraftState, extra: readonly string[]): Promise<readonly ToolRequest[]> {
    const argv = this.argvOf(state, extra);
    const prints = toolPrintsKind(argv);
    if (prints !== null) {
      return [{ tool: 'blast_formatter', stage: 'format', argv, directory: this.runDirectoryOf(state, []), prints }];
    }
    const archive = this.archiveInput(state);
    if (archive === undefined) throw new Error(MISSING_ARCHIVE[state.archiveSource]);
    const file = 'file' in archive ? archive.file : await this.deps.getSessionArchiveFile(archive.runId, archive.name);
    const directory = this.runDirectoryOf(state, [{ name: archive.name, file, origin: archive.origin }]);
    return [{ tool: 'blast_formatter', stage: 'format', argv, directory }];
  }

  // --- internals ----------------------------------------------------------------------------

  private argvOf(state: BlastFormatterDraftState, extra: readonly string[]): readonly string[] {
    const prints = toolPrintsKind(['blast_formatter', ...extra]);
    const archive = prints === null ? this.archiveInput(state) : undefined;
    return buildBlastFormatterArgv({
      ...(archive === undefined ? {} : { archiveName: archive.name }),
      outfmt: state.outfmt,
      outName: state.outName,
      values: state.values,
      ...(state.description === undefined ? {} : { description: state.description }),
      extra,
    });
  }

  /** The archive of the source chosen, if there is one. */
  private archiveInput(state: BlastFormatterDraftState): ArchiveInput | undefined {
    switch (state.archiveSource) {
      case 'file':
        return state.archiveFile === undefined
          ? undefined
          : { name: state.archiveFile.name, file: state.archiveFile.file, origin: { kind: 'picked' } };
      case 'paste':
        return state.archivePaste.trim() === ''
          ? undefined
          : { name: PASTED_ARCHIVE_NAME, file: new Blob([state.archivePaste], { type: 'text/plain' }), origin: { kind: 'pasted' } };
      case 'session':
        return state.sessionArchive === undefined
          ? undefined
          : {
              name: state.sessionArchive.name,
              runId: state.sessionArchive.runId,
              origin: { kind: 'archive', runNumber: state.sessionArchive.runNumber },
            };
    }
  }
}

/** What Run LOSAT says when the chosen source has no archive yet. */
const MISSING_ARCHIVE: Readonly<Record<ArchiveSourceKind, string>> = {
  file: 'Choose an archive file.',
  paste: 'Paste an archive.',
  session: 'Choose a search run.',
};
