// The files of reproduction (S15 instructions, items 3 and 5; design §12.3): settings files of
// the search form and of a run, "Edit Search" (a run's settings and Job Title put back in the
// search form, W4b decision 18), and the input FASTA of a run - the exact bytes that the engine
// searched, named as in the run's argv, so that its commands run as written. Every file is
// written through the Writer contract (export-writer.ts) in bounded blocks. Nothing here
// starts a search, and nothing reads the inputs' files again.
import { PASTED_NAMES } from '../domain/argv';
import type { InputRole } from '../domain/programs';
import { programById } from '../domain/programs';
import type { InputPart } from '../domain/reproduce';
import {
  parseSettings,
  settingsFileName,
  settingsOfRun,
  settingsText,
  SETTINGS_MAX_BYTES,
  SETTINGS_MIME,
  type AppliedSettings,
  type SearchSettings,
} from '../domain/settings-file';
import type { Downloader } from '../ports/download';
import type { RunView } from './coordinator';
import type { SearchDraft } from './draft';
import { BLOCK_CHARS, writeFile } from './export-writer';
import { Store } from './store';

/** The MIME type of a saved FASTA input (as the candidate tray's FASTA). */
const FASTA_MIME = 'text/plain';
/** Bytes per block of a saved input: the Writer's block size. */
export const INPUT_SLICE_BYTES = BLOCK_CHARS;

export interface FilesMessage {
  readonly kind: 'info' | 'error';
  readonly text: string;
}

export interface RunFilesState {
  /** What the latest settings file or "Edit Search" put in the search form, or why a file was refused. */
  readonly settings?: FilesMessage;
  /** The latest file saved from a run's reproduction panel, or why it was not saved. */
  readonly run?: FilesMessage & { readonly runId: string };
}

export interface RunFilesDeps {
  readonly draft: Pick<SearchDraft, 'settings' | 'applySettings' | 'state'>;
  readonly downloader: Pick<Downloader, 'open'>;
  /** The threads that the search form offers in this browser (domain/settings-file.ts `threadLimit`). */
  readonly maxThreads: () => number;
  /**
   * The run input of a role of a run loaded from a session file, rebuilt from the original FASTA
   * that was chosen again and matched in Run details (application/session.ts `attachedInput`),
   * with its SHA-256; undefined while none is attached.
   */
  readonly attachedInput?: (runId: string, role: InputRole) => Promise<{ readonly bytes: Uint8Array; readonly sha256: string } | undefined>;
}

export class RunFiles {
  readonly state = new Store<RunFilesState>({});
  /**
   * The record tables of the search form's sources, by their revision ID: a run input made of
   * such a revision is the whole source. The form makes a new revision for a source with
   * records left out (SearchDraft `revisionFor`), so a run's revision that is not here is a
   * selection of a source's records.
   */
  private readonly sources = new Map<string, InputPart>();

  constructor(private readonly deps: RunFilesDeps) {
    const learn = () => {
      const state = deps.draft.state.get();
      for (const role of ['query', 'subject'] as const) {
        for (const source of state[role].sources) {
          if (source.base === undefined || this.sources.has(source.base.revisionId)) continue;
          const name = source.origin === 'paste' ? PASTED_NAMES[role] : source.name;
          this.sources.set(source.base.revisionId, { origin: source.origin, name, records: source.base.records.length });
        }
      }
    };
    learn();
    deps.draft.state.subscribe(learn);
  }

  // --- settings files ---------------------------------------------------------------------

  /** Saves the search form's conditions as a settings file ("Save settings"). */
  async saveSettings(): Promise<void> {
    const settings = this.deps.draft.settings();
    try {
      await this.writeSettings(settings);
      this.setSettings({ kind: 'info', text: `Saved ${settingsFileName(settings.program)}: the program, the options and the threads.` });
    } catch (error) {
      this.setSettings({ kind: 'error', text: `The settings could not be saved: ${messageOf(error)}` });
    }
  }

  /**
   * Reads a settings file into the search form ("Load settings…"). A file that is not a
   * settings file of this LOSAT Web is refused and the form does not change.
   */
  async loadSettings(file: File): Promise<void> {
    const refuse = (reason: string) => this.setSettings({ kind: 'error', text: `${file.name} was not loaded: ${reason}` });
    // The size first, so that a large file is not read at all.
    if (file.size > SETTINGS_MAX_BYTES) {
      refuse(`it has ${file.size} bytes; a settings file has at most ${SETTINGS_MAX_BYTES} (1 MiB).`);
      return;
    }
    let bytes: Uint8Array;
    try {
      bytes = new Uint8Array(await file.arrayBuffer());
    } catch (error) {
      refuse(`it could not be read (${messageOf(error)}).`);
      return;
    }
    const read = parseSettings(bytes);
    if (!read.ok) {
      refuse(read.message);
      return;
    }
    const { settings } = read;
    const descriptor = programById(settings.program);
    if (descriptor.unavailable !== undefined) {
      refuse(`its program is ${descriptor.label}. ${descriptor.unavailable}`);
      return;
    }
    let applied: AppliedSettings;
    try {
      applied = await this.deps.draft.applySettings(settings, { maxThreads: this.deps.maxThreads() });
    } catch (error) {
      refuse(messageOf(error));
      return;
    }
    const words = settings.options.length;
    const threads = settings.threads === 'auto' ? 'Auto' : String(settings.threads);
    this.setSettings({
      kind: 'info',
      text:
        `Loaded ${file.name}: ${descriptor.label}, ${words === 0 ? 'the default options' : `${words} ${words === 1 ? 'word' : 'words'} of options`}, ` +
        `threads ${threads}. The inputs and the Job Title did not change.${notAppliedText(applied)}`,
    });
  }

  /** Saves the settings of a run, as its argv has them ("Save the settings of this run"). */
  async saveRunSettings(view: RunView): Promise<void> {
    const settings = settingsOfRun(view.snapshot);
    try {
      const name = settingsFileName(settings.program);
      await this.writeSettings(settings);
      this.setRun(view, { kind: 'info', text: `Saved ${name}: the settings of Run ${view.snapshot.number}.` });
    } catch (error) {
      this.setRun(view, { kind: 'error', text: `The settings could not be saved: ${messageOf(error)}` });
    }
  }

  /**
   * NCBI's "Edit Search": puts the run's settings and Job Title in the search form, as a
   * settings file would; the form keeps its own inputs and nothing is searched. Resolves
   * whether the form now has them (the caller then shows the Search tab).
   */
  async editSearch(view: RunView): Promise<boolean> {
    const { snapshot } = view;
    const descriptor = programById(snapshot.program);
    if (descriptor.unavailable !== undefined) {
      this.setRun(view, { kind: 'error', text: descriptor.unavailable });
      return false;
    }
    try {
      const applied = await this.deps.draft.applySettings(settingsOfRun(snapshot), {
        maxThreads: this.deps.maxThreads(),
        title: snapshot.title ?? '',
      });
      this.setSettings({
        kind: 'info',
        text: `The search form has the settings of Run ${snapshot.number}. The inputs are the form's own.${notAppliedText(applied)}`,
      });
      return true;
    } catch (error) {
      this.setRun(view, { kind: 'error', text: `The settings of Run ${snapshot.number} could not be put in the form: ${messageOf(error)}` });
      return false;
    }
  }

  // --- the input FASTA of a run -----------------------------------------------------------

  /** The name of the run's input in its argv: the file that its commands read. */
  inputName(view: RunView, role: InputRole): string {
    return view.snapshot[role].name;
  }

  /**
   * Where the run's input came from, part by part (domain/reproduce.ts `inputRelation`): a
   * whole source of the search form, or undefined for a selection of a source's records. A run
   * loaded from a session file has no revisions here; its sources are those that the file
   * recorded, with the records left out of each (code review L2). The file does not say which
   * source was pasted text, so each is named as the file that the input was given as.
   */
  inputParts(view: RunView, role: InputRole): ReadonlyArray<InputPart | undefined> {
    if (view.fromSession !== undefined) {
      return view.fromSession.inputs[role].sources.map(
        (source): InputPart => ({ origin: 'file', name: source.name, records: source.records, excluded: source.excluded.length }),
      );
    }
    return view.snapshot[role].revisionIds.map((revisionId) => this.sources.get(revisionId));
  }

  /** Saves the bytes that the engine searched for a role, under the argv's name, in bounded slices. */
  async saveInput(view: RunView, role: InputRole): Promise<void> {
    const name = this.inputName(view, role);
    try {
      const bytes = await this.inputBytes(view, role);
      await writeFile(this.deps.downloader, name, FASTA_MIME, async (writer) => {
        for (let at = 0; at < bytes.length; at += INPUT_SLICE_BYTES) await writer.bytes(bytes.subarray(at, at + INPUT_SLICE_BYTES));
      });
      this.setRun(view, { kind: 'info', text: `Saved ${name}: the ${role} that Run ${view.snapshot.number} searched.` });
    } catch (error) {
      this.setRun(view, { kind: 'error', text: `${name} was not saved: ${messageOf(error)}` });
    }
  }

  /**
   * The bytes that the engine searched for a role: the snapshot's, or for a run loaded from a
   * session file, the run input rebuilt from its attached original FASTA. A loaded run without
   * one is refused, never guessed; so is a rebuilt input whose SHA-256 is not the one that the
   * session file recorded (the original changed after it was attached; code review L1). The
   * rebuilt input is held whole in memory until it is written, as a search's snapshot holds the
   * bytes that it gives the engine.
   */
  private async inputBytes(view: RunView, role: InputRole): Promise<Uint8Array> {
    const bytes = view.snapshot[role].bytes;
    if (bytes !== undefined) return bytes;
    const attached = view.fromSession === undefined ? undefined : await this.deps.attachedInput?.(view.snapshot.runId, role);
    if (view.fromSession === undefined || attached === undefined) {
      throw new Error(
        `Run ${view.snapshot.number} was loaded from a session file; choose its original ${role} FASTA in Run details to save the input it searched.`,
      );
    }
    const recorded = view.fromSession.inputs[role].sha256;
    if (attached.sha256 !== recorded) {
      throw new Error(
        `the ${role} FASTA attached to Run ${view.snapshot.number} no longer makes the input that the run searched (SHA-256 ${attached.sha256}, ` +
          `not ${recorded}); the file changed after it was attached. Choose the original ${role} FASTA again in Run details.`,
      );
    }
    return attached.bytes;
  }

  // --- helpers ----------------------------------------------------------------------------

  private writeSettings(settings: SearchSettings): Promise<number> {
    return writeFile(this.deps.downloader, settingsFileName(settings.program), SETTINGS_MIME, async (writer) => {
      for (const part of settingsText(settings)) await writer.text(part);
    });
  }

  private setSettings(message: FilesMessage): void {
    const { run } = this.state.get();
    this.state.set(run === undefined ? { settings: message } : { settings: message, run });
  }

  private setRun(view: RunView, message: FilesMessage): void {
    const { settings } = this.state.get();
    const run = { ...message, runId: view.snapshot.runId };
    this.state.set(settings === undefined ? { run } : { settings, run });
  }
}

/** " Not applied: …." after a message, or nothing. */
function notAppliedText(applied: AppliedSettings): string {
  return applied.notApplied.length === 0 ? '' : ` Not applied: ${applied.notApplied.join('; ')}.`;
}

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
