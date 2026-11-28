// The parts of a tool tab's form draft that do not depend on the tool (WindowMasker's `ToolDraft`,
// Format Archive's `BlastFormatterDraft`, makeblastdb's `MakeblastdbDraft`; RG-F3p7 review R1 #3): the
// tool's `-xmlhelp`, the files of the run (refusing a second file under a name already there; a file
// chosen for an option, which the option then names), the check of the argv by
// `validateTool` (debounced, one at a time, an older answer dropped; it runs while the engine is busy,
// behind the run), the run directory (one builder for the check and the run), and the submit that
// waits for the check before it queues the runs on the shared Coordinator queue. A tool supplies its
// inputs and its argv through the abstract members.
import type { DirectoryFile, ToolDescription } from 'losat/domain';
import { isUnavailable, type ToolEngine, type ToolValidation } from 'losat/ports';
import { duplicateRunFileMessage, normalizeRunPath, pickedFilePath, RunDirectoryBuilder } from '../domain/blast-db';
import type { FieldValue } from '../domain/parameters';
import type { ToolFileOrigin } from '../domain/run';
import { listNameRefusal } from '../domain/list-flags';
import { fileListNames, toolCommandWords, toolValidationSummary } from '../domain/tool-argv';
import { TOOL_FILE_LIST_FLAGS, type ToolId } from '../domain/tools';
import { failureText } from './availability';
import type { EnqueueAllResult, ToolDirectoryFile, ToolRequest } from './coordinator';
import type { SubmitResult, Validation } from './draft';
import type { SessionFiles } from './session-files';
import { Store } from './store';

/** One file added to the run directory ("Run directory": `-conffile`, BLAST database files, `-ids`, …). */
export interface ToolExtraFile {
  readonly key: string;
  readonly name: string;
  readonly size: number;
  readonly file: File;
  /** Option flag that picked this file via `Choose file…` (`-ids`, `-exclude_ids`, `-mask_data`), if any. */
  readonly flag?: string;
  /** Where the file came from, for Run details (`made by run N`); a file that the user opened when absent. */
  readonly origin?: ToolFileOrigin;
}

/** The state that every tool draft holds. */
export interface ToolFormState {
  readonly runFiles: readonly ToolExtraFile[];
  readonly runFilesError?: string;
  readonly description?: ToolDescription;
  readonly descriptionError?: string;
  readonly descriptionUnavailable?: true;
  readonly validation: Validation;
  readonly submitting: boolean;
  readonly message?: { readonly kind: 'info' | 'error'; readonly text: string };
}

/** What Availability and Attention read of each tool tab's form. */
export interface ToolForm {
  readonly state: { subscribe(listener: () => void): () => void };
  /** Asks the engine again for what failed with an UnavailableError (`Availability.recheck`). */
  askAgain(): Promise<void>;
  /** The message of the first part of LOSAT that the form found unavailable, if any. */
  unavailable(): string | undefined;
  /** Whether the form holds inputs or files that a reload would lose. */
  holdsInputs(): boolean;
}

export interface ToolDraftCoreDeps {
  readonly engine: Pick<ToolEngine, 'describeTool' | 'validateTool'>;
  readonly enqueueToolAll: (requests: readonly ToolRequest[]) => Promise<EnqueueAllResult>;
  /** The holders of inflated session files, for a tool whose inputs are inflated into OPFS. */
  readonly sessionFiles?: SessionFiles;
  /** True while the Engine worker is running another job (so `validateTool` would wait behind it). */
  readonly isBusy?: () => boolean;
  readonly debounceMs?: number;
}

/** A change of a draft's state: a key set to `undefined` is removed. */
export type StatePatch<S> = { readonly [K in keyof S]?: S[K] | undefined };

/** The argv that a check sends, and the run directory that it sees. */
interface Check {
  readonly argv: readonly string[];
  readonly directory: readonly DirectoryFile[];
}

export abstract class ToolDraftCore<S extends ToolFormState> implements ToolForm {
  readonly state: Store<S>;
  protected nextKey = 1;
  private descriptionRequest: Promise<ToolDescription> | undefined;
  private validationTimer: ReturnType<typeof setTimeout> | undefined;
  private validationGeneration = 0;
  private validationTurn: Promise<void> = Promise.resolve();
  private validatedKey: string | undefined;
  private readonly asking = new Set<Promise<unknown>>();

  protected constructor(
    protected readonly tool: ToolId,
    private readonly core: ToolDraftCoreDeps,
    initial: S,
  ) {
    this.state = new Store<S>(initial);
  }

  // --- what a tool supplies -----------------------------------------------------------------

  /** The words of Additional options; throws the message of words that cannot be read. */
  protected abstract parseExtra(state: S): readonly string[];

  /**
   * The argv that the check sends (with placeholder names for inputs not chosen yet), or undefined
   * when the form has nothing to check yet (the check is then idle). Throws for words it cannot build.
   */
  protected abstract validationArgv(state: S, extra: readonly string[]): readonly string[] | undefined;

  /** The run's inputs that the run directory holds besides the files of the run; at the check, `file` may be a stand-in. */
  protected validationInputs?(state: S): readonly DirectoryFile[];

  /** The runs to queue for `state`; throws the message to show when the inputs are not ready. */
  protected abstract buildRequests(state: S, extra: readonly string[]): Promise<readonly ToolRequest[]>;

  /** The argv shown beside Run LOSAT. */
  abstract previewArgv(state?: S): readonly string[];

  /** Whether the form holds inputs or files that a reload would lose. */
  abstract holdsInputs(state?: S): boolean;

  /** The inflated session files that the draft's inputs hold (held while a submit builds its runs). */
  protected heldTokens?(state: S): readonly string[];

  /** The message of an input that could not be read because a part of LOSAT is unavailable. */
  protected inputsUnavailable?(state: S): string | undefined;

  /**
   * What else changes when files are added to the run directory (WindowMasker's `-in` name, makeblastdb's database
   * name); its `runFilesError` says why the files could not name what they hold, unless a name was refused already.
   */
  protected runFilesAdded?(state: S, runFiles: readonly ToolExtraFile[]): StatePatch<S>;

  /** The value of an option of the form (of the stage shown, for WindowMasker). */
  protected abstract optionValue(state: S, flag: string): FieldValue | undefined;

  /** Sets an option of the form; an empty value or `false` leaves it out. */
  abstract setField(flag: string, value: FieldValue): void;

  // --- the tool's description ---------------------------------------------------------------

  /** Loads the tool's `-xmlhelp` once, when its tab is first opened. */
  ensureLoaded(): void {
    if (this.state.get().description === undefined && this.descriptionRequest === undefined) {
      void this.loadDescription();
    }
  }

  /** Asks for `-xmlhelp` again and checks the argv again when LOSAT became unavailable (Availability.recheck). */
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

  unavailable(state: S = this.state.get()): string | undefined {
    if (state.descriptionUnavailable === true && state.descriptionError !== undefined) return state.descriptionError;
    if (state.validation.state === 'unavailable') return state.validation.message;
    return this.inputsUnavailable?.(state);
  }

  // --- the files of the run -------------------------------------------------------------------

  /** Adds files (or a folder picked with `webkitdirectory`) to the run directory; refuses a name already there. */
  addRunFiles(files: readonly File[], relativeNames?: readonly string[]): void {
    if (files.length === 0) return;
    const state = this.state.get();
    const builder = new RunDirectoryBuilder();
    for (const entry of state.runFiles) builder.add(entry.name, entry.file);
    const added: ToolExtraFile[] = [];
    let duplicateError: string | undefined;
    for (let i = 0; i < files.length; i++) {
      const file = files[i]!;
      const name = pickedFilePath(file, relativeNames?.[i]);
      if (name === undefined) continue;
      if (builder.has(name)) {
        duplicateError ??= duplicateRunFileMessage(name);
        continue;
      }
      try {
        builder.add(name, file);
      } catch (error) {
        duplicateError ??= errorMessage(error);
        continue;
      }
      added.push({ key: `rf-${this.nextKey++}`, name, size: file.size, file });
    }
    const runFiles = Object.freeze([...state.runFiles, ...added]);
    const found = this.runFilesAdded?.(state, runFiles);
    this.set({ ...found, ...common<S>({ runFiles, runFilesError: duplicateError ?? found?.runFilesError, message: undefined }) });
    this.scheduleValidation();
  }

  /** Removes a file of the run; a file chosen for an option (`-ids`) also leaves the option (its name leaves a list's). */
  removeRunFile(key: string): void {
    const state = this.state.get();
    const target = state.runFiles.find((entry) => entry.key === key);
    if (target === undefined) return;
    this.patch({
      runFiles: Object.freeze(state.runFiles.filter((entry) => entry.key !== key)),
      runFilesError: undefined,
      message: undefined,
    });
    this.scheduleValidation();
    if (target.flag === undefined) return;
    const value = this.optionValue(state, target.flag);
    if (TOOL_FILE_LIST_FLAGS.has(target.flag)) {
      const names = fileListNames(value);
      if (names.includes(target.name)) this.setField(target.flag, names.filter((name) => name !== target.name).join(','));
    } else if (value === target.name) {
      this.setField(target.flag, '');
    }
  }

  /**
   * Attaches a file chosen for a file-typed option (`-ids`, `-taxid_map`) and puts it in Files of the run under its
   * name (or `options.name`), which the option then names; undefined lets the option's files go. A file chosen again
   * replaces the option's file, except for an option that lists files (`-mask_data`, `TOOL_FILE_LIST_FLAGS`), whose list
   * it joins (the same file twice is one).
   */
  setOptionFile(flag: string, file: File | undefined, options: { readonly name?: string; readonly origin?: ToolFileOrigin } = {}): void {
    const state = this.state.get();
    if (file === undefined) {
      this.patch({ runFiles: Object.freeze(state.runFiles.filter((entry) => entry.flag !== flag)), runFilesError: undefined });
      this.setField(flag, '');
      return;
    }
    const given = options.name ?? file.name;
    const name = normalizeRunPath(given);
    if (name === undefined) {
      this.patch({ runFilesError: `Invalid file name: "${given}".` });
      return;
    }
    const refusal = listNameRefusal({ command: this.tool, flag }, name);
    if (refusal !== undefined) {
      this.patch({ runFilesError: refusal });
      return;
    }
    const list = TOOL_FILE_LIST_FLAGS.has(flag);
    const kept = list ? state.runFiles : state.runFiles.filter((entry) => entry.flag !== flag);
    const names = list ? fileListNames(this.optionValue(state, flag)) : [];
    const value = list ? [...names.filter((entry) => entry !== name), name].join(',') : name;
    const existing = kept.find((entry) => entry.name === name);
    if (existing !== undefined) {
      if (list && existing.flag === flag && sameRunFile(name, existing.file, file)) {
        this.patch({ runFilesError: undefined });
        this.setField(flag, value);
        return;
      }
      this.patch({ runFilesError: duplicateRunFileMessage(name) });
      return;
    }
    const entry: ToolExtraFile = {
      key: `rf-${this.nextKey++}`,
      name,
      size: file.size,
      file,
      flag,
      ...(options.origin === undefined ? {} : { origin: options.origin }),
    };
    this.patch({ runFiles: Object.freeze([...kept, entry]), runFilesError: undefined });
    this.setField(flag, value);
  }

  // --- the command line, the check and the submit --------------------------------------------

  previewCommandWords(state: S = this.state.get()): readonly string[] {
    return toolCommandWords(this.previewArgv(state));
  }

  /** Runs `validateTool` now (one at a time, so an older answer is dropped). */
  runValidationNow(): Promise<void> {
    const generation = ++this.validationGeneration;
    const stale = () => generation !== this.validationGeneration;
    const run = async () => {
      const state = this.state.get();
      if (stale()) return;
      let check: Check | undefined;
      try {
        check = this.check(state);
      } catch (error) {
        this.patch({ validation: { state: 'invalid', message: errorMessage(error) } });
        return;
      }
      if (check === undefined) {
        this.validatedKey = undefined;
        this.patch({ validation: { state: 'idle' } });
        return;
      }
      this.patch({ validation: { state: 'checking' } });
      let verdict: ToolValidation;
      try {
        verdict = await this.core.engine.validateTool(check.argv, { directory: check.directory });
      } catch (error) {
        if (stale()) return;
        const message = errorMessage(error);
        this.patch({
          validation: isUnavailable(error) ? { state: 'unavailable', message } : { state: 'invalid', message },
        });
        return;
      }
      if (stale()) return;
      this.validatedKey = checkKey(check);
      this.patch({
        validation: verdict.ok ? { state: 'ok' } : { state: 'invalid', message: verdict.message },
      });
    };
    const next = this.validationTurn.then(run, run);
    this.validationTurn = next.then(
      () => undefined,
      () => undefined,
    );
    return next;
  }

  /** Queues the form's runs on the shared Coordinator queue, after the check of the current form. */
  async submit(): Promise<SubmitResult> {
    if (this.state.get().submitting) return { ok: false, message: 'A submission is already in progress.' };
    this.patch({ submitting: true, message: undefined });
    try {
      const busy = this.core.isBusy?.() === true;
      if (!busy) await this.settleValidation();
      // Re-read the draft state after awaiting the check (same order as SearchDraft.prepareAndEnqueue),
      // so changes that arrived during the await (an inflate finishing, an option changing) are used.
      const state = this.state.get();
      // Hold every inflated session file before any await of buildRequests, so that Remove during the
      // submit cannot drop an OPFS session file out from under it.
      const holder = `submit:${this.nextKey++}`;
      const stateTokens = this.heldTokens?.(state) ?? [];
      if (stateTokens.length > 0) this.core.sessionFiles?.set(holder, stateTokens);
      try {
        let extra: readonly string[];
        try {
          extra = this.parseExtra(state);
        } catch (error) {
          const message = errorMessage(error);
          this.reportError(message);
          return { ok: false, message };
        }
        const validated = state.validation;
        // While the engine is busy the check may not have answered for this form: its refusal counts only when it was for it.
        if (validated.state === 'invalid' && (!busy || this.validatedKey === this.currentKey(state))) {
          const text = toolValidationSummary(validated.message);
          this.reportError(text);
          return { ok: false, message: text };
        }
        const requests = await this.buildRequests(state, extra);
        const requestTokens = requests.flatMap((request) => request.sessionTokens ?? []);
        if (requestTokens.length > 0) this.core.sessionFiles?.set(holder, new Set([...stateTokens, ...requestTokens]));
        const result = await this.core.enqueueToolAll(requests);
        if (!result.ok) {
          if (result.reason === 'counts') {
            this.reportError(result.message);
            return { ok: false, message: result.message };
          }
          const text = toolValidationSummary(result.message);
          this.patch({
            validation: { state: 'invalid', message: result.message },
            message: { kind: 'error', text },
          });
          return { ok: false, message: text };
        }
        return result;
      } finally {
        this.core.sessionFiles?.release(holder);
      }
    } catch (error) {
      const message = failureText(error);
      this.reportError(message);
      return { ok: false, message };
    } finally {
      this.patch({ submitting: false });
    }
  }

  // --- for the tools ------------------------------------------------------------------------

  /** The run directory of `inputs` and the files of the run; throws `duplicateRunFileMessage` when two files claim one name. */
  protected runDirectory<E extends DirectoryFile>(state: S, inputs: readonly E[], origin: (entry: ToolExtraFile) => E): readonly E[] {
    const builder = new RunDirectoryBuilder();
    const kept: E[] = [];
    for (const entry of [...inputs, ...state.runFiles.map(origin)]) {
      const existed = builder.has(entry.name);
      builder.add(entry.name, entry.file);
      if (!existed) kept.push(entry);
    }
    return Object.freeze(kept);
  }

  /** The run directory of a run: `inputs` and the files of the run, as picked files. */
  protected runDirectoryOf(state: S, inputs: readonly ToolDirectoryFile[]): readonly ToolDirectoryFile[] {
    return this.runDirectory(state, inputs, (entry) => ({ name: entry.name, file: entry.file, origin: entry.origin ?? { kind: 'picked' } }));
  }

  protected scheduleValidation(): void {
    if (this.validationTimer !== undefined) clearTimeout(this.validationTimer);
    const delay = this.core.debounceMs ?? 150;
    this.validationTimer = setTimeout(() => {
      this.validationTimer = undefined;
      void this.runValidationNow();
    }, delay);
  }

  protected reportError(text: string): void {
    this.patch({ message: { kind: 'error', text } });
  }

  /** Changes the state that every tool draft holds. */
  protected patch(change: StatePatch<ToolFormState>): void {
    this.set(common(change));
  }

  protected set(change: StatePatch<S>): void {
    const next = { ...this.state.get(), ...change } as Record<string, unknown>;
    for (const key of Object.keys(change)) {
      if ((change as Record<string, unknown>)[key] === undefined) delete next[key];
    }
    this.state.set(next as S);
  }

  // --- internals ----------------------------------------------------------------------------

  /** The argv and run directory that the check of `state` sends; undefined when there is nothing to check. */
  private check(state: S): Check | undefined {
    const argv = this.validationArgv(state, this.parseExtra(state));
    if (argv === undefined) return undefined;
    const inputs = this.validationInputs?.(state) ?? [];
    const directory = this.runDirectory(state, inputs, (entry): DirectoryFile => ({ name: entry.name, file: entry.file }));
    return { argv, directory: directory.map((entry) => ({ name: entry.name, file: entry.file })) };
  }

  private currentKey(state: S): string | undefined {
    try {
      const check = this.check(state);
      return check === undefined ? undefined : checkKey(check);
    } catch {
      return undefined;
    }
  }

  /** Waits until no check is scheduled or running (a check scheduled meanwhile runs now). */
  private async settleValidation(): Promise<void> {
    for (;;) {
      if (this.validationTimer !== undefined) {
        clearTimeout(this.validationTimer);
        this.validationTimer = undefined;
        await this.runValidationNow();
        continue;
      }
      const turn = this.validationTurn;
      await turn;
      if (this.validationTimer === undefined && this.validationTurn === turn) return;
    }
  }

  private loadDescription(force = false): Promise<ToolDescription | undefined> {
    if (!force && this.descriptionRequest !== undefined) return this.descriptionRequest;
    const request = this.core.engine.describeTool(this.tool);
    this.descriptionRequest = request;
    return request.then(
      (description) => {
        this.patch({ description, descriptionError: undefined, descriptionUnavailable: undefined });
        this.scheduleValidation();
        return description;
      },
      (error: unknown) => {
        if (this.descriptionRequest === request) this.descriptionRequest = undefined;
        this.patch({
          descriptionError: errorMessage(error),
          ...(isUnavailable(error) ? { descriptionUnavailable: true as const } : {}),
        });
        return undefined;
      },
    );
  }

  private trackAsking(promise: Promise<unknown>): void {
    this.asking.add(promise);
    void promise.finally(() => this.asking.delete(promise));
  }
}

/** A change of the common state as a change of a tool draft's state (whose type extends it). */
function common<S extends ToolFormState>(change: StatePatch<ToolFormState>): StatePatch<S> {
  return change as StatePatch<S>;
}

/** Whether two files under one name are one file (the same File, or a second pick of one file), as the run directory takes them. */
function sameRunFile(name: string, a: File, b: File): boolean {
  const builder = new RunDirectoryBuilder();
  builder.add(name, a);
  try {
    builder.add(name, b);
    return true;
  } catch {
    return false;
  }
}

/** The key of a check: its argv and the names of its run directory. */
function checkKey(check: Check): string {
  return JSON.stringify([check.argv, check.directory.map((entry) => entry.name)]);
}

export function errorMessage(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
