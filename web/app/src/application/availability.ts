// Whether LOSAT itself works (line WQA's sweep, D1 and D2, 2026-11-19). A part of LOSAT that could not be loaded (an
// engine file that could not be downloaded, a worker whose script could not start) or that ended for good (the Data
// worker) is no fault of the inputs or the options: the page says it once, with the remedy, which is to reload the
// page (UnavailableError). Before, each input card said "This input cannot be read: Failed to fetch" in the browser's
// words, Run LOSAT was refused for the pasted query, and the only sentence about the engine was in the closed
// Algorithm parameters. What the page says comes from what the search form found (`unavailableOf`), so it goes
// once the form reads its inputs and options again, and from the Data worker, which comes back only with a reload.
import { isUnavailable } from 'losat/ports';
import type { DraftState, SearchDraft } from './draft';
import { Store } from './store';
import type { BlastFormatterDraft, BlastFormatterDraftState } from './blast-formatter-draft';
import type { ToolDraft, ToolDraftState } from './tool-draft';

export interface AvailabilityState {
  /** What of LOSAT could not be loaded or stopped (an UnavailableError's message); undefined while LOSAT works. */
  readonly problem?: string;
}

export interface AvailabilityDeps {
  /**
   * The search form, whose describe, inputs, checks and validation say when a part of LOSAT could not be loaded, and
   * which asks the engine again (recheck).
   */
  readonly draft: { readonly state: Store<DraftState> } & Pick<SearchDraft, 'askAgain'>;
  /** The WindowMasker form, whose describeTool and validateTool can also observe an unavailable engine. */
  readonly toolDraft?: { readonly state: Store<ToolDraftState> } & Pick<ToolDraft, 'askAgain'>;
  /** The Format Archive form, whose describeTool and validateTool can also observe an unavailable engine. */
  readonly blastFormatterDraft?: { readonly state: Store<BlastFormatterDraftState> } & Pick<BlastFormatterDraft, 'askAgain'>;
  /** Settles when the Data worker could not start or ended: for good, as no call of it can succeed again. */
  readonly dataWorkerEnded: Promise<Error>;
}

export class Availability {
  readonly state = new Store<AvailabilityState>({});
  private ended: string | undefined;
  private draftWasUnavailable = false;

  constructor(private readonly deps: AvailabilityDeps) {
    deps.draft.state.subscribe(() => this.update());
    deps.toolDraft?.state.subscribe(() => this.update());
    deps.blastFormatterDraft?.state.subscribe(() => this.update());
    void deps.dataWorkerEnded.then((error) => {
      this.ended = error.message;
      this.update();
    });
    this.update();
  }

  /**
   * What of LOSAT is unavailable once the search form asked the engine again what it could not ask it: an engine file
   * that could not be downloaded may download now. The Data worker does not come back. An Example decides on it; it
   * was refused with "Reload the page", and asked nothing, also once the engine file downloaded again (the third
   * review of line WQA's D1 and D2, M1, 2026-11-19). It waits for the answers to what the form asked again, not for
   * the rest of the form's work (its fourth review, M1), also for those still on their way when nothing reads as
   * unavailable meanwhile: a second press was refused at once, or downloaded at once (its fifth review, M1).
   */
  async recheck(): Promise<string | undefined> {
    if (this.ended === undefined) {
      await Promise.all([this.deps.draft.askAgain(), this.deps.toolDraft?.askAgain(), this.deps.blastFormatterDraft?.askAgain()]);
    }
    return this.state.get().problem;
  }

  private update(): void {
    const draftProblem = unavailableOf(this.deps.draft.state.get());
    const toolState = this.deps.toolDraft?.state.get();
    const toolProblem = toolState === undefined ? undefined : unavailableOfTool(toolState);
    const bfState = this.deps.blastFormatterDraft?.state.get();
    const bfProblem = bfState === undefined ? undefined : unavailableOfBlastFormatter(bfState);
    if (this.draftWasUnavailable && draftProblem === undefined && this.ended === undefined) {
      if (toolProblem !== undefined) void this.deps.toolDraft?.askAgain();
      if (bfProblem !== undefined) void this.deps.blastFormatterDraft?.askAgain();
    }
    this.draftWasUnavailable = draftProblem !== undefined;
    const problem = this.ended ?? draftProblem ?? toolProblem ?? bfProblem;
    if (problem === this.state.get().problem) return;
    this.state.set(problem === undefined ? {} : { problem });
  }
}

export function unavailableOfBlastFormatter(state: BlastFormatterDraftState): string | undefined {
  if (state.descriptionUnavailable === true && state.descriptionError !== undefined) return state.descriptionError;
  if (state.validation.state === 'unavailable') return state.validation.message;
  return undefined;
}

export function unavailableOfTool(state: ToolDraftState): string | undefined {
  if (state.descriptionUnavailable === true && state.descriptionError !== undefined) return state.descriptionError;
  if (state.validation.state === 'unavailable') return state.validation.message;
  for (const file of state.fastaFiles) {
    if (file.status === 'failed' && file.unavailable === true && file.error !== undefined) return file.error;
  }
  return undefined;
}

/**
 * The message of the first part of LOSAT that the search form found unavailable: its describe, an input that could not
 * be read or checked, or its validation, failed with an UnavailableError; undefined when none did.
 */
export function unavailableOf(state: DraftState): string | undefined {
  if (state.descriptionUnavailable === true && state.descriptionError !== undefined) return state.descriptionError;
  if (state.validation.state === 'unavailable') return state.validation.message;
  for (const role of ['query', 'subject'] as const) {
    for (const source of state[role].sources) {
      if (source.status === 'failed' && source.unavailable === true && source.error !== undefined) return source.error;
      if (source.check?.state === 'error' && source.check.unavailable === true) return source.check.message;
    }
  }
  return undefined;
}

/** An error's message as the sentence that a refusal or a failed run shows: with the remedy when a part of LOSAT is unavailable. */
export function failureText(error: unknown): string {
  const message = error instanceof Error ? error.message : String(error);
  return isUnavailable(error) ? withRemedy(message) : message;
}

/**
 * The message of an UnavailableError and its remedy, as sentences; with `resultsLost`, the consequence that a reload
 * has when runs have results (ResumeNotice's words).
 */
export function withRemedy(message: string, resultsLost = false): string {
  const text = message.trimEnd();
  const sentence = /[.!?]["')\]]*$/.test(text) ? text : `${text}.`;
  return `${sentence} Reload the page${resultsLost ? '; the results of this session are then lost' : ''}.`;
}
