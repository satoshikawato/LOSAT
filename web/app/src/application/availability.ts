// Whether LOSAT itself works (line WQA's sweep, D1 and D2, 2026-11-19). A part of LOSAT that could not be loaded (an
// engine file that could not be downloaded, a worker whose script could not start) or that ended for good (the Data
// worker) is no fault of the inputs or the options: the page says it once, with the remedy, which is to reload the
// page (UnavailableError). Before, each input card said "This input cannot be read: Failed to fetch" in the browser's
// words, Run LOSAT was refused for the pasted query, and the only sentence about the engine was in the closed
// Algorithm parameters. What the page says comes from what the search form found (`unavailableOf`), so it goes
// once the form reads its inputs and options again, and from the Data worker, which comes back only with a reload.
import { isUnavailable } from 'losat/ports';
import type { DraftState } from './draft';
import { Store } from './store';

export interface AvailabilityState {
  /** What of LOSAT could not be loaded or stopped (an UnavailableError's message); undefined while LOSAT works. */
  readonly problem?: string;
}

export interface AvailabilityDeps {
  /** The search form, whose describe, inputs, checks and validation say when a part of LOSAT could not be loaded. */
  readonly draft: Store<DraftState>;
  /** Settles when the Data worker could not start or ended: for good, as no call of it can succeed again. */
  readonly dataWorkerEnded: Promise<Error>;
}

export class Availability {
  readonly state = new Store<AvailabilityState>({});
  private ended: string | undefined;

  constructor(private readonly deps: AvailabilityDeps) {
    deps.draft.subscribe(() => this.update());
    void deps.dataWorkerEnded.then((error) => {
      this.ended = error.message;
      this.update();
    });
    this.update();
  }

  private update(): void {
    const problem = this.ended ?? unavailableOf(this.deps.draft.get());
    if (problem === this.state.get().problem) return;
    this.state.set(problem === undefined ? {} : { problem });
  }
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
