// The words of a run that has not finished (RunStatus.vue). A run that prints the command line's help or version in place of a
// search (`RunSnapshot.prints`, line WQA's D14) runs no search, so its status says "Printing" and its warning "the run" (line
// WQA's D20 round 3, 2026-11-23).
import type { RunStatus } from '../domain/run';

const LABELS: Readonly<Record<RunStatus, string>> = {
  queued: 'Waiting',
  preparing: 'Preparing',
  running: 'Searching',
  finalizing: 'Organizing results',
  completed: 'Completed',
  cancelled: 'Cancelled',
  failed: 'Failed',
};

/** The phase of a run in the status block and for the screen reader; `printed`: the run prints help or version. */
export function phaseLabel(status: RunStatus, printed: boolean): string {
  return printed && status === 'running' ? 'Printing' : LABELS[status];
}

/** What closing or reloading the tab ends while the run has not finished. */
export function tabWarning(printed: boolean): string {
  return `Closing or reloading this tab ends the ${printed ? 'run' : 'search'}`;
}
