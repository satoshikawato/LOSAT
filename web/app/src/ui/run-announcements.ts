// How each run of the working session ends, told to a screen reader (WA2-AX-03; announcer.ts):
// a failed run with its reason as an alert, a cancelled run, and a completed run with its HSPs
// and subjects. The words are those of the results page ("Run 2 failed: … A failed run keeps no
// results."). A run whose results the page reads as it completes is told once they are read; the
// status of a run in progress is told by RunStatus, which shows it.
import { watch, type Ref } from 'vue';
import type { AppState } from '../application/coordinator';
import { noHspReason, unavailableMessage, type ResultsState } from '../application/results';
import { isTerminal, type RunStatus } from '../domain/run';
import { announce } from './announcer';
import { formatCounted } from './format';
import { NO_HSP_TEXT } from './no-hsp-text';

/** "Run 1 completed: 11 HSPs on 8 subjects." from the results read; the page's sentence when there is no HSP. */
function completedText(number: number, results: ResultsState): string {
  const table = results.loaded?.index.table;
  if (table === undefined) return `Run ${number} completed.`;
  if (table.count === 0) {
    const reason = noHspReason(results);
    return `Run ${number} completed: ${reason === undefined ? NO_HSP_TEXT['no-hits'] : NO_HSP_TEXT[reason]}`;
  }
  const subjects = new Set(table.sIdx).size;
  return `Run ${number} completed: ${formatCounted(table.count, 'HSP')} on ${formatCounted(subjects, 'subject')}.`;
}

export function useRunAnnouncements(app: Readonly<Ref<AppState>>, results: Readonly<Ref<ResultsState>>): void {
  const seen = new Map<string, RunStatus>();
  /** The runs of the page's start are not new: a run first seen later has just been added. */
  let started = false;
  /** The run that completed while the page shows it: told when its results are read. */
  let reading: { readonly runId: string; readonly number: number } | undefined;

  function settle(): void {
    if (reading === undefined) return;
    const state = results.value;
    if (state.runId !== reading.runId) {
      announce(`Run ${reading.number} completed.`);
    } else if (state.phase === 'ready') {
      announce(completedText(reading.number, state));
    } else if (state.phase === 'failed') {
      announce(`Run ${reading.number} completed.`);
      announce(state.message ?? '', 'alert');
    } else {
      return;
    }
    reading = undefined;
  }

  watch(
    () => app.value.runs,
    (runs) => {
      for (const run of runs) {
        const { runId, number } = run.snapshot;
        const before = seen.get(runId) ?? (started ? 'queued' : run.status);
        seen.set(runId, run.status);
        if (before === run.status || isTerminal(before) || !isTerminal(run.status)) continue;
        if (run.status === 'failed') announce(unavailableMessage(run), 'alert');
        else if (run.status === 'cancelled') announce(unavailableMessage(run));
        else if (results.value.runId === runId) {
          settle();
          reading = { runId, number };
          settle();
        } else announce(`Run ${number} completed.`);
      }
      started = true;
    },
    { immediate: true },
  );
  watch(() => [results.value.runId, results.value.phase], settle);
}
