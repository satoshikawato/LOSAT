<script setup lang="ts">
// The status of a run that has not finished, on its own results page (NCBI's "Status:
// Searching" page while a job runs): its phase and elapsed time (not progress by query, DW-5),
// what it waits for, and its Cancel. The run's record (engine path, threads, phases) is in Run
// details once it completes.
import { computed, onUnmounted, ref, watchEffect } from 'vue';
import type { Coordinator, RunView } from '../application/coordinator';
import { isTerminal, type RunStatus } from '../domain/run';
import { formatDuration } from './format';

const props = defineProps<{ coordinator: Coordinator; run: RunView; runs: readonly RunView[] }>();

const PHASE_LABELS: Readonly<Record<RunStatus, string>> = {
  queued: 'Waiting',
  preparing: 'Preparing',
  running: 'Searching',
  finalizing: 'Organizing results',
  completed: 'Completed',
  cancelled: 'Cancelled',
  failed: 'Failed',
};

const now = ref(Date.now());
let timer: ReturnType<typeof setInterval> | undefined;
const ticking = computed(() => props.run.record.startedAt !== undefined && !isTerminal(props.run.status));
watchEffect(() => {
  if (ticking.value && timer === undefined) {
    timer = setInterval(() => (now.value = Date.now()), 1000);
  } else if (!ticking.value && timer !== undefined) {
    clearInterval(timer);
    timer = undefined;
  }
});
onUnmounted(() => clearInterval(timer));

const elapsed = computed(() => {
  const { startedAt, endedAt } = props.run.record;
  return startedAt === undefined ? undefined : formatDuration((endedAt ?? now.value) - startedAt);
});

const cancellable = (run: RunView) => !isTerminal(run.status) && run.status !== 'finalizing';

/** The run that searches while this one waits (runs search one at a time). */
const ahead = computed(() => {
  if (props.run.status !== 'queued') return undefined;
  return props.runs.find((other) => !isTerminal(other.status) && other.status !== 'queued');
});
/** A run of a group shows the group's cancel while another of its runs can still be cancelled. */
const groupCancellable = computed(() => {
  const group = props.run.snapshot.group;
  if (group === undefined) return false;
  return props.runs.some((other) => other.snapshot.group?.groupId === group.groupId && other !== props.run && cancellable(other));
});

const phaseTimes = computed(() => {
  const { startedAt, phaseTimes: times } = props.run.record;
  if (startedAt === undefined || times === undefined) return '';
  return (['preparing', 'running', 'finalizing'] as const)
    .filter((phase) => times[phase] !== undefined)
    .map((phase) => `${PHASE_LABELS[phase]} at ${((times[phase]! - startedAt) / 1000).toFixed(1)} s`)
    .join(', ');
});
</script>

<template>
  <div class="run-status" :data-status="run.status" data-testid="run-status" :data-run="run.snapshot.number">
    <p class="run-status-line">
      <span class="muted">Status:</span>
      <span class="phase" data-testid="run-status-phase">{{ ahead ? `Waiting for run ${ahead.snapshot.number}` : PHASE_LABELS[run.status] }}</span>
      <span v-if="elapsed" class="elapsed muted" data-testid="run-status-elapsed">{{ elapsed }}</span>
      <span class="run-status-actions">
        <button v-if="groupCancellable" type="button" class="link" data-testid="run-status-cancel-group" @click="coordinator.cancelGroup(run.snapshot.group!.groupId)">
          Cancel the group
        </button>
        <button v-if="cancellable(run)" type="button" data-testid="run-status-cancel" @click="coordinator.cancel(run.snapshot.runId)">Cancel</button>
      </span>
    </p>
    <p v-if="phaseTimes" class="muted small" data-testid="run-status-phases">{{ phaseTimes }}</p>
    <!-- The one consequence that matters while a search runs (stop-describing-ui rule 3). -->
    <p class="muted small run-status-warning" data-testid="run-status-tab-warning">
      <svg aria-hidden="true" width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M12 9v4m0 4h.01M10.3 3.9 1.8 18a2 2 0 0 0 1.7 3h17a2 2 0 0 0 1.7-3L13.7 3.9a2 2 0 0 0-3.4 0z" /></svg>
      Closing or reloading this tab ends the search
    </p>
  </div>
</template>
