<script setup lang="ts">
// The runs of this working session: their status, the phase and elapsed time of the
// active one (not progress by query, DW-5), and what the engine runtime did (RunRecord).
import { computed, onUnmounted, ref, watchEffect } from 'vue';
import type { Coordinator, RunView } from '../application/coordinator';
import { isTerminal, type RunStatus } from '../domain/run';
import { programById } from '../domain/programs';
import { formatBytes, formatDuration } from './format';

const props = defineProps<{ coordinator: Coordinator; runs: readonly RunView[] }>();
defineEmits<{ 'open-results': [runId: string] }>();

/**
 * The runs in the order that matters while working (S12's screen review L-e): the running
 * one first, then those waiting in queue order, then the finished ones, newest first.
 */
const ordered = computed(() => {
  const active = props.runs.filter((run) => !isTerminal(run.status) && run.status !== 'queued');
  const waiting = props.runs.filter((run) => run.status === 'queued');
  const finished = props.runs.filter((run) => isTerminal(run.status)).reverse();
  return [...active, ...waiting, ...finished];
});

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
const anyActive = computed(() => props.runs.some((run) => run.record.startedAt !== undefined && !isTerminal(run.status)));
watchEffect(() => {
  if (anyActive.value && timer === undefined) {
    timer = setInterval(() => (now.value = Date.now()), 1000);
  } else if (!anyActive.value && timer !== undefined) {
    clearInterval(timer);
    timer = undefined;
  }
});
onUnmounted(() => clearInterval(timer));

/** The options of the run's argv: the words after the program and the two inputs (plan §5.3). */
function options(run: RunView): readonly string[] {
  const words = run.snapshot.argv.slice(5);
  return words.length === 0 ? ['defaults'] : words;
}

const cancellable = (run: RunView) => !isTerminal(run.status) && run.status !== 'finalizing';

function elapsed(run: RunView): string | undefined {
  const { startedAt, endedAt } = run.record;
  if (startedAt === undefined) return undefined;
  return formatDuration((endedAt ?? now.value) - startedAt);
}

/** The first listed run of a group that can be cancelled shows the group's cancel button. */
function groupCancellable(run: RunView): boolean {
  const group = run.snapshot.group;
  if (group === undefined) return false;
  return ordered.value.find((other) => other.snapshot.group?.groupId === group.groupId && cancellable(other)) === run;
}

function phaseTimes(run: RunView): string {
  const { startedAt, phaseTimes: times } = run.record;
  if (startedAt === undefined || times === undefined) return '';
  return (['preparing', 'running', 'finalizing'] as const)
    .filter((phase) => times[phase] !== undefined)
    .map((phase) => `${PHASE_LABELS[phase]} at ${((times[phase]! - startedAt) / 1000).toFixed(1)} s`)
    .join(', ');
}
</script>

<template>
  <section id="queue" class="queue-panel">
    <h2>Queue</h2>
    <p v-if="runs.length === 0" class="muted">No runs yet.</p>
    <ol class="queue" data-testid="queue">
      <li
        v-for="run in ordered"
        :key="run.snapshot.runId"
        class="run"
        :data-status="run.status"
        :data-testid="`run-${run.snapshot.number}`"
      >
        <!-- The run's state wraps beside its button, so that the button shares a line with it
             (S13 screen review L3: "Open results" took a line of its own). -->
        <div class="run-line">
          <div class="run-state">
            <strong>Run {{ run.snapshot.number }}</strong>
            <span>{{ programById(run.snapshot.program).label }}</span>
            <span class="status" :data-status="run.status" :data-testid="`run-${run.snapshot.number}-status`">{{
              run.status
            }}</span>
            <span v-if="!isTerminal(run.status) && run.status !== 'queued'" class="phase" :data-testid="`run-${run.snapshot.number}-phase`">
              {{ PHASE_LABELS[run.status] }}
            </span>
            <span v-if="elapsed(run)" class="elapsed muted" :data-testid="`run-${run.snapshot.number}-elapsed`">
              {{ elapsed(run) }}
            </span>
          </div>
          <button
            v-if="cancellable(run)"
            type="button"
            :data-testid="`run-${run.snapshot.number}-cancel`"
            @click="coordinator.cancel(run.snapshot.runId)"
          >
            Cancel
          </button>
          <button
            v-if="run.status === 'completed'"
            type="button"
            class="run-open"
            :data-testid="`run-${run.snapshot.number}-open`"
            @click="$emit('open-results', run.snapshot.runId)"
          >
            Open results
          </button>
        </div>
        <div v-if="run.snapshot.title" class="run-title" :data-testid="`run-${run.snapshot.number}-title`">
          {{ run.snapshot.title }}
        </div>
        <div class="run-inputs muted">
          {{ run.snapshot.query.name }} vs {{ run.snapshot.subject.name }}
          <template v-if="run.snapshot.group">
            · group run {{ run.snapshot.group.position }} of {{ run.snapshot.group.size }}
          </template>
        </div>
        <div class="run-inputs muted" :data-testid="`run-${run.snapshot.number}-options`">
          Options:
          <template v-for="(word, i) in options(run)" :key="i"><span class="word">{{ word }}</span>{{ ' ' }}</template>
        </div>
        <button
          v-if="groupCancellable(run)"
          type="button"
          class="link"
          :data-testid="`run-${run.snapshot.number}-cancel-group`"
          @click="coordinator.cancelGroup(run.snapshot.group!.groupId)"
        >
          Cancel the group
        </button>
        <p v-if="run.record.error" class="error" :data-testid="`run-${run.snapshot.number}-error`">{{ run.record.error }}</p>
        <details v-if="run.record.startedAt !== undefined" class="diagnostics">
          <summary>Details</summary>
          <dl :data-testid="`run-${run.snapshot.number}-details`">
            <template v-if="run.record.runtimePath">
              <dt>Engine</dt>
              <dd data-detail="path">
                {{ run.record.runtimePath }}, {{ run.record.threads }}
                {{ run.record.threads === 1 ? 'thread' : 'threads' }}
                (requested: {{ run.snapshot.requestedThreads === 'auto' ? 'Auto' : run.snapshot.requestedThreads }})
              </dd>
            </template>
            <template v-if="run.record.fallbackReason">
              <dt>Serial because</dt>
              <dd data-detail="fallback">{{ run.record.fallbackReason }}</dd>
            </template>
            <template v-if="run.record.subjectRetained !== undefined">
              <dt>Subject</dt>
              <dd data-detail="subject-retained">
                {{ run.record.subjectRetained ? 'kept from the previous search' : 'read for this search' }}
              </dd>
            </template>
            <template v-if="run.record.engineBuild">
              <dt>Build</dt>
              <dd>{{ run.record.engineBuild }}</dd>
            </template>
            <template v-if="run.record.runtimeGeneration !== undefined">
              <dt>Runtime</dt>
              <dd>generation {{ run.record.runtimeGeneration }}</dd>
            </template>
            <template v-if="run.record.memory">
              <dt>Engine memory</dt>
              <dd>
                {{ formatBytes(run.record.memory.linearBytesAfter) }} after the search ({{ run.record.memory.instanceRuns }}
                {{ run.record.memory.instanceRuns === 1 ? 'search' : 'searches' }} on this instance)
              </dd>
            </template>
            <template v-if="phaseTimes(run)">
              <dt>Phases</dt>
              <dd>{{ phaseTimes(run) }}</dd>
            </template>
          </dl>
        </details>
      </li>
    </ol>
  </section>
</template>
