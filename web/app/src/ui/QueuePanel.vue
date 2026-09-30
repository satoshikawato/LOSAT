<script setup lang="ts">
import type { Coordinator, RunView } from '../application/coordinator';
import { isTerminal } from '../domain/run';
import { programById } from '../domain/programs';

defineProps<{ coordinator: Coordinator; runs: readonly RunView[] }>();

const cancellable = (run: RunView) => !isTerminal(run.status) && run.status !== 'finalizing';
</script>

<template>
  <h2>Queue</h2>
  <p v-if="runs.length === 0">No runs yet.</p>
  <ol class="queue" data-testid="queue">
    <li v-for="run in runs" :key="run.snapshot.runId" :data-testid="`run-${run.snapshot.number}`">
      <span>Run {{ run.snapshot.number }} · {{ programById(run.snapshot.program).label }}</span>
      <span class="status" :data-testid="`run-${run.snapshot.number}-status`">{{ run.status }}</span>
      <button v-if="cancellable(run)" @click="coordinator.cancel(run.snapshot.runId)">Cancel</button>
      <span v-if="run.record.error" class="error">{{ run.record.error }}</span>
    </li>
  </ol>
</template>
