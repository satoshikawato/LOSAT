<script setup lang="ts">
// The results screen (plan §5.7, design §3.4): Run -> Query -> Subject list -> HSP, with
// the HSP's alignment and the pair's dot plot, the run's details and its outputs. Every
// part follows one selection, held by the HSP's identity (application/results.ts).
import { computed, ref, watch } from 'vue';
import type { Coordinator, RunView } from '../application/coordinator';
import type { ResultsBrowser } from '../application/results';
import { programById } from '../domain/programs';
import { useStore } from './useStore';
import DotPlot from './DotPlot.vue';
import HspDetail from './HspDetail.vue';
import HspTable from './HspTable.vue';
import OutputsView from './OutputsView.vue';
import QueryPicker from './QueryPicker.vue';
import ResultFilters from './ResultFilters.vue';
import ResultNotices from './ResultNotices.vue';
import RunDetails from './RunDetails.vue';
import SubjectTable from './SubjectTable.vue';
import VerificationBadge from './VerificationBadge.vue';

const props = defineProps<{ coordinator: Coordinator; results: ResultsBrowser; runs: readonly RunView[] }>();
const state = useStore(props.results.state);
const view = ref<'hits' | 'details' | 'outputs'>('hits');
const pane = ref<'alignment' | 'dotplot'>('alignment');

/** Every run of the working session, newest first: a run without results says why. */
const choices = computed(() => [...props.runs].reverse());
const selected = computed(() => props.runs.find((run) => run.snapshot.runId === state.value.runId));
const loaded = computed(() => state.value.loaded);

// With no run chosen, show the newest completed one.
watch(
  () => props.runs.filter((run) => run.status === 'completed').at(-1)?.snapshot.runId,
  (runId) => {
    if (runId !== undefined && state.value.runId === undefined) void props.results.open(runId);
  },
  { immediate: true },
);

function runLabel(run: RunView): string {
  const { number, program, query, subject } = run.snapshot;
  const status = run.status === 'completed' ? '' : ` (${run.status})`;
  return `Run ${number} · ${programById(program).label} · ${query.name} vs ${subject.name}${status}`;
}

/** The options of the run's argv: the words after the program and the two inputs (plan §5.3). */
function options(run: RunView): string {
  const words = run.snapshot.argv.slice(5);
  return words.length === 0 ? 'defaults' : words.join(' ');
}
</script>

<template>
  <section class="results" data-testid="results-panel">
    <h2>Results</h2>
    <p v-if="runs.length === 0" class="muted" data-testid="results-empty">No runs yet. Searches appear here when they complete.</p>
    <template v-else>
      <div class="results-run">
        <label class="results-run-select">
          <span>Run</span>
          <select
            :value="state.runId ?? ''"
            data-testid="results-run"
            @change="results.open(($event.target as HTMLSelectElement).value)"
          >
            <option value="" disabled>Choose a run</option>
            <option v-for="run in choices" :key="run.snapshot.runId" :value="run.snapshot.runId">{{ runLabel(run) }}</option>
          </select>
        </label>
        <VerificationBadge v-if="loaded" :badge="loaded.badge" compact @details="view = 'details'" />
      </div>
      <p v-if="selected" class="run-inputs muted" data-testid="results-run-options">
        Options: <span class="word">{{ options(selected) }}</span>
        <template v-if="selected.snapshot.group">
          · group run {{ selected.snapshot.group.position }} of {{ selected.snapshot.group.size }}
        </template>
      </p>

      <p v-if="state.phase === 'loading'" class="muted" data-testid="results-status" data-phase="loading">Reading the results…</p>
      <p
        v-else-if="state.phase === 'unavailable' || state.phase === 'failed'"
        class="results-message"
        :class="{ error: selected?.status === 'failed' || state.phase === 'failed' }"
        data-testid="results-status"
        :data-phase="state.phase"
        :data-run-status="selected?.status"
      >
        {{ state.message }}
      </p>

      <template v-if="state.phase === 'ready' && loaded">
        <nav class="tabs" aria-label="Results view">
          <button :aria-pressed="view === 'hits'" data-testid="results-view-hits" @click="view = 'hits'">Hits</button>
          <button :aria-pressed="view === 'details'" data-testid="results-view-details" @click="view = 'details'">Run details</button>
          <button :aria-pressed="view === 'outputs'" data-testid="results-view-outputs" @click="view = 'outputs'">Outputs</button>
        </nav>

        <div v-show="view === 'hits'" class="hits-view" data-testid="results-hits" :data-run="loaded.run.snapshot.number">
          <div class="hits-top">
            <QueryPicker :results="results" :state="state" />
            <ResultFilters :results="results" :filters="state.filters" />
          </div>
          <ResultNotices :state="state" @clear="results.setFilters({ queriesWithHitsOnly: state.filters.queriesWithHitsOnly ?? false })" />
          <SubjectTable v-if="state.subjects.length > 0" :results="results" :state="state" />
          <template v-if="state.hsps.length > 0">
            <HspTable :results="results" :state="state" />
            <nav class="tabs" aria-label="Selected HSP">
              <button :aria-pressed="pane === 'alignment'" data-testid="pane-alignment" @click="pane = 'alignment'">Alignment</button>
              <button :aria-pressed="pane === 'dotplot'" data-testid="pane-dotplot" @click="pane = 'dotplot'">Dot plot</button>
            </nav>
            <HspDetail v-if="pane === 'alignment'" :state="state" />
            <DotPlot v-else :results="results" :state="state" />
          </template>
        </div>
        <RunDetails v-if="view === 'details'" :run="loaded.run" :loaded="loaded" />
        <OutputsView v-if="view === 'outputs'" :coordinator="coordinator" :run="loaded.run" />
      </template>
    </template>
  </section>
</template>
