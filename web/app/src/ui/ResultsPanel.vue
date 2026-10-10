<script setup lang="ts">
// The results screen (plan §5.7, design §3.4), in the order and words of NCBI BLAST's results
// page (S13b, docs/web/ncbi_ui_mapping.md §2): the run's header block with "Filter Results"
// beside it, "Results for" (the query) and the notices, then the tabs Descriptions, Graphic
// Summary, Alignments and Dot Plot, and LOSAT's Run details and Outputs. Every tab follows one
// selection, held by the HSP's identity (application/results.ts). The Descriptions, the
// Alignments and the dot plot's popup add HSPs to the candidate tray (S14); a short line
// confirms each addition.
import { computed, nextTick, onUnmounted, ref, watch } from 'vue';
import type { CandidateTray } from '../application/candidates';
import type { Coordinator, RunView } from '../application/coordinator';
import type { HspId, ResultsBrowser } from '../application/results';
import type { ResultExporter } from '../application/result-export';
import type { Session } from '../application/session';
import { programById, residueUnit } from '../domain/programs';
import { useStore } from './useStore';
import AlignmentsView from './AlignmentsView.vue';
import CommandText from './CommandText.vue';
import DotPlot from './DotPlot.vue';
import { formatCount, formatCounted } from './format';
import GraphicSummary from './GraphicSummary.vue';
import HspTable from './HspTable.vue';
import OutputsView from './OutputsView.vue';
import QueryPicker from './QueryPicker.vue';
import ResultFilters from './ResultFilters.vue';
import ResultNotices from './ResultNotices.vue';
import RunDetails from './RunDetails.vue';
import SubjectTable from './SubjectTable.vue';
import VerificationBadge from './VerificationBadge.vue';

export type ResultsView = 'hits' | 'graphic' | 'alignment' | 'dotplot' | 'details' | 'outputs';

const props = defineProps<{
  coordinator: Coordinator;
  results: ResultsBrowser;
  candidates: CandidateTray;
  exporter: ResultExporter;
  runs: readonly RunView[];
  /** Session files: a loaded run's origin and the re-attachment of its original FASTA in Run details. */
  session?: Session;
}>();
/** The tab shown. The main view keeps it, so that another run, or the results shown again, open on the same tab. */
const view = defineModel<ResultsView>('view', { default: 'hits' });
const state = useStore(props.results.state);
const trayState = useStore(props.candidates.state);
/** The keys of the HSPs in the tray: the Alignments and the dot plot say "In candidates" for them. */
const inTray = computed<ReadonlySet<string>>(() => new Set(trayState.value.candidates.map((candidate) => candidate.key)));
const heading = ref<HTMLElement>();
const alignments = ref<InstanceType<typeof AlignmentsView>>();

/** NCBI's tabs, then LOSAT's. The test IDs are W4's (the hits view, the HSP's panes, the run's views). */
const TABS: readonly { readonly view: ResultsView; readonly label: string; readonly testid: string }[] = [
  { view: 'hits', label: 'Descriptions', testid: 'results-view-hits' },
  { view: 'graphic', label: 'Graphic Summary', testid: 'results-view-graphic' },
  { view: 'alignment', label: 'Alignments', testid: 'pane-alignment' },
  { view: 'dotplot', label: 'Dot Plot', testid: 'pane-dotplot' },
  { view: 'details', label: 'Run details', testid: 'results-view-details' },
  { view: 'outputs', label: 'Outputs', testid: 'results-view-outputs' },
];
const HITS_VIEWS: readonly ResultsView[] = ['hits', 'graphic', 'alignment', 'dotplot'];

/**
 * Brings the panel's heading into view and moves the focus to it (a run opened from the
 * queue). The scroll is instant: Firefox ended a smooth one early while the run was read.
 */
function showHeading(): void {
  heading.value?.scrollIntoView({ block: 'start' });
  heading.value?.focus({ preventScroll: true });
}
/** Shows an HSP's Range in the Alignments with the focus on it ("Show in results" of a candidate). */
async function showHsp(id: HspId): Promise<void> {
  view.value = 'alignment';
  await nextTick();
  await alignments.value?.reveal(id, true);
}
defineExpose({ showHeading, showHsp });

// --- adding to the candidate tray ----------------------------------------------------------------

/** The confirmation of the last addition (or why it was refused), shown for a few seconds. */
const added = ref<{ readonly text: string; readonly error: boolean }>();
const CONFIRMATION_MS = 4000;
let confirmationTimer: ReturnType<typeof setTimeout> | undefined;
onUnmounted(() => clearTimeout(confirmationTimer));

function addCandidates(ids: readonly HspId[]): void {
  let text: string;
  let error = false;
  try {
    const result = props.candidates.add(props.results.candidateSources(ids));
    if (!result.ok) {
      text = result.message;
      error = true;
    } else if (result.added === 0) {
      text = result.already === 1 ? 'This HSP is already in Candidates.' : 'These HSPs are already in Candidates.';
    } else {
      text = `${formatCounted(result.added, 'HSP')} added to Candidates.`;
      if (result.already > 0) text += ` ${formatCount(result.already)} ${result.already === 1 ? 'was' : 'were'} already there.`;
    }
  } catch (failure) {
    text = `The HSPs could not be added: ${failure instanceof Error ? failure.message : String(failure)}`;
    error = true;
  }
  added.value = { text, error };
  clearTimeout(confirmationTimer);
  confirmationTimer = setTimeout(() => (added.value = undefined), CONFIRMATION_MS);
}

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

/** The program and the task that the search used (the argv's, or the engine's default once the run is read). */
const programText = computed(() => {
  const run = selected.value;
  if (run === undefined) return '';
  const task = loaded.value?.run.snapshot.runId === run.snapshot.runId ? loaded.value.task : undefined;
  const label = programById(run.snapshot.program).label;
  return task === undefined ? label : `${label} (task ${task})`;
});
const units = computed(() => {
  const program = selected.value === undefined ? undefined : programById(selected.value.snapshot.program);
  return program === undefined ? { query: '', subject: '' } : { query: residueUnit(program.query), subject: residueUnit(program.subject) };
});
const plural = (count: number, one: string) => `${formatCount(count)} ${count === 1 ? one : `${one}s`}`;
/** "Results for" (the query list) is for runs of more than one query, as NCBI's. */
const multiQuery = computed(() => (loaded.value?.run.snapshot.query.records.length ?? 0) > 1);

/** "Show alignment" of the Graphic Summary and the Dot Plot: the Alignments tab, with the HSP's Range in view. */
async function toAlignments(id: HspId): Promise<void> {
  view.value = 'alignment';
  await nextTick();
  alignments.value?.reveal(id, true);
}
</script>

<template>
  <section class="results" data-testid="results-panel">
    <h2 ref="heading" tabindex="-1" data-testid="results-heading">Results</h2>
    <p v-if="runs.length === 0" class="muted" data-testid="results-empty">No runs yet. Searches appear here when they complete.</p>
    <template v-else>
      <div class="results-top">
        <dl class="results-summary" data-testid="results-summary">
          <template v-if="selected?.snapshot.title">
            <dt>Job Title</dt>
            <dd class="results-title" data-testid="results-job-title">{{ selected.snapshot.title }}</dd>
          </template>
          <dt>Run</dt>
          <dd class="results-run">
            <label class="results-run-select">
              <span class="visually-hidden">Run</span>
              <select
                :value="state.runId ?? ''"
                data-testid="results-run"
                @change="results.open(($event.target as HTMLSelectElement).value)"
              >
                <option value="" disabled>Choose a run</option>
                <option v-for="run in choices" :key="run.snapshot.runId" :value="run.snapshot.runId">{{ runLabel(run) }}</option>
              </select>
            </label>
            <button v-if="state.phase === 'ready'" type="button" class="link" data-testid="results-download-all" @click="view = 'outputs'">
              Download All
            </button>
          </dd>
          <template v-if="selected">
            <dt>Program</dt>
            <dd class="results-program">
              <span data-testid="results-program">{{ programText }}</span>
              <VerificationBadge v-if="loaded" :badge="loaded.badge" compact @details="view = 'details'" />
            </dd>
            <dt>Options</dt>
            <dd class="run-inputs" data-testid="results-run-options">
              <CommandText :text="options(selected)" />
              <template v-if="selected.snapshot.group">
                · group run {{ selected.snapshot.group.position }} of {{ selected.snapshot.group.size }}
              </template>
            </dd>
            <template v-if="selected.fromSession">
              <dt>Loaded</dt>
              <dd class="run-inputs" data-testid="results-origin">
                from {{ selected.fromSession.fileName }}, run {{ selected.fromSession.number }} there (not searched again)
              </dd>
            </template>
            <template v-if="selected.snapshot.query.records.length === 1">
              <dt>Query ID</dt>
              <dd data-testid="results-query-id">{{ selected.snapshot.query.records[0]!.id }}</dd>
              <dt>Query Length</dt>
              <dd data-testid="results-query-length">{{ formatCount(selected.snapshot.query.records[0]!.length) }} {{ units.query }}</dd>
            </template>
            <template v-if="selected.snapshot.subject.records.length === 1">
              <dt>Subject ID</dt>
              <dd data-testid="results-subject-id">{{ selected.snapshot.subject.records[0]!.id }}</dd>
              <dt>Subject Length</dt>
              <dd data-testid="results-subject-length">{{ formatCount(selected.snapshot.subject.records[0]!.length) }} {{ units.subject }}</dd>
            </template>
            <template v-else>
              <dt>Subjects</dt>
              <dd data-testid="results-subjects">
                {{ selected.snapshot.subject.name }}, {{ plural(selected.snapshot.subject.records.length, 'record') }}
              </dd>
            </template>
          </template>
        </dl>
        <ResultFilters v-if="state.phase === 'ready' && loaded" :results="results" :filters="state.filters" />
      </div>

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
        <QueryPicker v-if="multiQuery" :results="results" :state="state" />
        <ResultNotices :state="state" @clear="results.setFilters({ queriesWithHitsOnly: state.filters.queriesWithHitsOnly ?? false })" />

        <nav class="tabs results-tabs" aria-label="Results view">
          <button
            v-for="tab in TABS"
            :key="tab.view"
            type="button"
            :aria-pressed="view === tab.view"
            :data-testid="tab.testid"
            @click="view = tab.view"
          >
            {{ tab.label }}
          </button>
        </nav>

        <div v-show="HITS_VIEWS.includes(view)" class="hits-view" data-testid="results-hits" :data-run="loaded.run.snapshot.number">
          <SubjectTable v-if="view === 'hits' && state.subjects.length > 0" :results="results" :state="state" @add-candidates="addCandidates" />
          <GraphicSummary
            v-else-if="view === 'graphic' && state.subjects.length > 0"
            :results="results"
            :state="state"
            @show-alignment="toAlignments"
          />
          <AlignmentsView
            v-else-if="view === 'alignment' && state.hsps.length > 0"
            ref="alignments"
            :results="results"
            :state="state"
            :in-tray="inTray"
            @descriptions="view = 'hits'"
            @add-candidates="addCandidates"
          />
          <template v-else-if="view === 'dotplot' && state.hsps.length > 0">
            <DotPlot :results="results" :state="state" :in-tray="inTray" @show-alignment="toAlignments" @add-candidates="addCandidates" />
            <HspTable :results="results" :state="state" />
          </template>
        </div>
        <RunDetails v-if="view === 'details'" :run="loaded.run" :loaded="loaded" :session="session" />
        <OutputsView v-if="view === 'outputs'" :coordinator="coordinator" :run="loaded.run" :exporter="exporter" :state="state" />
      </template>
    </template>
    <!-- A live region that stays in the page, so that each confirmation is announced. -->
    <div class="confirmation-host" role="status">
      <p v-if="added" class="confirmation" :class="{ error: added.error }" data-testid="candidates-added" :data-error="added.error">
        {{ added.text }}
      </p>
    </div>
  </section>
</template>
