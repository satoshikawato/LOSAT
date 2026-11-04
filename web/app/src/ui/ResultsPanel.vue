<script setup lang="ts">
// The results screen (plan §5.7, design §3.4), in the order and words of NCBI BLAST's results
// page (S13b, docs/web/ncbi_ui_mapping.md §2): the run's header block with "Filter Results"
// beside it, "Results for" (the query) and the notices, then the hits as NCBI's classic page
// shows them, on one page from top to bottom: the Graphic Summary, the Descriptions and the
// Alignments of the selected subject (Owner, 2026-10-10 and 2026-11-02: not behind tabs). The
// Dot Plot and LOSAT's Run details and Outputs are tabs beside the hits. A run that has not
// finished shows its status here (RunStatus) until it completes. Every view follows one
// selection, held by the HSP's identity (application/results.ts).
import { computed, nextTick, onMounted, onUnmounted, ref, useId, watch } from 'vue';
import type { Coordinator, RunView } from '../application/coordinator';
import { OutputSaver } from '../application/outputs';
import { noHspReason, type HspId, type ResultsBrowser } from '../application/results';
import { recordName } from '../domain/dataset';
import { CUSTOM_OUTPUT, customOutfmt, type RunOutput } from '../domain/output-format';
import { programById, residueUnit } from '../domain/programs';
import { isTerminal } from '../domain/run';
import { useStore } from './useStore';
import AlignmentsView from './AlignmentsView.vue';
import CommandText from './CommandText.vue';
import DotPlot from './DotPlot.vue';
import { formatCount } from './format';
import GraphicSummary from './GraphicSummary.vue';
import HspTable from './HspTable.vue';
import LosatLoader from './LosatLoader.vue';
import { NO_HSP_TEXT } from './no-hsp-text';
import OutputsView from './OutputsView.vue';
import QueryPicker from './QueryPicker.vue';
import ResultFilters from './ResultFilters.vue';
import ResultNotices from './ResultNotices.vue';
import RunDetails from './RunDetails.vue';
import RunStatus from './RunStatus.vue';
import SubjectTable from './SubjectTable.vue';
import VerificationBadge from './VerificationBadge.vue';

export type ResultsView = 'hits' | 'dotplot' | 'details' | 'outputs';

const props = defineProps<{ coordinator: Coordinator; results: ResultsBrowser; runs: readonly RunView[] }>();
/** "Run a search" of the empty results page: the Search tab. */
const emit = defineEmits<{ search: [] }>();
/** The tab shown. The main view keeps it, so that another run, or the results shown again, open on the same tab. */
const view = defineModel<ResultsView>('view', { default: 'hits' });
const state = useStore(props.results.state);
const heading = ref<HTMLElement>();
const alignments = ref<InstanceType<typeof AlignmentsView>>();

/**
 * The hits (Graphic Summary, Descriptions and Alignments on one page), NCBI's Dot Plot, then
 * LOSAT's views. The test IDs are W4's (the hits view, the dot plot's pane, the run's views).
 */
const TABS: readonly { readonly view: ResultsView; readonly label: string; readonly testid: string }[] = [
  { view: 'hits', label: 'Hits', testid: 'results-view-hits' },
  { view: 'dotplot', label: 'Dot Plot', testid: 'pane-dotplot' },
  { view: 'details', label: 'Run details', testid: 'results-view-details' },
  { view: 'outputs', label: 'Outputs', testid: 'results-view-outputs' },
];

const emptyPlot = computed(() => (state.value.hsps.length > 0 ? undefined : noHspReason(state.value)));

/**
 * "Download All" lists the outputs of the run, as NCBI's lists its formats, by BLAST+'s names of
 * them (its -outfmt help), the output of the run's own -outfmt with its value (W10); the one
 * chosen is saved as one file (WA-RS-05). The list opens in the page's flow, under the run's line,
 * so that it fits a phone's width; an output that is not saved says why there, until the next
 * download or another run (WA2-BR-01).
 */
const FORMAT_NAMES = { 0: 'Pairwise', 6: 'Tabular', 7: 'Tabular with comment lines' } as const;
const downloads = computed(() => {
  const items: { readonly format: RunOutput; readonly label: string; readonly testid: string }[] = ([0, 6, 7] as const).map((format) => ({
    format,
    label: `${FORMAT_NAMES[format]} (outfmt ${format})`,
    testid: `results-download-${format}`,
  }));
  const custom = customOutfmt(loaded.value?.run.snapshot.argv ?? []);
  if (custom !== undefined) {
    const name = FORMAT_NAMES[Number(custom.split(' ', 1)[0]) as 6 | 7];
    items.push({ format: CUSTOM_OUTPUT, label: `${name} (outfmt ${custom})`, testid: 'results-download-custom' });
  }
  return items;
});
const downloadsOpen = ref(false);
const downloadsId = useId();
const downloadsButton = ref<HTMLButtonElement>();
const downloadsList = ref<HTMLElement>();
const saver = new OutputSaver(props.coordinator);
const downloadProblem = useStore(saver.problem);

function download(format: RunOutput): void {
  downloadsOpen.value = false;
  // The chosen item leaves the page with the list; the focus stays on "Download All".
  downloadsButton.value?.focus();
  const run = loaded.value?.run;
  if (run !== undefined) void saver.save(run, format);
}

/** Escape closes the list, and the focus goes back to its button (Safari does not focus a clicked button). */
function onDownloadsKey(event: KeyboardEvent): void {
  if (event.key !== 'Escape' || !downloadsOpen.value || event.defaultPrevented) return;
  event.preventDefault();
  downloadsOpen.value = false;
  downloadsButton.value?.focus();
}

/** A press outside the list and its button closes the list. */
function onPointerDown(event: PointerEvent): void {
  const target = event.target as Node | null;
  if (!downloadsOpen.value || downloadsList.value?.contains(target) || downloadsButton.value?.contains(target)) return;
  downloadsOpen.value = false;
}
onMounted(() => {
  document.addEventListener('pointerdown', onPointerDown);
  document.addEventListener('keydown', onDownloadsKey);
});
onUnmounted(() => {
  document.removeEventListener('pointerdown', onPointerDown);
  document.removeEventListener('keydown', onDownloadsKey);
});

/**
 * Brings the panel's heading into view and moves the focus to it (a search just started).
 * The scroll is instant: Firefox ended a smooth one early while the run was read.
 */
function showHeading(): void {
  heading.value?.scrollIntoView({ block: 'start' });
  heading.value?.focus({ preventScroll: true });
}
defineExpose({ showHeading });

/** Every run of the working session, newest first: a run without results says why. */
const choices = computed(() => [...props.runs].reverse());
const selected = computed(() => props.runs.find((run) => run.snapshot.runId === state.value.runId));
const loaded = computed(() => state.value.loaded);
/** The selected run waits or searches: its status is shown in place of its results. */
const active = computed(() => selected.value !== undefined && !isTerminal(selected.value.status));

// Another run, or the same run read again, starts with the list closed and without the reason of a download.
watch([() => state.value.runId, () => state.value.phase], () => {
  downloadsOpen.value = false;
  saver.reset();
});

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
  const status = run.status === 'completed' ? '' : ` (${run.status === 'queued' ? 'waiting' : run.status})`;
  return `Run ${number} · ${programById(program).label} · ${query.name} vs ${subject.name}${status}`;
}

/**
 * The options of the run's argv: the words after the program and the two inputs (plan §5.3),
 * which CommandText quotes as on a command line, so that they can be typed again ("Additional options").
 */
function options(run: RunView): readonly string[] {
  return run.snapshot.argv.slice(5);
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

/** "Show alignment" of the Graphic Summary and the Dot Plot: the Alignments, with the HSP's Range in view. */
async function toAlignments(id: HspId): Promise<void> {
  view.value = 'hits';
  await nextTick();
  alignments.value?.reveal(id, true);
}
</script>

<template>
  <section class="results" data-testid="results-panel">
    <h2 ref="heading" tabindex="-1" data-testid="results-heading">Results</h2>
    <p v-if="runs.length === 0" data-testid="results-empty">
      <button type="button" class="primary-action" data-testid="results-run-a-search" @click="emit('search')">Run a search</button>
    </p>
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
                <option
                  v-for="run in choices"
                  :key="run.snapshot.runId"
                  :value="run.snapshot.runId"
                  :data-testid="`run-${run.snapshot.number}-option`"
                  :data-status="run.status"
                >
                  {{ runLabel(run) }}
                </option>
              </select>
            </label>
            <button
              v-if="state.phase === 'ready'"
              ref="downloadsButton"
              type="button"
              class="link"
              data-testid="results-download-all"
              :aria-expanded="downloadsOpen"
              :aria-controls="downloadsId"
              @click="downloadsOpen = !downloadsOpen"
            >
              Download All
            </button>
            <ul v-if="state.phase === 'ready' && downloadsOpen" :id="downloadsId" ref="downloadsList" class="download-list" data-testid="results-download-menu">
              <li v-for="item in downloads" :key="item.format">
                <button type="button" :data-testid="item.testid" @click="download(item.format)">{{ item.label }}</button>
              </li>
            </ul>
            <p v-if="state.phase === 'ready' && downloadProblem" class="results-message error" role="alert" data-testid="results-download-problem">
              {{ downloadProblem }}
            </p>
          </dd>
          <template v-if="selected">
            <dt>Mode</dt>
            <dd class="results-program">
              <span data-testid="results-program">{{ programText }}</span>
              <VerificationBadge v-if="loaded" :badge="loaded.badge" compact @details="view = 'details'" />
            </dd>
            <dt>Options</dt>
            <dd class="run-inputs" data-testid="results-run-options">
              <CommandText v-if="options(selected).length > 0" :words="options(selected)" :takes-value="state.grammar?.takesValue" />
              <template v-else>defaults</template>
              <template v-if="selected.snapshot.group">
                · group run {{ selected.snapshot.group.position }} of {{ selected.snapshot.group.size }}
              </template>
            </dd>
            <template v-if="selected.snapshot.query.records.length === 1">
              <dt>Query ID</dt>
              <dd data-testid="results-query-id">{{ recordName(selected.snapshot.query.records[0]!.id, 'query', 0) }}</dd>
              <dt>Query Length</dt>
              <dd data-testid="results-query-length">{{ formatCount(selected.snapshot.query.records[0]!.length) }} {{ units.query }}</dd>
            </template>
            <template v-if="selected.snapshot.subject.records.length === 1">
              <dt>Subject ID</dt>
              <dd data-testid="results-subject-id">{{ recordName(selected.snapshot.subject.records[0]!.id, 'subject', 0) }}</dd>
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

      <p v-if="state.phase === 'loading'" class="muted" data-testid="results-status" data-phase="loading"><LosatLoader />Reading the results…</p>
      <RunStatus v-else-if="state.phase === 'unavailable' && active && selected" :coordinator="coordinator" :run="selected" :runs="runs" />
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

        <!-- NCBI's classic results page: the three sections on one page, in its order. The Alignments
             show the selected subject; a row of the Descriptions or a bar of the Graphic Summary selects it. -->
        <div v-show="view === 'hits'" class="hits-view" data-testid="results-hits" :data-run="loaded.run.snapshot.number">
          <template v-if="state.subjects.length > 0">
            <section class="hits-section" data-testid="section-graphic">
              <h3 class="section-heading">Graphic Summary</h3>
              <GraphicSummary :results="results" :state="state" @show-alignment="toAlignments" />
            </section>
            <section class="hits-section" data-testid="section-descriptions">
              <h3 class="section-heading">Descriptions</h3>
              <SubjectTable :results="results" :state="state" />
            </section>
            <section v-if="state.hsps.length > 0" class="hits-section" data-testid="section-alignments">
              <h3 class="section-heading">Alignments</h3>
              <AlignmentsView ref="alignments" :results="results" :state="state" />
            </section>
          </template>
        </div>
        <div v-if="view === 'dotplot' && state.hsps.length > 0" class="hits-view results-group" data-testid="results-dotplot">
          <DotPlot :results="results" :state="state" @show-alignment="toAlignments" />
          <HspTable :results="results" :state="state" />
        </div>
        <p v-else-if="view === 'dotplot' && emptyPlot" class="results-message" data-testid="results-dotplot-empty" :data-kind="emptyPlot">
          {{ NO_HSP_TEXT[emptyPlot] }}
        </p>
        <div v-if="view === 'details'" class="results-group">
          <RunDetails :run="loaded.run" :loaded="loaded" />
        </div>
        <div v-if="view === 'outputs'" class="results-group">
          <OutputsView :coordinator="coordinator" :run="loaded.run" />
        </div>
      </template>
    </template>
  </section>
</template>
