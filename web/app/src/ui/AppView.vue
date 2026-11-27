<script setup lang="ts">
import { computed, nextTick, onMounted, onUnmounted, ref, shallowRef, useId, watchEffect, type Ref } from 'vue';
import { holdsResults, type Attention } from '../application/attention';
import { withRemedy, type Availability } from '../application/availability';
import type { BlastFormatterDraft } from '../application/blast-formatter-draft';
import type { Coordinator, MadeDatabase, MaskDataHandover, SessionArchiveEntry, SessionCountsEntry } from '../application/coordinator';
import type { SearchDraft } from '../application/draft';
import type { CountsLoader } from '../application/counts';
import type { ExampleLoader } from '../application/examples';
import type { MakeblastdbDraft } from '../application/makeblastdb-draft';
import type { ResultsBrowser } from '../application/results';
import type { ToolDraft } from '../application/tool-draft';
import { useStore } from './useStore';
import FormatArchivePanel from './FormatArchivePanel.vue';
import LiveAnnouncer from './LiveAnnouncer.vue';
import LosatLogo from './LosatLogo.vue';
import MakeblastdbPanel from './MakeblastdbPanel.vue';
import { holdWhilePressed } from './press-hold';
import ResultsPanel, { type ResultsView } from './ResultsPanel.vue';
import ResumeNotice from './ResumeNotice.vue';
import SearchPanel from './SearchPanel.vue';
import { useRunAnnouncements } from './run-announcements';
import StorageStatus from './StorageStatus.vue';
import { tabNavigation, type MainTab } from './tab-history';
import TabList from './TabList.vue';
import WindowMaskerPanel from './WindowMaskerPanel.vue';

const props = defineProps<{
  coordinator: Coordinator;
  draft: SearchDraft;
  toolDraft: ToolDraft;
  blastFormatterDraft: BlastFormatterDraft;
  makeblastdbDraft: MakeblastdbDraft;
  examples: ExampleLoader;
  counts: CountsLoader;
  results: ResultsBrowser;
  attention: Attention;
  availability: Availability;
  usesFakeEngine: boolean;
}>();
const state = useStore(props.coordinator.state);
const attentionState = useStore(props.attention.state);
const availabilityState = useStore(props.availability.state);
/**
 * What of LOSAT could not be loaded or stopped, and the remedy, once for the page (line WQA's sweep, D1 and D2): a
 * reload also loses the results of this session, which the page then says.
 */
const unavailable = computed(() => {
  const problem = availabilityState.value.problem;
  return problem === undefined ? undefined : withRemedy(problem, state.value.runs.some(holdsResults));
});
/** The tab shown. Choosing one goes through the page's history, so Back on the results returns to the form (tab-history.ts). */
const shown = ref<MainTab>('search');
const navigation = tabNavigation(window.history);
/** How far each form tab was scrolled, and what had the focus, when it was left: Back returns there, as between pages. */
const forms = {
  search: { scroll: 0, focus: null as Element | null },
  windowmasker: { scroll: 0, focus: null as Element | null },
  makeblastdb: { scroll: 0, focus: null as Element | null },
  blast_formatter: { scroll: 0, focus: null as Element | null },
};
function leave(): void {
  if (shown.value === 'results') return;
  forms[shown.value].scroll = window.scrollY;
  forms[shown.value].focus = document.activeElement;
}
const tab = computed<MainTab>({
  get: () => shown.value,
  set: (value) => {
    if (value !== shown.value) leave();
    shown.value = value;
    navigation.choose(value);
  },
});
// Load tool `-xmlhelp` descriptions (and saved counts) lazily on the first visit to each tool tab
// so startup on Search does not issue extra worker calls.
watchEffect(() => {
  if (shown.value === 'windowmasker') {
    props.toolDraft.ensureLoaded();
    props.counts.refreshList();
  } else if (shown.value === 'makeblastdb') {
    props.makeblastdbDraft.ensureLoaded();
  } else if (shown.value === 'blast_formatter') {
    props.blastFormatterDraft.ensureLoaded();
  }
});
const searchPanel = ref<InstanceType<typeof SearchPanel>>();
const resultsPanel = ref<InstanceType<typeof ResultsPanel>>();
const mainTabs = ref<{ focus: (options?: FocusOptions) => void }>();
/** Search, Results, WindowMasker, makeblastdb, and Format Archive: tabs, each with its panel (WA2-AX-12). */
const ids = useId();
const MAIN_TABS = [
  { value: 'search', label: 'Search', testid: 'tab-search', id: `${ids}-search-tab`, panel: `${ids}-search` },
  { value: 'results', label: 'Results', testid: 'tab-results', id: `${ids}-results-tab`, panel: `${ids}-results` },
  { value: 'windowmasker', label: 'WindowMasker', testid: 'tab-windowmasker', id: `${ids}-wm-tab`, panel: `${ids}-wm` },
  { value: 'makeblastdb', label: 'makeblastdb', testid: 'tab-makeblastdb', id: `${ids}-mk-tab`, panel: `${ids}-mk` },
  { value: 'blast_formatter', label: 'Format Archive', testid: 'tab-blast-formatter', id: `${ids}-bf-tab`, panel: `${ids}-bf` },
] as const;
useRunAnnouncements(state, useStore(props.results.state));
/** The results tab's view, kept while the search tab is shown and when another run opens. */
const resultsView = ref<ResultsView>('hits');

/**
 * Opens a run's results: "Run LOSAT" goes straight to the results page of the run it started, as
 * NCBI BLAST's "BLAST" button does (Owner, 2026-11-02); the page shows the run's status until it
 * completes, then its results. The results panel comes into view with the focus on its heading.
 */
async function openResults(runId: string): Promise<void> {
  tab.value = 'results';
  void props.results.open(runId);
  await nextTick();
  resultsPanel.value?.showHeading();
}

/**
 * "Dismiss" leaves the page with its notice: the focus goes to the tab shown, the next control
 * (WA2-AX-04), and not to the page's body, from where Tab would start again at the top.
 */
async function dismissResume(): Promise<void> {
  props.attention.dismissResume();
  await nextTick();
  const active = document.activeElement;
  if (active === null || active === document.body || !active.isConnected) mainTabs.value?.focus();
}

/** "Run a search" of the empty Results tab leaves with the tab: the focus goes to the Search tab (WA2-AX-04). */
async function showSearch(): Promise<void> {
  tab.value = 'search';
  await nextTick();
  mainTabs.value?.focus();
}

/** "Make counts…" in BLASTN's repeats row: switches to WindowMasker on Make counts. */
async function showMakeCounts(): Promise<void> {
  props.toolDraft.setStage('mk_counts');
  tab.value = 'windowmasker';
  await nextTick();
  window.scrollTo(0, 0);
}

/** "Use in BLASTN" on a counts output: selects the file in BLASTN's repeats row and focuses it (DESIGN §3). */
async function useCountsInBlastn(file: File): Promise<void> {
  props.draft.useCounts(file);
  tab.value = 'search';
  await nextTick();
  await searchPanel.value?.focusRepeats();
}

/** "Use in BLASTN" after saving counts for a taxonomy ID: selects `-window_masker_taxid N` in BLASTN (DESIGN §3). */
async function useTaxidInBlastn(taxid: number): Promise<void> {
  props.draft.useTaxid(taxid);
  tab.value = 'search';
  await nextTick();
  await searchPanel.value?.focusRepeats();
}

/** "Mask with these counts" on a counts output: switches to WindowMasker Mask with -ustat preselected (DESIGN §3). */
async function maskWithCounts(entry: SessionCountsEntry): Promise<void> {
  props.toolDraft.useCountsForMask(entry.runId, entry.runNumber, entry.key, entry.name, entry.size);
  tab.value = 'windowmasker';
  await nextTick();
  window.scrollTo(0, 0);
}

/** "Format archive…" on a completed `-outfmt 11` search output: opens it in Format Archive. */
async function formatArchiveFromRun(entry: SessionArchiveEntry): Promise<void> {
  props.blastFormatterDraft.useSessionArchive(entry);
  tab.value = 'blast_formatter';
  await nextTick();
  window.scrollTo(0, 0);
}

/**
 * "Use as Database (-db)" on a makeblastdb run's database: Search with the run's files as the database folder, the
 * database selected, in BLASTN for a nucleotide database and BLASTP for a protein one (line WW, WW4b2).
 */
async function useAsDatabase(made: MadeDatabase): Promise<void> {
  props.draft.setProgram(made.kind === 'nucleotide' ? 'blastn' : 'blastp');
  tab.value = 'search';
  await props.draft.setSubjectDatabaseFiles(made.files, made.names, { made: true });
  await nextTick();
  searchPanel.value?.focusDatabase();
}

/** "Use as Repeats database (-filtering_db)" on a makeblastdb run's nucleotide database: BLASTN's repeats database. */
async function useAsRepeatsDatabase(made: MadeDatabase): Promise<void> {
  props.draft.setProgram('blastn');
  tab.value = 'search';
  await props.draft.setFilteringDbFiles(made.files, made.names, { made: true });
  await nextTick();
  await searchPanel.value?.focusFilteringDb();
}

/** "Use as -mask_data in makeblastdb" on a masking run's mask data: the file joins makeblastdb's `-mask_data`. */
async function useAsMaskData(data: MaskDataHandover): Promise<void> {
  props.makeblastdbDraft.addMaskData(data.file, data.name, data.origin);
  tab.value = 'makeblastdb';
  await nextTick();
  window.scrollTo(0, 0);
}

/**
 * Back or Forward moved between the tabs' history entries (tab-history.ts). The form comes back where
 * it was scrolled, with the focus where it was (it stays mounted, so "Run LOSAT" has it again after
 * Back); a focus that is not there to give back goes to the form's tab, without scrolling, when it
 * was on the other tab or went with the results. The results come back from their heading, as
 * "Run LOSAT" shows them: their panel mounts again and its views read their content afterwards
 * (outputs, alignments), so a position further down could not be met (review of 6fff89265).
 */
async function arrived(): Promise<void> {
  const next = navigation.arrived();
  if (next === shown.value) return;
  leave();
  shown.value = next;
  await nextTick();
  if (next === 'results') {
    resultsPanel.value?.showHeading();
    return;
  }
  const saved = forms[next];
  window.scrollTo(0, saved.scroll);
  const panelId = MAIN_TABS.find((main) => main.value === next)!.panel;
  const panel = document.getElementById(panelId);
  const active = document.activeElement;
  const onTab = MAIN_TABS.some((main) => main.id === active?.id);
  if (saved.focus instanceof HTMLElement && saved.focus.isConnected && panel?.contains(saved.focus)) {
    saved.focus.focus({ preventScroll: true });
  } else if (active === null || active === document.body || !active.isConnected || onTab) {
    mainTabs.value?.focus({ preventScroll: true });
  }
}
const onPopState = () => void arrived();

/**
 * From a press until its click, the main area keeps its height, and what the page says above the tabs stays as it was,
 * so the page does not move under the press (press-hold.ts).
 */
const mainArea = ref<HTMLElement>();
const pressing = ref(false);
let stopHold: () => void = () => undefined;
/** `source` as it was when a press started, until its click; as it is between presses. */
function steady<T>(source: () => T): Readonly<Ref<T>> {
  const shown = shallowRef(source());
  watchEffect(() => {
    const now = source();
    if (!pressing.value) shown.value = now;
  });
  return shown;
}
const unavailableShown = steady(() => unavailable.value);
const resumeShown = steady(() => attentionState.value.resume);

// A file dropped outside an input's drop zone would make the browser open it in place of
// the application (and end the searches of this tab).
const keep = (event: DragEvent) => {
  if (event.dataTransfer?.types.includes('Files')) event.preventDefault();
};
onMounted(() => {
  window.addEventListener('dragover', keep);
  window.addEventListener('drop', keep);
  window.addEventListener('popstate', onPopState);
  stopHold = holdWhilePressed(
    () => mainArea.value,
    (now) => {
      pressing.value = now;
    },
  );
});
onUnmounted(() => {
  window.removeEventListener('dragover', keep);
  window.removeEventListener('drop', keep);
  window.removeEventListener('popstate', onPopState);
  stopHold();
});
</script>

<template>
  <header class="app-header">
    <h1><LosatLogo /><span class="visually-hidden">LOSAT</span></h1>
    <span class="tagline">A miniaturized reimplementation of NCBI BLAST algorithm</span>
  </header>
  <p v-if="usesFakeEngine" class="banner" data-testid="fake-engine-banner">
    Development build: a fake engine is in use. Outputs are not search results.
  </p>
  <p v-if="unavailableShown" class="error unavailable" role="alert" data-testid="losat-unavailable">{{ unavailableShown }}</p>
  <ResumeNotice v-if="resumeShown" :report="resumeShown" @dismiss="dismissResume" />
  <!-- The tabs stay in the page's navigation landmark, which they were before they were tabs (WA2-AX-12). -->
  <nav aria-label="Main">
    <TabList ref="mainTabs" v-model="tab" class="main-tabs" label="Main" :tabs="MAIN_TABS" />
  </nav>
  <!-- data-losat: whether LOSAT works now, also while a press keeps the alert above as it was (the E2E tests wait for it). -->
  <main ref="mainArea" class="layout" :data-losat="unavailable === undefined ? 'available' : 'unavailable'">
    <section class="primary">
      <!-- The search form stays mounted, so the next search keeps its edits while results are viewed. -->
      <SearchPanel
        v-show="tab === 'search'"
        ref="searchPanel"
        :id="MAIN_TABS[0].panel"
        role="tabpanel"
        :aria-labelledby="MAIN_TABS[0].id"
        :draft="draft"
        :examples="examples"
        :counts="counts"
        :attention="attention"
        :attention-state="attentionState"
        @submitted="openResults"
        @make-counts="showMakeCounts"
      />
      <ResultsPanel
        v-if="tab === 'results'"
        :id="MAIN_TABS[1].panel"
        role="tabpanel"
        :aria-labelledby="MAIN_TABS[1].id"
        ref="resultsPanel"
        v-model:view="resultsView"
        :coordinator="coordinator"
        :results="results"
        :runs="state.runs"
        :counts="counts"
        @search="showSearch"
        @use-in-blastn="useCountsInBlastn"
        @use-taxid-in-blastn="useTaxidInBlastn"
        @mask-with-counts="maskWithCounts"
        @format-archive="formatArchiveFromRun"
        @use-as-database="useAsDatabase"
        @use-as-repeats-database="useAsRepeatsDatabase"
        @use-as-mask-data="useAsMaskData"
      />
      <WindowMaskerPanel
        v-show="tab === 'windowmasker'"
        :id="MAIN_TABS[2].panel"
        role="tabpanel"
        :aria-labelledby="MAIN_TABS[2].id"
        :draft="toolDraft"
        :coordinator="coordinator"
        :counts="counts"
        @submitted="openResults"
      />
      <MakeblastdbPanel
        v-show="tab === 'makeblastdb'"
        :id="MAIN_TABS[3].panel"
        role="tabpanel"
        :aria-labelledby="MAIN_TABS[3].id"
        :draft="makeblastdbDraft"
        :coordinator="coordinator"
        @submitted="openResults"
      />
      <FormatArchivePanel
        v-show="tab === 'blast_formatter'"
        :id="MAIN_TABS[4].panel"
        role="tabpanel"
        :aria-labelledby="MAIN_TABS[4].id"
        :draft="blastFormatterDraft"
        :coordinator="coordinator"
        @submitted="openResults"
      />
    </section>
    <!-- Under the page: the storage's warnings, when there are any (design §3.4). -->
    <aside class="notes">
      <StorageStatus :storage="state.storage" />
    </aside>
  </main>
  <LiveAnnouncer />
</template>
