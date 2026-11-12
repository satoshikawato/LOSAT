<script setup lang="ts">
import { nextTick, onMounted, onUnmounted, ref, useId } from 'vue';
import type { Attention } from '../application/attention';
import type { Coordinator } from '../application/coordinator';
import type { SearchDraft } from '../application/draft';
import type { CountsLoader } from '../application/counts';
import type { ExampleLoader } from '../application/examples';
import type { ResultsBrowser } from '../application/results';
import { useStore } from './useStore';
import LiveAnnouncer from './LiveAnnouncer.vue';
import LosatLogo from './LosatLogo.vue';
import ResultsPanel, { type ResultsView } from './ResultsPanel.vue';
import ResumeNotice from './ResumeNotice.vue';
import SearchPanel from './SearchPanel.vue';
import { useRunAnnouncements } from './run-announcements';
import StorageStatus from './StorageStatus.vue';
import TabList from './TabList.vue';

const props = defineProps<{
  coordinator: Coordinator;
  draft: SearchDraft;
  examples: ExampleLoader;
  counts: CountsLoader;
  results: ResultsBrowser;
  attention: Attention;
  usesFakeEngine: boolean;
}>();
const state = useStore(props.coordinator.state);
const attentionState = useStore(props.attention.state);
const tab = ref<'search' | 'results'>('search');
const resultsPanel = ref<InstanceType<typeof ResultsPanel>>();
const mainTabs = ref<{ focus: () => void }>();
/** Search and Results: tabs, each with its panel (WA2-AX-12). */
const ids = useId();
const MAIN_TABS = [
  { value: 'search', label: 'Search', testid: 'tab-search', id: `${ids}-search-tab`, panel: `${ids}-search` },
  { value: 'results', label: 'Results', testid: 'tab-results', id: `${ids}-results-tab`, panel: `${ids}-results` },
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

// A file dropped outside an input's drop zone would make the browser open it in place of
// the application (and end the searches of this tab).
const keep = (event: DragEvent) => {
  if (event.dataTransfer?.types.includes('Files')) event.preventDefault();
};
onMounted(() => {
  window.addEventListener('dragover', keep);
  window.addEventListener('drop', keep);
});
onUnmounted(() => {
  window.removeEventListener('dragover', keep);
  window.removeEventListener('drop', keep);
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
  <ResumeNotice v-if="attentionState.resume" :report="attentionState.resume" @dismiss="dismissResume" />
  <!-- The tabs stay in the page's navigation landmark, which they were before they were tabs (WA2-AX-12). -->
  <nav aria-label="Main">
    <TabList ref="mainTabs" v-model="tab" class="main-tabs" label="Main" :tabs="MAIN_TABS" />
  </nav>
  <main class="layout">
    <section class="primary">
      <!-- The search form stays mounted, so the next search keeps its edits while results are viewed. -->
      <SearchPanel
        v-show="tab === 'search'"
        :id="MAIN_TABS[0].panel"
        role="tabpanel"
        :aria-labelledby="MAIN_TABS[0].id"
        :draft="draft"
        :examples="examples"
        :counts="counts"
        :attention="attention"
        :attention-state="attentionState"
        @submitted="openResults"
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
        @search="showSearch"
      />
    </section>
    <!-- Under the page: the storage's warnings, when there are any (design §3.4). -->
    <aside class="notes">
      <StorageStatus :storage="state.storage" />
    </aside>
  </main>
  <LiveAnnouncer />
</template>
