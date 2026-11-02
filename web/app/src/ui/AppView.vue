<script setup lang="ts">
import { nextTick, onMounted, onUnmounted, ref } from 'vue';
import type { Attention } from '../application/attention';
import type { Coordinator } from '../application/coordinator';
import type { SearchDraft } from '../application/draft';
import type { ResultsBrowser } from '../application/results';
import { useStore } from './useStore';
import ResultsPanel, { type ResultsView } from './ResultsPanel.vue';
import ResumeNotice from './ResumeNotice.vue';
import SearchPanel from './SearchPanel.vue';
import StorageStatus from './StorageStatus.vue';

const props = defineProps<{
  coordinator: Coordinator;
  draft: SearchDraft;
  results: ResultsBrowser;
  attention: Attention;
  usesFakeEngine: boolean;
}>();
const state = useStore(props.coordinator.state);
const attentionState = useStore(props.attention.state);
const tab = ref<'search' | 'results'>('search');
const resultsPanel = ref<InstanceType<typeof ResultsPanel>>();
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
    <h1>LOSAT Web</h1>
    <span class="local-note">
      <svg aria-hidden="true" width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><rect x="4" y="11" width="16" height="10" rx="2" /><path d="M8 11V7a4 4 0 0 1 8 0v4" /></svg>
      Nothing leaves this browser
    </span>
  </header>
  <p v-if="usesFakeEngine" class="banner" data-testid="fake-engine-banner">
    Development build: a fake engine is in use. Outputs are not search results.
  </p>
  <ResumeNotice v-if="attentionState.resume" :report="attentionState.resume" @dismiss="attention.dismissResume()" />
  <nav class="tabs main-tabs" aria-label="Main">
    <button :aria-pressed="tab === 'search'" data-testid="tab-search" @click="tab = 'search'">Search</button>
    <button :aria-pressed="tab === 'results'" data-testid="tab-results" @click="tab = 'results'">Results</button>
  </nav>
  <main class="layout">
    <section class="primary">
      <!-- The search form stays mounted, so the next search keeps its edits while results are viewed. -->
      <SearchPanel v-show="tab === 'search'" :draft="draft" :attention="attention" :attention-state="attentionState" @submitted="openResults" />
      <ResultsPanel
        v-if="tab === 'results'"
        ref="resultsPanel"
        v-model:view="resultsView"
        :coordinator="coordinator"
        :results="results"
        :runs="state.runs"
        @search="tab = 'search'"
      />
    </section>
    <!-- Under the page: the storage's warnings, when there are any (design §3.4). -->
    <aside class="notes">
      <StorageStatus :storage="state.storage" />
    </aside>
  </main>
</template>
