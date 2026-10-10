<script setup lang="ts">
import { nextTick, onMounted, onUnmounted, ref } from 'vue';
import type { Attention } from '../application/attention';
import type { CandidateTray } from '../application/candidates';
import type { Coordinator } from '../application/coordinator';
import type { SearchDraft } from '../application/draft';
import type { HspId, ResultsBrowser } from '../application/results';
import type { Session } from '../application/session';
import { formatCount } from './format';
import { useStore } from './useStore';
import AttentionPanel from './AttentionPanel.vue';
import CandidatesPanel from './CandidatesPanel.vue';
import QueuePanel from './QueuePanel.vue';
import ResultsPanel, { type ResultsView } from './ResultsPanel.vue';
import ResumeNotice from './ResumeNotice.vue';
import SearchPanel from './SearchPanel.vue';
import SessionPanel from './SessionPanel.vue';
import StorageStatus from './StorageStatus.vue';

const props = defineProps<{
  coordinator: Coordinator;
  draft: SearchDraft;
  results: ResultsBrowser;
  candidates: CandidateTray;
  attention: Attention;
  usesFakeEngine: boolean;
  session: Session;
}>();
const state = useStore(props.coordinator.state);
const attentionState = useStore(props.attention.state);
const trayState = useStore(props.candidates.state);
const tab = ref<'search' | 'results' | 'candidates'>('search');
const resultsPanel = ref<InstanceType<typeof ResultsPanel>>();
/** The results tab's view, kept while the search tab is shown and when another run opens. */
const resultsView = ref<ResultsView>('hits');

/**
 * Opens a run's results from the queue (S12's screen review L7), and brings the results
 * panel into view with the focus on its heading: on a narrow screen the queue is below the
 * results (S13 screen review M3).
 */
async function openResults(runId: string): Promise<void> {
  tab.value = 'results';
  void props.results.open(runId);
  await nextTick();
  resultsPanel.value?.showHeading();
}

/**
 * "Show in results" of a candidate (REQ-14): its run, query, subject and HSP are selected (the
 * view filters that hide it are cleared), and the Alignments show its Range with the focus on it.
 */
async function showCandidate(id: HspId): Promise<void> {
  const revealing = props.results.reveal(id);
  resultsView.value = 'alignment';
  tab.value = 'results';
  const shown = await revealing;
  await nextTick();
  if (shown) resultsPanel.value?.showHsp(id);
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
    <span class="local-note">Local processing: your sequences stay in this browser.</span>
  </header>
  <p v-if="usesFakeEngine" class="banner" data-testid="fake-engine-banner">
    Development build: a fake engine is in use. Outputs are not search results.
  </p>
  <ResumeNotice v-if="attentionState.resume" :report="attentionState.resume" @dismiss="attention.dismissResume()" />
  <nav class="tabs main-tabs" aria-label="Main">
    <button :aria-pressed="tab === 'search'" data-testid="tab-search" @click="tab = 'search'">Search</button>
    <button :aria-pressed="tab === 'results'" data-testid="tab-results" @click="tab = 'results'">Results</button>
    <button :aria-pressed="tab === 'candidates'" data-testid="tab-candidates" @click="tab = 'candidates'">
      Candidates <span class="tab-count" data-testid="tab-candidates-count">{{ formatCount(trayState.candidates.length) }}</span>
    </button>
  </nav>
  <main class="layout">
    <section class="primary">
      <!-- The search form stays mounted, so the next job keeps its edits while results are viewed. -->
      <SearchPanel v-show="tab === 'search'" :draft="draft" :runs="state.runs" />
      <ResultsPanel
        v-if="tab === 'results'"
        ref="resultsPanel"
        v-model:view="resultsView"
        :coordinator="coordinator"
        :results="results"
        :candidates="candidates"
        :runs="state.runs"
        :session="session"
      />
      <!-- Kept mounted, as the search form: the extraction's choices and the order shown stay while other tabs are viewed. -->
      <CandidatesPanel v-show="tab === 'candidates'" :candidates="candidates" :state="trayState" :runs="state.runs" @reveal="showCandidate" />
    </section>
    <aside class="secondary">
      <QueuePanel :coordinator="coordinator" :runs="state.runs" @open-results="openResults" />
      <SessionPanel :session="session" :runs="state.runs" />
      <AttentionPanel :attention="attention" :state="attentionState" />
      <StorageStatus :storage="state.storage" />
    </aside>
  </main>
</template>
