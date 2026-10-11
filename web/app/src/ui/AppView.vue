<script setup lang="ts">
import { nextTick, onMounted, onUnmounted, ref } from 'vue';
import type { Analytics } from '../application/analytics';
import type { Attention } from '../application/attention';
import type { CandidateTray } from '../application/candidates';
import type { Coordinator } from '../application/coordinator';
import type { SearchDraft } from '../application/draft';
import type { HspId, ResultsBrowser } from '../application/results';
import type { ResultExporter } from '../application/result-export';
import type { Session } from '../application/session';
import type { RunFiles } from '../application/run-files';
import type { Examples } from '../application/examples';
import type { GettingStarted as GettingStartedState } from '../application/onboarding';
import type { SelfTest } from '../application/self-test';
import type { ReleaseManager } from '../application/release';
import { formatCount } from './format';
import { useStore } from './useStore';
import AnalyticsConsent from './AnalyticsConsent.vue';
import AnalyticsSetting from './AnalyticsSetting.vue';
import AttentionPanel from './AttentionPanel.vue';
import CandidatesPanel from './CandidatesPanel.vue';
import GettingStarted from './GettingStarted.vue';
import QueuePanel from './QueuePanel.vue';
import OfflineNotice from './OfflineNotice.vue';
import OfflineStatus from './OfflineStatus.vue';
import ResultsPanel, { type ResultsView } from './ResultsPanel.vue';
import ResumeNotice from './ResumeNotice.vue';
import SearchPanel from './SearchPanel.vue';
import SessionPanel from './SessionPanel.vue';
import StorageStatus from './StorageStatus.vue';
import UpdateNotice from './UpdateNotice.vue';

const props = defineProps<{
  coordinator: Coordinator;
  draft: SearchDraft;
  results: ResultsBrowser;
  candidates: CandidateTray;
  exporter: ResultExporter;
  attention: Attention;
  usesFakeEngine: boolean;
  session: Session;
  runFiles: RunFiles;
  examples: Examples;
  selfTest: SelfTest;
  gettingStarted: GettingStartedState;
  release: ReleaseManager;
  analytics: Analytics;
}>();
const state = useStore(props.coordinator.state);
const attentionState = useStore(props.attention.state);
const trayState = useStore(props.candidates.state);
const gettingStartedState = useStore(props.gettingStarted.state);
const releaseState = useStore(props.release.state);
const tab = ref<'search' | 'results' | 'candidates'>('search');
const resultsPanel = ref<InstanceType<typeof ResultsPanel>>();
const searchPanel = ref<InstanceType<typeof SearchPanel>>();
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
  resultsView.value = 'hits';
  tab.value = 'results';
  const shown = await revealing;
  await nextTick();
  if (shown) resultsPanel.value?.showHsp(id);
}

/** "Edit Search" put a run's settings in the search form: the form is shown at its settings line, which says so. */
async function showEditedSearch(): Promise<void> {
  tab.value = 'search';
  await nextTick();
  searchPanel.value?.showSettings();
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
    <!-- At the right, together: a fixed-size offline status never moves the link or the page. -->
    <div class="header-tools">
      <button
        type="button"
        class="link header-link"
        :aria-expanded="gettingStartedState.open"
        aria-controls="getting-started"
        data-testid="getting-started-toggle"
        @click="gettingStarted.toggle()"
      >
        Getting started
      </button>
      <OfflineStatus :release="release" :state="releaseState" />
    </div>
  </header>
  <p v-if="usesFakeEngine" class="banner" data-testid="fake-engine-banner">
    Development build: a fake engine is in use. Outputs are not search results.
  </p>
  <UpdateNotice :release="release" :state="releaseState" :runs="state.runs" :candidates="trayState.candidates" />
  <OfflineNotice :release="release" :state="releaseState" />
  <!-- Kept mounted, so that the header's aria-controls names it and the self-test's section keeps its state. -->
  <GettingStarted v-show="gettingStartedState.open" :getting-started="gettingStarted" :self-test="selfTest" />
  <ResumeNotice v-if="attentionState.resume" :report="attentionState.resume" @dismiss="attention.dismissResume()" />
  <AnalyticsConsent :analytics="analytics" />
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
      <SearchPanel v-show="tab === 'search'" ref="searchPanel" :draft="draft" :runs="state.runs" :run-files="runFiles" :examples="examples" />
      <ResultsPanel
        v-if="tab === 'results'"
        ref="resultsPanel"
        v-model:view="resultsView"
        :coordinator="coordinator"
        :results="results"
        :candidates="candidates"
        :exporter="exporter"
        :runs="state.runs"
        :session="session"
        :run-files="runFiles"
        @edit-search="showEditedSearch"
      />
      <!-- Kept mounted, as the search form: the extraction's choices and the order shown stay while other tabs are viewed. -->
      <CandidatesPanel v-show="tab === 'candidates'" :candidates="candidates" :state="trayState" :runs="state.runs" @reveal="showCandidate" />
    </section>
    <aside class="secondary">
      <QueuePanel :coordinator="coordinator" :runs="state.runs" @open-results="openResults" />
      <SessionPanel :session="session" :runs="state.runs" :candidates="trayState.candidates" />
      <AttentionPanel :attention="attention" :state="attentionState" />
      <StorageStatus :storage="state.storage" />
    </aside>
  </main>
  <AnalyticsSetting :analytics="analytics" />
</template>
