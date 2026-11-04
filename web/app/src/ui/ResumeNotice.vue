<script setup lang="ts">
import { computed } from 'vue';
import type { ResumeReport } from '../application/attention';
import { formatDuration, formatRuns } from './format';
import LosatLoader from './LosatLoader.vue';

const props = defineProps<{ report: ResumeReport }>();
defineEmits<{ dismiss: [] }>();

/** What became of the stored results that the Data worker could not read back (WA2-SC-05); empty if none. */
const lost = computed(() => {
  const { removed, unreadable } = props.report;
  const sentences = [
    ...(removed.length > 0 ? [`The browser removed the stored results of ${formatRuns(removed)}.`] : []),
    ...(unreadable.length > 0 ? [`The stored results of ${formatRuns(unreadable)} could not be read.`] : []),
  ];
  if (sentences.length === 0) return '';
  return `${sentences.join(' ')} Run ${removed.length + unreadable.length === 1 ? 'it' : 'them'} again.`;
});
</script>

<template>
  <div class="resume-notice" role="status" data-testid="resume-notice">
    <p>
      This tab was in the background for {{ formatDuration(report.hiddenMs) }} while searches were in progress.
    </p>
    <ul>
      <li v-for="run in report.runs" :key="run.number" :data-testid="`resume-run-${run.number}`">
        Run {{ run.number }} was {{ run.before }} when the tab was hidden{{ run.now === run.before ? ' and still is' : `; it is ${run.now} now` }}.
        <template v-if="run.now === 'preparing' || run.now === 'running'">
          If it does not end, the browser may have stopped it in the background: cancel it and run it again.
        </template>
      </li>
    </ul>
    <p v-if="report.dataWorker === 'checking'" class="muted" data-testid="resume-data"><LosatLoader />Checking the stored results…</p>
    <p v-else-if="report.dataWorker === 'responding' && lost !== ''" class="error" data-testid="resume-data">{{ lost }}</p>
    <p v-else-if="report.dataWorker === 'responding'" class="muted" data-testid="resume-data">The stored results are available.</p>
    <p v-else class="error" data-testid="resume-data">
      The part of LOSAT that keeps the results did not answer. If results cannot be opened, reload the page; the
      results of this session are then lost.
    </p>
    <button type="button" data-testid="resume-dismiss" @click="$emit('dismiss')">Dismiss</button>
  </div>
</template>
