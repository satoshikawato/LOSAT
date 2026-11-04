<script setup lang="ts">
import { computed } from 'vue';
import type { ResumeReport } from '../application/attention';
import { formatDuration, lostResultsText } from './format';
import LosatLoader from './LosatLoader.vue';

const props = defineProps<{ report: ResumeReport }>();
defineEmits<{ dismiss: [] }>();

/** The runs whose stored results the Data worker could not read back (WA2-SC-05); empty if none. */
const lost = computed(() => lostResultsText(props.report.removed, props.report.unreadable));
/** Until the Data worker has answered and the check of the stored results has ended (it has no time limit). */
const checking = computed(() => props.report.dataWorker === 'checking' || props.report.stored === 'checking');
/** The Data worker did not answer, or refused the check of the stored results. */
const silent = computed(() => props.report.dataWorker === 'not-responding' || props.report.stored === 'failed');
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
    <p v-if="silent" class="error" data-testid="resume-data">
      The part of LOSAT that keeps the results did not answer. If results cannot be opened, reload the page; the
      results of this session are then lost.
    </p>
    <p v-else-if="checking" class="muted" data-testid="resume-data"><LosatLoader />Checking the stored results…</p>
    <p v-else-if="lost !== ''" class="error" data-testid="resume-data">{{ lost }}</p>
    <p v-else class="muted" data-testid="resume-data">The stored results are available.</p>
    <button type="button" data-testid="resume-dismiss" @click="$emit('dismiss')">Dismiss</button>
  </div>
</template>
