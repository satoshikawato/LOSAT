<script setup lang="ts">
// What a switch of release does to this tab's work (plan §5.9: download and check -> what the
// switch does to running and unsaved work -> switch). Counted from the app's state as it is
// now, so the list shown next to "Switch ... and reload" includes the runs that finished and
// the candidates and notes added while the release was downloaded.
import { computed } from 'vue';
import type { RunView } from '../application/coordinator';
import { switchImpact } from '../application/release';
import { formatCount, formatCounted } from './format';

const props = defineProps<{
  latest: string;
  runs: readonly RunView[];
  /** The tray's candidates (application/candidates.ts). */
  candidates: ReadonlyArray<{ readonly note: string }>;
}>();
const impact = computed(() => switchImpact(props.runs, props.candidates));
</script>

<template>
  <p>Switching to release <strong>{{ latest }}</strong> reloads this tab:</p>
  <ul data-testid="update-impact">
    <li v-if="impact.running + impact.queued > 0" class="warning" data-testid="update-impact-runs">
      The runs in progress are stopped and leave no results ({{ formatCount(impact.running) }} running,
      {{ formatCount(impact.queued) }} queued).
    </li>
    <li v-if="impact.completed + impact.candidates > 0" data-testid="update-impact-saved">
      {{ formatCounted(impact.completed, 'completed run') }}, {{ formatCounted(impact.candidates, 'candidate') }} and
      {{ formatCounted(impact.notes, 'note') }} are not kept by the reload, and LOSAT Web does not restore them by itself. To keep
      them, use <strong>Save session</strong> in the Session panel first, and open the saved file after the switch.
    </li>
    <li v-else data-testid="update-impact-saved">No completed run or candidate is lost.</li>
    <li>The search form starts empty: its inputs and settings are not kept.</li>
  </ul>
</template>
