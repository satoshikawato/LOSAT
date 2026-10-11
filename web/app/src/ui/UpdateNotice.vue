<script setup lang="ts">
// The explicit update (plan §5.9, design §14): the site serves another release -> the user
// reviews what a reload does to this tab's work -> the release is downloaded and checked ->
// the user switches, and the page reloads. Nothing changes before the user chooses, and a
// failure leaves this tab on its release. Release IDs come from the site and are shown as text.
import { computed } from 'vue';
import type { RunView } from '../application/coordinator';
import { switchImpact, type ReleaseManager, type ReleaseState } from '../application/release';
import { formatBytes, formatCount, formatCounted } from './format';

const props = defineProps<{
  release: ReleaseManager;
  state: ReleaseState;
  runs: readonly RunView[];
  /** The tray's candidates (application/candidates.ts). */
  candidates: ReadonlyArray<{ readonly note: string }>;
}>();
const update = computed(() => props.state.update);
const shown = computed(() => update.value.kind !== 'none' && update.value.kind !== 'dismissed');
const latest = computed(() => ('latest' in update.value ? update.value.latest : ''));
const impact = computed(() => switchImpact(props.runs, props.candidates));
const direct = computed(() => props.state.offline.kind === 'unsupported');
</script>

<template>
  <p v-if="state.elsewhere" class="notice" role="status" data-testid="update-other-tab">
    LOSAT Web was switched to release {{ state.elsewhere }} in another tab. This tab keeps release {{ state.release }} until it is
    reloaded.
  </p>
  <section v-if="shown" class="update-notice" aria-label="Update" data-testid="update-notice" :data-state="update.kind">
    <template v-if="update.kind === 'available'">
      <p>
        This site now serves LOSAT Web release <strong>{{ latest }}</strong>; this tab uses release {{ state.release }}. Nothing changes
        until you choose to switch.
      </p>
      <div class="update-actions">
        <button type="button" data-testid="update-review" @click="release.review()">Review update</button>
        <button type="button" data-testid="update-not-now" @click="release.notNow()">Not now</button>
      </div>
    </template>

    <template v-else-if="update.kind === 'reviewing'">
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
      <p v-if="direct" class="hint">This browser cannot store releases, so the new release is not checked before the switch.</p>
      <div class="update-actions">
        <button v-if="direct" type="button" data-testid="update-switch" @click="release.switchNow()">
          Switch to release {{ latest }} and reload
        </button>
        <button v-else type="button" data-testid="update-download" @click="release.download()">Download and check</button>
        <button type="button" data-testid="update-cancel" @click="release.cancel()">Cancel</button>
      </div>
    </template>

    <template v-else-if="update.kind === 'downloading'">
      <p data-testid="update-progress">
        Downloading and checking release {{ latest }}<template v-if="update.progress">:
          {{ formatCount(update.progress.files) }} of {{ formatCount(update.progress.filesTotal) }} files
          ({{ formatBytes(update.progress.bytes) }} of {{ formatBytes(update.progress.bytesTotal) }})</template>…
      </p>
      <progress
        :value="update.progress?.bytes"
        :max="update.progress?.bytesTotal || undefined"
        :aria-label="`Download of release ${latest}`"
      ></progress>
    </template>

    <template v-else-if="update.kind === 'downloaded'">
      <p data-testid="update-downloaded">Release {{ latest }} is downloaded and every file matched its checksum.</p>
      <div class="update-actions">
        <button type="button" data-testid="update-switch" @click="release.switchNow()">Switch to release {{ latest }} and reload</button>
        <button type="button" data-testid="update-cancel" @click="release.cancel()">Cancel</button>
      </div>
    </template>

    <template v-else-if="update.kind === 'failed'">
      <p class="error" data-testid="update-failed">
        The update could not be downloaded completely ({{ update.reason }}). LOSAT Web stays on release {{ state.release }}.
      </p>
      <div class="update-actions">
        <button type="button" data-testid="update-retry" @click="release.download()">Try again</button>
        <button type="button" data-testid="update-cancel" @click="release.cancel()">Cancel</button>
      </div>
    </template>

    <template v-else-if="update.kind === 'not-kept'">
      <p class="error" data-testid="update-not-kept">
        The switch cannot be kept: the browser did not store the choice of release ({{ update.reason }}). LOSAT Web stays on release
        {{ state.release }}.
      </p>
      <div class="update-actions">
        <button type="button" data-testid="update-cancel" @click="release.cancel()">Cancel</button>
      </div>
    </template>

    <template v-else-if="update.kind === 'switched'">
      <p data-testid="update-switched">
        LOSAT Web uses release {{ latest }} from the next reload. This tab keeps release {{ state.release }} until it is reloaded.
      </p>
      <div class="update-actions">
        <button type="button" data-testid="update-reload" @click="release.reloadNow()">Reload now</button>
      </div>
    </template>
  </section>
</template>
