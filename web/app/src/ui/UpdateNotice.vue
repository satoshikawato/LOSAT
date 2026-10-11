<script setup lang="ts">
// The explicit update (plan §5.9, design §14): the site serves another release -> the user
// reviews what a reload does to this tab's work -> the release is downloaded and checked ->
// the screen says again what the reload does, counted now -> the user switches, and the page
// reloads. Nothing changes before the user chooses, and a failure leaves this tab on its
// release; the notice then takes the error style (FAILURES). Release IDs come from the site and
// are shown as text.
import { computed } from 'vue';
import type { RunView } from '../application/coordinator';
import { switchImpact, type ReleaseManager, type ReleaseState } from '../application/release';
import { formatBytesUnbroken, formatCount, formatCounted } from './format';
import UpdateImpact from './UpdateImpact.vue';

const props = defineProps<{
  release: ReleaseManager;
  state: ReleaseState;
  runs: readonly RunView[];
  /** The tray's candidates (application/candidates.ts). */
  candidates: ReadonlyArray<{ readonly note: string }>;
}>();
const update = computed(() => props.state.update);
const shown = computed(() => update.value.kind !== 'none' && update.value.kind !== 'dismissed');
/** The states in which the update failed or was not made: the notice is an error, not information (screen review 1 L5). */
const FAILURES: ReadonlySet<string> = new Set(['failed', 'removed', 'not-kept']);
const latest = computed(() => ('latest' in update.value ? update.value.latest : ''));
/** What a reload of this tab drops (the leave guard covers only runs in progress). */
const impact = computed(() => switchImpact(props.runs, props.candidates));
/** Only a browser without service workers (or where they are blocked) switches without a download and check first. */
const direct = computed(() => {
  // Read again when the offline state changes: the manager learns of a block at its first registration.
  void props.state.offline;
  return props.release.switchesDirectly();
});
</script>

<template>
  <p v-if="state.elsewhere" class="notice" role="status" data-testid="update-other-tab">
    LOSAT Web was switched to release {{ state.elsewhere }} in another tab. This tab keeps release {{ state.release }} until it is
    reloaded.
    <template v-if="impact.completed + impact.candidates > 0">
      <br />
      <span data-testid="update-other-tab-loss">
        Reloading does not keep this tab's {{ formatCounted(impact.completed, 'completed run') }} and
        {{ formatCounted(impact.candidates, 'candidate') }}: use <strong>Save session</strong> in the Session panel first.
      </span>
    </template>
  </p>
  <section
    v-if="shown"
    class="update-notice"
    :class="{ failure: FAILURES.has(update.kind) }"
    aria-label="Update"
    data-testid="update-notice"
    :data-state="update.kind"
  >
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
      <UpdateImpact :latest="latest" :runs="runs" :candidates="candidates" />
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
          ({{ formatBytesUnbroken(update.progress.bytes) }} of {{ formatBytesUnbroken(update.progress.bytesTotal) }})</template>…
      </p>
      <progress
        :value="update.progress?.bytes"
        :max="update.progress?.bytesTotal || undefined"
        :aria-label="`Download of release ${latest}`"
      ></progress>
    </template>

    <template v-else-if="update.kind === 'downloaded'">
      <p data-testid="update-downloaded">Release {{ latest }} is downloaded and every file matched its checksum.</p>
      <UpdateImpact :latest="latest" :runs="runs" :candidates="candidates" />
      <div class="update-actions">
        <button type="button" data-testid="update-switch" @click="release.switchNow()">Switch to release {{ latest }} and reload</button>
        <button type="button" data-testid="update-cancel" @click="release.cancel()">Cancel</button>
      </div>
    </template>

    <template v-else-if="update.kind === 'failed'">
      <p class="error" data-testid="update-failed">
        Release {{ latest }} could not be downloaded and checked: {{ update.reason }}. LOSAT Web stays on release {{ state.release }};
        nothing was changed.
      </p>
      <div class="update-actions">
        <button type="button" data-testid="update-retry" @click="release.download()">Try again</button>
        <button type="button" data-testid="update-cancel" @click="release.cancel()">Cancel</button>
      </div>
    </template>

    <template v-else-if="update.kind === 'removed'">
      <p class="error" data-testid="update-removed">
        Release {{ latest }} is no longer stored in this browser, so LOSAT Web did not switch to it and stays on release
        {{ state.release }}. Download and check it again to switch.
      </p>
      <div class="update-actions">
        <button type="button" data-testid="update-retry" @click="release.download()">Download and check again</button>
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
