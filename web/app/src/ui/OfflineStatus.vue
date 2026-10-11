<script setup lang="ts">
// The release this tab runs and its offline use, in the header (plan §5.9): "Offline ready"
// only once the site's service worker has stored and checked every file of the release, and
// while the pinned release (which an offline start opens) is this tab's.
// The status keeps one size in every state, at every width: a fixed width (at most the
// header's) and one line (two on a phone: the release, then the state), so that a change of
// state never moves the page under the user's pointer. A longer text ends in an ellipsis; the
// whole text is the status's title.
import { computed } from 'vue';
import type { OfflineState, ReleaseManager, ReleaseState } from '../application/release';
import { formatBytes, formatCount } from './format';

const props = defineProps<{ release: ReleaseManager; state: ReleaseState }>();
const offline = computed(() => props.state.offline);

const text = computed(() => describe(offline.value));

/** The line's text, and the whole text for the title (the files of a preparation as well). */
function describe(state: OfflineState): { line: string; whole: string } {
  switch (state.kind) {
    case 'ready':
      return same('Offline ready');
    case 'stored':
      return same(
        state.pin === undefined
          ? 'Not offline ready: an offline start needs this browser to keep the choice of release, and it did not.'
          : `Stored for offline use; this browser now starts release ${state.pin}.`,
      );
    case 'preparing': {
      const p = state.progress;
      if (p === undefined) return same('Preparing offline use…');
      return {
        line: `Preparing offline use: ${formatBytes(p.bytes)} of ${formatBytes(p.bytesTotal)}…`,
        whole: `Preparing offline use: ${formatCount(p.files)} of ${formatCount(p.filesTotal)} files (${formatBytes(p.bytes)} of ${formatBytes(p.bytesTotal)})…`,
      };
    }
    case 'failed':
      return same(`Offline use could not be prepared (${state.reason}).`);
    case 'unsupported':
      return same(state.reason === undefined ? 'Offline use is not available in this browser.' : `Offline use is not available (${state.reason}).`);
    case 'waiting':
      return same('Offline use is prepared once the page has loaded.');
  }
}

function same(line: string): { line: string; whole: string } {
  return { line, whole: line };
}
</script>

<template>
  <span class="offline-status" data-testid="offline-status" :data-state="offline.kind" :title="`Release ${state.release}. ${text.whole}`">
    <span class="offline-release">Release <span class="release-id" data-testid="release-id">{{ state.release }}</span></span>
    <span class="offline-state">
      <span
        class="offline-message"
        :class="{ ok: offline.kind === 'ready', error: offline.kind === 'failed' }"
        :data-testid="`offline-${offline.kind}`"
        >{{ text.line }}</span
      >
      <button v-if="offline.kind === 'failed'" type="button" class="link" data-testid="offline-retry" @click="release.retryOffline()">
        Try again
      </button>
    </span>
  </span>
</template>
