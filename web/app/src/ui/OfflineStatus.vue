<script setup lang="ts">
// The release this tab runs and its offline use, in the header (plan §5.9): "Offline ready"
// only once the site's service worker has stored and checked every file of the release, and
// while the pinned release (which an offline start opens) is this tab's.
// The status is one line of one fixed size in every state, at every width, so that a change of
// state never moves the page under the user's pointer: a short form of the state, which is
// never cut, and the release ID, which ends in an ellipsis where the box is too narrow. The
// whole text is the status's title (the ID is whole in the page's text as well); a failure's
// reason and "Try again" are in the notice under the header (OfflineNotice.vue).
import { computed } from 'vue';
import type { OfflineState, ReleaseManager, ReleaseState } from '../application/release';
import { formatBytesUnbroken, formatCount } from './format';

const props = defineProps<{ release: ReleaseManager; state: ReleaseState }>();
const offline = computed(() => props.state.offline);
const text = computed(() => describe(offline.value));

/** The short form on the line, and the whole text for the title. */
function describe(state: OfflineState): { line: string; whole: string } {
  switch (state.kind) {
    case 'ready':
      return { line: 'Offline ready', whole: 'Offline ready: this release starts without the network.' };
    case 'stored':
      return {
        line: 'Not offline ready',
        whole:
          state.pin === undefined
            ? 'Stored, but an offline start needs this browser to keep the choice of release, and it did not.'
            : `Stored; this browser now starts release ${state.pin}, chosen in another tab.`,
      };
    case 'preparing': {
      const p = state.progress;
      if (p === undefined) return { line: 'Preparing offline use…', whole: 'Preparing offline use…' };
      const bytes = `${formatBytesUnbroken(p.bytes)} of ${formatBytesUnbroken(p.bytesTotal)}`;
      return {
        line: `Offline: ${bytes}…`,
        whole: `Preparing offline use: ${formatCount(p.files)} of ${formatCount(p.filesTotal)} files (${bytes})…`,
      };
    }
    case 'failed':
      return { line: 'Offline use failed', whole: `Offline use could not be prepared (${state.reason}).` };
    case 'unsupported':
      return {
        line: 'No offline use',
        whole: state.reason === undefined ? 'Offline use is not available in this browser.' : `Offline use is not available (${state.reason}).`,
      };
    case 'waiting':
      return { line: 'Offline: waiting', whole: 'Offline use is prepared once the page has loaded.' };
  }
}
</script>

<template>
  <span class="offline-status" data-testid="offline-status" :data-state="offline.kind" :title="`Release ${state.release}. ${text.whole}`">
    <span class="offline-release">Release <span class="release-id" data-testid="release-id">{{ state.release }}</span></span>
    <span class="offline-message" :class="{ ok: offline.kind === 'ready', error: offline.kind === 'failed' }" data-testid="offline-message">{{
      text.line
    }}</span>
  </span>
</template>
