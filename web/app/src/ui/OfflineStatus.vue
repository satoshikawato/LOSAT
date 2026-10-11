<script setup lang="ts">
// The release this tab runs and its offline use, in the header (plan §5.9): "Offline ready"
// only once the site's service worker has stored and checked every file of the release, and
// while the pinned release (which an offline start opens) is this tab's.
import type { ReleaseManager, ReleaseState } from '../application/release';
import { formatBytes, formatCount } from './format';

defineProps<{ release: ReleaseManager; state: ReleaseState }>();
</script>

<template>
  <span class="offline-status" data-testid="offline-status" :data-state="state.offline.kind">
    <span>Release <span class="release-id" data-testid="release-id">{{ state.release }}</span></span>
    <span v-if="state.offline.kind === 'ready'" class="ok" data-testid="offline-ready">Offline ready</span>
    <span v-else-if="state.offline.kind === 'stored'" data-testid="offline-stored">
      {{
        state.offline.pin === undefined
          ? 'Not offline ready: an offline start needs this browser to keep the choice of release, and it did not.'
          : `Stored for offline use; this browser now starts release ${state.offline.pin}.`
      }}
    </span>
    <span v-else-if="state.offline.kind === 'preparing'" data-testid="offline-preparing">
      Preparing offline use<template v-if="state.offline.progress">:
        {{ formatCount(state.offline.progress.files) }} of {{ formatCount(state.offline.progress.filesTotal) }} files
        ({{ formatBytes(state.offline.progress.bytes) }} of {{ formatBytes(state.offline.progress.bytesTotal) }})</template>…
    </span>
    <span v-else-if="state.offline.kind === 'failed'" data-testid="offline-failed">
      <span class="error">Offline use could not be prepared ({{ state.offline.reason }}).</span>
      <button type="button" class="link" data-testid="offline-retry" @click="release.retryOffline()">Try again</button>
    </span>
    <span v-else-if="state.offline.kind === 'unsupported'" data-testid="offline-unsupported">
      {{ state.offline.reason === undefined ? 'Offline use is not available in this browser.' : `Offline use is not available (${state.offline.reason}).` }}
    </span>
    <span v-else data-testid="offline-waiting">Offline use is prepared once the page has loaded.</span>
  </span>
</template>
