<script setup lang="ts">
// What the header's short offline status cannot say (OfflineStatus.vue): why offline use could
// not be prepared, with "Try again", and that this browser cannot start offline because it did
// not keep the choice of release. Shown under the header only in those states, which are rare
// and come at most once per preparation; texts from the worker are shown as text.
import type { ReleaseManager, ReleaseState } from '../application/release';

defineProps<{ release: ReleaseManager; state: ReleaseState }>();
</script>

<template>
  <section v-if="state.offline.kind === 'failed'" class="notice offline-notice" aria-label="Offline use" data-testid="offline-notice">
    <p data-testid="offline-failed">
      Offline use could not be prepared ({{ state.offline.reason }}). LOSAT Web works online as before.
    </p>
    <button type="button" data-testid="offline-retry" @click="release.retryOffline()">Try again</button>
  </section>
  <section
    v-else-if="state.offline.kind === 'stored' && state.offline.pin === undefined"
    class="notice offline-notice"
    aria-label="Offline use"
    data-testid="offline-notice"
  >
    <p data-testid="offline-not-pinned">
      LOSAT Web cannot start offline in this browser: an offline start needs the browser to keep the choice of release (site
      storage), and it did not.
    </p>
  </section>
</template>
