<script setup lang="ts">
import { computed } from 'vue';
import type { StorageInfo } from '../ports/data';

const props = defineProps<{ storage: StorageInfo | undefined }>();

const removed = computed(() =>
  props.storage?.cleanup.state === 'done' ? props.storage.cleanup.removedSessions : undefined,
);
</script>

<template>
  <section
    v-if="storage"
    class="storage"
    data-testid="storage-status"
    :data-backend="storage.backend"
    :data-cleanup="storage.cleanup.state"
    :data-removed-sessions="removed"
  >
    <!-- Nothing in the normal case (stop-describing-ui rule 4): the warnings only. -->
    <!-- One sentence in LOSAT's words; the browser's error, for developers, is in data-error (line WQA's sweep, D5). -->
    <p v-if="storage.fallbackReason" class="warning" data-testid="storage-memory-note" :data-error="storage.fallbackError">
      Results are kept in this tab's memory only and are lost when it is closed or reloaded, because {{ storage.fallbackReason }}.
    </p>
    <p v-if="storage.cleanup.state === 'unavailable'" class="note" data-testid="storage-cleanup-note" :data-error="storage.cleanup.error">
      Temporary data of closed tabs is not removed: {{ storage.cleanup.reason }}.
    </p>
  </section>
</template>
