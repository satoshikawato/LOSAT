<script setup lang="ts">
import { computed } from 'vue';
import { MEMORY_RESULTS_CAPACITY, type StorageInfo } from '../ports/data';

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
    <!-- One sentence in LOSAT's words; the browser's error, for developers, is in data-error (line WQA's sweep, D5). It says
         what differs from a tab that stores files: memory, and its budget. A reload loses results there too, since LOSAT
         Web keeps nothing for the next visit (REQ-20), so that is no consequence of the reason (review of D5, M1). -->
    <p v-if="storage.fallbackReason" class="warning" data-testid="storage-memory-note" :data-error="storage.fallbackError">
      Results are kept in this tab's memory (at most {{ MEMORY_RESULTS_CAPACITY }}), because {{ storage.fallbackReason }}.
    </p>
    <p v-if="storage.cleanup.state === 'unavailable'" class="note" data-testid="storage-cleanup-note" :data-error="storage.cleanup.error">
      Temporary data of closed tabs is not removed: {{ storage.cleanup.reason }}.
    </p>
  </section>
</template>
