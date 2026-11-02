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
    <p v-if="storage.fallbackReason" class="warning" data-testid="storage-memory-note">
      Results are kept in memory only and are lost when this tab is closed or reloaded ({{ storage.fallbackReason }}).
    </p>
    <p v-if="storage.cleanup.state === 'unavailable'" class="note">
      Temporary data of closed tabs is not removed: {{ storage.cleanup.reason }}.
    </p>
  </section>
</template>
