<script setup lang="ts">
import { computed } from 'vue';
import type { StorageInfo } from '../ports/data';
import { formatBytes } from './format';

const props = defineProps<{ storage: StorageInfo | undefined }>();

const usage = computed(() => {
  const storage = props.storage;
  if (storage === undefined) return '';
  const place = storage.backend === 'opfs' ? 'Browser storage (OPFS)' : 'Memory';
  const site =
    storage.estimate === undefined
      ? ''
      : `; this site uses about ${formatBytes(storage.estimate.usage)} of ` +
        `${formatBytes(storage.estimate.quota)} (browser estimate)`;
  return `${place}: ${formatBytes(storage.sessionBytes)} used by this tab${site}.`;
});
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
    <h2>Temporary storage</h2>
    <p data-testid="storage-usage">{{ usage }}</p>
    <p v-if="storage.fallbackReason" class="note">
      Results are kept in memory because OPFS cannot be used: {{ storage.fallbackReason }}.
    </p>
    <p v-if="storage.cleanup.state === 'unavailable'" class="note">
      Temporary data left by closed tabs is not removed automatically: {{ storage.cleanup.reason }}.
    </p>
    <p v-else-if="removed" class="note">
      Removed the temporary data left by {{ removed }} closed {{ removed === 1 ? 'tab' : 'tabs' }}.
    </p>
  </section>
</template>
