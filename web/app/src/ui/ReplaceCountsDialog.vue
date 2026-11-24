<script setup lang="ts">
// Confirmation dialog before replacing an existing saved WindowMasker counts copy for a taxonomy ID
// (WW3b DESIGN §3): native <dialog> opened with showModal(), Escape cancels, and closing returns
// focus to the Save button in ToolRunView.
import { computed, nextTick, onMounted, ref, useId, watch } from 'vue';
import { COUNTS_HOST } from '../domain/counts';
import { dateOf } from '../domain/download-sentences';
import type { CountsCopy } from '../ports/counts';
import { formatBytes } from './format';

const props = defineProps<{
  open: boolean;
  taxid: number;
  existing: CountsCopy | undefined;
  newFileName: string;
  runNumber: number;
}>();

const emit = defineEmits<{
  confirm: [];
  cancel: [];
}>();

const dialog = ref<HTMLDialogElement | null>(null);
const titleId = useId();
const descId = useId();

const question = computed(() => {
  const copy = props.existing;
  if (copy === undefined) return '';
  let source: string;
  let date: string;
  if (copy.made !== undefined) {
    source = `made in this browser by run ${copy.made.run}`;
    date = new Date(copy.made.madeAt).toISOString().slice(0, 10);
  } else {
    try {
      source = copy.url === '' ? COUNTS_HOST : new URL(copy.url).host;
    } catch {
      source = COUNTS_HOST;
    }
    date = dateOf(copy.lastModified) ?? new Date(copy.savedAt).toISOString().slice(0, 10);
  }
  const size = copy.complete ? formatBytes(copy.size) : `${formatBytes(copy.savedBytes)} of ${formatBytes(copy.size)}`;
  return `Counts for ${props.taxid} are saved (${source}, ${date}, ${size}). Replace them with ${props.newFileName} of run ${props.runNumber}?`;
});

function syncOpen(open: boolean): void {
  const el = dialog.value;
  if (el === null) return;
  if (open) {
    if (!el.open) el.showModal();
  } else if (el.open) {
    el.close();
  }
}

watch(
  () => props.open,
  (open) => {
    if (!open) {
      // Close synchronously so the page is no longer inert before the parent restores focus.
      syncOpen(false);
    } else {
      void nextTick(() => syncOpen(props.open));
    }
  },
  { flush: 'sync' },
);

onMounted(() => syncOpen(props.open));

function onCancel(event: Event): void {
  event.preventDefault();
  emit('cancel');
}
</script>

<template>
  <dialog
    ref="dialog"
    class="replace-dialog"
    :aria-labelledby="titleId"
    :aria-describedby="descId"
    data-testid="replace-counts-dialog"
    @cancel="onCancel"
  >
    <h3 :id="titleId">Replace saved counts?</h3>
    <p :id="descId" data-testid="replace-counts-message">{{ question }}</p>
    <div class="replace-dialog-actions">
      <button type="button" class="primary-action" data-testid="replace-counts-confirm" @click="emit('confirm')">
        Replace
      </button>
      <button type="button" data-testid="replace-counts-cancel" @click="emit('cancel')">
        Cancel
      </button>
    </div>
  </dialog>
</template>
