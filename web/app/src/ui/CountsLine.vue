<script setup lang="ts">
// The line under the taxonomy ID box of "Species-specific repeats filter for:" (line WW2,
// docs/web/fetch_design.md §10), in the Example line's sentences: what pressing "Download" sends
// ("Downloads wmasker.obinary for 9606 from ftp.ncbi.nlm.nih.gov"), the download's progress with
// Cancel, a download that stopped with Continue and Delete, or the file saved in this browser with
// Delete. Nothing is sent before "Download" is pressed. The rest waits behind "?"
// (stop-describing-ui). A download that ends or fails is told to a screen reader, and a button that
// leaves with the focus on it hands the focus to the taxonomy ID box (WA2-AX-04).
import { computed, onBeforeUnmount, onBeforeUpdate, onUpdated, ref, watch } from 'vue';
import type { CountsLoader } from '../application/counts';
import { COUNTS_FILE, COUNTS_HOST } from '../domain/counts';
import { dateOf } from '../domain/download-sentences';
import { announce } from './announcer';
import { formatBytes } from './format';
import HelpTip from './HelpTip.vue';
import LosatLoader from './LosatLoader.vue';
import { useStore } from './useStore';

const props = defineProps<{ counts: CountsLoader; taxid: number; focusId: string }>();
const state = useStore(props.counts.state);
const root = ref<HTMLElement | null>(null);

watch(
  () => props.taxid,
  (taxid) => props.counts.refresh(taxid),
  { immediate: true },
);

const copy = computed(() => state.value.saved[props.taxid] ?? undefined);
/** The download of this taxonomy ID, or the failure of it. */
const progress = computed(() => (state.value.load?.taxid === props.taxid ? state.value.load.progress : undefined));
const busy = computed(() => progress.value !== undefined && progress.value.state !== 'failed');
const failure = computed(() => (progress.value?.state === 'failed' ? progress.value.message : undefined));
const note = computed(() => (state.value.note?.taxid === props.taxid ? state.value.note.text : undefined));
/** A download that stopped after some bytes was saved. */
const partial = computed(() => (copy.value !== undefined && !copy.value.complete && copy.value.savedBytes > 0 ? copy.value : undefined));
const saved = computed(() => (copy.value?.complete === true ? copy.value : undefined));

const help = computed(
  () =>
    `Download asks ${COUNTS_HOST} only for the WindowMasker counts of this taxonomy ID and keeps the file in this browser; ` +
    'searches read it from there without a request.',
);

const sentence = computed(() => {
  const now = progress.value;
  switch (now?.state) {
    case 'reading':
      return 'Opening the saved counts…';
    case 'downloading':
      if (now.bytes === 0 && now.total === undefined) return `Downloading from ${COUNTS_HOST}…`;
      return `Downloading ${formatBytes(now.bytes)}${now.total === undefined ? '' : ` of ${formatBytes(now.total)}`} from ${COUNTS_HOST}…`;
  }
  if (saved.value) {
    const date = dateOf(saved.value.lastModified) ?? `saved ${new Date(saved.value.savedAt).toISOString().slice(0, 10)}`;
    // A browser that keeps nothing across sessions: the file is this tab's.
    return `${state.value.inTab[props.taxid] ? 'Kept in this tab only' : 'Saved in this browser'} (${date}, ${formatBytes(saved.value.size)})`;
  }
  if (partial.value) return `Downloaded ${formatBytes(partial.value.savedBytes)} of ${formatBytes(partial.value.size)}`;
  return `Downloads ${COUNTS_FILE} for ${props.taxid} from ${COUNTS_HOST}`;
});

function press(): void {
  if (!busy.value) props.counts.download(props.taxid);
}

// A download runs while its line is shown, as an Example set's load does while its mode is: another
// taxonomy ID, the row unchecked or a counts file chosen stops it (what came stays; Continue asks for
// the rest), so that no download goes on without its progress and Cancel.
watch(
  () => props.taxid,
  (_now, before) => {
    if (state.value.load?.taxid === before) props.counts.cancel();
  },
);
onBeforeUnmount(() => {
  if (state.value.load?.taxid === props.taxid) props.counts.cancel();
});

watch(failure, (message) => {
  if (message !== undefined) announce(message, 'alert');
});
// A download that ends with the whole file saved is told; the saved counts found at the start are not.
watch(busy, (now, before) => {
  if (before && !now && saved.value) announce(`Counts for ${props.taxid} ${state.value.inTab[props.taxid] ? 'kept in this tab' : 'saved in this browser'}.`);
});

/** Cancel, Download / Continue and Delete come and go with the state: a focused one that leaves hands the focus to the box. */
let focusedInside = false;
onBeforeUpdate(() => {
  focusedInside = root.value !== null && root.value.contains(document.activeElement);
});
onUpdated(() => {
  const active = document.activeElement;
  if (focusedInside && (active === null || active === document.body)) document.getElementById(props.focusId)?.focus();
  focusedInside = false;
});
</script>

<template>
  <div ref="root" class="counts-line" data-testid="counts-line">
    <button
      v-if="!busy && !saved && !partial"
      type="button"
      data-testid="counts-download"
      @click="press"
    >
      Download
    </button>
    <button v-else-if="!busy && partial" type="button" data-testid="counts-continue" @click="press">Continue</button>
    <span class="hint" data-testid="counts-note"><LosatLoader v-if="busy" />{{ sentence }}</span>
    <button v-if="busy" type="button" class="link" data-testid="counts-cancel" @click="counts.cancel()">Cancel</button>
    <button v-else-if="saved || partial" type="button" class="link" data-testid="counts-delete" @click="void counts.delete(taxid)">Delete</button>
    <HelpTip :id="`counts-${taxid}`" label="About the counts download" :text="help" />
    <p v-if="failure" class="error" data-testid="counts-error">{{ failure }}</p>
    <p v-if="note" class="muted" data-testid="counts-note-text">{{ note }}</p>
  </div>
</template>
