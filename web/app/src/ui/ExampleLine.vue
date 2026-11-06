<script setup lang="ts">
// "Example" under the mode tabs (line WF-2, docs/web/fetch_design.md §9): a button per set of the
// mode, each with its in-place sentence: what pressing it sends ("Downloads 8.0 MB from
// ftp.ncbi.nlm.nih.gov"), or the copy saved in this browser with "Save to folder…" (or
// "Download", line WF-3) and Delete, or the download's progress with Cancel. The rest waits
// behind "?" (stop-describing-ui). A set that loads, or that cannot, is told to a screen reader
// (WA2-AX-03), and a button that cannot be pressed for a while is marked aria-disabled, so that it
// keeps the focus (WA2-AX-04).
import { computed, nextTick, watch } from 'vue';
import type { ExampleLoader } from '../application/examples';
import { exampleById, exampleHost, examplesFor, type ExampleId, type ExampleSet } from '../domain/examples';
import type { ProgramId } from '../domain/programs';
import { announce } from './announcer';
import { formatBytes } from './format';
import HelpTip from './HelpTip.vue';
import LosatLoader from './LosatLoader.vue';
import { useStore } from './useStore';

const props = defineProps<{ examples: ExampleLoader; program: ProgramId }>();
const state = useStore(props.examples.state);
const sets = computed(() => examplesFor(props.program));
const load = computed(() => state.value.load);
/** A note on a set of this mode. */
const note = computed(() => (sets.value.some((set) => set.id === state.value.note?.id) ? state.value.note : undefined));

watch(
  () => props.program,
  () => props.examples.refresh(),
  { immediate: true },
);

const help = computed(
  () =>
    'Example puts a query in the query box and a public set in the subject list. It asks only for the file named, ' +
    'from the site named, and keeps the file in this browser: the next time, it is read from there without a request. ' +
    (props.examples.copyTo === undefined
      ? ''
      : `${props.examples.copyTo === 'folder' ? 'Save to folder…' : 'Download'} writes that file outside this browser's storage, ` +
        'which clearing the site data empties; "Open FASTA files…" opens it again. ') +
    sets.value.map((set) => `${set.label}: ${set.provenance}`).join(' '),
);

/** The sentence beside a set's button. */
function sentence(set: ExampleSet): string {
  const progress = load.value?.id === set.id ? load.value.progress : undefined;
  switch (progress?.state) {
    case 'reading':
      return 'Opening the saved copy…';
    case 'downloading':
      if (progress.bytes === 0 && progress.total === undefined) return `Downloading from ${exampleHost(set)}…`;
      return `Downloading ${formatBytes(progress.bytes)}${progress.total === undefined ? '' : ` of ${formatBytes(progress.total)}`} from ${exampleHost(set)}…`;
    case 'saving':
      return 'Saving in this browser…';
  }
  const copy = state.value.saved[set.id];
  if (copy) {
    const version = copy.version ?? `saved ${new Date(copy.savedAt).toISOString().slice(0, 10)}`;
    return `Saved in this browser (${version}, ${formatBytes(copy.bytes)})`;
  }
  return `Downloads ${formatBytes(set.downloadBytes)} from ${exampleHost(set)}`;
}

const busy = (set: ExampleSet) => {
  const progress = load.value?.id === set.id ? load.value.progress.state : undefined;
  return progress !== undefined && progress !== 'failed';
};

/** The set being loaded, until its load ends: loaded (its files are in the inputs), failed or cancelled. */
let loading: ExampleId | undefined;
let cancelled = false;
watch(load, (now) => {
  if (now !== undefined && now.progress.state !== 'failed') {
    if (now.id !== loading) cancelled = false;
    loading = now.id;
    return;
  }
  const id = loading;
  loading = undefined;
  if (id === undefined) return;
  if (now?.progress.state === 'failed') announce(now.progress.message, 'alert');
  else if (!cancelled) announce(`Example ${exampleById(id).label} loaded.${note.value?.id === id ? ` ${note.value.text}` : ''}`);
  cancelled = false;
});

const buttons: Partial<Record<ExampleId, HTMLButtonElement>> = {};
function keepButton(id: ExampleId, element: unknown): void {
  if (element instanceof HTMLButtonElement) buttons[id] = element;
  else delete buttons[id];
}
function press(set: ExampleSet): void {
  if (!busy(set)) props.examples.load(set.id);
}
/** Cancel leaves with the load: the focus goes back to the set's button. */
async function cancel(set: ExampleSet): Promise<void> {
  cancelled = true;
  props.examples.cancel();
  await nextTick();
  buttons[set.id]?.focus();
}
function save(set: ExampleSet): void {
  if (state.value.copying !== set.id) props.examples.saveCopy(set.id);
}
function remove(set: ExampleSet): void {
  if (state.value.copying !== set.id) void props.examples.delete(set.id);
}
</script>

<template>
  <div v-if="sets.length > 0" class="example-line" data-testid="example-line">
    <span class="example-title">Example</span>
    <ul class="example-sets">
      <li v-for="set in sets" :key="set.id" :data-testid="`example-${set.id}-row`">
        <button
          :ref="(element) => keepButton(set.id, element)"
          type="button"
          :aria-disabled="busy(set) || undefined"
          :data-testid="`example-${set.id}`"
          @click="press(set)"
        >
          {{ set.label }}
        </button>
        <span class="hint" :data-testid="`example-${set.id}-note`"><LosatLoader v-if="busy(set)" />{{ sentence(set) }}</span>
        <button
          v-if="busy(set)"
          type="button"
          class="link"
          :data-testid="`example-${set.id}-cancel`"
          @click="cancel(set)"
        >
          Cancel
        </button>
        <template v-else-if="state.saved[set.id]">
          <button
            v-if="examples.copyTo"
            type="button"
            class="link"
            :aria-disabled="state.copying === set.id || undefined"
            :data-testid="`example-${set.id}-save`"
            @click="save(set)"
          >
            {{ examples.copyTo === 'folder' ? 'Save to folder…' : 'Download' }}
          </button>
          <button
            type="button"
            class="link"
            :aria-disabled="state.copying === set.id || undefined"
            :data-testid="`example-${set.id}-delete`"
            @click="remove(set)"
          >
            Delete
          </button>
        </template>
      </li>
    </ul>
    <HelpTip id="example" label="About the examples" :text="help" />
    <p v-if="load?.progress.state === 'failed'" class="error" data-testid="example-error">
      {{ load.progress.message }}
    </p>
    <p v-if="note" class="muted" data-testid="example-note">{{ note.text }}</p>
  </div>
</template>
