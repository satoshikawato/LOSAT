<script setup lang="ts">
// "Example" under the mode tabs (line WF-2, docs/web/fetch_design.md §9): a button per set of the
// mode, each with its in-place sentence: what pressing it sends ("Downloads 8.0 MB from
// ftp.ncbi.nlm.nih.gov"), or the copy saved in this browser with "Save to folder…" (or
// "Download", line WF-3) and Delete, or the download's progress with Cancel. The rest waits
// behind "?" (stop-describing-ui).
import { computed, watch } from 'vue';
import type { ExampleLoader } from '../application/examples';
import { exampleHost, examplesFor, type ExampleSet } from '../domain/examples';
import type { ProgramId } from '../domain/programs';
import { formatBytes } from './format';
import HelpTip from './HelpTip.vue';
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
</script>

<template>
  <div v-if="sets.length > 0" class="example-line" data-testid="example-line">
    <span class="example-title">Example</span>
    <ul class="example-sets">
      <li v-for="set in sets" :key="set.id" :data-testid="`example-${set.id}-row`">
        <button type="button" :disabled="busy(set)" :data-testid="`example-${set.id}`" @click="examples.load(set.id)">
          {{ set.label }}
        </button>
        <span class="hint" :data-testid="`example-${set.id}-note`">{{ sentence(set) }}</span>
        <button
          v-if="busy(set)"
          type="button"
          class="link"
          :data-testid="`example-${set.id}-cancel`"
          @click="examples.cancel()"
        >
          Cancel
        </button>
        <template v-else-if="state.saved[set.id]">
          <button
            v-if="examples.copyTo"
            type="button"
            class="link"
            :disabled="state.copying === set.id"
            :data-testid="`example-${set.id}-save`"
            @click="examples.saveCopy(set.id)"
          >
            {{ examples.copyTo === 'folder' ? 'Save to folder…' : 'Download' }}
          </button>
          <button type="button" class="link" :data-testid="`example-${set.id}-delete`" @click="examples.delete(set.id)">
            Delete
          </button>
        </template>
      </li>
    </ul>
    <HelpTip id="example" label="About the examples" :text="help" />
    <p v-if="load?.progress.state === 'failed'" class="error" role="alert" data-testid="example-error">
      {{ load.progress.message }}
    </p>
    <p v-if="note" class="muted" data-testid="example-note">{{ note.text }}</p>
  </div>
</template>
