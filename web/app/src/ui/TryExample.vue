<script setup lang="ts">
// "Try an example" (S16) at the end of the program's sentence under the tabs: an example fills
// both inputs and belongs to the program, so it sits with the program rather than with one of
// the inputs. It asks before it replaces anything that the user entered (inputs, the program's
// settings, the Job Title), naming what, and never starts a search:
// the line it shows says to choose "Run LOSAT". The program's sentence is the slot.
import { computed, nextTick, ref, watch } from 'vue';
import type { DraftState } from '../application/draft';
import { holdsExample, type Examples, type Replaced } from '../application/examples';
import { programById, type ProgramId } from '../domain/programs';
import { useStore } from './useStore';
import './onboarding.css';

const props = defineProps<{ examples: Examples; state: DraftState }>();
const examplesState = useStore(props.examples.state);
const unavailable = computed(() => props.examples.unavailable(props.state.program));
const loading = computed(() => examplesState.value.loading !== undefined);
/** The program whose example waits for the user's answer, and what it would replace. */
const confirming = ref<{ readonly program: ProgramId; readonly replaced: Replaced }>();
/** "your inputs, your BLASTN settings ... and your Job Title". */
const replacedText = computed(() => {
  const asked = confirming.value;
  if (asked === undefined) return '';
  const parts = [
    ...(asked.replaced.inputs ? ['your inputs'] : []),
    ...(asked.replaced.settings ? [`your ${programById(asked.program).label} settings (they return to their defaults)`] : []),
    ...(asked.replaced.title ? ['your Job Title'] : []),
  ];
  return parts.length > 1 ? `${parts.slice(0, -1).join(', ')} and ${parts.at(-1)!}` : (parts[0] ?? '');
});
const replaceButton = ref<HTMLButtonElement>();
/** "Example loaded", while the inputs are still the example's files. */
const loaded = computed(() => {
  const example = examplesState.value.loaded;
  return example !== undefined && holdsExample(props.state, example.program) ? example.text : undefined;
});

// A question about another program's example does not outlive the change of program.
watch(
  () => props.state.program,
  () => (confirming.value = undefined),
);

async function tryExample(): Promise<void> {
  if (unavailable.value !== undefined) return;
  const program = props.state.program;
  const replaced = props.examples.replaces(program);
  if (replaced !== undefined) {
    confirming.value = { program, replaced };
    await nextTick();
    replaceButton.value?.focus();
    return;
  }
  await props.examples.load(props.state.program);
}

function replace(): void {
  const asked = confirming.value;
  confirming.value = undefined;
  if (asked !== undefined) void props.examples.load(asked.program);
}
</script>

<template>
  <div class="program-line">
    <slot />
    <button
      type="button"
      class="link try-example"
      data-testid="try-example"
      :disabled="unavailable !== undefined || loading"
      :aria-describedby="unavailable !== undefined ? 'try-example-reason' : undefined"
      @click="tryExample"
    >
      Try an example
    </button>
    <span v-if="unavailable !== undefined" id="try-example-reason" class="visually-hidden">
      No example can run: {{ unavailable }}
    </span>
  </div>
  <div v-if="confirming" class="notice example-confirm" role="group" aria-labelledby="example-confirm-text" data-testid="try-example-question">
    <p id="example-confirm-text">
      Replace the query and subject in the form with the {{ programById(confirming.program).label }} example? It removes
      {{ replacedText }} from the form. Runs in the queue keep theirs.
    </p>
    <div class="example-confirm-actions">
      <button ref="replaceButton" type="button" data-testid="try-example-confirm" @click="replace">Replace with the example</button>
      <button type="button" data-testid="try-example-cancel" @click="confirming = undefined">Keep what I entered</button>
    </div>
  </div>
  <p v-if="loading" class="hint" role="status">Loading the example…</p>
  <p v-else-if="loaded" class="note example-loaded" role="status" data-testid="example-loaded">{{ loaded }}</p>
  <p v-if="examplesState.error && !loading" class="error" role="alert" data-testid="example-error">{{ examplesState.error }}</p>
</template>
