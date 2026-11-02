<script setup lang="ts">
// The search screen in the order of NCBI BLAST's two-sequence page (docs/web/ncbi_ui_mapping.md
// §1): the program tabs and sentence, "Enter Query Sequence", "Enter Subject Sequence",
// "Program Selection", the "Run LOSAT" button (NCBI's BLAST button) with the line that says
// what it searches, and "Algorithm parameters", which opens and closes, with a second button
// below it. LOSAT's own controls (threads, keep the screen on) are under the first button, with
// the engine's refusal of the options when there is one. "Run LOSAT" starts the search and the
// page goes to its results (`submitted`), as NCBI's BLAST button does. Explanatory copy follows
// stop-describing-ui (Owner, 2026-11-02): NCBI's sentences stay, LOSAT's explanations go.
import { computed, ref, watch } from 'vue';
import type { Attention, AttentionState } from '../application/attention';
import type { SearchDraft } from '../application/draft';
import { effectiveValue, formParameters, placementFlags, sectionsAt } from '../domain/parameters';
import { PROGRAMS, programById, searchSummary, type InputRole } from '../domain/programs';
import { useStore } from './useStore';
import InputPanel from './InputPanel.vue';
import ParameterField from './ParameterField.vue';
import ParameterForm from './ParameterForm.vue';

const props = defineProps<{ draft: SearchDraft; attention: Attention; attentionState: AttentionState }>();
/** A search was started: the ID of its first run (a group of separate searches starts several). */
const emit = defineEmits<{ submitted: [runId: string] }>();
const state = useStore(props.draft.state);
const program = computed(() => programById(state.value.program));
const options = computed(() => state.value.description?.parameters);
const values = computed(() => state.value.values[state.value.program]);
const roles: readonly InputRole[] = ['query', 'subject'];

/** Threads offered besides Auto: up to the logical processors (at most 32, Auto's maximum). */
const threadChoices = computed(() => {
  const hardware = Number.isInteger(navigator.hardwareConcurrency) ? navigator.hardwareConcurrency : 4;
  return Array.from({ length: Math.max(1, Math.min(32, hardware)) }, (_, i) => i + 1);
});

const runCount = computed(() => {
  void state.value; // recompute with the draft
  return props.draft.runCount();
});
/**
 * The button runs LOSAT (the Owner's words, 2026-10-10: not "BLAST", which names NCBI's
 * service), and says how many runs separate searches make.
 */
const buttonLabel = computed(() => {
  if (state.value.submitting) return 'Starting…';
  return runCount.value > 1 && program.value.unavailable === undefined ? `Run LOSAT (${runCount.value} runs)` : 'Run LOSAT';
});
/**
 * Whether the engine's describe of the program has come. It comes after the page or the
 * program changes, and adds the task's help and the task to the line next to the button,
 * which moves the controls below them; the E2E tests wait for it before they click.
 */
const describeState = computed(() => {
  if (state.value.description !== undefined) return 'ready';
  return state.value.descriptionError === undefined ? 'loading' : 'failed';
});
const taskFields = computed(() => sectionsAt(program.value, options.value, 'program').flatMap((section) => section.fields));
const summaryLine = computed(() => searchSummary(program.value, effectiveValue('-task', values.value, options.value)));
/** The Algorithm parameters that the argv writes: each differs from the engine's default. */
const changed = computed(() => formParameters(program.value, values.value, options.value, 'algorithm').length);
/**
 * "Algorithm parameters" starts closed, as on NCBI's page, unless one of its values is
 * written to the argv; a program whose values differ from the defaults opens it when it is
 * chosen, so that no written value is hidden from view. It never closes by itself.
 */
const open = ref(changed.value > 0);
watch(
  () => state.value.program,
  () => {
    if (changed.value > 0) open.value = true;
  },
);
/** Where the latest search was asked for: its message is shown next to that button, while it is shown. */
const submittedFrom = ref<'top' | 'bottom'>('top');
const messageBelow = computed(() => submittedFrom.value === 'bottom' && open.value);

async function submit(from: 'top' | 'bottom'): Promise<void> {
  submittedFrom.value = from;
  const result = await props.draft.submit();
  const first = result.ok ? result.runIds?.[0] : undefined;
  if (first !== undefined) emit('submitted', first);
}

function restoreDefaults(): void {
  props.draft.resetFields(placementFlags(program.value, 'algorithm'));
}

function onThreads(event: Event): void {
  const value = (event.target as HTMLSelectElement).value;
  props.draft.setThreads(value === 'auto' ? 'auto' : Number(value));
}
</script>

<template>
  <div class="search-panel" data-testid="search-panel" :data-describe="describeState">
    <h2 class="visually-hidden">Search</h2>
    <fieldset class="program-tabs" data-testid="program-tabs">
      <legend class="visually-hidden">Program</legend>
      <label
        v-for="p in PROGRAMS"
        :key="p.id"
        class="program-tab"
        :class="{ selected: state.program === p.id, unavailable: p.unavailable !== undefined }"
      >
        <input
          type="radio"
          name="program"
          :value="p.id"
          :checked="state.program === p.id"
          :data-testid="`program-${p.id}`"
          @change="draft.setProgram(p.id)"
        />
        {{ p.label }}
      </label>
    </fieldset>
    <p class="program-summary" data-testid="program-summary">{{ program.summary }}</p>
    <p v-if="program.unavailable" class="notice" data-testid="program-unavailable">{{ program.unavailable }}</p>

    <InputPanel v-for="role in roles" :key="role" :draft="draft" :state="state" :role="role" />

    <fieldset v-if="taskFields.length > 0" class="search-block" data-testid="program-selection">
      <legend>Program Selection</legend>
      <ParameterField v-for="field in taskFields" :key="field.flag" :draft="draft" :state="state" :field="field" />
    </fieldset>

    <div class="blast-row">
      <button
        class="primary-action"
        data-testid="run-losat"
        :disabled="state.submitting || program.unavailable !== undefined"
        @click="submit('top')"
      >
        {{ buttonLabel }}
      </button>
      <p v-if="!program.unavailable" class="search-summary-line" data-testid="search-summary-line">{{ summaryLine }}</p>
    </div>

    <div class="blast-notes">
      <div class="run-options">
        <label>
          Threads
          <select data-testid="threads" :value="String(state.threads)" @change="onThreads">
            <option value="auto">Auto</option>
            <option v-for="n in threadChoices" :key="n" :value="String(n)">{{ n }}</option>
          </select>
        </label>
        <label v-if="attentionState.wakeLock !== 'unsupported'" class="flag-label" :data-wake-lock="attentionState.wakeLock">
          <input
            type="checkbox"
            :checked="attentionState.keepAwake"
            data-testid="keep-awake"
            @change="attention.setKeepAwake(($event.target as HTMLInputElement).checked)"
          />
          Keep screen on while searching
        </label>
        <span v-if="attentionState.wakeLock === 'failed'" class="warning" data-testid="wake-lock-failed">
          The browser did not keep the screen on: {{ attentionState.wakeLockError }}
        </span>
      </div>
      <!-- The engine's verdict on the options: only a refusal is shown (stop-describing-ui rule 7). -->
      <p class="validation" :data-state="state.validation.state" data-testid="argv-validation" aria-live="polite">
        <span v-if="state.validation.state === 'invalid'" class="error">{{ state.validation.message }}</span>
      </p>
      <p v-if="state.message && !messageBelow" role="status" class="error" data-testid="search-message">
        {{ state.message.text }}
      </p>
    </div>

    <section v-if="!program.unavailable" class="algorithm-parameters" data-testid="algorithm-parameters">
      <h3 class="parameters-heading">
        <button
          type="button"
          class="parameters-bar"
          :aria-expanded="open"
          aria-controls="algorithm-parameters-body"
          data-testid="algorithm-parameters-toggle"
          @click="open = !open"
        >
          <span class="bar-mark" aria-hidden="true">{{ open ? '−' : '+' }}</span>
          Algorithm parameters
          <span v-if="changed > 0" class="bar-count" data-testid="algorithm-parameters-changed">({{ changed }} changed)</span>
        </button>
      </h3>
      <div v-if="open" id="algorithm-parameters-body" class="parameters-body">
        <div class="parameters-tools">
          <p class="parameters-note">
            Parameter values that differ from the default are highlighted in yellow and marked with
            <span aria-hidden="true">♦</span><span class="visually-hidden">a diamond</span> sign
          </p>
          <button type="button" data-testid="restore-defaults" @click="restoreDefaults">Restore default search parameters</button>
        </div>
        <ParameterForm :draft="draft" :state="state" />
        <div class="blast-row">
          <button
            class="primary-action"
            data-testid="run-losat-bottom"
            :disabled="state.submitting || program.unavailable !== undefined"
            @click="submit('bottom')"
          >
            {{ buttonLabel }}
          </button>
          <p class="search-summary-line">{{ summaryLine }}</p>
        </div>
        <p v-if="state.message && messageBelow" role="status" class="error" data-testid="search-message">
          {{ state.message.text }}
        </p>
      </div>
    </section>
  </div>
</template>
