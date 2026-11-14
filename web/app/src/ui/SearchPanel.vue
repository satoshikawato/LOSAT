<script setup lang="ts">
// The search screen in the order of NCBI BLAST's two-sequence page (docs/web/ncbi_ui_mapping.md
// §1): the mode tabs and sentence, "Enter Query Sequence", "Enter Subject Sequence",
// "Task Selection" (NCBI's "Program Selection"), the "Run LOSAT" button (NCBI's BLAST button) with the line that says
// what it searches, and "Algorithm parameters", which opens and closes, with a second button
// below it. LOSAT's own controls (threads, keep searching while away) are under the first button, with
// the engine's refusal of the options when there is one. "Run LOSAT" starts the search and the
// page goes to its results (`submitted`), as NCBI's BLAST button does. Explanatory copy follows
// stop-describing-ui (Owner, 2026-11-02): NCBI's sentences stay, LOSAT's explanations go.
import { computed, onUnmounted, ref, watch } from 'vue';
import type { InputRole } from 'losat/domain';
import type { Attention, AttentionState } from '../application/attention';
import type { SearchDraft } from '../application/draft';
import type { CountsLoader } from '../application/counts';
import type { ExampleLoader } from '../application/examples';
import { changedFields, effectiveValue, placementFlags, sectionsAt } from '../domain/parameters';
import { programById, PROGRAMS, searchSummaryParts } from '../domain/programs';
import { announce } from './announcer';
import { REFUSAL_ID } from './refusal';
import { useStore } from './useStore';
import ExampleLine from './ExampleLine.vue';
import HelpTip from './HelpTip.vue';
import InputPanel from './InputPanel.vue';
import LosatLoader from './LosatLoader.vue';
import ParameterField from './ParameterField.vue';
import ParameterForm from './ParameterForm.vue';

const props = defineProps<{
  draft: SearchDraft;
  examples: ExampleLoader;
  counts: CountsLoader;
  attention: Attention;
  attentionState: AttentionState;
}>();
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
  return runCount.value > 1 ? `Run LOSAT (${runCount.value} runs)` : 'Run LOSAT';
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
/**
 * "Loading the engine…" (DW-43) once the describe has waited 0.4 s: the first describe of a
 * mode fetches its engine, which takes seconds; the next ones end at once and show nothing.
 */
const engineWaitShown = ref(false);
let engineWait: ReturnType<typeof setTimeout> | undefined;
watch(
  describeState,
  (now) => {
    clearTimeout(engineWait);
    engineWaitShown.value = false;
    if (now === 'loading') engineWait = setTimeout(() => (engineWaitShown.value = true), 400);
  },
  { immediate: true },
);
onUnmounted(() => clearTimeout(engineWait));
const taskFields = computed(() => sectionsAt(program.value, options.value, 'program').flatMap((section) => section.fields));
const summaryParts = computed(() => searchSummaryParts(program.value, effectiveValue('-task', values.value, options.value)));
/**
 * The Algorithm parameters that differ from the engine's defaults ("♦"); the words of
 * "Additional options" count as one.
 */
const changed = computed(
  () =>
    changedFields(program.value, values.value, options.value, 'algorithm') +
    (state.value.extra[state.value.program].trim() === '' ? 0 : 1),
);
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

/**
 * "Run LOSAT" is marked aria-disabled, not disabled, while the search starts: the focus stays on
 * it (Chromium drops the focus of a button that becomes disabled), and a press then does nothing
 * (WA2-AX-04). Why a search cannot start is told as an alert (WA2-AX-03).
 */
async function submit(from: 'top' | 'bottom'): Promise<void> {
  if (state.value.submitting) return;
  submittedFrom.value = from;
  const result = await props.draft.submit();
  const first = result.ok ? result.runIds?.[0] : undefined;
  if (first !== undefined) emit('submitted', first);
}

watch(
  () => state.value.message,
  (message) => {
    if (message !== undefined) announce(message.text, message.kind === 'error' ? 'alert' : 'status');
  },
);

function restoreDefaults(): void {
  props.draft.resetFields(placementFlags(program.value, 'algorithm'));
  props.draft.setExtra('');
  // The file of a -conffile that Additional options gave goes with them (line RG E1w-3b).
  props.draft.setConffile(undefined);
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
      <legend class="visually-hidden">Mode</legend>
      <label
        v-for="p in PROGRAMS"
        :key="p.id"
        class="program-tab"
        :class="{ selected: state.program === p.id }"
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
    <ExampleLine :examples="examples" :program="state.program" />

    <InputPanel v-for="role in roles" :key="role" :draft="draft" :state="state" :role="role" />

    <fieldset v-if="taskFields.length > 0" class="search-block" data-testid="program-selection">
      <legend>Task Selection</legend>
      <ParameterField v-for="field in taskFields" :key="field.flag" :draft="draft" :state="state" :field="field" />
    </fieldset>

    <div class="blast-row">
      <button
        class="primary-action"
        data-testid="run-losat"
        :aria-disabled="state.submitting || undefined"
        @click="submit('top')"
      >
        <LosatLoader v-if="state.submitting" />{{ buttonLabel }}
      </button>
      <p class="search-summary-line" data-testid="search-summary-line">
        <template v-for="(part, index) in summaryParts" :key="index"><strong v-if="part.strong">{{ part.text }}</strong><template v-else>{{ part.text }}</template></template>
      </p>
      <p v-if="engineWaitShown" class="muted" role="status" data-testid="engine-loading"><LosatLoader />Loading the engine…</p>
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
        <!-- The screen wake lock, on unless turned off; the choice is shown on handheld devices, whose
             screens turn off by themselves and end the search (docs/web/ncbi_ui_mapping.md §4, DW-41). -->
        <template v-if="attentionState.wakeLock !== 'unsupported' && attentionState.handheld">
          <label class="flag-label" :data-wake-lock="attentionState.wakeLock">
            <input
              type="checkbox"
              :checked="attentionState.keepAwake"
              data-testid="keep-awake"
              @change="attention.setKeepAwake(($event.target as HTMLInputElement).checked)"
            />
            Keep searching while you're away
          </label>
          <HelpTip
            id="keep-awake"
            label="About keeping the search going"
            text="Keeps the screen on until the search ends. A phone or tablet ends the search when its screen turns off; switching to another app can end it too."
          />
          <span v-if="attentionState.wakeLock === 'failed'" class="warning" data-testid="wake-lock-failed">
            The browser did not keep the screen on: {{ attentionState.wakeLockError }}
          </span>
        </template>
      </div>
      <!-- The engine's verdict on the options: only a refusal is shown (stop-describing-ui rule 7). The
           fields that it names point to it (refusal.ts, WA2-AX-05). -->
      <p class="validation" :data-state="state.validation.state" data-testid="argv-validation" aria-live="polite">
        <span v-if="state.validation.state === 'invalid'" :id="REFUSAL_ID" class="error">{{ state.validation.message }}</span>
      </p>
      <p v-if="state.message && !messageBelow" class="error" data-testid="search-message">
        {{ state.message.text }}
      </p>
    </div>

    <section class="algorithm-parameters" data-testid="algorithm-parameters">
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
        <ParameterForm :draft="draft" :state="state" :counts="counts" />
        <div class="blast-row">
          <button
            class="primary-action"
            data-testid="run-losat-bottom"
            :aria-disabled="state.submitting || undefined"
            @click="submit('bottom')"
          >
            <LosatLoader v-if="state.submitting" />{{ buttonLabel }}
          </button>
          <p class="search-summary-line">
            <template v-for="(part, index) in summaryParts" :key="index"><strong v-if="part.strong">{{ part.text }}</strong><template v-else>{{ part.text }}</template></template>
          </p>
        </div>
        <p v-if="state.message && messageBelow" class="error" data-testid="search-message">
          {{ state.message.text }}
        </p>
      </div>
    </section>
  </div>
</template>
