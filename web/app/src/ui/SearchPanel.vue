<script setup lang="ts">
// The search screen in the order of NCBI BLAST's two-sequence page (docs/web/ncbi_ui_mapping.md
// §1): the program tabs and sentence, "Enter Query Sequence", "Enter Subject Sequence",
// "Program Selection", the "Run LOSAT" button (NCBI's BLAST button) with the line that says
// what it searches, and "Algorithm parameters", which opens and closes, with a second button
// below it. LOSAT's own lines (threads, readiness, the engine's check of the options, the
// queue) are under the first button.
import { computed, ref, watch } from 'vue';
import type { RunView } from '../application/coordinator';
import type { SearchDraft } from '../application/draft';
import { isTerminal } from '../domain/run';
import { effectiveValue, formParameters, placementFlags, sectionsAt } from '../domain/parameters';
import { PROGRAMS, programById, searchSummary, type InputRole } from '../domain/programs';
import { useStore } from './useStore';
import InputPanel from './InputPanel.vue';
import ParameterField from './ParameterField.vue';
import ParameterForm from './ParameterForm.vue';

const props = defineProps<{ draft: SearchDraft; runs: readonly RunView[] }>();
const state = useStore(props.draft.state);
const program = computed(() => programById(state.value.program));
const options = computed(() => state.value.description?.parameters);
const values = computed(() => state.value.values[state.value.program]);
const roles: readonly InputRole[] = ['query', 'subject'];

/** Threads offered besides Auto: up to the logical processors (at most 16). */
const threadChoices = computed(() => {
  const hardware = Number.isInteger(navigator.hardwareConcurrency) ? navigator.hardwareConcurrency : 4;
  return Array.from({ length: Math.max(1, Math.min(16, hardware)) }, (_, i) => i + 1);
});

const waiting = computed(() => props.runs.filter((run) => !isTerminal(run.status)).length);
const active = computed(() => props.runs.find((run) => !isTerminal(run.status) && run.status !== 'queued'));
const runCount = computed(() => {
  void state.value; // recompute with the draft
  return props.draft.runCount();
});
const readiness = computed(() => {
  void state.value;
  return props.draft.readiness();
});
/**
 * The button runs LOSAT (the Owner's words, 2026-10-10: not "BLAST", which names NCBI's
 * service), and says when the search waits in the queue and how many runs it makes.
 */
const buttonLabel = computed(() => {
  const notes = [];
  if (runCount.value > 1 && program.value.unavailable === undefined) notes.push(`${runCount.value} runs`);
  if (waiting.value > 0) notes.push('add to queue');
  return notes.length === 0 ? 'Run LOSAT' : `Run LOSAT (${notes.join(', ')})`;
});
const taskFields = computed(() => sectionsAt(program.value, options.value, 'program').flatMap((section) => section.fields));
const summaryLine = computed(() => searchSummary(program.value, effectiveValue('-task', values.value, options.value)));
/** One line about the queue next to the button; on narrow screens the queue is far below. */
const queueSummary = computed(() => {
  if (waiting.value === 0) return undefined;
  const queued = waiting.value - (active.value === undefined ? 0 : 1);
  const parts = [];
  if (active.value !== undefined) parts.push(`Run ${active.value.snapshot.number} is running`);
  if (queued > 0) parts.push(`${queued} waiting`);
  return parts.join(', ');
});

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

function submit(from: 'top' | 'bottom'): void {
  submittedFrom.value = from;
  void props.draft.submit();
}

function restoreDefaults(): void {
  props.draft.resetFields(placementFlags(program.value, 'algorithm'));
}

function showQueue(): void {
  document.getElementById('queue')?.scrollIntoView({ behavior: 'smooth', block: 'start' });
}

function onThreads(event: Event): void {
  const value = (event.target as HTMLSelectElement).value;
  props.draft.setThreads(value === 'auto' ? 'auto' : Number(value));
}
</script>

<template>
  <div class="search-panel">
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
        data-testid="add-to-queue"
        :disabled="state.submitting || program.unavailable !== undefined"
        @click="submit('top')"
      >
        {{ buttonLabel }}
      </button>
      <p v-if="program.unavailable" class="hint">{{ program.label }} is not available yet.</p>
      <p v-else class="search-summary-line" data-testid="search-summary-line">{{ summaryLine }}</p>
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
        <span class="hint">Auto runs small searches on one thread and larger ones on up to four.</span>
      </div>
      <p v-if="readiness" class="notice" data-testid="draft-readiness">{{ readiness }}</p>
      <p class="validation" :data-state="state.validation.state" data-testid="argv-validation" aria-live="polite">
        <template v-if="state.validation.state === 'invalid'">
          <span class="error">The engine refuses these options: {{ state.validation.message }}</span>
        </template>
        <template v-else-if="state.validation.state === 'checking'">Checking the options with the engine…</template>
        <template v-else-if="state.validation.state === 'ok'">The engine accepts these options.</template>
      </p>
      <p v-if="queueSummary" class="hint" data-testid="queue-summary">
        {{ queueSummary }}.
        <button type="button" class="link" @click="showQueue">Show the queue</button>
      </p>
      <p
        v-if="state.message && !messageBelow"
        role="status"
        :class="state.message.kind === 'error' ? 'error' : 'note'"
        data-testid="search-message"
      >
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
            data-testid="add-to-queue-bottom"
            :disabled="state.submitting || program.unavailable !== undefined"
            @click="submit('bottom')"
          >
            {{ buttonLabel }}
          </button>
          <p class="search-summary-line">{{ summaryLine }}</p>
        </div>
        <p
          v-if="state.message && messageBelow"
          role="status"
          :class="state.message.kind === 'error' ? 'error' : 'note'"
          data-testid="search-message"
        >
          {{ state.message.text }}
        </p>
      </div>
    </section>
  </div>
</template>
