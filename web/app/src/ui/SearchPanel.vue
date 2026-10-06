<script setup lang="ts">
import { computed } from 'vue';
import type { RunView } from '../application/coordinator';
import type { SearchDraft } from '../application/draft';
import { isTerminal } from '../domain/run';
import { PROGRAMS, programById, type InputRole } from '../domain/programs';
import { useStore } from './useStore';
import InputPanel from './InputPanel.vue';
import ParameterForm from './ParameterForm.vue';

const props = defineProps<{ draft: SearchDraft; runs: readonly RunView[] }>();
const state = useStore(props.draft.state);
const program = computed(() => programById(state.value.program));
const roles: readonly InputRole[] = ['query', 'subject'];

/** Threads offered besides Auto: up to the logical processors (at most 16). */
const threadChoices = computed(() => {
  const hardware = Number.isInteger(navigator.hardwareConcurrency) ? navigator.hardwareConcurrency : 4;
  return Array.from({ length: Math.max(1, Math.min(16, hardware)) }, (_, i) => i + 1);
});

const waiting = computed(() => props.runs.filter((run) => !isTerminal(run.status)).length);
const runCount = computed(() => {
  void state.value; // recompute with the draft
  return props.draft.runCount();
});
const buttonLabel = computed(() => {
  const runs = runCount.value > 1 ? ` (${runCount.value} runs)` : '';
  return waiting.value === 0 ? `Run locally${runs}` : `Add to queue${runs}`;
});

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
    <p class="program-summary">{{ program.summary }}</p>
    <p v-if="program.unavailable" class="notice" data-testid="program-unavailable">{{ program.unavailable }}</p>

    <div class="inputs">
      <InputPanel v-for="role in roles" :key="role" :draft="draft" :state="state" :role="role" />
    </div>

    <details v-if="!program.unavailable" class="parameters" open>
      <summary>Program parameters</summary>
      <ParameterForm :draft="draft" :state="state" />
    </details>

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

    <p class="validation" :data-state="state.validation.state" data-testid="argv-validation">
      <template v-if="state.validation.state === 'invalid'">
        <span class="error">The engine refuses these options: {{ state.validation.message }}</span>
      </template>
      <template v-else-if="state.validation.state === 'checking'">Checking the options with the engine…</template>
      <template v-else-if="state.validation.state === 'ok'">The engine accepts these options.</template>
    </p>

    <div class="actions">
      <button
        class="primary-action"
        data-testid="add-to-queue"
        :disabled="state.submitting || program.unavailable !== undefined"
        @click="draft.submit()"
      >
        {{ buttonLabel }}
      </button>
      <span v-if="waiting > 0" class="hint">{{ waiting }} {{ waiting === 1 ? 'run' : 'runs' }} in the queue.</span>
    </div>
    <p
      v-if="state.message"
      role="status"
      :class="state.message.kind === 'error' ? 'error' : 'note'"
      data-testid="search-message"
    >
      {{ state.message.text }}
    </p>
  </div>
</template>
