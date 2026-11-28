<script setup lang="ts">
// The WindowMasker main tab (WW3b DESIGN §1-§2):
// Stage switcher (Make counts | Mask | Convert counts), NCBI help sentence, Input / Counts / Output
// blocks, Run LOSAT with live command preview, in-place validateTool refusals, and collapsible
// "WindowMasker options" showing every stage option from `-xmlhelp` plus "Files of the run" and
// "Additional options".
import { computed, ref, watch } from 'vue';
import type { CountsLoader } from '../application/counts';
import type { Coordinator } from '../application/coordinator';
import type { ToolDraft } from '../application/tool-draft';
import { changedToolOptionsCount, toolValidationSummary } from '../domain/tool-argv';
import {
  stageBlockParameter,
  toolTakesValue,
  toToolParameterField,
  WINDOWMASKER_STAGES,
  windowMaskerStageSummary,
} from '../domain/tools';
import { announce } from './announcer';
import { windowMaskerFormModel, WM_REFUSAL_ID } from './form-model';
import OptionName from './OptionName.vue';
import ParameterField from './ParameterField.vue';
import ParameterForm from './ParameterForm.vue';
import ToolInputs from './ToolInputs.vue';
import ToolRunBar from './ToolRunBar.vue';
import ToolRunFiles from './ToolRunFiles.vue';
import { useStore } from './useStore';

const props = defineProps<{
  draft: ToolDraft;
  coordinator: Coordinator;
  counts: CountsLoader;
}>();

const emit = defineEmits<{ submitted: [runId: string] }>();

const state = useStore(props.draft.state);
const appState = useStore(props.coordinator.state);
const countsState = useStore(props.counts.state);
const formModel = computed(() => windowMaskerFormModel(props.draft, state.value));

const describeState = computed(() => {
  if (state.value.description !== undefined) return 'ready';
  return state.value.descriptionError === undefined ? 'loading' : 'failed';
});

const stageSummary = computed(() => windowMaskerStageSummary(state.value.stage, state.value.description));

const sessionCounts = computed(() => {
  void appState.value.runs;
  return props.coordinator.sessionMadeCounts();
});

const formatField = computed(() => {
  const flag = state.value.stage === 'mask' ? '-outfmt' : '-sformat';
  const param = stageBlockParameter(state.value.description, state.value.stage, flag);
  return param === undefined ? undefined : toToolParameterField(param);
});

const commandWords = computed(() => props.draft.previewCommandWords(state.value));
const takesValue = computed(() => toolTakesValue(state.value.description));

const runCount = computed(() => props.draft.runCount(state.value));

const buttonLabel = computed(() => {
  if (state.value.submitting) return 'Starting…';
  return runCount.value > 1 ? `Run LOSAT (${runCount.value} runs)` : 'Run LOSAT';
});

const changed = computed(() =>
  changedToolOptionsCount(
    state.value.stage,
    state.value.values[state.value.stage],
    state.value.extra[state.value.stage],
    state.value.description,
  ),
);

const open = ref(changed.value > 0);
watch(
  () => state.value.stage,
  () => {
    if (changed.value > 0) open.value = true;
  },
);

const submittedFrom = ref<'top' | 'bottom'>('top');
const messageBelow = computed(() => submittedFrom.value === 'bottom' && open.value);

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

</script>

<template>
  <div class="search-panel wm-panel" data-testid="wm-panel" :data-describe="describeState">
    <h2 class="visually-hidden">WindowMasker</h2>

    <fieldset class="program-tabs" data-testid="wm-stage-tabs">
      <legend class="visually-hidden">WindowMasker stage</legend>
      <label
        v-for="stage in WINDOWMASKER_STAGES"
        :key="stage.id"
        class="program-tab"
        :class="{ selected: state.stage === stage.id }"
      >
        <input
          type="radio"
          name="wm-stage"
          :value="stage.id"
          :checked="state.stage === stage.id"
          :data-testid="`wm-stage-${stage.id}`"
          @change="draft.setStage(stage.id)"
        />
        {{ stage.label }}
      </label>
    </fieldset>
    <p class="program-summary" data-testid="wm-stage-summary">{{ stageSummary }}</p>

    <ToolInputs
      :draft="draft"
      :state="state"
      :session-counts="sessionCounts"
      :saved-counts="countsState.list"
      :model="formModel"
    />

    <fieldset class="search-block" data-testid="wm-output-block">
      <legend>Output</legend>
      <div class="param-row">
        <label class="param-label" for="wm-out">
          <OptionName class="flag" name="-out" />
        </label>
        <div class="param-control">
          <input
            id="wm-out"
            type="text"
            autocomplete="off"
            spellcheck="false"
            placeholder="- (standard output)"
            :value="state.outNames[state.stage]"
            data-testid="wm-out"
            @input="draft.setOutName(($event.target as HTMLInputElement).value)"
          />
          <span class="param-hint">Output file name (type <code>-</code> for standard output)</span>
        </div>
      </div>
      <ParameterField v-if="formatField" :model="formModel" :field="formatField" />
    </fieldset>

    <ToolRunBar
      :submitting="state.submitting"
      :button-label="buttonLabel"
      button-test-id="wm-run-losat"
      :command-words="commandWords"
      :takes-value="takesValue"
      command-test-id="wm-command"
      @run="submit('top')"
    />

    <div class="blast-notes">
      <p class="validation" :data-state="state.validation.state" data-testid="wm-argv-validation" aria-live="polite">
        <span v-if="state.validation.state === 'invalid'" :id="WM_REFUSAL_ID" class="error">{{ toolValidationSummary(state.validation.message) }}</span>
      </p>
      <p v-if="state.message && !messageBelow" class="error" data-testid="wm-message">{{ state.message.text }}</p>
    </div>

    <section class="algorithm-parameters" data-testid="wm-options">
      <h3 class="parameters-heading">
        <button
          type="button"
          class="parameters-bar"
          :aria-expanded="open"
          aria-controls="wm-options-body"
          data-testid="wm-options-toggle"
          @click="open = !open"
        >
          <span class="bar-mark" aria-hidden="true">{{ open ? '−' : '+' }}</span>
          WindowMasker options
          <span v-if="changed > 0" class="bar-count" data-testid="wm-options-changed">({{ changed }} changed)</span>
        </button>
      </h3>
      <div v-if="open" id="wm-options-body" class="parameters-body">
        <div v-if="!state.descriptionUnavailable" class="parameters-tools">
          <p class="parameters-note">
            Parameter values that differ from the default are highlighted in yellow and marked with
            <span aria-hidden="true">♦</span><span class="visually-hidden">a diamond</span> sign
          </p>
          <button type="button" data-testid="wm-restore-defaults" @click="draft.resetOptions()">
            Restore default options
          </button>
        </div>
        <ParameterForm :model="formModel" test-id="wm-parameter-form">
          <template #before-extra>
            <div class="search-block parameter-section" data-testid="wm-run-files">
              <div class="section-legend">Files of the run</div>
              <ToolRunFiles :draft="draft" :files="state.runFiles" :error="state.runFilesError" prefix="wm"
                >Files opened by name (<OptionName name="-ids" />, <OptionName name="-exclude_ids" />, <OptionName name="-conffile" />, <OptionName name="-fa_list" />, BLAST database files)</ToolRunFiles
              >
            </div>
          </template>
        </ParameterForm>
        <ToolRunBar
          :submitting="state.submitting"
          :button-label="buttonLabel"
          button-test-id="wm-run-losat-bottom"
          :command-words="commandWords"
          :takes-value="takesValue"
          @run="submit('bottom')"
        />
        <p v-if="state.message && messageBelow" class="error" data-testid="wm-message">{{ state.message.text }}</p>
      </div>
    </section>
  </div>
</template>
