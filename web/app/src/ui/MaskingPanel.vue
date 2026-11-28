<script setup lang="ts">
// The Masking main tab (WW4c): dustmasker, segmasker and convert2blastmask in one tab. A tool switcher
// (DustMasker | SegMasker | convert2blastmask) with the tool's NCBI sentence from `-xmlhelp`, the Input
// block (`-infmt`; FASTA paste/files or the `-in` name of a BLAST database, as WindowMasker's), the
// options that the command line requires (convert2blastmask's `-masking_algorithm` and
// `-masking_options`), the Output block (`-out`, `-outfmt`), Run LOSAT with the command line,
// validateTool's refusals in place, and the collapsible "<Tool> options" with every other option of
// `-xmlhelp`, "Files of the run" and "Additional options".
import { computed, ref, watch } from 'vue';
import type { MaskingDraft } from '../application/masking-draft';
import { changedMaskingOptionsCount, toolValidationSummary } from '../domain/tool-argv';
import {
  MASKING_TOOLS,
  maskingBlockParameter,
  maskingRequiredParameters,
  maskingToolSummary,
  toolLabel,
  toolTakesValue,
  toToolParameterField,
} from '../domain/tools';
import { announce } from './announcer';
import { maskingFormModel, MS_REFUSAL_ID } from './form-model';
import OptionName from './OptionName.vue';
import ParameterField from './ParameterField.vue';
import ParameterForm from './ParameterForm.vue';
import ToolFastaInputs from './ToolFastaInputs.vue';
import ToolRunBar from './ToolRunBar.vue';
import ToolRunFiles from './ToolRunFiles.vue';
import { useStore } from './useStore';

const props = defineProps<{ draft: MaskingDraft }>();

const emit = defineEmits<{ submitted: [runId: string] }>();

const state = useStore(props.draft.state);
const formModel = computed(() => maskingFormModel(props.draft, state.value));

const describeState = computed(() => {
  if (state.value.description !== undefined) return 'ready';
  return state.value.descriptionError === undefined ? 'loading' : 'failed';
});

const summary = computed(() => maskingToolSummary(state.value.tool, state.value.description));
const label = computed(() => toolLabel(state.value.tool));

const blockField = (flag: string) =>
  computed(() => {
    const param = maskingBlockParameter(state.value.description, flag);
    return param === undefined ? undefined : toToolParameterField(param);
  });
const infmtField = blockField('-infmt');
const outfmtField = blockField('-outfmt');
const requiredFields = computed(() => maskingRequiredParameters(state.value.description, state.value.tool).map(toToolParameterField));

const commandWords = computed(() => props.draft.previewCommandWords(state.value));
const takesValue = computed(() => toolTakesValue(state.value.description));

const runCount = computed(() => props.draft.runCount(state.value));
const buttonLabel = computed(() => {
  if (state.value.submitting) return 'Starting…';
  return runCount.value > 1 ? `Run LOSAT (${runCount.value} runs)` : 'Run LOSAT';
});

const changed = computed(() =>
  changedMaskingOptionsCount(state.value.tool, state.value.values[state.value.tool], state.value.extra[state.value.tool], state.value.description),
);

const open = ref(changed.value > 0);
watch(
  () => state.value.tool,
  () => {
    if (changed.value > 0) open.value = true;
  },
);

const dragging = ref(false);

function onDrop(event: DragEvent): void {
  dragging.value = false;
  const files = event.dataTransfer?.files;
  if (files !== undefined && files.length > 0) props.draft.dropInputFiles([...files]);
}

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
  <div class="search-panel wm-panel ms-panel" data-testid="ms-panel" :data-describe="describeState" :data-tool="state.tool">
    <h2 class="visually-hidden">Masking</h2>

    <fieldset class="program-tabs" data-testid="ms-tool-tabs">
      <legend class="visually-hidden">Masking tool</legend>
      <label v-for="tool in MASKING_TOOLS" :key="tool.id" class="program-tab" :class="{ selected: state.tool === tool.id }">
        <input
          type="radio"
          name="ms-tool"
          :value="tool.id"
          :checked="state.tool === tool.id"
          :data-testid="`ms-tool-${tool.id}`"
          @change="draft.setTool(tool.id)"
        />
        {{ tool.label }}
      </label>
    </fieldset>
    <p class="program-summary" data-testid="ms-tool-summary">{{ summary }}</p>

    <fieldset
      class="search-block input-panel"
      :class="{ dragging }"
      data-testid="ms-input-block"
      @dragover.prevent="dragging = true"
      @dragleave="dragging = false"
      @drop.prevent="onDrop"
    >
      <legend>Input</legend>
      <ParameterField v-if="infmtField" :model="formModel" :field="infmtField" />
      <ToolFastaInputs :draft="draft" :state="state" prefix="ms" />
    </fieldset>

    <fieldset v-if="requiredFields.length > 0" class="search-block" data-testid="ms-required-block">
      <legend>Masking algorithm</legend>
      <ParameterField v-for="field in requiredFields" :key="field.flag" :model="formModel" :field="field" />
    </fieldset>

    <fieldset class="search-block" data-testid="ms-output-block">
      <legend>Output</legend>
      <div class="param-row">
        <label class="param-label" for="ms-out">
          <OptionName class="flag" name="-out" />
        </label>
        <div class="param-control">
          <input
            id="ms-out"
            type="text"
            autocomplete="off"
            spellcheck="false"
            placeholder="- (standard output)"
            :value="state.outNames[state.tool]"
            data-testid="ms-out"
            @input="draft.setOutName(($event.target as HTMLInputElement).value)"
          />
          <span class="param-hint">Output file name (leave empty for standard output)</span>
        </div>
      </div>
      <ParameterField v-if="outfmtField" :model="formModel" :field="outfmtField" />
    </fieldset>

    <ToolRunBar
      :submitting="state.submitting"
      :button-label="buttonLabel"
      button-test-id="ms-run-losat"
      :command-words="commandWords"
      :takes-value="takesValue"
      command-test-id="ms-command"
      @run="submit('top')"
    />

    <div class="blast-notes">
      <p class="validation" :data-state="state.validation.state" data-testid="ms-argv-validation" aria-live="polite">
        <span v-if="state.validation.state === 'invalid'" :id="MS_REFUSAL_ID" class="error">{{ toolValidationSummary(state.validation.message) }}</span>
      </p>
      <p v-if="state.message && !messageBelow" class="error" data-testid="ms-message">{{ state.message.text }}</p>
    </div>

    <section class="algorithm-parameters" data-testid="ms-options">
      <h3 class="parameters-heading">
        <button
          type="button"
          class="parameters-bar"
          :aria-expanded="open"
          aria-controls="ms-options-body"
          data-testid="ms-options-toggle"
          @click="open = !open"
        >
          <span class="bar-mark" aria-hidden="true">{{ open ? '−' : '+' }}</span>
          {{ label }} options
          <span v-if="changed > 0" class="bar-count" data-testid="ms-options-changed">({{ changed }} changed)</span>
        </button>
      </h3>
      <div v-if="open" id="ms-options-body" class="parameters-body">
        <div v-if="!state.descriptionUnavailable" class="parameters-tools">
          <p class="parameters-note">
            Parameter values that differ from the default are highlighted in yellow and marked with
            <span aria-hidden="true">♦</span><span class="visually-hidden">a diamond</span> sign
          </p>
          <button type="button" data-testid="ms-restore-defaults" @click="draft.resetOptions()">
            Restore default options
          </button>
        </div>
        <ParameterForm :model="formModel" test-id="ms-parameter-form">
          <template #before-extra>
            <div class="search-block parameter-section" data-testid="ms-run-files">
              <div class="section-legend">Files of the run</div>
              <ToolRunFiles :draft="draft" :files="state.runFiles" :error="state.runFilesError" prefix="ms"
                >Files opened by name (<OptionName name="-conffile" />, BLAST database files)</ToolRunFiles
              >
            </div>
          </template>
        </ParameterForm>
        <ToolRunBar
          :submitting="state.submitting"
          :button-label="buttonLabel"
          button-test-id="ms-run-losat-bottom"
          :command-words="commandWords"
          :takes-value="takesValue"
          @run="submit('bottom')"
        />
        <p v-if="state.message && messageBelow" class="error" data-testid="ms-message">{{ state.message.text }}</p>
      </div>
    </section>
  </div>
</template>
