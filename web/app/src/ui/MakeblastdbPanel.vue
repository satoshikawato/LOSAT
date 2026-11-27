<script setup lang="ts">
// The makeblastdb main tab (line WW, WW4b2; design DESIGN_ww4b.md §3): builds a BLAST database in the
// browser with every option of makeblastdb's `-xmlhelp`. The Input block (`-input_type`, the input that
// `-in` names: FASTA pasted or opened, ASN.1 files, or a BLAST database opened as files of the run, and
// `-dbtype`), the Output block (`-out`), the command beside Run LOSAT (ToolRunBar, the argv that the engine
// runs), validateTool's refusals in place, and "makeblastdb options" in the groups of NCBI's `-help` with
// Files of the run and Additional options.
import { computed, ref, watch } from 'vue';
import type { Coordinator } from '../application/coordinator';
import type { MakeblastdbDraft } from '../application/makeblastdb-draft';
import type { ParameterField as Field } from '../domain/programs';
import { changedMakeblastdbOptionsCount, toolValidationSummary } from '../domain/tool-argv';
import { toolTakesValue, toToolParameterField } from '../domain/tools';
import { announce } from './announcer';
import { formatBytes } from './format';
import { makeblastdbFormModel, MK_REFUSAL_ID } from './form-model';
import OptionName from './OptionName.vue';
import ParameterField from './ParameterField.vue';
import ParameterForm from './ParameterForm.vue';
import ToolRunBar from './ToolRunBar.vue';
import ToolRunFiles from './ToolRunFiles.vue';
import { useStore } from './useStore';

const props = defineProps<{
  draft: MakeblastdbDraft;
  coordinator: Coordinator;
}>();

const emit = defineEmits<{ submitted: [runId: string] }>();

const state = useStore(props.draft.state);
const formModel = computed(() => makeblastdbFormModel(props.draft, state.value));

const describeState = computed(() => {
  if (state.value.description !== undefined) return 'ready';
  return state.value.descriptionError === undefined ? 'loading' : 'failed';
});

const summary = computed(() => state.value.description?.description || 'Application to create BLAST databases');

/** A field of the Input and Output blocks, from `-xmlhelp`; `-dbtype` is a pair of radio buttons (it has no default). */
function blockField(flag: string): Field | undefined {
  const param = state.value.description?.parameters.find((entry) => entry.flag === flag);
  if (param === undefined) return undefined;
  const field = toToolParameterField(param);
  return flag === '-dbtype' ? { ...field, kind: 'radio' } : field;
}
const inputTypeField = computed(() => blockField('-input_type'));
const dbtypeField = computed(() => blockField('-dbtype'));
const outField = computed(() => blockField('-out'));

const inputType = computed(() => props.draft.inputType(state.value));

const dragging = ref(false);
const inputFilesInput = ref<HTMLInputElement | null>(null);
const dbFilesInput = ref<HTMLInputElement | null>(null);
const dbFolderInput = ref<HTMLInputElement | null>(null);

function picked(event: Event): File[] {
  const input = event.target as HTMLInputElement;
  const files = input.files !== null ? [...input.files] : [];
  input.value = '';
  return files;
}

function onInputFiles(event: Event): void {
  props.draft.addInputFiles(picked(event));
}

function onDbFiles(event: Event): void {
  const files = picked(event);
  if (files.length > 0) props.draft.addRunFiles(files);
}

function onDrop(event: DragEvent): void {
  dragging.value = false;
  const files = event.dataTransfer?.files;
  if (files === undefined || files.length === 0) return;
  if (inputType.value === 'blastdb') props.draft.addRunFiles([...files]);
  else props.draft.addInputFiles([...files]);
}

const commandWords = computed(() => props.draft.previewCommandWords(state.value));
const takesValue = computed(() => toolTakesValue(state.value.description));
const buttonLabel = computed(() => (state.value.submitting ? 'Starting…' : 'Run LOSAT'));

const changed = computed(() => changedMakeblastdbOptionsCount(state.value.values, state.value.extra, state.value.description));
const open = ref(changed.value > 0);
// An option set from outside the tab ("Use as -mask_data in makeblastdb") opens the options, where it is shown.
watch(changed, (count, before) => {
  if (count > before) open.value = true;
});
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
  <div class="search-panel wm-panel" data-testid="mk-panel" :data-describe="describeState">
    <h2 class="visually-hidden">makeblastdb</h2>

    <p class="program-summary" data-testid="mk-summary">{{ summary }}</p>

    <fieldset
      class="search-block input-panel"
      :class="{ dragging }"
      data-testid="mk-input-block"
      @dragover.prevent="dragging = true"
      @dragleave="dragging = false"
      @drop.prevent="onDrop"
    >
      <legend>Input</legend>
      <ParameterField v-if="inputTypeField" :model="formModel" :field="inputTypeField" />

      <template v-if="inputType === 'blastdb'">
        <div class="param-row">
          <label class="param-label" for="mk-db-in">
            Database name <OptionName class="flag" name="-in" />
          </label>
          <div class="param-control wm-in-row">
            <input
              id="mk-db-in"
              type="text"
              autocomplete="off"
              spellcheck="false"
              placeholder="db/mydb"
              :value="state.dbIn"
              data-testid="mk-db-in"
              @input="draft.setDbIn(($event.target as HTMLInputElement).value)"
            />
            <button type="button" data-testid="mk-db-open-files" @click="dbFilesInput?.click()">Open files…</button>
            <button type="button" data-testid="mk-db-open-folder" @click="dbFolderInput?.click()">Open folder…</button>
            <input
              ref="dbFilesInput"
              class="visually-hidden"
              type="file"
              multiple
              tabindex="-1"
              aria-hidden="true"
              aria-label="Open database files"
              data-testid="mk-db-files"
              @change="onDbFiles"
            />
            <input
              ref="dbFolderInput"
              class="visually-hidden"
              type="file"
              multiple
              webkitdirectory
              tabindex="-1"
              aria-hidden="true"
              aria-label="Open database folder"
              data-testid="mk-db-folder"
              @change="onDbFiles"
            />
          </div>
        </div>
        <p v-if="state.runFiles.length > 0" class="muted small wm-in-count" data-testid="mk-db-files-count">
          {{ state.runFiles.length }} {{ state.runFiles.length === 1 ? 'file' : 'files' }} in the run directory.
        </p>
        <p v-if="state.runFilesError" class="error small">{{ state.runFilesError }}</p>
      </template>

      <template v-else>
        <div v-if="inputType === 'fasta'" class="param-row">
          <label class="param-label" for="mk-fasta-input">
            FASTA sequences <OptionName class="flag" name="-in" />
          </label>
          <div class="param-control fill">
            <textarea
              id="mk-fasta-input"
              rows="5"
              spellcheck="false"
              placeholder="Enter FASTA sequences"
              :value="state.fastaPaste"
              data-testid="mk-fasta-input"
              @input="draft.setFastaPaste(($event.target as HTMLTextAreaElement).value)"
            />
            <div v-if="state.fastaPaste !== ''" class="wm-paste-tools">
              <button type="button" class="link" data-testid="mk-fasta-clear" @click="draft.setFastaPaste('')">Clear</button>
            </div>
          </div>
        </div>

        <div class="param-row">
          <span class="param-label">
            <template v-if="inputType === 'fasta'">Or, upload files</template>
            <template v-else>Input files <OptionName class="flag" name="-in" /></template>
          </span>
          <div class="param-control dropzone">
            <button type="button" data-testid="mk-open-input" @click="inputFilesInput?.click()">
              {{ inputType === 'fasta' ? 'Open FASTA files…' : 'Open files…' }}
            </button>
            <input
              ref="inputFilesInput"
              class="visually-hidden"
              type="file"
              multiple
              tabindex="-1"
              aria-hidden="true"
              aria-label="Open the files of -in"
              data-testid="mk-input-files"
              @change="onInputFiles"
            />
            <span class="muted">{{ inputType === 'fasta' ? 'Drop FASTA files here (.fa, .fna, .faa, .gz, .bz2, .zst)' : 'Drop ASN.1 files here' }}</span>
          </div>
        </div>

        <ul v-if="state.inputFiles.length > 0" class="sources" data-testid="mk-input-list">
          <li v-for="(entry, index) in state.inputFiles" :key="entry.key" class="source">
            <div class="source-head">
              <strong class="source-name" :data-testid="`mk-input-name-${index}`">{{ entry.name }}</strong>
              <span class="muted">{{ formatBytes(entry.size) }}</span>
              <button type="button" class="link" :data-testid="`mk-input-remove-${index}`" @click="draft.removeInputFile(entry.key)">
                Remove
              </button>
            </div>
          </li>
        </ul>
      </template>

      <ParameterField v-if="dbtypeField" :model="formModel" :field="dbtypeField" />
    </fieldset>

    <fieldset class="search-block" data-testid="mk-output-block">
      <legend>Output</legend>
      <ParameterField v-if="outField" :model="formModel" :field="outField" />
    </fieldset>

    <ToolRunBar
      :submitting="state.submitting"
      :button-label="buttonLabel"
      button-test-id="mk-run-losat"
      :command-words="commandWords"
      :takes-value="takesValue"
      command-test-id="mk-command"
      @run="submit('top')"
    />

    <div class="blast-notes">
      <p class="validation" :data-state="state.validation.state" data-testid="mk-argv-validation" aria-live="polite">
        <span v-if="state.validation.state === 'invalid'" :id="MK_REFUSAL_ID" class="error">{{ toolValidationSummary(state.validation.message) }}</span>
      </p>
      <p v-if="state.message && !messageBelow" class="error" data-testid="mk-message">{{ state.message.text }}</p>
    </div>

    <section class="algorithm-parameters" data-testid="mk-options">
      <h3 class="parameters-heading">
        <button
          type="button"
          class="parameters-bar"
          :aria-expanded="open"
          aria-controls="mk-options-body"
          data-testid="mk-options-toggle"
          @click="open = !open"
        >
          <span class="bar-mark" aria-hidden="true">{{ open ? '−' : '+' }}</span>
          makeblastdb options
          <span v-if="changed > 0" class="bar-count" data-testid="mk-options-changed">({{ changed }} changed)</span>
        </button>
      </h3>
      <div v-if="open" id="mk-options-body" class="parameters-body">
        <div v-if="!state.descriptionUnavailable" class="parameters-tools">
          <p class="parameters-note">
            Parameter values that differ from the default are highlighted in yellow and marked with
            <span aria-hidden="true">♦</span><span class="visually-hidden">a diamond</span> sign
          </p>
          <button type="button" data-testid="mk-restore-defaults" @click="draft.resetOptions()">
            Restore default options
          </button>
        </div>
        <ParameterForm :model="formModel" test-id="mk-parameter-form">
          <template #before-extra>
            <div class="search-block parameter-section" data-testid="mk-run-files">
              <div class="section-legend">Files of the run</div>
              <ToolRunFiles :draft="draft" :files="state.runFiles" :error="state.runFilesError" prefix="mk"
                >Files opened by name (<OptionName name="-mask_data" />, <OptionName name="-taxid_map" />, <OptionName name="-conffile" />, BLAST database files)</ToolRunFiles
              >
            </div>
          </template>
        </ParameterForm>
        <ToolRunBar
          :submitting="state.submitting"
          :button-label="buttonLabel"
          button-test-id="mk-run-losat-bottom"
          :command-words="commandWords"
          :takes-value="takesValue"
          @run="submit('bottom')"
        />
        <p v-if="state.message && messageBelow" class="error" data-testid="mk-message-bottom">{{ state.message.text }}</p>
      </div>
    </section>
  </div>
</template>
