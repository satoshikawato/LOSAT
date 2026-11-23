<script setup lang="ts">
// Input and Counts blocks of the WindowMasker tab (WW3b DESIGN §2.1):
// - Make counts & Mask: -infmt (fasta, blastdb, seqids), FASTA paste/files/dropzone with
//   Combined (combined_input.fa) or Separate runs, or an editable -in name with Open files /
//   Open folder (added into the run directory) for blastdb, seqids or -fa_list.
// - Convert counts: counts file picker (session-made counts or Choose file…).
// - Mask: Counts (-ustat) picker (session-made counts or Choose file…).
import { computed, ref } from 'vue';
import type { SessionCountsEntry } from '../application/coordinator';
import type { ToolCountsChoice, ToolDraft, ToolDraftState } from '../application/tool-draft';
import { stageBlockParameter, toToolParameterField } from '../domain/tools';
import { formatBytes } from './format';
import type { FormModel } from './form-model';
import LosatLoader from './LosatLoader.vue';
import OptionName from './OptionName.vue';
import ParameterField from './ParameterField.vue';

const props = defineProps<{
  draft: ToolDraft;
  state: ToolDraftState;
  sessionCounts: readonly SessionCountsEntry[];
  model: FormModel;
}>();

const dragging = ref(false);
const fastaInput = ref<HTMLInputElement | null>(null);
const inFilesInput = ref<HTMLInputElement | null>(null);
const inFolderInput = ref<HTMLInputElement | null>(null);
const ustatInput = ref<HTMLInputElement | null>(null);
const convertInput = ref<HTMLInputElement | null>(null);

const infmtField = computed(() => {
  const param = stageBlockParameter(props.state.description, props.state.stage, '-infmt');
  return param === undefined ? undefined : toToolParameterField(param);
});
const usesFasta = computed(() => props.draft.usesFastaBox(props.state));
const fastaCount = computed(() => (props.state.fastaPaste.trim() !== '' ? 1 : 0) + props.state.fastaFiles.length);

function onFastaFiles(event: Event): void {
  const input = event.target as HTMLInputElement;
  if (input.files !== null && input.files.length > 0) props.draft.addFastaFiles([...input.files]);
  input.value = '';
}

function onDrop(event: DragEvent): void {
  dragging.value = false;
  const files = event.dataTransfer?.files;
  if (files !== undefined && files.length > 0) {
    if (usesFasta.value) props.draft.addFastaFiles([...files]);
    else props.draft.addRunFiles([...files]);
  }
}

function onInFiles(event: Event): void {
  const input = event.target as HTMLInputElement;
  if (input.files !== null && input.files.length > 0) props.draft.addRunFiles([...input.files]);
  input.value = '';
}

function onUstatFile(event: Event): void {
  const input = event.target as HTMLInputElement;
  const file = input.files?.[0];
  if (file !== undefined) props.draft.setUstatFile(file);
  input.value = '';
}

function onConvertFile(event: Event): void {
  const input = event.target as HTMLInputElement;
  const file = input.files?.[0];
  if (file !== undefined) props.draft.setConvertFile(file);
  input.value = '';
}

function selectValue(choice: ToolCountsChoice): string {
  return choice.kind === 'session' ? `session:${choice.runId}:${choice.key}` : '';
}

function onSelectCounts(target: 'ustat' | 'convert', value: string): void {
  if (value === '') {
    if (target === 'ustat') props.draft.clearUstat();
    else props.draft.clearConvert();
    return;
  }
  const match = props.sessionCounts.find((entry) => `session:${entry.runId}:${entry.key}` === value);
  if (match === undefined) return;
  if (target === 'ustat') {
    props.draft.setUstatSession(match.runId, match.runNumber, match.key, match.name, match.size);
  } else {
    props.draft.setConvertSession(match.runId, match.runNumber, match.key, match.name, match.size);
  }
}
</script>

<template>
  <fieldset
    class="search-block input-panel"
    :class="{ dragging }"
    data-testid="wm-input-block"
    @dragover.prevent="dragging = true"
    @dragleave="dragging = false"
    @drop.prevent="onDrop"
  >
    <legend>Input</legend>

    <template v-if="state.stage === 'convert'">
      <div class="param-row">
        <label class="param-label" for="wm-convert-select">
          Counts file <OptionName class="flag" name="-in" />
        </label>
        <div class="param-control wm-counts-pick">
          <span v-if="sessionCounts.length > 0" class="select-clip">
            <select
              id="wm-convert-select"
              :value="selectValue(state.convertSource)"
              data-testid="wm-convert-select"
              @change="onSelectCounts('convert', ($event.target as HTMLSelectElement).value)"
            >
              <option value="">Choose a counts file from this session…</option>
              <option v-for="entry in sessionCounts" :key="`${entry.runId}:${entry.key}`" :value="`session:${entry.runId}:${entry.key}`">
                Run {{ entry.runNumber }}: {{ entry.name }} ({{ formatBytes(entry.size) }})
              </option>
            </select>
          </span>
          <button type="button" data-testid="wm-convert-choose" @click="convertInput?.click()">Choose file…</button>
          <input
            ref="convertInput"
            class="visually-hidden"
            type="file"
            tabindex="-1"
            aria-hidden="true"
            aria-label="Choose counts file to convert"
            data-testid="wm-convert-file"
            @change="onConvertFile"
          />
          <span v-if="state.convertSource.kind !== 'none'" class="file-name" data-testid="wm-convert-chosen">
            <span>{{ state.convertSource.name }} ({{ formatBytes(state.convertSource.size) }})<template v-if="state.convertSource.kind === 'session'">, from run {{ state.convertSource.runNumber }}</template></span>
            <button type="button" class="link" data-testid="wm-convert-clear" @click="draft.clearConvert()">Clear</button>
          </span>
        </div>
      </div>
    </template>

    <template v-else>
      <ParameterField v-if="infmtField" :model="model" :field="infmtField" />

      <template v-if="usesFasta">
        <div class="param-row">
          <label class="param-label" for="wm-fasta-input">
            FASTA sequences <OptionName class="flag" name="-in" />
          </label>
          <div class="param-control fill">
            <textarea
              id="wm-fasta-input"
              rows="5"
              spellcheck="false"
              placeholder="Enter FASTA sequences"
              :value="state.fastaPaste"
              data-testid="wm-fasta-input"
              @input="draft.setFastaPaste(($event.target as HTMLTextAreaElement).value)"
            />
            <div v-if="state.fastaPaste !== ''" class="wm-paste-tools">
              <button type="button" class="link" data-testid="wm-fasta-clear" @click="draft.setFastaPaste('')">Clear</button>
            </div>
          </div>
        </div>

        <div class="param-row">
          <span class="param-label">Or, upload files</span>
          <div class="param-control dropzone">
            <button type="button" data-testid="wm-open-fasta" @click="fastaInput?.click()">Open FASTA files…</button>
            <input
              ref="fastaInput"
              class="visually-hidden"
              type="file"
              multiple
              tabindex="-1"
              aria-hidden="true"
              aria-label="Open FASTA files"
              data-testid="wm-fasta-files"
              @change="onFastaFiles"
            />
            <span class="muted">Drop FASTA files here (.fa, .fna, .gz, .bz2, .zst)</span>
          </div>
        </div>

        <ul v-if="state.fastaFiles.length > 0" class="sources" data-testid="wm-fasta-list">
          <li v-for="(entry, index) in state.fastaFiles" :key="entry.key" class="source" :data-status="entry.status">
            <div class="source-head">
              <strong class="source-name" :data-testid="`wm-fasta-name-${index}`">{{ entry.name }}</strong>
              <span class="muted">{{ formatBytes(entry.size) }}</span>
              <span v-if="entry.inflatedFrom" class="muted" :data-testid="`wm-fasta-inflated-${index}`">from {{ entry.inflatedFrom }}</span>
              <span v-if="entry.status === 'inflating'" class="muted"><LosatLoader />Reading…</span>
              <button type="button" class="link" :data-testid="`wm-fasta-remove-${index}`" @click="draft.removeFastaFile(entry.key)">
                Remove
              </button>
            </div>
            <p v-if="entry.status === 'failed' && entry.error" class="error">{{ entry.error }}</p>
          </li>
        </ul>

        <fieldset v-if="fastaCount >= 2" class="mode" data-testid="wm-fasta-mode">
          <legend>Multiple FASTA inputs</legend>
          <label>
            <input
              type="radio"
              name="wm-fasta-mode"
              value="combined"
              :checked="state.fastaMode === 'combined'"
              data-testid="wm-fasta-mode-combined"
              @change="draft.setFastaMode('combined')"
            />
            Combined: one run of all (<code>combined_input.fa</code>)
          </label>
          <label>
            <input
              type="radio"
              name="wm-fasta-mode"
              value="separate"
              :checked="state.fastaMode === 'separate'"
              data-testid="wm-fasta-mode-separate"
              @change="draft.setFastaMode('separate')"
            />
            Separate: one run for each
          </label>
        </fieldset>
      </template>

      <template v-else>
        <div class="param-row">
          <label class="param-label" for="wm-in-name">
            Input name <OptionName class="flag" name="-in" />
          </label>
          <div class="param-control wm-in-row">
            <input
              id="wm-in-name"
              type="text"
              autocomplete="off"
              spellcheck="false"
              placeholder="db/wm"
              :value="state.blastDbIn"
              data-testid="wm-in-name"
              @input="draft.setBlastDbIn(($event.target as HTMLInputElement).value)"
            />
            <button type="button" data-testid="wm-in-open-files" @click="inFilesInput?.click()">Open files…</button>
            <button type="button" data-testid="wm-in-open-folder" @click="inFolderInput?.click()">Open folder…</button>
            <input
              ref="inFilesInput"
              class="visually-hidden"
              type="file"
              multiple
              tabindex="-1"
              aria-hidden="true"
              aria-label="Open input files"
              data-testid="wm-in-files"
              @change="onInFiles"
            />
            <input
              ref="inFolderInput"
              class="visually-hidden"
              type="file"
              multiple
              webkitdirectory
              tabindex="-1"
              aria-hidden="true"
              aria-label="Open input folder"
              data-testid="wm-in-folder"
              @change="onInFiles"
            />
          </div>
        </div>
        <p v-if="state.runFiles.length > 0" class="muted small wm-in-count" data-testid="wm-in-files-count">
          {{ state.runFiles.length }} {{ state.runFiles.length === 1 ? 'file' : 'files' }} in the run directory.
        </p>
        <p v-if="state.runFilesError" class="error small">{{ state.runFilesError }}</p>
      </template>
    </template>
  </fieldset>

  <fieldset v-if="state.stage === 'mask'" class="search-block" data-testid="wm-counts-block">
    <legend>Counts</legend>
    <div class="param-row">
      <label class="param-label" for="wm-ustat-select">
        Counts file <OptionName class="flag" name="-ustat" />
      </label>
      <div class="param-control wm-counts-pick">
        <span v-if="sessionCounts.length > 0" class="select-clip">
          <select
            id="wm-ustat-select"
            :value="selectValue(state.ustatSource)"
            data-testid="wm-ustat-select"
            @change="onSelectCounts('ustat', ($event.target as HTMLSelectElement).value)"
          >
            <option value="">Choose a counts file from this session…</option>
            <option v-for="entry in sessionCounts" :key="`${entry.runId}:${entry.key}`" :value="`session:${entry.runId}:${entry.key}`">
              Run {{ entry.runNumber }}: {{ entry.name }} ({{ formatBytes(entry.size) }})
            </option>
          </select>
        </span>
        <button type="button" data-testid="wm-ustat-choose" @click="ustatInput?.click()">Choose file…</button>
        <input
          ref="ustatInput"
          class="visually-hidden"
          type="file"
          tabindex="-1"
          aria-hidden="true"
          aria-label="Choose counts file for -ustat"
          data-testid="wm-ustat-file"
          @change="onUstatFile"
        />
        <span v-if="state.ustatSource.kind !== 'none'" class="file-name" data-testid="wm-ustat-chosen">
          <span>{{ state.ustatSource.name }} ({{ formatBytes(state.ustatSource.size) }})<template v-if="state.ustatSource.kind === 'session'">, from run {{ state.ustatSource.runNumber }}</template></span>
          <button type="button" class="link" data-testid="wm-ustat-clear" @click="draft.clearUstat()">Clear</button>
        </span>
      </div>
    </div>
  </fieldset>
</template>
