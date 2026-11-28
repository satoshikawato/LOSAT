<script setup lang="ts">
// The `-in` of a tool tab's Input block, one component for WindowMasker and Masking (WW4c): FASTA
// paste, Open FASTA files… (and the files dropped on the block, `.gz`, `.bz2`, `.zst` inflated), the
// list of the files with Remove, and Combined (`combined_input.fa`) or Separate runs for several
// inputs; or, for an input that is not FASTA (`-infmt blastdb`, `-fa_list`), the editable `-in` name
// with Open files… / Open folder… (added into the run directory). The prefix names the ids and test ids.
import { computed, ref } from 'vue';
import type { FastaInputForm, FastaInputState } from '../application/fasta-tool-draft';
import { formatBytes } from './format';
import LosatLoader from './LosatLoader.vue';
import OptionName from './OptionName.vue';

const props = defineProps<{
  draft: FastaInputForm;
  state: FastaInputState;
  /** The prefix of the ids and test ids (`wm`, `ms`). */
  prefix: string;
}>();

const fastaInput = ref<HTMLInputElement | null>(null);
const inFilesInput = ref<HTMLInputElement | null>(null);
const inFolderInput = ref<HTMLInputElement | null>(null);

const usesFasta = computed(() => props.draft.usesFastaBox(props.state));
const fastaCount = computed(() => (props.state.fastaPaste.trim() !== '' ? 1 : 0) + props.state.fastaFiles.length);

function readPickedFiles(event: Event): File[] {
  const input = event.target as HTMLInputElement;
  const files = input.files !== null ? [...input.files] : [];
  input.value = '';
  return files;
}

function onFastaFiles(event: Event): void {
  const files = readPickedFiles(event);
  if (files.length > 0) props.draft.addFastaFiles(files);
}

function onInFiles(event: Event): void {
  const files = readPickedFiles(event);
  if (files.length > 0) props.draft.addRunFiles(files);
}
</script>

<template>
  <template v-if="usesFasta">
    <div class="param-row">
      <label class="param-label" :for="`${prefix}-fasta-input`">
        FASTA sequences <OptionName class="flag" name="-in" />
      </label>
      <div class="param-control fill">
        <textarea
          :id="`${prefix}-fasta-input`"
          rows="5"
          spellcheck="false"
          placeholder="Enter FASTA sequences"
          :value="state.fastaPaste"
          :data-testid="`${prefix}-fasta-input`"
          @input="draft.setFastaPaste(($event.target as HTMLTextAreaElement).value)"
        />
        <div v-if="state.fastaPaste !== ''" class="wm-paste-tools">
          <button type="button" class="link" :data-testid="`${prefix}-fasta-clear`" @click="draft.setFastaPaste('')">Clear</button>
        </div>
      </div>
    </div>

    <div class="param-row">
      <span class="param-label">Or, upload files</span>
      <div class="param-control dropzone">
        <button type="button" :data-testid="`${prefix}-open-fasta`" @click="fastaInput?.click()">Open FASTA files…</button>
        <input
          ref="fastaInput"
          class="visually-hidden"
          type="file"
          multiple
          tabindex="-1"
          aria-hidden="true"
          aria-label="Open FASTA files"
          :data-testid="`${prefix}-fasta-files`"
          @change="onFastaFiles"
        />
        <span class="muted">Drop FASTA files here (.fa, .fna, .gz, .bz2, .zst)</span>
      </div>
    </div>

    <ul v-if="state.fastaFiles.length > 0" class="sources" :data-testid="`${prefix}-fasta-list`">
      <li v-for="(entry, index) in state.fastaFiles" :key="entry.key" class="source" :data-status="entry.status">
        <div class="source-head">
          <strong class="source-name" :data-testid="`${prefix}-fasta-name-${index}`">{{ entry.name }}</strong>
          <span class="muted">{{ formatBytes(entry.size) }}</span>
          <span v-if="entry.inflatedFrom" class="muted" :data-testid="`${prefix}-fasta-inflated-${index}`">from {{ entry.inflatedFrom }}</span>
          <span v-if="entry.status === 'inflating'" class="muted"><LosatLoader />Reading…</span>
          <button type="button" class="link" :data-testid="`${prefix}-fasta-remove-${index}`" @click="draft.removeFastaFile(entry.key)">
            Remove
          </button>
        </div>
        <p v-if="entry.status === 'failed' && entry.error" class="error">{{ entry.error }}</p>
      </li>
    </ul>

    <fieldset v-if="fastaCount >= 2" class="mode" :data-testid="`${prefix}-fasta-mode`">
      <legend>Multiple FASTA inputs</legend>
      <label>
        <input
          type="radio"
          :name="`${prefix}-fasta-mode`"
          value="combined"
          :checked="state.fastaMode === 'combined'"
          :data-testid="`${prefix}-fasta-mode-combined`"
          @change="draft.setFastaMode('combined')"
        />
        Combined: one run of all (<code>combined_input.fa</code>)
      </label>
      <label>
        <input
          type="radio"
          :name="`${prefix}-fasta-mode`"
          value="separate"
          :checked="state.fastaMode === 'separate'"
          :data-testid="`${prefix}-fasta-mode-separate`"
          @change="draft.setFastaMode('separate')"
        />
        Separate: one run for each
      </label>
    </fieldset>
  </template>

  <template v-else>
    <div class="param-row">
      <label class="param-label" :for="`${prefix}-in-name`">
        Input name <OptionName class="flag" name="-in" />
      </label>
      <div class="param-control wm-in-row">
        <input
          :id="`${prefix}-in-name`"
          type="text"
          autocomplete="off"
          spellcheck="false"
          placeholder="db/wm"
          :value="state.blastDbIn"
          :data-testid="`${prefix}-in-name`"
          @input="draft.setBlastDbIn(($event.target as HTMLInputElement).value)"
        />
        <button type="button" :data-testid="`${prefix}-in-open-files`" @click="inFilesInput?.click()">Open files…</button>
        <button type="button" :data-testid="`${prefix}-in-open-folder`" @click="inFolderInput?.click()">Open folder…</button>
        <input
          ref="inFilesInput"
          class="visually-hidden"
          type="file"
          multiple
          tabindex="-1"
          aria-hidden="true"
          aria-label="Open input files"
          :data-testid="`${prefix}-in-files`"
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
          :data-testid="`${prefix}-in-folder`"
          @change="onInFiles"
        />
      </div>
    </div>
    <p v-if="state.runFiles.length > 0" class="muted small wm-in-count" :data-testid="`${prefix}-in-files-count`">
      {{ state.runFiles.length }} {{ state.runFiles.length === 1 ? 'file' : 'files' }} in the run directory.
    </p>
    <p v-if="state.runFilesError" class="error small">{{ state.runFilesError }}</p>
  </template>
</template>
