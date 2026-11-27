<script setup lang="ts">
// The Format Archive (`blast_formatter`) main tab (RG-F3p7):
// Re-formats `-outfmt 11` BLAST archives from a local file, pasted ASN.1 text, or a completed
// `-outfmt 11` search output of this session, with optional run directory files (local BLAST
// database folder for `-db` archives, or `-conffile`), `-outfmt` / `-out` controls, live command
// preview via ToolRunBar, in-place validateTool refusals, and collapsible "Formatting options"
// driven by `-xmlhelp`.
import { computed, nextTick, ref, watch } from 'vue';
import type { ArchiveSourceKind, BlastFormatterDraft } from '../application/blast-formatter-draft';
import type { Coordinator } from '../application/coordinator';
import { changedBlastFormatterOptionsCount, toolValidationSummary } from '../domain/tool-argv';
import {
  BLAST_FORMATTER_OUTFMT_CHOICES,
  toolParameterHelpText,
  toolTakesValue,
} from '../domain/tools';
import { announce } from './announcer';
import { formatBytes, optionHelp } from './format';
import { BF_REFUSAL_ID, blastFormatterFormModel } from './form-model';
import HelpTip from './HelpTip.vue';
import OptionName from './OptionName.vue';
import ParameterForm from './ParameterForm.vue';
import ToolRunBar from './ToolRunBar.vue';
import { useStore } from './useStore';

const props = defineProps<{
  draft: BlastFormatterDraft;
  coordinator: Coordinator;
}>();

const emit = defineEmits<{ submitted: [runId: string] }>();

const state = useStore(props.draft.state);
const appState = useStore(props.coordinator.state);
const formModel = computed(() => blastFormatterFormModel(props.draft, state.value));

const archiveFileInput = ref<HTMLInputElement | null>(null);
const runFilesInput = ref<HTMLInputElement | null>(null);
const runFolderInput = ref<HTMLInputElement | null>(null);
const customOutfmtInput = ref<HTMLInputElement | null>(null);
const dragActive = ref(false);

const OTHER_OUTFMT = 'other';
const presetValues = new Set(BLAST_FORMATTER_OUTFMT_CHOICES.map((choice) => choice.value));
const customOutfmt = ref(!presetValues.has(state.value.outfmt.trim()));

watch(
  () => state.value.outfmt,
  (val) => {
    const trimmed = val.trim();
    if (trimmed !== '' && !presetValues.has(trimmed)) customOutfmt.value = true;
  },
);

const outfmtSelectValue = computed(() => {
  if (customOutfmt.value) return OTHER_OUTFMT;
  const trimmed = state.value.outfmt.trim();
  return presetValues.has(trimmed) ? trimmed : OTHER_OUTFMT;
});

function onOutfmtSelect(next: string): void {
  if (next === OTHER_OUTFMT) {
    customOutfmt.value = true;
    void nextTick(() => customOutfmtInput.value?.focus());
    return;
  }
  customOutfmt.value = false;
  props.draft.setOutfmt(next);
}

const describeState = computed(() => {
  if (state.value.description !== undefined) return 'ready';
  return state.value.descriptionError === undefined ? 'loading' : 'failed';
});

const summary = computed(
  () =>
    state.value.description?.description ||
    'Re-format a BLAST archive (-outfmt 11) into another output format without searching again',
);

const sessionArchives = computed(() => {
  void appState.value.runs;
  return props.coordinator.sessionOutfmt11Archives();
});

watch(
  [() => state.value.archiveSource, sessionArchives],
  ([source, list]) => {
    if (source === 'session' && state.value.sessionArchive === undefined && list.length > 0) {
      props.draft.selectSessionArchive(list[0]);
    }
  },
  { immediate: true },
);

const paramHelp = (flag: string): string => {
  const param = state.value.description?.parameters.find((p) => p.flag === flag);
  return param === undefined ? '' : optionHelp(toolParameterHelpText(param));
};
const archiveHelp = computed(() => paramHelp('-archive'));
const outfmtHelp = computed(() => paramHelp('-outfmt'));
const outHelp = computed(() => paramHelp('-out'));

const commandWords = computed(() => props.draft.previewCommandWords(state.value));
const takesValue = computed(() => toolTakesValue(state.value.description));
const buttonLabel = computed(() => (state.value.submitting ? 'Starting…' : 'Run LOSAT'));

const changed = computed(() =>
  changedBlastFormatterOptionsCount(state.value.values, state.value.extra, state.value.description),
);
const open = ref(changed.value > 0);
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

const ARCHIVE_SOURCES: readonly { readonly id: ArchiveSourceKind; readonly label: string }[] = [
  { id: 'file', label: 'Archive file' },
  { id: 'paste', label: 'Paste archive' },
  { id: 'session', label: 'Session outfmt 11' },
];

function onArchiveFile(event: Event): void {
  const input = event.target as HTMLInputElement;
  const file = input.files?.[0];
  if (file !== undefined) props.draft.setArchiveFile(file);
  input.value = '';
}

function onArchiveDrop(event: DragEvent): void {
  dragActive.value = false;
  const file = event.dataTransfer?.files?.[0];
  if (file !== undefined) props.draft.setArchiveFile(file);
}

function onSessionChange(event: Event): void {
  const runId = (event.target as HTMLSelectElement).value;
  const found = sessionArchives.value.find((entry) => entry.runId === runId);
  props.draft.selectSessionArchive(found);
}

function onRunFiles(event: Event): void {
  const input = event.target as HTMLInputElement;
  if (input.files !== null && input.files.length > 0) props.draft.addRunFiles([...input.files]);
  input.value = '';
}
</script>

<template>
  <div class="search-panel wm-panel" data-testid="bf-panel" :data-describe="describeState">
    <h2 class="visually-hidden">Format Archive</h2>

    <p class="program-summary" data-testid="bf-summary">{{ summary }}</p>

    <fieldset class="search-block" data-testid="bf-archive-block">
      <legend>BLAST archive (<OptionName name="-archive" />)</legend>

      <div class="param-row">
        <span class="param-label">Source</span>
        <div class="param-control">
          <div class="control-tip">
            <div class="bf-archive-mode" role="radiogroup" aria-label="Archive source">
              <label v-for="src in ARCHIVE_SOURCES" :key="src.id" class="bf-radio-option">
                <input
                  type="radio"
                  name="bf-archive-source"
                  :value="src.id"
                  :checked="state.archiveSource === src.id"
                  :data-testid="`bf-source-${src.id}`"
                  @change="draft.setArchiveSource(src.id)"
                />
                {{ src.label }}
                <span v-if="src.id === 'session' && sessionArchives.length > 0" class="muted">({{ sessionArchives.length }})</span>
              </label>
            </div>
            <HelpTip v-if="archiveHelp" id="bf-archive-help" :text="archiveHelp" label="Help for -archive" />
          </div>
        </div>
      </div>

      <div v-if="state.archiveSource === 'file'" class="param-row">
        <span class="param-label">File</span>
        <div class="param-control">
          <div
            class="bf-dropzone"
            :class="{ active: dragActive }"
            data-testid="bf-archive-dropzone"
            @dragover.prevent="dragActive = true"
            @dragleave="dragActive = false"
            @drop.prevent="onArchiveDrop"
          >
            <button type="button" data-testid="bf-choose-archive" @click="archiveFileInput?.click()">
              Choose archive file…
            </button>
            <input
              ref="archiveFileInput"
              class="visually-hidden"
              type="file"
              accept=".asn,.txt,text/plain"
              tabindex="-1"
              aria-hidden="true"
              aria-label="Choose BLAST archive file"
              data-testid="bf-archive-file-input"
              @change="onArchiveFile"
            />
            <span v-if="state.archiveFile" class="bf-chosen-file" data-testid="bf-archive-chosen">
              <code>{{ state.archiveFile.name }}</code>
              <span class="muted">({{ formatBytes(state.archiveFile.size) }})</span>
              <button type="button" class="link" data-testid="bf-archive-clear" @click="draft.setArchiveFile(undefined)">
                Remove
              </button>
            </span>
            <span v-else class="param-hint">Choose or drop a BLAST archive file written with <OptionName name="-outfmt" /> 11</span>
          </div>
        </div>
      </div>

      <div v-else-if="state.archiveSource === 'paste'" class="param-row">
        <label class="param-label" for="bf-archive-paste">ASN.1 text</label>
        <div class="param-control">
          <textarea
            id="bf-archive-paste"
            class="bf-paste-area"
            rows="6"
            spellcheck="false"
            placeholder="Blast4-archive ::= { ... }"
            :value="state.archivePaste"
            data-testid="bf-archive-paste"
            @input="draft.setArchivePaste(($event.target as HTMLTextAreaElement).value)"
          />
        </div>
      </div>

      <div v-else class="param-row">
        <label class="param-label" for="bf-session-select">Search run</label>
        <div class="param-control">
          <select
            v-if="sessionArchives.length > 0"
            id="bf-session-select"
            :value="state.sessionArchive?.runId ?? ''"
            data-testid="bf-session-select"
            @change="onSessionChange"
          >
            <option v-for="entry in sessionArchives" :key="entry.runId" :value="entry.runId">
              Run {{ entry.runNumber }} · {{ entry.program }} ({{ entry.name }})
            </option>
          </select>
          <span v-else class="param-hint" data-testid="bf-session-empty">
            No completed <OptionName name="-outfmt" /> 11 searches in this session yet. Run a search with
            <code>-outfmt 11</code> in Additional options, or choose an archive file.
          </span>
        </div>
      </div>

      <div class="param-row" data-testid="bf-run-files">
        <span class="param-label">Database / run files</span>
        <div class="param-control wm-run-files-controls">
          <button type="button" data-testid="bf-open-run-files" @click="runFilesInput?.click()">Open files…</button>
          <button type="button" data-testid="bf-open-run-folder" @click="runFolderInput?.click()">Open folder…</button>
          <input
            ref="runFilesInput"
            class="visually-hidden"
            type="file"
            multiple
            tabindex="-1"
            aria-hidden="true"
            aria-label="Open database or run files"
            data-testid="bf-run-files-input"
            @change="onRunFiles"
          />
          <input
            ref="runFolderInput"
            class="visually-hidden"
            type="file"
            multiple
            webkitdirectory
            tabindex="-1"
            aria-hidden="true"
            aria-label="Open database folder for the run directory"
            data-testid="bf-run-folder-input"
            @change="onRunFiles"
          />
          <span class="param-hint">Optional: local BLAST database files referenced by a <OptionName name="-db" /> archive, or <OptionName name="-conffile" /></span>
        </div>
      </div>
      <p v-if="state.runFilesError" class="error small" role="alert" data-testid="bf-run-files-error">{{ state.runFilesError }}</p>
      <ul v-if="state.runFiles.length > 0" class="wm-file-list" data-testid="bf-run-files-list">
        <li v-for="entry in state.runFiles" :key="entry.key" class="wm-file-item">
          <code>{{ entry.name }}</code>
          <span class="muted">({{ formatBytes(entry.size) }})</span>
          <button type="button" class="link" @click="draft.removeRunFile(entry.key)">Remove</button>
        </li>
      </ul>
    </fieldset>

    <fieldset class="search-block" data-testid="bf-output-block">
      <legend>Output</legend>
      <div class="param-row">
        <label class="param-label" for="bf-outfmt">
          Output format <OptionName class="flag" name="-outfmt" />
        </label>
        <div class="param-control">
          <div v-if="!customOutfmt" class="control-tip bf-outfmt-row">
            <select
              id="bf-outfmt"
              :value="outfmtSelectValue"
              data-testid="bf-outfmt-select"
              @change="onOutfmtSelect(($event.target as HTMLSelectElement).value)"
            >
              <option v-for="choice in BLAST_FORMATTER_OUTFMT_CHOICES" :key="choice.value" :value="choice.value">
                {{ choice.label }}
              </option>
              <option :value="OTHER_OUTFMT">Other / custom format…</option>
            </select>
            <HelpTip v-if="outfmtHelp" id="bf-outfmt-help" :text="outfmtHelp" label="Help for -outfmt" />
          </div>
          <div v-else class="bf-outfmt-custom-wrap">
            <select
              id="bf-outfmt"
              :value="outfmtSelectValue"
              data-testid="bf-outfmt-select"
              @change="onOutfmtSelect(($event.target as HTMLSelectElement).value)"
            >
              <option v-for="choice in BLAST_FORMATTER_OUTFMT_CHOICES" :key="choice.value" :value="choice.value">
                {{ choice.label }}
              </option>
              <option :value="OTHER_OUTFMT">Other / custom format…</option>
            </select>
            <div class="control-tip bf-outfmt-row">
              <input
                ref="customOutfmtInput"
                type="text"
                autocomplete="off"
                spellcheck="false"
                placeholder="e.g. 6 qseqid sseqid pident"
                aria-label="Custom -outfmt value"
                :value="state.outfmt"
                data-testid="bf-outfmt-custom"
                @input="draft.setOutfmt(($event.target as HTMLInputElement).value)"
              />
              <HelpTip v-if="outfmtHelp" id="bf-outfmt-help" :text="outfmtHelp" label="Help for -outfmt" />
            </div>
          </div>
        </div>
      </div>
      <div class="param-row">
        <label class="param-label" for="bf-out">
          <OptionName class="flag" name="-out" />
        </label>
        <div class="param-control">
          <input
            id="bf-out"
            type="text"
            autocomplete="off"
            spellcheck="false"
            placeholder="- (standard output)"
            :value="state.outName"
            data-testid="bf-out"
            @input="draft.setOutName(($event.target as HTMLInputElement).value)"
          />
          <span class="hint-line">
            <span class="param-hint">Output file name (leave blank or type <code>-</code> for standard output)</span>
            <HelpTip v-if="outHelp" id="bf-out-help" :text="outHelp" label="Help for -out" />
          </span>
        </div>
      </div>
    </fieldset>

    <ToolRunBar
      :submitting="state.submitting"
      :button-label="buttonLabel"
      button-test-id="bf-run-losat"
      :command-words="commandWords"
      :takes-value="takesValue"
      command-test-id="bf-command"
      @run="submit('top')"
    />

    <div class="blast-notes">
      <p class="validation" :data-state="state.validation.state" data-testid="bf-argv-validation" aria-live="polite">
        <span v-if="state.validation.state === 'invalid'" :id="BF_REFUSAL_ID" class="error">{{ toolValidationSummary(state.validation.message) }}</span>
      </p>
      <p v-if="state.message && !messageBelow" class="error" data-testid="bf-message">{{ state.message.text }}</p>
    </div>

    <section class="algorithm-parameters" data-testid="bf-options">
      <h3 class="parameters-heading">
        <button
          type="button"
          class="parameters-bar"
          :aria-expanded="open"
          aria-controls="bf-options-body"
          data-testid="bf-options-toggle"
          @click="open = !open"
        >
          <span class="bar-mark" aria-hidden="true">{{ open ? '−' : '+' }}</span>
          Formatting options
          <span v-if="changed > 0" class="bar-count" data-testid="bf-options-changed">({{ changed }} changed)</span>
        </button>
      </h3>
      <div v-if="open" id="bf-options-body" class="parameters-body">
        <div v-if="!state.descriptionUnavailable" class="parameters-tools">
          <p class="parameters-note">
            Parameter values that differ from the default are highlighted in yellow and marked with
            <span aria-hidden="true">♦</span><span class="visually-hidden">a diamond</span> sign
          </p>
          <button type="button" data-testid="bf-restore-defaults" @click="draft.resetOptions()">
            Restore default options
          </button>
        </div>
        <ParameterForm :model="formModel" test-id="bf-parameter-form" />
        <ToolRunBar
          :submitting="state.submitting"
          :button-label="buttonLabel"
          button-test-id="bf-run-losat-bottom"
          :command-words="commandWords"
          :takes-value="takesValue"
          @run="submit('bottom')"
        />
        <p v-if="state.message && messageBelow" class="error" data-testid="bf-message-bottom">{{ state.message.text }}</p>
      </div>
    </section>
  </div>
</template>
