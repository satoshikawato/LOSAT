<script setup lang="ts">
// The files of a tool run's working directory ("Run directory"), one component for every tool tab
// (WindowMasker, Format Archive; RG-F3p7 review R1 #3): Open files…, Open folder…, the refusal of a
// second file under a name already there, and the list with Remove. The slot names what the tool
// opens by name.
import { ref } from 'vue';
import type { ToolExtraFile } from '../application/tool-draft-core';
import { formatBytes } from './format';

const props = defineProps<{
  draft: { addRunFiles(files: readonly File[]): void; removeRunFile(key: string): void };
  files: readonly ToolExtraFile[];
  error?: string | undefined;
  /** The prefix of the test ids (`wm`, `bf`). */
  prefix: string;
}>();

const filesInput = ref<HTMLInputElement | null>(null);
const folderInput = ref<HTMLInputElement | null>(null);

function onFiles(event: Event): void {
  const input = event.target as HTMLInputElement;
  if (input.files !== null && input.files.length > 0) props.draft.addRunFiles([...input.files]);
  input.value = '';
}
</script>

<template>
  <div class="param-row">
    <span class="param-label">Run directory</span>
    <div class="param-control wm-run-files-controls">
      <button type="button" :data-testid="`${prefix}-open-run-files`" @click="filesInput?.click()">Open files…</button>
      <button type="button" :data-testid="`${prefix}-open-run-folder`" @click="folderInput?.click()">Open folder…</button>
      <input
        ref="filesInput"
        class="visually-hidden"
        type="file"
        multiple
        tabindex="-1"
        aria-hidden="true"
        aria-label="Open files for the run directory"
        :data-testid="`${prefix}-run-files-input`"
        @change="onFiles"
      />
      <input
        ref="folderInput"
        class="visually-hidden"
        type="file"
        multiple
        webkitdirectory
        tabindex="-1"
        aria-hidden="true"
        aria-label="Open folder for the run directory"
        :data-testid="`${prefix}-run-folder-input`"
        @change="onFiles"
      />
      <span class="param-hint"><slot /></span>
    </div>
  </div>
  <p v-if="error" class="error small wm-run-files-error" role="alert" :data-testid="`${prefix}-run-files-error`">
    {{ error }}
  </p>
  <ul v-if="files.length > 0" class="wm-file-list" :data-testid="`${prefix}-run-files-list`">
    <li v-for="(entry, index) in files" :key="entry.key" class="wm-file-item">
      <span class="source-name" :data-testid="`${prefix}-run-file-name-${index}`">{{ entry.name }}</span>
      <span class="muted">{{ formatBytes(entry.size) }}</span>
      <button type="button" class="link" :data-testid="`${prefix}-run-file-remove-${index}`" @click="draft.removeRunFile(entry.key)">
        Remove
      </button>
    </li>
  </ul>
</template>
