<script setup lang="ts">
// Shared counts-file picker row used by WindowMasker Convert (-in) and Mask (-ustat).
import { ref } from 'vue';
import type { SessionCountsEntry } from '../application/coordinator';
import type { ToolCountsChoice } from '../application/tool-draft';
import { COUNTS_HOST, countsName } from '../domain/counts';
import type { SavedCountsEntry } from '../ports/counts';
import { formatBytes } from './format';

const props = defineProps<{
  prefix: 'convert' | 'ustat';
  source: ToolCountsChoice;
  sessionCounts: readonly SessionCountsEntry[];
  savedCounts: readonly SavedCountsEntry[];
  fileAriaLabel: string;
}>();

const emit = defineEmits<{
  'select-session': [entry: SessionCountsEntry];
  'select-saved': [entry: SavedCountsEntry];
  'choose-file': [file: File];
  clear: [];
}>();

const fileInput = ref<HTMLInputElement | null>(null);

function selectValue(choice: ToolCountsChoice): string {
  if (choice.kind === 'session') return `session:${choice.runId}:${choice.key}`;
  if (choice.kind === 'saved') return `saved:${choice.taxid}`;
  return '';
}

function onSelectChange(value: string): void {
  if (value === '') {
    emit('clear');
    return;
  }
  const sessionMatch = props.sessionCounts.find((entry) => `session:${entry.runId}:${entry.key}` === value);
  if (sessionMatch !== undefined) {
    emit('select-session', sessionMatch);
    return;
  }
  const savedMatch = props.savedCounts.find((entry) => `saved:${entry.taxid}` === value);
  if (savedMatch !== undefined) emit('select-saved', savedMatch);
}

function onFileChange(event: Event): void {
  const target = event.target as HTMLInputElement;
  const file = target.files?.[0];
  if (file !== undefined) emit('choose-file', file);
  target.value = '';
}
</script>

<template>
  <div class="param-control wm-counts-pick">
    <span v-if="sessionCounts.length > 0 || savedCounts.length > 0" class="select-clip">
      <select
        :id="`wm-${prefix}-select`"
        :value="selectValue(source)"
        :data-testid="`wm-${prefix}-select`"
        @change="onSelectChange(($event.target as HTMLSelectElement).value)"
      >
        <option value="">Choose a counts file…</option>
        <option v-for="entry in sessionCounts" :key="`session:${entry.runId}:${entry.key}`" :value="`session:${entry.runId}:${entry.key}`">
          Run {{ entry.runNumber }}: {{ entry.name }} ({{ formatBytes(entry.size) }})
        </option>
        <option v-for="entry in savedCounts" :key="`saved:${entry.taxid}`" :value="`saved:${entry.taxid}`">
          {{ countsName(entry.taxid) }} ({{ entry.origin === 'made' ? 'made in this browser' : COUNTS_HOST }}, {{ entry.date }}, {{ formatBytes(entry.size) }})
        </option>
      </select>
    </span>
    <button type="button" :data-testid="`wm-${prefix}-choose`" @click="fileInput?.click()">Choose file…</button>
    <input
      ref="fileInput"
      class="visually-hidden"
      type="file"
      tabindex="-1"
      aria-hidden="true"
      :aria-label="fileAriaLabel"
      :data-testid="`wm-${prefix}-file`"
      @change="onFileChange"
    />
    <span v-if="source.kind !== 'none'" class="file-name" :data-testid="`wm-${prefix}-chosen`">
      <span>
        {{ source.name }} ({{ formatBytes(source.size) }})<template v-if="source.kind === 'session'">, from run {{ source.runNumber }}</template><template v-else-if="source.kind === 'saved'">, saved for taxonomy ID {{ source.taxid }}</template>
      </span>
      <button type="button" class="link" :data-testid="`wm-${prefix}-clear`" @click="emit('clear')">Clear</button>
    </span>
  </div>
</template>
