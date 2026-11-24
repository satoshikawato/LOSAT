<script setup lang="ts">
// Shared counts-file picker row used by WindowMasker Convert (-in) and Mask (-ustat).
import { ref } from 'vue';
import type { SessionCountsEntry } from '../application/coordinator';
import type { ToolCountsChoice } from '../application/tool-draft';
import { formatBytes } from './format';

const props = defineProps<{
  prefix: 'convert' | 'ustat';
  source: ToolCountsChoice;
  sessionCounts: readonly SessionCountsEntry[];
  fileAriaLabel: string;
}>();

const emit = defineEmits<{
  'select-session': [entry: SessionCountsEntry];
  'choose-file': [file: File];
  clear: [];
}>();

const fileInput = ref<HTMLInputElement | null>(null);

function selectValue(choice: ToolCountsChoice): string {
  return choice.kind === 'session' ? `session:${choice.runId}:${choice.key}` : '';
}

function onSelectChange(value: string): void {
  if (value === '') {
    emit('clear');
    return;
  }
  const match = props.sessionCounts.find((entry) => `session:${entry.runId}:${entry.key}` === value);
  if (match !== undefined) emit('select-session', match);
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
    <span v-if="sessionCounts.length > 0" class="select-clip">
      <select
        :id="`wm-${prefix}-select`"
        :value="selectValue(source)"
        :data-testid="`wm-${prefix}-select`"
        @change="onSelectChange(($event.target as HTMLSelectElement).value)"
      >
        <option value="">Choose a counts file from this session…</option>
        <option v-for="entry in sessionCounts" :key="`${entry.runId}:${entry.key}`" :value="`session:${entry.runId}:${entry.key}`">
          Run {{ entry.runNumber }}: {{ entry.name }} ({{ formatBytes(entry.size) }})
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
      <span>{{ source.name }} ({{ formatBytes(source.size) }})<template v-if="source.kind === 'session'">, from run {{ source.runNumber }}</template></span>
      <button type="button" class="link" :data-testid="`wm-${prefix}-clear`" @click="emit('clear')">Clear</button>
    </span>
  </div>
</template>
