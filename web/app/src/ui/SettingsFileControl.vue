<script setup lang="ts">
// The search form's settings file (S15 instructions, item 3): "Save settings" writes the form's
// program, options and threads (never the inputs, names or title); "Load settings…" reads one
// back into the form. The line below says what was loaded, what was not applied, or why a file
// was refused; "Edit Search" of the results reports here too, since it fills this form.
import { ref } from 'vue';
import type { RunFiles } from '../application/run-files';
import { useStore } from './useStore';

const props = defineProps<{ runFiles: RunFiles; disabled: boolean }>();
const state = useStore(props.runFiles.state);
const root = ref<HTMLElement>();
const fileInput = ref<HTMLInputElement>();
const busy = ref(false);

/** Brings the control and its message into view ("Edit Search" opens the form at it). */
function show(): void {
  root.value?.scrollIntoView({ block: 'center' });
}
defineExpose({ show });

async function act(work: () => Promise<unknown>): Promise<void> {
  busy.value = true;
  try {
    await work();
  } finally {
    busy.value = false;
  }
}

function onFile(event: Event): void {
  const input = event.target as HTMLInputElement;
  const file = input.files?.[0];
  input.value = '';
  if (file !== undefined) void act(() => props.runFiles.loadSettings(file));
}
</script>

<template>
  <div ref="root" class="settings-file" data-testid="settings-file">
    <span class="settings-file-label">Search settings</span>
    <button type="button" data-testid="settings-save" :disabled="disabled || busy" @click="act(() => runFiles.saveSettings())">
      Save settings
    </button>
    <button type="button" data-testid="settings-load-button" :disabled="busy" @click="fileInput?.click()">Load settings…</button>
    <input
      ref="fileInput"
      class="visually-hidden"
      type="file"
      accept=".json,application/json"
      tabindex="-1"
      aria-label="Load a settings file"
      data-testid="settings-load"
      @change="onFile"
    />
    <span class="hint">The program, the options and the threads, without the inputs (a LOSAT Web file).</span>
  </div>
  <p
    v-if="state.settings"
    role="status"
    class="settings-message"
    :class="state.settings.kind === 'error' ? 'error' : 'note'"
    data-testid="settings-message"
    :data-kind="state.settings.kind"
  >
    {{ state.settings.text }}
  </p>
</template>
