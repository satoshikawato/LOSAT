<script setup lang="ts">
// The files that a run made besides its outputs (line RG E1w-3b; docs/web/abi_v2.md §5, stream 8):
// the file of `-logfile NAME`, as the engine wrote it. Each file shows its name, its size and its
// text (a log, read once as UTF-8: bytes that are not UTF-8 show as U+FFFD) as the Outputs show an
// output (OutputText), and saves its bytes unchanged under the run's name for its files and its own
// name (domain/output-format.ts, runFileName). The Outputs of a completed run and the status of a
// failed run show them; each run has its own (keyed by the run), so that the reason of a save that
// failed stays with its run.
import { computed } from 'vue';
import type { Coordinator, RunView } from '../application/coordinator';
import { OutputSaver } from '../application/outputs';
import { formatBytes } from './format';
import OutputText from './OutputText.vue';
import { useStore } from './useStore';

const props = defineProps<{ coordinator: Coordinator; run: RunView }>();
const files = computed(() => props.run.record.runFiles ?? []);
/** The text of each file, decoded once for the run. */
const texts = computed(() => {
  const decoder = new TextDecoder();
  return files.value.map((file) => decoder.decode(file.bytes));
});
const saver = new OutputSaver(props.coordinator);
const problem = useStore(saver.problem);
</script>

<template>
  <div class="run-files" data-testid="run-files">
    <section v-for="(file, index) in files" :key="`${run.snapshot.runId}:${index}`" class="run-file" :data-testid="`run-file-${index}`">
      <p class="run-file-name">
        <span class="option-name">-logfile</span> <span data-testid="run-file-name">{{ file.name }}</span><span class="muted small">, {{ formatBytes(file.bytes.length) }}</span>
      </p>
      <OutputText v-if="file.bytes.length > 0" :text="texts[index]!" :label="file.name" data-testid="run-file-text" />
      <p v-else class="muted" data-testid="run-file-empty">The file is empty.</p>
      <button type="button" :data-testid="`run-file-download-${index}`" @click="saver.saveRunFile(run, index)">Download {{ file.name }}</button>
    </section>
    <p v-if="problem" class="results-message error" role="alert" data-testid="run-files-problem">{{ problem }}</p>
  </div>
</template>
