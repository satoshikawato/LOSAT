<script setup lang="ts">
// The files that a run made besides its outputs (line RG E1w-3b; docs/web/abi_v2.md §5, stream 8):
// the file of `-logfile NAME`, as the engine wrote it. Each file shows its name, its size and its
// text (a log, read as UTF-8: bytes that are not UTF-8 show as U+FFFD), and saves its bytes unchanged under
// the run's name for its files and its own name (domain/output-format.ts, runFileName). The Outputs
// of a completed run and the status of a failed run show them.
import { computed } from 'vue';
import type { Coordinator, RunView } from '../application/coordinator';
import { OutputSaver } from '../application/outputs';
import { formatBytes } from './format';
import { vScrollFocus } from './scroll-focus';
import { useStore } from './useStore';

const props = defineProps<{ coordinator: Coordinator; run: RunView }>();
const files = computed(() => props.run.record.runFiles ?? []);
const saver = new OutputSaver(props.coordinator);
const problem = useStore(saver.problem);
const decoder = new TextDecoder();
</script>

<template>
  <div class="run-files" data-testid="run-files">
    <section v-for="(file, index) in files" :key="`${run.snapshot.runId}:${index}`" class="run-file" :data-testid="`run-file-${index}`">
      <p class="run-file-name">
        <span class="option-name">-logfile</span> <span data-testid="run-file-name">{{ file.name }}</span><span class="muted small">, {{ formatBytes(file.bytes.length) }}</span>
      </p>
      <pre v-if="file.bytes.length > 0" v-scroll-focus="file.name" class="output" data-testid="run-file-text">{{ decoder.decode(file.bytes) }}</pre>
      <p v-else class="muted" data-testid="run-file-empty">The file is empty.</p>
      <button type="button" :data-testid="`run-file-download-${index}`" @click="saver.saveRunFile(run, index)">Download {{ file.name }}</button>
    </section>
    <p v-if="problem" class="results-message error" role="alert" data-testid="run-files-problem">{{ problem }}</p>
  </div>
</template>
