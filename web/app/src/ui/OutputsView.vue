<script setup lang="ts">
// The compatibility outputs of a run, byte for byte as the engine wrote them, with the CLI
// command of each format and their export (plan §5.8; the view filters never apply here).
// LOSAT Web's own formats follow in a section of their own (ExportFormats.vue, S15).
import { ref, watch } from 'vue';
import type { Coordinator, RunView } from '../application/coordinator';
import type { ResultExporter } from '../application/result-export';
import type { ResultsState } from '../application/results';
import { toShellCommand } from '../domain/argv';
import { OUTPUT_FORMATS, type OutputFormat } from '../domain/output-format';
import ExportFormats from './ExportFormats.vue';

const props = defineProps<{ coordinator: Coordinator; run: RunView; exporter: ResultExporter; state: ResultsState }>();
const format = ref<OutputFormat>(6);
const text = ref('');
/** Which run and format `text` shows ("<number>:<format>"), once it is read. */
const shown = ref('');

watch(
  [() => props.run, format],
  async ([run, fmt]) => {
    shown.value = '';
    const value = await props.coordinator.readOutput(run.snapshot.runId, fmt);
    // A later selection may have been read first.
    if (run !== props.run || fmt !== format.value) return;
    text.value = value;
    shown.value = `${run.snapshot.number}:${fmt}`;
  },
  { immediate: true },
);
</script>

<template>
  <div class="outputs-view" data-testid="results-outputs">
    <p class="muted">The outputs as the engine wrote them. The view filters do not change them.</p>
    <nav class="tabs" aria-label="Output format">
      <button
        v-for="f in OUTPUT_FORMATS"
        :key="f"
        :aria-pressed="format === f"
        :data-testid="`format-${f}`"
        @click="format = f"
      >
        outfmt {{ f }}
      </button>
    </nav>
    <p class="command">
      Command: <code data-testid="result-command">{{ toShellCommand(run.snapshot.argv, format) }}</code>
    </p>
    <pre class="output" data-testid="result-output" :data-shown="shown">{{ text }}</pre>
    <button data-testid="export-output" @click="coordinator.exportOutput(run.snapshot.runId, format)">
      Export outfmt {{ format }}
    </button>
    <ExportFormats :exporter="exporter" :state="state" />
  </div>
</template>
