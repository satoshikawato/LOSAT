<script setup lang="ts">
import { computed, ref, watch } from 'vue';
import type { Coordinator, RunView } from '../application/coordinator';
import { toShellCommand } from '../domain/argv';
import { programById } from '../domain/programs';
import { OUTPUT_FORMATS, type OutputFormat } from '../ports/engine';

const props = defineProps<{ coordinator: Coordinator; runs: readonly RunView[] }>();
const completed = computed(() => props.runs.filter((run) => run.status === 'completed'));
const selectedId = ref<string>();
const format = ref<OutputFormat>(6);
const text = ref('');

const selected = computed(
  () => completed.value.find((run) => run.snapshot.runId === selectedId.value) ?? completed.value.at(-1),
);

watch(
  [selected, format],
  async ([run, fmt]) => {
    text.value = run === undefined ? '' : await props.coordinator.readOutput(run.snapshot.runId, fmt);
  },
  { immediate: true },
);
</script>

<template>
  <h2>Results</h2>
  <p v-if="selected === undefined">No completed runs yet.</p>
  <template v-else>
    <label>
      Run
      <select v-model="selectedId" data-testid="result-run">
        <option v-for="run in completed" :key="run.snapshot.runId" :value="run.snapshot.runId">
          Run {{ run.snapshot.number }} · {{ programById(run.snapshot.program).label }}
        </option>
      </select>
    </label>
    <p class="command">
      Command: <code data-testid="result-command">{{ toShellCommand(selected.snapshot.argv) }}</code>
    </p>
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
    <pre class="output" data-testid="result-output">{{ text }}</pre>
    <button data-testid="export-output" @click="coordinator.exportOutput(selected.snapshot.runId, format)">
      Export outfmt {{ format }}
    </button>
  </template>
</template>
