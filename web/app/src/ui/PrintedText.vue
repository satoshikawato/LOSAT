<script setup lang="ts">
// The text that a run printed in place of a search (line WQA's D14, 2026-11-20; docs/web/abi_v2.md §5,
// stream 11): the command line's help (`-h`, `-help`, `-help-full`, `--help`, the XML of `-xmlhelp`) or
// version (`-version` and the full versions), byte for byte as the engine wrote it, shown as the
// Outputs show an output (OutputText, the bytes read once as UTF-8), with the command that prints it,
// one line that says so, and its export under the run's name (domain/output-format.ts,
// printedFileName). Each run has its own (keyed by the run), so that the reason of a save that failed
// stays with its run.
import { computed } from 'vue';
import type { Coordinator, RunView } from '../application/coordinator';
import { OutputSaver } from '../application/outputs';
import { userArgv } from '../domain/argv';
import CommandText from './CommandText.vue';
import OutputText from './OutputText.vue';
import { useStore } from './useStore';

const props = defineProps<{ coordinator: Coordinator; run: RunView }>();
const kind = computed(() => props.run.snapshot.prints ?? 'help');
const text = computed(() => new TextDecoder().decode(props.run.record.printed ?? new Uint8Array(0)));
const saver = new OutputSaver(props.coordinator);
const problem = useStore(saver.problem);
</script>

<template>
  <div class="printed-text" data-testid="printed-text" :data-kind="kind">
    <p class="command">
      Command: <code data-testid="result-command"><CommandText :words="['losat', ...userArgv(run.snapshot)]" /></code>
    </p>
    <p class="muted" data-testid="printed-note">The command above prints this text and runs no search.</p>
    <OutputText :text="text" :label="`The command line's ${kind}`" data-testid="printed-output" />
    <p v-if="problem" class="results-message error" role="alert" data-testid="printed-problem">{{ problem }}</p>
    <button type="button" data-testid="export-printed" @click="saver.savePrinted(run)">Export {{ kind }}</button>
  </div>
</template>
