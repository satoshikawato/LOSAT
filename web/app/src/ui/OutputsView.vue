<script setup lang="ts">
// The compatibility outputs of a run, byte for byte as the engine wrote them, with the CLI
// command of each format and their export (plan §5.8; the view filters never apply here): the
// formats of every run, and the output of the run's own -outfmt (W10, W11). A binary output
// (outfmt 9, and the ZIP archive of the files of 13 and 14) is exported, not shown. An output that cannot be read, or an export that is not
// saved, says why (application/outputs.ts).
import { computed, onMounted, onUnmounted, ref, useId, watch } from 'vue';
import type { Coordinator, RunView } from '../application/coordinator';
import { OutputReader } from '../application/outputs';
import { toShellCommand } from '../domain/argv';
import { outputChunks } from '../domain/output-chunks';
import { CUSTOM_OUTPUT, isBinaryOutput, outfmtOf, runOutputs, type RunOutput } from '../domain/output-format';
import { vScrollFocus } from './scroll-focus';
import TabList from './TabList.vue';
import { useStore } from './useStore';

const props = defineProps<{ coordinator: Coordinator; run: RunView }>();
const chosen = ref<RunOutput>(6);
const outputs = computed(() => runOutputs(props.run.snapshot.argv));
/** The output shown: the one chosen, or outfmt 6 for a run that does not write it. */
const format = computed(() => (outputs.value.includes(chosen.value) ? chosen.value : 6));
const binary = computed(() => isBinaryOutput(props.run.snapshot.argv, format.value));
/** The formats: tabs of one panel, the output shown (WA2-AX-12). */
const ids = useId();
const shown = computed({ get: () => format.value, set: (value: RunOutput) => (chosen.value = value) });
const tabs = computed(() =>
  outputs.value.map((f) => ({
    value: f,
    label: `outfmt ${outfmtOf(props.run.snapshot.argv, f)}`,
    testid: f === CUSTOM_OUTPUT ? 'format-custom' : `format-${f}`,
    id: `${ids}-tab-${f}`,
    panel: `${ids}-panel`,
    ...(f === CUSTOM_OUTPUT ? { class: 'custom-output' } : {}),
  })),
);
const reader = new OutputReader(props.coordinator);
const output = useStore(reader.state);
const chunks = computed(() => outputChunks(output.value.text));
/** A new output has new blocks: a block remembers the height it was last laid out with. */
const generation = ref(0);
const view = ref<HTMLElement>();
let marked = false;

/**
 * The blocks that a selection over more than one text spans are laid out while it lasts: a browser
 * may leave text that it did not lay out out of the selection (Firefox 151: 1 % of the output on
 * select all). A selection inside one text needs none: its block is laid out already.
 */
function onSelectionChange(): void {
  const selection = document.getSelection();
  const ranges: Range[] = [];
  for (let i = 0; i < (selection?.rangeCount ?? 0); i++) {
    const range = selection!.getRangeAt(i);
    const oneText = range.startContainer === range.endContainer && range.startContainer.nodeType === Node.TEXT_NODE;
    if (!range.collapsed && !oneText) ranges.push(range);
  }
  if (ranges.length === 0 && !marked) return;
  marked = false;
  for (const block of view.value?.children ?? []) {
    const selected = ranges.some((range) => range.intersectsNode(block));
    block.classList.toggle('output-chunk-selected', selected);
    marked ||= selected;
  }
}
onMounted(() => document.addEventListener('selectionchange', onSelectionChange));
// The blocks of a new output may come inside a selection, which then does not change.
watch(generation, onSelectionChange, { flush: 'post' });
onUnmounted(() => document.removeEventListener('selectionchange', onSelectionChange));

// The reader keeps the latest choice: a later one may be read first.
watch([() => props.run, format], ([run, fmt]) => void reader.show(run, fmt, !isBinaryOutput(run.snapshot.argv, fmt)), { immediate: true });
watch(() => output.value.text, () => generation.value++);
</script>

<template>
  <div class="outputs-view" data-testid="results-outputs">
    <TabList v-model="shown" label="Output format" :tabs="tabs" />
    <div :id="`${ids}-panel`" role="tabpanel" :aria-labelledby="`${ids}-tab-${format}`" class="outputs-panel">
      <p class="command">
        Command: <code data-testid="result-command">{{ toShellCommand(run.snapshot.argv, format) }}</code>
      </p>
      <!-- A Tab stop that the arrow keys scroll, also in WebKit (WA2-AX-08). -->
      <pre v-if="!output.failed && !binary" ref="view" v-scroll-focus="`outfmt ${outfmtOf(run.snapshot.argv, format)} output`" class="output" data-testid="result-output" :data-shown="output.shown"><span v-for="(chunk, index) in chunks" :key="`${generation}:${index}`" class="output-chunk" :style="{ '--output-chunk-lines': chunk.lines }">{{ chunk.text }}</span></pre>
      <!-- Next to Export, whose failure it also tells; where the text would be when the text could not be read. -->
      <p v-if="output.problem" class="results-message error" role="alert" data-testid="outputs-problem" :data-failed="output.failed">{{ output.problem }}</p>
      <button data-testid="export-output" @click="reader.export(run, format)">
        Export outfmt {{ outfmtOf(run.snapshot.argv, format) }}
      </button>
    </div>
  </div>
</template>
