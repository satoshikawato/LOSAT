<script setup lang="ts">
// The compatibility outputs of a run, byte for byte as the engine wrote them, with the CLI
// command of each format and their export (plan §5.8; the view filters never apply here).
import { computed, onMounted, onUnmounted, ref, watch } from 'vue';
import type { Coordinator, RunView } from '../application/coordinator';
import { toShellCommand } from '../domain/argv';
import { outputChunks } from '../domain/output-chunks';
import { OUTPUT_FORMATS, type OutputFormat } from '../domain/output-format';

const props = defineProps<{ coordinator: Coordinator; run: RunView }>();
const format = ref<OutputFormat>(6);
const text = ref('');
/** Which run and format `text` shows ("<number>:<format>"), once it is read. */
const shown = ref('');
const chunks = computed(() => outputChunks(text.value));
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

watch(
  [() => props.run, format],
  async ([run, fmt]) => {
    shown.value = '';
    const value = await props.coordinator.readOutput(run.snapshot.runId, fmt);
    // A later selection may have been read first.
    if (run !== props.run || fmt !== format.value) return;
    text.value = value;
    generation.value++;
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
    <pre ref="view" class="output" data-testid="result-output" :data-shown="shown"><span v-for="(chunk, index) in chunks" :key="`${generation}:${index}`" class="output-chunk" :style="{ '--output-chunk-lines': chunk.lines }">{{ chunk.text }}</span></pre>
    <button data-testid="export-output" @click="coordinator.exportOutput(run.snapshot.runId, format)">
      Export outfmt {{ format }}
    </button>
  </div>
</template>
