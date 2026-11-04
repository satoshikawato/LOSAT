<script setup lang="ts">
// The compatibility outputs of a run, byte for byte as the engine wrote them, with the CLI
// command of each format and their export (plan §5.8; the view filters never apply here). An
// output that cannot be read, or an export that is not saved, says why (application/outputs.ts).
import { computed, onMounted, onUnmounted, ref, watch } from 'vue';
import type { Coordinator, RunView } from '../application/coordinator';
import { OutputReader } from '../application/outputs';
import { toShellCommand } from '../domain/argv';
import { outputChunks } from '../domain/output-chunks';
import { OUTPUT_FORMATS, type OutputFormat } from '../domain/output-format';
import { useStore } from './useStore';

const props = defineProps<{ coordinator: Coordinator; run: RunView }>();
const format = ref<OutputFormat>(6);
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
watch([() => props.run, format], ([run, fmt]) => void reader.show(run, fmt), { immediate: true });
watch(() => output.value.text, () => generation.value++);
</script>

<template>
  <div class="outputs-view" data-testid="results-outputs">
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
    <pre v-if="!output.failed" ref="view" class="output" data-testid="result-output" :data-shown="output.shown"><span v-for="(chunk, index) in chunks" :key="`${generation}:${index}`" class="output-chunk" :style="{ '--output-chunk-lines': chunk.lines }">{{ chunk.text }}</span></pre>
    <!-- Next to Export, whose failure it also tells; where the text would be when the text could not be read. -->
    <p v-if="output.problem" class="results-message error" role="alert" data-testid="outputs-problem" :data-failed="output.failed">{{ output.problem }}</p>
    <button data-testid="export-output" @click="reader.export(run, format)">
      Export outfmt {{ format }}
    </button>
  </div>
</template>
