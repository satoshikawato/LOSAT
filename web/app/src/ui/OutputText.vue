<script setup lang="ts">
// A text as the Outputs show it (domain/output-chunks.ts): blocks of lines that the browser lays out
// only near the view, so that a large output does not keep the shaped text of all its lines, in a
// box that is a Tab stop which the arrow keys scroll, also in WebKit (WA2-AX-08). The Outputs view
// shows an output with it, and RunFiles a file of the run (line RG E1w-3b). The attributes given to
// it (a test ID) go to its box.
import { computed, onMounted, onUnmounted, ref, watch } from 'vue';
import { outputChunks } from '../domain/output-chunks';
import { vScrollFocus } from './scroll-focus';

const props = defineProps<{ text: string; label: string }>();
const chunks = computed(() => outputChunks(props.text));
/** A new text has new blocks: a block remembers the height it was last laid out with. */
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
// The blocks of a new text may come inside a selection, which then does not change.
watch(generation, onSelectionChange, { flush: 'post' });
onUnmounted(() => document.removeEventListener('selectionchange', onSelectionChange));
watch(() => props.text, () => generation.value++);
</script>

<template>
  <pre ref="view" v-scroll-focus="label" class="output"><span v-for="(chunk, index) in chunks" :key="`${generation}:${index}`" class="output-chunk" :style="{ '--output-chunk-lines': chunk.lines }">{{ chunk.text }}</span></pre>
</template>
