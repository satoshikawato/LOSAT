<script setup lang="ts">
// A list that draws only the rows in view (fixed row height), so lists of a run with many
// queries, subjects or HSPs stay responsive (S12's note: 100,000 records). The parent
// draws each row with the `row` slot, keeps the selected row in view with `reveal` (also when it
// is drawn again), and
// shows the first rows again after a sort with `orderKey`. It reports the width of its scroll
// bar (`gutter`), so that a table's header can keep its columns over the rows' columns.
import { computed, onMounted, onUnmounted, ref, watch } from 'vue';

const props = withDefaults(
  defineProps<{
    count: number;
    rowPx: number;
    /** Rows shown at most before the list scrolls. */
    maxRows?: number;
    /** Position of the selected row, brought into view when the selection or the count changes. */
    reveal?: number | undefined;
    /** Identity of the selected row: a sort moves the row without selecting another, and does not scroll to it. */
    revealKey?: string | undefined;
    /** The order of the rows (the sort): when it changes, the list shows its first rows. */
    orderKey?: string | undefined;
    label: string;
    testid: string;
  }>(),
  { maxRows: 10, reveal: undefined, revealKey: undefined, orderKey: undefined },
);
const emit = defineEmits<{ gutter: [px: number] }>();

const SPARE_ROWS = 6;
const viewport = ref<HTMLElement>();
const scrollTop = ref(0);
const heightPx = computed(() => Math.max(1, Math.min(props.count, props.maxRows)) * props.rowPx);
const first = computed(() => Math.max(0, Math.floor(scrollTop.value / props.rowPx) - SPARE_ROWS));
const last = computed(() => Math.min(props.count, first.value + props.maxRows + 2 * SPARE_ROWS));
const positions = computed(() => Array.from({ length: Math.max(0, last.value - first.value) }, (_, i) => first.value + i));

/** Scrolls the list so that the selected row is in view. */
function revealSelected(): void {
  const element = viewport.value;
  const position = props.reveal;
  if (position === undefined || position < 0 || element === undefined) return;
  const top = position * props.rowPx;
  if (top < element.scrollTop) element.scrollTop = top;
  else if (top + props.rowPx > element.scrollTop + heightPx.value) element.scrollTop = top + props.rowPx - heightPx.value;
  scrollTop.value = element.scrollTop;
}
watch(() => [props.revealKey, props.count] as const, revealSelected, { flush: 'post' });
watch(
  () => props.orderKey,
  () => {
    if (viewport.value === undefined) return;
    viewport.value.scrollTop = 0;
    scrollTop.value = 0;
  },
  { flush: 'post' },
);
watch(
  () => props.count,
  () => {
    if (viewport.value !== undefined && viewport.value.scrollTop > props.count * props.rowPx) {
      viewport.value.scrollTop = 0;
      scrollTop.value = 0;
    }
  },
);

/** The width of the vertical scroll bar, which appears with more rows than `maxRows` (none where scroll bars overlay). */
let gutter: number | undefined;
function measureGutter(): void {
  const element = viewport.value;
  if (element === undefined) return;
  const style = getComputedStyle(element);
  const borders = parseFloat(style.borderLeftWidth) + parseFloat(style.borderRightWidth);
  const width = Math.max(0, element.offsetWidth - element.clientWidth - borders);
  if (width !== gutter) {
    gutter = width;
    emit('gutter', width);
  }
}
watch(() => props.count, measureGutter, { flush: 'post' });
let observer: ResizeObserver | undefined;
onMounted(() => {
  observer = new ResizeObserver(measureGutter);
  if (viewport.value !== undefined) observer.observe(viewport.value);
  measureGutter();
  // A list drawn again (the results' tabs mount their tables) shows its selected row.
  revealSelected();
});
onUnmounted(() => observer?.disconnect());
</script>

<template>
  <div
    ref="viewport"
    class="virtual-rows"
    role="list"
    :aria-label="label"
    :data-testid="testid"
    :data-count="count"
    :style="{ height: `${heightPx}px` }"
    @scroll="scrollTop = ($event.target as HTMLElement).scrollTop"
  >
    <div class="virtual-spacer" :style="{ height: `${count * rowPx}px` }">
      <div
        v-for="position in positions"
        :key="position"
        role="listitem"
        class="virtual-row"
        :style="{ top: `${position * rowPx}px`, height: `${rowPx}px` }"
      >
        <slot name="row" :position="position" />
      </div>
    </div>
  </div>
</template>
