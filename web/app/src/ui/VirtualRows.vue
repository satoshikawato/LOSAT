<script setup lang="ts">
// A list that draws only the rows in view (fixed row height), so lists of a run with many
// queries, subjects or HSPs stay responsive (S12's note: 100,000 records). The parent
// draws each row with the `row` slot and keeps the selected row in view with `reveal`.
import { computed, ref, watch } from 'vue';

const props = withDefaults(
  defineProps<{
    count: number;
    rowPx: number;
    /** Rows shown at most before the list scrolls. */
    maxRows?: number;
    /** Position of the row to bring into view when it changes. */
    reveal?: number | undefined;
    label: string;
    testid: string;
  }>(),
  { maxRows: 10, reveal: undefined },
);

const SPARE_ROWS = 6;
const viewport = ref<HTMLElement>();
const scrollTop = ref(0);
const heightPx = computed(() => Math.max(1, Math.min(props.count, props.maxRows)) * props.rowPx);
const first = computed(() => Math.max(0, Math.floor(scrollTop.value / props.rowPx) - SPARE_ROWS));
const last = computed(() => Math.min(props.count, first.value + props.maxRows + 2 * SPARE_ROWS));
const positions = computed(() => Array.from({ length: Math.max(0, last.value - first.value) }, (_, i) => first.value + i));

watch(
  () => [props.reveal, props.count] as const,
  ([position]) => {
    const element = viewport.value;
    if (position === undefined || position < 0 || element === undefined) return;
    const top = position * props.rowPx;
    if (top < element.scrollTop) element.scrollTop = top;
    else if (top + props.rowPx > element.scrollTop + heightPx.value) element.scrollTop = top + props.rowPx - heightPx.value;
    scrollTop.value = element.scrollTop;
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
