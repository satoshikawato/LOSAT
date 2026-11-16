<script setup lang="ts">
// A list that draws only the rows in view (fixed row height), so lists of a run with many
// queries, subjects or HSPs stay responsive (S12's note: 100,000 records). The parent
// draws each row with the `row` slot, keeps the selected row in view with `reveal` (also when it
// is drawn again, and when the list is shown again after `display: none`), and
// shows the first rows again after a sort with `orderKey`. It reports the width of its scroll
// bar (`gutter`), so that a table's header can keep its columns over the rows' columns.
//
// The rows are one Tab stop (a roving tabindex, WA2-AX-01): the slot gives each row its
// `tabindex`, 0 for the active row (the row last focused, else the selected row, else the first)
// and -1 for the others. Up and Down, Home and End, Page Up and Page Down move the focus between
// rows, scrolling rows into view and drawing them; Enter and Space click a row that is not a
// button. The active row and the focused row are drawn wherever the list is scrolled, so that
// Tab comes back to the row and the focus never falls to the page (WA2-SC-02). The scroll box
// itself is not a Tab stop (Firefox makes a box that scrolls one), as the keys on the rows scroll it.
//
// `kind` 'list' is a list of items (aria-setsize and aria-posinset give each item's place among
// all, WA2-SC-03). `kind` 'grid' is the row group of a grid that the parent owns: the slot draws
// each `role="row"` with its aria-rowindex, and the parent gives the grid its aria-rowcount
// (WA2-AX-02).
import { computed, nextTick, onMounted, onUnmounted, ref, watch } from 'vue';

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
    /** A list of items, or the rows of the parent's grid. */
    kind?: 'list' | 'grid';
    /** The list's name (a grid has its name on the parent's grid). */
    label?: string | undefined;
    testid: string;
  }>(),
  { maxRows: 10, reveal: undefined, revealKey: undefined, orderKey: undefined, kind: 'list', label: undefined },
);
const emit = defineEmits<{ gutter: [px: number] }>();
defineSlots<{ row(props: { position: number; tabindex: 0 | -1 }): unknown }>();

const SPARE_ROWS = 6;
const viewport = ref<HTMLElement>();
const scrollTop = ref(0);
const heightPx = computed(() => Math.max(1, Math.min(props.count, props.maxRows)) * props.rowPx);
const first = computed(() => Math.max(0, Math.floor(scrollTop.value / props.rowPx) - SPARE_ROWS));
const last = computed(() => Math.min(props.count, first.value + props.maxRows + 2 * SPARE_ROWS));

/** The row that Tab stops at: the row last focused, else the selected row, else the first. */
const active = ref<number>();
/** The row that has the focus. */
const focused = ref<number>();
const inRange = (position: number | undefined): position is number => position !== undefined && position >= 0 && position < props.count;
const activePosition = computed(() => {
  if (inRange(active.value)) return active.value;
  if (inRange(props.reveal)) return props.reveal;
  return props.count > 0 ? 0 : undefined;
});
const positions = computed(() => {
  const drawn = Array.from({ length: Math.max(0, last.value - first.value) }, (_, i) => first.value + i);
  for (const kept of new Set([activePosition.value, focused.value])) {
    if (inRange(kept) && (kept < first.value || kept >= last.value)) drawn.push(kept);
  }
  return drawn.sort((a, b) => a - b);
});

/**
 * Scrolls the list so that the row at `position` is in view. False when the list is not shown
 * (`display: none`), where scrolling has no effect.
 */
function revealPosition(position: number): boolean {
  const element = viewport.value;
  if (element === undefined || element.clientHeight === 0) return false;
  const top = position * props.rowPx;
  if (top < element.scrollTop) element.scrollTop = top;
  else if (top + props.rowPx > element.scrollTop + heightPx.value) element.scrollTop = top + props.rowPx - heightPx.value;
  scrollTop.value = element.scrollTop;
  return true;
}
/** The selected row could not be brought into view while the list was not shown (WA2-SC-01). */
let revealPending = false;
/** Scrolls the list so that the selected row is in view. */
function revealSelected(): void {
  const position = props.reveal;
  if (position === undefined || position < 0) return;
  revealPending = !revealPosition(position);
}
watch(() => [props.revealKey, props.count] as const, revealSelected, { flush: 'post' });
// Another selection, order or set of rows (a view filter): Tab stops at the selected row again.
watch(
  () => [props.revealKey, props.orderKey, props.count] as const,
  () => {
    active.value = undefined;
    if (viewport.value === undefined || !viewport.value.contains(document.activeElement)) focused.value = undefined;
  },
);
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

/** The position of the row that holds `target`, if it is in this list. */
function positionOf(target: EventTarget | null): number | undefined {
  if (!(target instanceof Element) || viewport.value === undefined || !viewport.value.contains(target)) return undefined;
  const holder = target.closest<HTMLElement>('[data-position]');
  return holder === null ? undefined : Number(holder.dataset['position']);
}
/** The element that the slot drew for the row at `position` (when it is drawn). */
function rowElement(position: number): HTMLElement | undefined {
  const holder = viewport.value?.querySelector<HTMLElement>(`:scope > .virtual-spacer > [data-position="${position}"]`);
  const row = holder?.firstElementChild;
  return row instanceof HTMLElement ? row : undefined;
}
function onFocusin(event: FocusEvent): void {
  const position = positionOf(event.target);
  if (position === undefined) return;
  focused.value = position;
  active.value = position;
}
function onFocusout(event: FocusEvent): void {
  if (positionOf(event.relatedTarget) === undefined) focused.value = undefined;
}
async function focusRow(position: number): Promise<void> {
  active.value = position;
  revealPosition(position);
  await nextTick();
  const row = rowElement(position);
  if (row === undefined) return;
  row.focus({ preventScroll: true });
  // The list is scrolled already: this brings the row into the window (the page above the list may
  // have moved, as the Alignments do to show a chosen HSP's Range), never the table sideways.
  row.scrollIntoView({ block: 'nearest', inline: 'nearest' });
}
function onKeydown(event: KeyboardEvent): void {
  if (event.altKey || event.ctrlKey || event.metaKey || event.shiftKey) return;
  const from = positionOf(event.target);
  // Keys on a row only, not on a control inside one.
  if (from === undefined || event.target !== rowElement(from)) return;
  const row = event.target as HTMLElement;
  const page = Math.max(1, Math.min(props.count, props.maxRows) - 1);
  let to: number;
  switch (event.key) {
    case 'ArrowDown':
      to = from + 1;
      break;
    case 'ArrowUp':
      to = from - 1;
      break;
    case 'PageDown':
      to = from + page;
      break;
    case 'PageUp':
      to = from - page;
      break;
    case 'Home':
      to = 0;
      break;
    case 'End':
      to = props.count - 1;
      break;
    case 'Enter':
    case ' ':
      // A button clicks itself.
      if (row instanceof HTMLButtonElement) return;
      event.preventDefault();
      row.click();
      return;
    default:
      return;
  }
  event.preventDefault();
  void focusRow(Math.max(0, Math.min(props.count - 1, to)));
}

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
/**
 * Whether the list was shown at the last resize. A list in a tab that `v-show` hid is drawn
 * again with its selected row in view when the selection changed meanwhile (WA2-SC-01), else
 * where it was.
 */
let shown = false;
function onResize(): void {
  measureGutter();
  const element = viewport.value;
  if (element === undefined) return;
  const visible = element.clientHeight > 0;
  if (visible && !shown) {
    if (revealPending) revealSelected();
    else {
      element.scrollTop = scrollTop.value;
      scrollTop.value = element.scrollTop;
    }
  }
  shown = visible;
}
let observer: ResizeObserver | undefined;
onMounted(() => {
  observer = new ResizeObserver(onResize);
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
    :role="kind === 'grid' ? 'rowgroup' : 'list'"
    :aria-label="kind === 'grid' ? undefined : label"
    :data-testid="testid"
    :data-count="count"
    tabindex="-1"
    :style="{ height: `${heightPx}px` }"
    @scroll="scrollTop = ($event.target as HTMLElement).scrollTop"
    @focusin="onFocusin"
    @focusout="onFocusout"
    @keydown="onKeydown"
  >
    <div class="virtual-spacer" :style="{ height: `${count * rowPx}px` }">
      <div
        v-for="position in positions"
        :key="position"
        :role="kind === 'grid' ? undefined : 'listitem'"
        :aria-setsize="kind === 'grid' ? undefined : count"
        :aria-posinset="kind === 'grid' ? undefined : position + 1"
        class="virtual-row"
        :data-position="position"
        :style="{ top: `${position * rowPx}px`, height: `${rowPx}px` }"
      >
        <slot name="row" :position="position" :tabindex="position === activePosition ? 0 : -1" />
      </div>
    </div>
  </div>
</template>
