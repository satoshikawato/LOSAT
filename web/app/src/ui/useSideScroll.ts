// Whether a table is wider than its box (it scrolls sideways), kept up to date as the box
// or the table's columns change size (S13 screen review M2: the tables say so while they do).
import { onMounted, onUnmounted, ref, type Ref } from 'vue';

export function useSideScroll(scroller: Ref<HTMLElement | undefined>): Ref<boolean> {
  const overflows = ref(false);
  const update = () => {
    const element = scroller.value;
    overflows.value = element !== undefined && element.scrollWidth > element.clientWidth + 1;
  };
  let observer: ResizeObserver | undefined;
  onMounted(() => {
    observer = new ResizeObserver(update);
    const element = scroller.value;
    if (element !== undefined) {
      observer.observe(element);
      // The header's width follows the columns' widths (for example the frames column).
      if (element.firstElementChild !== null) observer.observe(element.firstElementChild);
    }
    update();
  });
  onUnmounted(() => observer?.disconnect());
  return overflows;
}

/**
 * Focuses a row that the mouse or a finger presses without scrolling its table: the browsers
 * scrolled a table that is wider than its box sideways to the pressed row's button (Firefox on a
 * phone, W4 screen review middle 2). A row reached with the keyboard still scrolls into view.
 */
export function focusPressed(event: MouseEvent): void {
  event.preventDefault();
  (event.currentTarget as HTMLElement).focus({ preventScroll: true });
}
