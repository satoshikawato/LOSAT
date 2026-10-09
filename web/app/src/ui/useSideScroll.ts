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
