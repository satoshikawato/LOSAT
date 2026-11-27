// What an element keeps of its height while the page changes what it shows (min-height): the main area from a press until
// its click (press-hold.ts), and the inputs' blocks and cards while their inputs are read again (holdHeightWhile).
import { watch } from 'vue';

/** The height that `element` has now, without its padding and border unless its box includes them (what min-height sets). */
export function contentHeight(element: HTMLElement): number {
  const style = getComputedStyle(element);
  const edges =
    style.boxSizing === 'border-box'
      ? 0
      : [style.paddingTop, style.paddingBottom, style.borderTopWidth, style.borderBottomWidth].reduce((sum, value) => sum + (parseFloat(value) || 0), 0);
  return element.getBoundingClientRect().height - edges;
}

/**
 * From when `holding()` becomes true until it is false again, `element()` keeps at least the height that it had just before
 * (line WQA's D25, 2026-11-27). The component's watcher runs before the component shows the change, so the height is the
 * one on the screen before it; the hold ends with the change that shows its end, so no frame shows the element without it.
 * What the element then shows takes the place that it had, and it may grow.
 */
export function holdHeightWhile(element: () => HTMLElement | undefined, holding: () => boolean): void {
  let held: HTMLElement | undefined;
  watch(holding, (now) => {
    if (now) {
      const target = element();
      if (held !== undefined || target === undefined) return;
      target.style.minHeight = `${contentHeight(target)}px`;
      held = target;
    } else if (held !== undefined) {
      held.style.removeProperty('min-height');
      held = undefined;
    }
  });
}
