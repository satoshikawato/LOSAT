// A text box that scrolls is a Tab stop, so that the keyboard alone can read its text further
// (WA2-AX-08). Chromium and Firefox let the focus reach a scrolling box of their own accord;
// WebKit does not. While its text is larger than the box (it scrolls down or sideways), the box
// has tabindex 0 and is a region named by the directive's value, as "outfmt 0 output"; a box that
// does not scroll is not a stop. The size is looked at again when the box changes size (one
// observer for every box: it also looks at a box when it starts observing it) or the component
// that shows it is updated (its text changed).
import type { Directive } from 'vue';

const names = new WeakMap<HTMLElement, string>();

function mark(element: HTMLElement): void {
  const scrolls = element.scrollHeight > element.clientHeight + 1 || element.scrollWidth > element.clientWidth + 1;
  const name = names.get(element) ?? '';
  if (scrolls) {
    if (element.getAttribute('tabindex') !== '0') element.setAttribute('tabindex', '0');
    if (element.getAttribute('role') !== 'region') element.setAttribute('role', 'region');
    if (element.getAttribute('aria-label') !== name) element.setAttribute('aria-label', name);
  } else {
    element.removeAttribute('tabindex');
    element.removeAttribute('role');
    element.removeAttribute('aria-label');
  }
}

let observer: ResizeObserver | undefined;
const boxes = () => (observer ??= new ResizeObserver((entries) => entries.forEach((entry) => mark(entry.target as HTMLElement))));

export const vScrollFocus: Directive<HTMLElement, string> = {
  mounted(element, binding) {
    names.set(element, binding.value);
    boxes().observe(element);
  },
  updated(element, binding) {
    names.set(element, binding.value);
    mark(element);
  },
  unmounted(element) {
    boxes().unobserve(element);
  },
};
