// The height that an element keeps while the page changes what it shows (ui/height-hold.ts; line WQA's D25, 2026-11-27):
// an input's block and card read again keep their height until the input is read.
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { nextTick, ref } from 'vue';
import { contentHeight, holdHeightWhile } from '../../src/ui/height-hold';

/** An element of `height` px whose height the test changes as its content would, and what the hold set on it. */
function element(height: number) {
  const styles = new Map<string, string>();
  let now = height;
  const target = {
    style: {
      set minHeight(value: string) {
        styles.set('min-height', value);
      },
      removeProperty: (name: string) => styles.delete(name),
    },
    getBoundingClientRect: () => ({ height: now }),
  } as unknown as HTMLElement;
  return {
    target,
    held: () => styles.get('min-height'),
    show: (next: number) => {
      now = next;
    },
  };
}

beforeEach(() => {
  vi.stubGlobal('getComputedStyle', () => ({ boxSizing: 'border-box' }));
});
afterEach(() => {
  vi.unstubAllGlobals();
});

describe('a height held while the page changes what it shows', () => {
  it('keeps the height that the element had before, until the hold ends', async () => {
    const e = element(171);
    const reading = ref(false);
    holdHeightWhile(() => e.target, () => reading.value);
    await nextTick();
    expect(e.held()).toBeUndefined();
    reading.value = true;
    await nextTick();
    expect(e.held()).toBe('171px');
    // What it shows while held is shorter: the height stays the one from before.
    e.show(57);
    reading.value = false;
    await nextTick();
    expect(e.held()).toBeUndefined();
    reading.value = true;
    await nextTick();
    expect(e.held()).toBe('57px');
  });

  it('holds nothing without an element', async () => {
    const reading = ref(false);
    holdHeightWhile(() => undefined, () => reading.value);
    reading.value = true;
    await nextTick();
    reading.value = false;
    await nextTick();
  });

  it('measures a box without its padding and border unless the box includes them', () => {
    vi.stubGlobal('getComputedStyle', () => ({
      boxSizing: 'content-box',
      paddingTop: '8px',
      paddingBottom: '10px',
      borderTopWidth: '1px',
      borderBottomWidth: '1px',
    }));
    expect(contentHeight(element(171).target)).toBe(151);
    vi.stubGlobal('getComputedStyle', () => ({ boxSizing: 'border-box', paddingTop: '8px' }));
    expect(contentHeight(element(171).target)).toBe(171);
  });
});
