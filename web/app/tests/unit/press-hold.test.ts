// The hold of a press (ui/press-hold.ts; line WQA, W-2 and W-3, and review r2 of W-2, M1): it starts with a press of
// the main button and ends with its click, and every way that a press ends without one ends it too.
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { holdWhilePressed } from '../../src/ui/press-hold';

/** A page with one held element of 900 px: its document and window as event targets, what the hold did to both. */
function page() {
  const view = new EventTarget();
  const document = Object.assign(new EventTarget(), { defaultView: view }) as unknown as Document;
  const styles = new Map<string, string>();
  const element = {
    style: {
      set minHeight(value: string) {
        styles.set('min-height', value);
      },
      removeProperty: (name: string) => styles.delete(name),
    },
    getBoundingClientRect: () => ({ height: 900 }),
  } as unknown as HTMLElement;
  const said: boolean[] = [];
  const stop = holdWhilePressed(
    () => element,
    (now) => said.push(now),
    document,
  );
  const fire = (target: EventTarget, type: string, fields: object = {}) => target.dispatchEvent(Object.assign(new Event(type), fields));
  const mouse = { pointerId: 1, isPrimary: true, pointerType: 'mouse' };
  return {
    view,
    document,
    stop,
    fire,
    said,
    held: () => styles.get('min-height'),
    press: (fields: object = {}) => fire(document, 'pointerdown', { ...mouse, button: 0, buttons: 1, ...fields }),
    release: (fields: object = {}) => fire(document, 'pointerup', { ...mouse, button: 0, buttons: 0, ...fields }),
    move: (buttons: number) => fire(document, 'pointermove', { ...mouse, button: -1, buttons }),
    click: () => fire(view, 'click'),
  };
}

beforeEach(() => {
  vi.useFakeTimers();
  vi.stubGlobal('getComputedStyle', () => ({ boxSizing: 'border-box' }));
});
afterEach(() => {
  vi.useRealTimers();
  vi.unstubAllGlobals();
});

describe('a press holds until its click', () => {
  it('keeps the height that the element had at the press, and says when the press starts and ends', () => {
    const p = page();
    p.press();
    expect(p.held()).toBe('900px');
    expect(p.said).toEqual([true]);
    p.release();
    expect(p.held()).toBe('900px');
    p.click();
    expect(p.held()).toBeUndefined();
    expect(p.said).toEqual([true, false]);
  });

  it('holds while the pointer moves with its button down', () => {
    const p = page();
    p.press();
    p.move(1);
    expect(p.held()).toBe('900px');
  });

  it('holds nothing for a press of another button: a right press makes no click', () => {
    const p = page();
    p.press({ button: 2, buttons: 2 });
    expect(p.held()).toBeUndefined();
    expect(p.said).toEqual([]);
  });

  it('ends when the page stops the hold', () => {
    const p = page();
    p.press();
    p.stop();
    expect(p.held()).toBeUndefined();
    expect(p.said).toEqual([true, false]);
  });
});

describe('a press that ends without a click ends the hold', () => {
  it('when the browser cancels it', () => {
    const p = page();
    p.press();
    p.fire(p.document, 'pointercancel', { pointerId: 1 });
    expect(p.held()).toBeUndefined();
    expect(p.said).toEqual([true, false]);
  });

  it('when a key is pressed after the release', () => {
    const p = page();
    p.press();
    p.release();
    p.fire(p.document, 'keydown', { key: 'Tab' });
    expect(p.held()).toBeUndefined();
    expect(p.said).toEqual([true, false]);
  });

  it('when the pointer moves with no button down after the release', () => {
    const p = page();
    p.press();
    p.release();
    p.move(0);
    expect(p.held()).toBeUndefined();
    expect(p.said).toEqual([true, false]);
  });

  it('a second after the release, when no click came', () => {
    const p = page();
    p.press();
    p.release();
    vi.advanceTimersByTime(999);
    expect(p.held()).toBe('900px');
    vi.advanceTimersByTime(1);
    expect(p.held()).toBeUndefined();
    expect(p.said).toEqual([true, false]);
  });

  it('when the pointer moves with no button down after a release that the page did not see', () => {
    const p = page();
    p.press();
    p.move(0);
    expect(p.held()).toBeUndefined();
    expect(p.said).toEqual([true, false]);
  });

  it('when the window loses the focus during the press', () => {
    const p = page();
    p.press();
    p.fire(p.view, 'blur');
    expect(p.held()).toBeUndefined();
    expect(p.said).toEqual([true, false]);
    // The release that comes back later ends nothing more.
    p.release();
    p.click();
    expect(p.said).toEqual([true, false]);
  });
});
