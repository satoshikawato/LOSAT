// The hold of a press (ui/press-hold.ts; line WQA, W-2 and W-3, and review r2 of W-2, M1): it starts with a press of
// the main button and ends with its click, and every way that a press ends without one ends it too.
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { holdWhilePressed, pressHeld } from '../../src/ui/press-hold';

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

/**
 * A page with a held element (the main area) and a button in it, laid out at `top` px from the page's top: the window's
 * scroll, its height and the page's, a MutationObserver and frames that the test runs, what the hold set on the styles.
 */
function laidOut(options: { top: number; scroll: number; window: number; height: number; position?: string; resize?: string }) {
  const style = () => {
    const values = new Map<string, string>();
    return {
      values,
      style: {
        setProperty: (name: string, value: string) => values.set(name, value),
        getPropertyValue: (name: string) => values.get(name) ?? '',
        removeProperty: (name: string) => values.delete(name),
        set minHeight(value: string) {
          values.set('min-height', value);
        },
      },
    };
  };
  const page = { ...style(), scrollTop: options.scroll, scrollHeight: options.height, clientHeight: options.window, parentElement: null };
  const body = { ...style(), scrollTop: 0, parentElement: page };
  const main = {
    ...style(),
    scrollTop: 0,
    parentElement: body,
    getBoundingClientRect: () => ({ top: -page.scrollTop, height: page.scrollHeight }),
    contains: (other: unknown) => other === button,
  };
  const at = { top: options.top };
  const button = {
    ...style(),
    offsetParent: main,
    isConnected: true,
    scrollTop: 0,
    parentElement: main,
    parentNode: main,
    getClientRects: () => [{}],
    getBoundingClientRect: () => ({ top: at.top - page.scrollTop, height: 40 }),
  };
  const styles = new Map<unknown, object>([
    [button, { position: options.position ?? 'static', resize: options.resize ?? 'none', width: '120px', height: '40px' }],
    [main, { position: 'static', boxSizing: 'border-box', paddingTop: '0px' }],
    [body, { position: 'static' }],
  ]);
  vi.stubGlobal('getComputedStyle', (element: unknown) => styles.get(element) ?? { boxSizing: 'border-box' });
  const observers: Array<{ callback: () => void; observing: boolean }> = [];
  const frames = new Map<number, () => void>();
  let frameId = 0;
  const scrolls: number[] = [];
  const view = Object.assign(new EventTarget(), {
    MutationObserver: class {
      entry: { callback: () => void; observing: boolean };
      constructor(callback: () => void) {
        this.entry = { callback, observing: false };
        observers.push(this.entry);
      }
      observe() {
        this.entry.observing = true;
      }
      disconnect() {
        this.entry.observing = false;
      }
    },
    requestAnimationFrame: (callback: () => void) => {
      frames.set(++frameId, callback);
      return frameId;
    },
    cancelAnimationFrame: (id: number) => frames.delete(id),
    scrollBy: ({ top }: { top: number }) => {
      const before = page.scrollTop;
      page.scrollTop = Math.max(0, Math.min(page.scrollHeight - page.clientHeight, page.scrollTop + top));
      if (page.scrollTop !== before) scrolls.push(page.scrollTop - before);
    },
  });
  const document = Object.assign(new EventTarget(), { defaultView: view, documentElement: page, scrollingElement: page, body });
  const stop = holdWhilePressed(() => main as unknown as HTMLElement, () => undefined, document as unknown as Document);
  const fire = (target: EventTarget, type: string, fields: object = {}) => target.dispatchEvent(Object.assign(new Event(type), fields));
  const mouse = { pointerId: 1, isPrimary: true, pointerType: 'mouse' };
  return {
    page,
    body,
    main,
    button,
    scrolls,
    stop,
    document,
    /** The layout moves the button to `top`, and the document says that it changed. */
    move: (top: number) => {
      at.top = top;
      for (const observer of observers) if (observer.observing) observer.callback();
    },
    /** The layout moves the button to `top` with no change of the document; the next frame comes. */
    moveInFrame: (top: number) => {
      at.top = top;
      const now = [...frames.values()];
      frames.clear();
      for (const frame of now) frame();
    },
    frames: () => frames.size,
    /** A press on the button (the event's target is the button, as the browser's would be). */
    press: () => {
      const event = Object.assign(new Event('pointerdown'), { ...mouse, button: 0, buttons: 1 });
      Object.defineProperty(event, 'target', { value: button });
      document.dispatchEvent(event);
    },
    pressKey: () => fire(document, 'keydown', { key: ' ' }),
    release: () => fire(document, 'pointerup', { ...mouse, button: 0, buttons: 0 }),
    click: () => fire(view, 'click'),
  };
}

describe('the pressed control stays under the pointer until its click (line WQA, D25 and D27)', () => {
  it('scrolls the page by as much as the layout moves the pressed element, and only until its click', () => {
    const p = laidOut({ top: 1500, scroll: 1000, window: 700, height: 3000 });
    p.press();
    p.move(1534);
    expect(p.scrolls).toEqual([34]);
    p.move(1400);
    expect(p.scrolls).toEqual([34, -134]);
    p.release();
    p.click();
    p.move(1450);
    p.moveInFrame(1480);
    expect(p.scrolls).toEqual([34, -134]);
    expect(p.frames()).toBe(0);
  });

  it('looks again at each frame, for a layout that changes with no change of the document', () => {
    const p = laidOut({ top: 1500, scroll: 1000, window: 700, height: 3000 });
    p.press();
    p.moveInFrame(1520);
    expect(p.scrolls).toEqual([20]);
    expect(p.frames()).toBe(1);
  });

  it('follows a scroll during the press: the element moves with it, and the page does not undo it', () => {
    const p = laidOut({ top: 1500, scroll: 1000, window: 700, height: 3000 });
    p.press();
    p.page.scrollTop = 1100;
    p.move(1500);
    p.moveInFrame(1500);
    expect(p.scrolls).toEqual([]);
  });

  it('makes the held element longer where the page cannot scroll so far down, until the click', () => {
    const p = laidOut({ top: 500, scroll: 0, window: 700, height: 710 });
    p.press();
    expect(p.main.values.get('min-height')).toBe('710px');
    p.move(534);
    expect(p.main.values.get('min-height')).toBe('734px');
    p.page.scrollHeight = 734;
    // The scroll that the room makes possible comes at the next look.
    p.moveInFrame(534);
    expect(p.scrolls).toEqual([10, 24]);
    p.release();
    p.click();
    expect(p.main.values.has('min-height')).toBe(false);
  });

  it('puts space above what the held element shows where the page cannot scroll so far up, until the click', () => {
    const p = laidOut({ top: 500, scroll: 10, window: 700, height: 3000 });
    p.press();
    // Above the button, 34 px went: the page scrolls up 10 px, to its top, and 24 px of space take the button down.
    p.move(466);
    expect(p.scrolls).toEqual([-10]);
    expect(p.main.values.get('padding-top')).toBe('24px');
    p.release();
    p.click();
    expect(p.main.values.has('padding-top')).toBe(false);
  });

  it('keeps the size of the pressed element and turns the browsers’ scroll anchoring off, until the click', () => {
    const p = laidOut({ top: 500, scroll: 0, window: 700, height: 3000 });
    p.button.values.set('min-width', '2em');
    p.press();
    expect(p.button.values.get('min-width')).toBe('120px');
    expect(p.button.values.get('min-height')).toBe('40px');
    expect(p.body.values.get('overflow-anchor')).toBe('none');
    p.release();
    p.click();
    expect(p.button.values.get('min-width')).toBe('2em');
    expect(p.button.values.has('min-height')).toBe(false);
    expect(p.body.values.has('overflow-anchor')).toBe(false);
  });

  it('says that a press holds from the press until its click', () => {
    const p = laidOut({ top: 500, scroll: 0, window: 700, height: 3000 });
    expect(pressHeld(p.document as unknown as Document)).toBe(false);
    p.press();
    expect(pressHeld(p.document as unknown as Document)).toBe(true);
    p.release();
    p.click();
    expect(pressHeld(p.document as unknown as Document)).toBe(false);
  });

  it('does not move the page for a press with the keys', () => {
    const p = laidOut({ top: 1500, scroll: 1000, window: 700, height: 3000 });
    p.pressKey();
    p.move(1534);
    expect(p.scrolls).toEqual([]);
  });

  for (const [what, options] of [
    ['fixed to the window', { position: 'fixed' }],
    ['fixed to a box that scrolls', { position: 'sticky' }],
    ['that the press may resize', { resize: 'both' }],
  ] as const) {
    it(`keeps nothing of an element ${what}`, () => {
      const p = laidOut({ top: 1500, scroll: 1000, window: 700, height: 3000, ...options });
      p.press();
      p.move(1534);
      expect(p.scrolls).toEqual([]);
      expect(p.button.values.has('min-width')).toBe(false);
      expect(p.frames()).toBe(0);
    });
  }
});
