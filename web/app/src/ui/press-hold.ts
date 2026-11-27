// A press keeps the page's height until its click (line WQA, W-2, 2026-11-20). A click goes where the
// press and the release both were. When the page became shorter between them than where it was
// scrolled to (a run that ends, its results read, a cancelled or failed run's message in place of
// its status, another run shown, a list closed by the press), the browser moved the page up by as
// much, and the release came on another control: the press was no click. So from a press until its
// click, the element (the page's main area) keeps at least the height that the page had at the
// press; a press never makes it shorter.
//
// What the page says above its tabs is not in the main area: when it went (or came) between a press
// and its release, everything below it moved by as much, with the same end (W-3, 2026-11-21: a mode
// pressed as LOSAT was asked again was not chosen). So the page also learns when a press starts and
// ends (`pressed`), and keeps what it says there as it was until then (AppView).
//
// Only a press of the main button (a mouse's left one, a touch, a pen's tip) can be a click, so only
// it holds (a right press's release may go to its context menu). The hold ends on events: as the click's dispatch starts (its target is chosen by then, and
// what the click itself does is not held), or when the page learns that the press ended without one:
// the browser cancels it, or the pointer is released and then moves or a key is pressed. A release
// that the page never sees (another window taken during the press) ends it as well: the page sees the
// pointer move with no button down, or the window lose the focus (review r2 of W-2, M1: until then,
// such a press held until the next one). A timer only ends a hold whose click never came after its
// release (a touch that drags, a long touch): a tap's click may come in a later task than its release.
//
// What came, grew, shrank or went above the pressed control between its press and its release moved the control itself,
// with the same end (line WQA's D25 and D27, 2026-11-27): the inputs read again under -parse_deflines ("Reading the
// records…" for a moment in place of each card, then "Checking…" in place of a card's warning) moved Run LOSAT up by
// 35 to 175 px and back, and so did the warning that the query card's check adds as it ends; the label "Run LOSAT
// (2 runs)" lost its count while the inputs were read, and the button became narrower. WebKit has no scroll anchoring,
// and Chromium's and Firefox's do not act at the page's top, nor for a change between what they keep in place and the
// control. So from a press until its click, while the pointer is over the pressed element, the page keeps that element
// where it was in the window and at least as large: when the page's layout moves its top, the page scrolls by as much.
// Its top, as what it shows below its top stays under the pointer when the element itself grows or shrinks at its end (a
// warning, a section that comes in it); the top of its box in the layout, as a pressed button is scaled about its middle.
// Where the page cannot scroll that far (its end, its top), the held element makes the room until the click: it is
// longer, or has space above what it shows; what came back above takes that space back first. The page looks at the
// element after each change of the document, before the next event (the release is one), and at each frame (a layout
// that changes with no change of the document). The browsers' own scroll anchoring is off meanwhile, as it would move
// the page a second time. The keeping ends as the pointer leaves the element (a move, a scroll, a release elsewhere),
// where a release would not click it. The page moves for nothing else: not for a press with the keys (the focused
// control takes it wherever it is), nor between presses, nor for what the click does, nor against a scroll during the
// press (a scroll does not change where the element is in the page), nor for an element fixed to the window or to a
// box that scrolls (position fixed or sticky), which the page's scroll does not move, nor for one that the press may
// resize (a text box's corner). A scroll of the page's own (AlignmentsView keeping the Range being read in place) waits
// while a press keeps a control that is in the page (`pressHeld`), and only then: a press on the space between things (a drag that
// selects text) keeps its element too, but no click is at stake.

/** The height that `element` has now, without its padding and border unless its box includes them (what min-height sets). */
function contentHeight(element: HTMLElement): number {
  const style = getComputedStyle(element);
  const edges =
    style.boxSizing === 'border-box'
      ? 0
      : [style.paddingTop, style.paddingBottom, style.borderTopWidth, style.borderBottomWidth].reduce((sum, value) => sum + (parseFloat(value) || 0), 0);
  return element.getBoundingClientRect().height - edges;
}

/** How long after its release a press whose click has not come still holds. */
const UNCLICKED_MS = 1000;

/** What a press can click: the elements that take the focus, or say by their role that a press on them acts. */
const CONTROLS = 'a[href], button, input, select, textarea, summary, label, [tabindex], [role="button"], [role="tab"], [role="row"], [role="option"], [role="menuitem"], [role="link"]';

/** The pages in which a press keeps the pressed element in its place now. */
const keeping = new WeakSet<Document>();

/**
 * Whether a press keeps a pressed control in its place in `page` now: from the press until its click, while the pointer
 * is over the control and the control is in the page.
 */
export function pressHeld(page: Document = document): boolean {
  return keeping.has(page);
}

/** A point in the window. */
interface Point {
  readonly x: number;
  readonly y: number;
}

/**
 * The box of `element` in the window, as shown (`shown`) and as it is in the layout (`laid`): without the scale of its
 * transform, which a pressed button has about its middle (scale(0.97)).
 */
function boxesOf(element: HTMLElement): { shown: DOMRect; laid: { top: number; bottom: number; left: number; right: number } } {
  const shown = element.getBoundingClientRect();
  const scale = /^matrix\(([^,]+),[^,]+,[^,]+,([^,]+),/.exec(getComputedStyle(element).transform ?? '');
  const width = shown.width / (Math.abs(Number(scale?.[1])) || 1);
  const height = shown.height / (Math.abs(Number(scale?.[2])) || 1);
  const x = shown.left + shown.width / 2;
  const y = shown.top + shown.height / 2;
  return { shown, laid: { top: y - height / 2, bottom: y + height / 2, left: x - width / 2, right: x + width / 2 } };
}

/**
 * Where the top of `element` is in the page: the top of its box in the layout, from the page's top, with the scroll of
 * the window and of every box around it added back, so that only a change of the page's layout moves it, not a scroll
 * nor a scale. Undefined when the element is not shown: gone from the page, or hidden.
 */
function placeOf(element: HTMLElement): number | undefined {
  if (!element.isConnected || element.getClientRects().length === 0) return undefined;
  let place = boxesOf(element).laid.top;
  for (let around = element.parentElement; around !== null; around = around.parentElement) place += around.scrollTop;
  return place;
}

/** Whether `at` is on `element`, as it is shown or as it is in the layout. */
function under(element: HTMLElement, at: Point): boolean {
  const { shown, laid } = boxesOf(element);
  const inside = (box: { top: number; bottom: number; left: number; right: number }): boolean =>
    at.x >= box.left && at.x <= box.right && at.y >= box.top && at.y <= box.bottom;
  return inside(shown) || inside(laid);
}

/**
 * The element that a press is on, to keep: the press's target, or the nearest element around it that is laid out as a
 * box (a button's icon is in the button). None for the page itself; for an element that the press may resize (a text
 * box's corner), whose size and place the press changes; or for one fixed to the window or to a box that scrolls, which
 * the page's scroll does not move.
 */
function pressedElement(event: Event, page: Document): HTMLElement | undefined {
  let node = event.target as Node | null;
  while (node && !('offsetParent' in node)) node = node.parentNode;
  const element = node as HTMLElement | null;
  if (!element || element === page.body || element === page.documentElement) return undefined;
  if (getComputedStyle(element).resize !== 'none') return undefined;
  for (let at: HTMLElement | null = element; at !== null; at = at.parentElement) {
    const position = getComputedStyle(at).position;
    if (position === 'fixed' || position === 'sticky') return undefined;
  }
  return element;
}

/**
 * Holds the height of `element()` from each press on `target` until its click, keeps the pressed element where it was in
 * the window and as large, and says to `pressed` when a press starts (true) and ends (false); returns the function that
 * stops.
 */
export function holdWhilePressed(
  element: () => HTMLElement | undefined,
  pressed: (now: boolean) => void = () => undefined,
  target: Document = document,
): () => void {
  const view = target.defaultView;
  /** The pointers down. */
  const down = new Set<number>();
  let held: HTMLElement | undefined;
  let pressing = false;
  let unclicked: ReturnType<typeof setTimeout> | undefined;
  /**
   * The element pressed and kept, where its top is in the page (`placeOf`) when it is in its place in the window,
   * scrolled as now, the pointer that pressed it, and where that pointer is in the window.
   */
  let kept: { element: HTMLElement; place: number; pointer: number; at: Point } | undefined;
  /** The styles that the hold set, besides the held element's height, with the values that they had before. */
  const styles = new Map<CSSStyleDeclaration, Map<string, string>>();
  const set = (style: CSSStyleDeclaration, name: string, value: string): void => {
    const before = styles.get(style) ?? new Map<string, string>();
    styles.set(style, before);
    if (!before.has(name)) before.set(name, style.getPropertyValue(name));
    style.setProperty(name, value);
  };
  /** Gives the styles that the hold set back their values from before: those of `style`, or all. */
  const restore = (style?: CSSStyleDeclaration): void => {
    for (const [each, before] of styles) {
      if (style !== undefined && each !== style) continue;
      for (const [name, value] of before) {
        if (value === '') each.removeProperty(name);
        else each.setProperty(name, value);
      }
      styles.delete(each);
    }
  };
  /** The space that the hold put above what the held element shows, and the padding that it had there before. */
  let room = 0;
  let padding = 0;
  let frame: number | undefined;
  /**
   * Takes the pressed element back to its place in the window when the page's layout moved it: the page scrolls by as
   * much. Where it cannot scroll that far, the held element makes the room: longer by what is missing for a scroll down;
   * for a scroll up past the page's top, space above what it shows, which takes the element down (one inside it). What
   * comes back above the element takes that space back first.
   */
  const keep = (): boolean => {
    if (kept === undefined || view === null) return false;
    // An element that left the page (a placeholder that its content took the place of) is no longer pressed: its place is
    // gone, and nothing is left to keep.
    if (!kept.element.isConnected) {
      unkeep();
      return false;
    }
    const now = placeOf(kept.element);
    if (now === undefined) return true;
    let moved = now - kept.place;
    // Less than the layout's unit (1/64 px) is no move of the layout, only the arithmetic of a scale.
    if (Math.abs(moved) < 1 / 64) return true;
    const page = target.scrollingElement ?? target.documentElement;
    // The held element when space above what it shows takes the pressed element down: the pressed element is in it.
    const around = held !== undefined && held !== kept.element && held.contains(kept.element) ? held : undefined;
    if (around !== undefined && moved > 0 && room > 0) {
      const taken = Math.min(room, moved);
      room -= taken;
      set(around.style, 'padding-top', `${padding + room}px`);
      moved -= taken;
    }
    let added = 0;
    if (held !== undefined && moved > 0) {
      // How far the page reaches: its scroll height, or where it fits the window (which that height then is), its own.
      const length = page.scrollHeight > page.clientHeight ? page.scrollHeight : target.documentElement.getBoundingClientRect().height;
      const missing = moved - (length - page.clientHeight - page.scrollTop);
      if (missing > 0) held.style.minHeight = `${contentHeight(held) + missing}px`;
    } else if (around !== undefined && moved < 0) {
      added = Math.max(0, -moved - page.scrollTop);
      if (added > 0) {
        if (room === 0) padding = parseFloat(getComputedStyle(around).paddingTop) || 0;
        room += added;
        set(around.style, 'padding-top', `${padding + room}px`);
      }
    }
    if (moved + added !== 0) {
      const scrolled = page.scrollTop;
      view.scrollBy({ top: moved + added, behavior: 'instant' });
      // What the scroll did not take back (the page could not scroll so far) is taken at the next look.
      kept.place += page.scrollTop - scrolled;
    }
    // The styles set here changed the document: that is no change to look at again.
    observer?.takeRecords();
    return true;
  };
  const observer = view?.MutationObserver === undefined ? undefined : new view.MutationObserver(keep);
  const eachFrame = (): void => {
    frame = undefined;
    if (keep()) frame = view?.requestAnimationFrame(eachFrame);
  };
  /** Stops keeping the pressed element in its place: its size and the browsers' scroll anchoring are as before. */
  const unkeep = (): void => {
    if (kept === undefined) return;
    restore(kept.element.style);
    if (target.body) restore(target.body.style);
    kept = undefined;
    keeping.delete(target);
    observer?.disconnect();
    if (frame !== undefined) view?.cancelAnimationFrame(frame);
    frame = undefined;
  };
  /** The pointer that pressed the kept element is at `at`: the element stays in its place while the pointer is on it. */
  const follow = (at?: Point): void => {
    if (kept === undefined) return;
    if (at !== undefined) kept.at = at;
    if (keep() && !under(kept.element, kept.at)) unkeep();
  };
  const release = (): void => {
    clearTimeout(unclicked);
    unclicked = undefined;
    unkeep();
    held?.style.removeProperty('min-height');
    held = undefined;
    restore();
    room = 0;
    if (!pressing) return;
    pressing = false;
    pressed(false);
  };
  const onDown = (event: PointerEvent): void => {
    if (event.button !== 0) return;
    // A primary pointer starts a new press: a pointer whose release never came is forgotten.
    if (event.isPrimary) down.clear();
    down.add(event.pointerId);
    clearTimeout(unclicked);
    unclicked = undefined;
    if (!pressing) {
      pressing = true;
      pressed(true);
    }
    const next = element();
    if (next !== undefined) {
      if (held !== undefined && held !== next) held.style.removeProperty('min-height');
      // The height that the element has now (a hold of an earlier press included), without its padding and border.
      next.style.minHeight = `${contentHeight(next)}px`;
      held = next;
    }
    // The element pressed (the latest press's, of several pointers) keeps its place and its size from here, while the
    // pointer is on it: the size that its box has in the layout (what min-width and min-height set; its pressed look,
    // scaled, is not in it).
    unkeep();
    const on = pressedElement(event, target);
    const at = { x: event.clientX, y: event.clientY };
    const place = on === undefined || !under(on, at) ? undefined : placeOf(on);
    if (on === undefined || place === undefined || observer === undefined) return;
    kept = { element: on, place, pointer: event.pointerId, at };
    // Only a press on a control can be a click that the page's moving loses: a press on the space between things
    // (a drag that selects text) keeps its element in its place, and the page's other scrolls go on.
    if (on.closest?.(CONTROLS) != null) keeping.add(target);
    if (on !== held) {
      const size = getComputedStyle(on);
      if (size.width.endsWith('px')) set(on.style, 'min-width', size.width);
      if (size.height.endsWith('px')) set(on.style, 'min-height', size.height);
    }
    if (target.body) set(target.body.style, 'overflow-anchor', 'none');
    observer.observe(target.documentElement, { subtree: true, childList: true, attributes: true, characterData: true });
    frame = view?.requestAnimationFrame(eachFrame);
  };
  const onUp = (event: PointerEvent): void => {
    down.delete(event.pointerId);
    // Released elsewhere than on the kept element, the press does not click it.
    if (event.pointerId === kept?.pointer) follow({ x: event.clientX, y: event.clientY });
    if (down.size === 0 && pressing) unclicked = setTimeout(release, UNCLICKED_MS);
  };
  const onCancel = (event: PointerEvent): void => {
    down.delete(event.pointerId);
    if (down.size === 0) release();
  };
  /** The click of the press is dispatched (the window's listener in the capture phase comes first). */
  const onClick = (): void => {
    if (down.size === 0) release();
  };
  /** A press released without a click (also a release that the page did not see): a move with no button down, or a key. */
  const onMove = (event: PointerEvent): void => {
    if (event.buttons === 0) down.delete(event.pointerId);
    if (pressing && down.size === 0) release();
    else if (event.pointerId === kept?.pointer) follow({ x: event.clientX, y: event.clientY });
  };
  /** The page or a box in it scrolled: the kept element may no longer be under the pointer. */
  const onScroll = (): void => follow();
  const onKey = (): void => {
    if (down.size === 0) release();
  };
  /** The window loses the focus: the page no longer sees the pointer, nor its release. */
  const onBlur = (event: Event): void => {
    if (event.target !== view) return;
    down.clear();
    release();
  };
  const capture = { capture: true, passive: true };
  target.addEventListener('pointerdown', onDown, capture);
  target.addEventListener('pointerup', onUp, capture);
  target.addEventListener('pointercancel', onCancel, capture);
  target.addEventListener('pointermove', onMove, capture);
  target.addEventListener('keydown', onKey, capture);
  target.addEventListener('scroll', onScroll, capture);
  view?.addEventListener('click', onClick, capture);
  view?.addEventListener('blur', onBlur, { passive: true });
  return () => {
    target.removeEventListener('pointerdown', onDown, capture);
    target.removeEventListener('pointerup', onUp, capture);
    target.removeEventListener('pointercancel', onCancel, capture);
    target.removeEventListener('pointermove', onMove, capture);
    target.removeEventListener('keydown', onKey, capture);
    target.removeEventListener('scroll', onScroll, capture);
    view?.removeEventListener('click', onClick, capture);
    view?.removeEventListener('blur', onBlur);
    down.clear();
    release();
  };
}
