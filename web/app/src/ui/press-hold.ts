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

/** How long after its release a press whose click has not come still holds. */
const UNCLICKED_MS = 1000;

/**
 * Holds the height of `element()` from each press on `target` until its click, and says to `pressed` when a press
 * starts (true) and ends (false); returns the function that stops.
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
  const release = (): void => {
    clearTimeout(unclicked);
    unclicked = undefined;
    held?.style.removeProperty('min-height');
    held = undefined;
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
    if (next === undefined) return;
    if (held !== undefined && held !== next) held.style.removeProperty('min-height');
    // The height that the element has now (a hold of an earlier press included), without its padding and border.
    const style = getComputedStyle(next);
    const edges =
      style.boxSizing === 'border-box'
        ? 0
        : [style.paddingTop, style.paddingBottom, style.borderTopWidth, style.borderBottomWidth].reduce((sum, value) => sum + (parseFloat(value) || 0), 0);
    next.style.minHeight = `${next.getBoundingClientRect().height - edges}px`;
    held = next;
  };
  const onUp = (event: PointerEvent): void => {
    down.delete(event.pointerId);
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
  };
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
  view?.addEventListener('click', onClick, capture);
  view?.addEventListener('blur', onBlur, { passive: true });
  return () => {
    target.removeEventListener('pointerdown', onDown, capture);
    target.removeEventListener('pointerup', onUp, capture);
    target.removeEventListener('pointercancel', onCancel, capture);
    target.removeEventListener('pointermove', onMove, capture);
    target.removeEventListener('keydown', onKey, capture);
    view?.removeEventListener('click', onClick, capture);
    view?.removeEventListener('blur', onBlur);
    down.clear();
    release();
  };
}
