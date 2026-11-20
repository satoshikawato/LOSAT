// A press keeps the page's height until its click (line WQA, W-2, 2026-11-20). A click goes where the
// press and the release both were. When the page became shorter between them than where it was
// scrolled to (a run that ends, its results read, a cancelled or failed run's message in place of
// its status, another run shown, a list closed by the press), the browser moved the page up by as
// much, and the release came on another control: the press was no click. So from a press until its
// click, the element (the page's main area) keeps at least the height that the page had at the
// press; a press never makes it shorter. A press that ends without a click (cancelled by the
// browser, or released elsewhere and followed by a move or a key) ends the hold as well.
//
// The hold ends on events, not on a timer: a tap's click may come in a later task than its release.

/** Holds the height of `element()` from each press on `target` until its click; returns the function that stops. */
export function holdWhilePressed(element: () => HTMLElement | undefined, target: Document = document): () => void {
  const view = target.defaultView;
  /** The pointers down. */
  const down = new Set<number>();
  let held: HTMLElement | undefined;
  const release = (): void => {
    held?.style.removeProperty('min-height');
    held = undefined;
  };
  const onDown = (event: PointerEvent): void => {
    // A primary pointer starts a new press: a pointer whose release never came is forgotten.
    if (event.isPrimary) down.clear();
    down.add(event.pointerId);
    const next = element();
    if (next === undefined) return;
    if (held !== undefined && held !== next) release();
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
  };
  const onCancel = (event: PointerEvent): void => {
    down.delete(event.pointerId);
    if (down.size === 0) release();
  };
  /** The click of the press has been dispatched (the window's listener in the bubble phase comes last). */
  const onClick = (): void => {
    if (down.size === 0) release();
  };
  /** A press released without a click: the pointer moves on with no button down, or a key is pressed. */
  const onMove = (event: PointerEvent): void => {
    if (held !== undefined && down.size === 0 && event.buttons === 0) release();
  };
  const onKey = (): void => {
    if (down.size === 0) release();
  };
  const capture = { capture: true, passive: true };
  target.addEventListener('pointerdown', onDown, capture);
  target.addEventListener('pointerup', onUp, capture);
  target.addEventListener('pointercancel', onCancel, capture);
  target.addEventListener('pointermove', onMove, capture);
  target.addEventListener('keydown', onKey, capture);
  view?.addEventListener('click', onClick, { passive: true });
  return () => {
    target.removeEventListener('pointerdown', onDown, capture);
    target.removeEventListener('pointerup', onUp, capture);
    target.removeEventListener('pointercancel', onCancel, capture);
    target.removeEventListener('pointermove', onMove, capture);
    target.removeEventListener('keydown', onKey, capture);
    view?.removeEventListener('click', onClick);
    down.clear();
    release();
  };
}
