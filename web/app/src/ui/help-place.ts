// Where a shown help (the "?" of HelpTip.vue) goes across the page: it floats under its control, at most as wide as
// the page less a margin at each side, and moves left when it would end past the page's right edge, but never so far
// that it starts nearer the left edge than the margin.

/** The space the help keeps from each edge of the page, in pixels. */
export const HELP_MARGIN = 16;

/** The widest a help may be in a page `page` pixels wide. */
export function helpRoom(page: number): number {
  return page - 2 * HELP_MARGIN;
}

/**
 * The pixels a help moves left. `left` and `right` are its box as measured, with the `shifted` pixels it has moved
 * already, so a help that moved for a narrow window comes back when the window widens. `page` is the page's width:
 * the window's less a scroll bar.
 */
export function helpShift(left: number, right: number, shifted: number, page: number): number {
  const over = Math.ceil(right + shifted - (page - HELP_MARGIN));
  const most = Math.floor(left + shifted - HELP_MARGIN);
  return Math.min(Math.max(0, over), Math.max(0, most));
}
