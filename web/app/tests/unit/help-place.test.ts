// Where a shown help goes across the page (ui/help-place.ts, HelpTip.vue). D32's review r2, L2 (2026-11-28): a help
// opened in a wide window kept its place when the window narrowed and the page scrolled sideways; it is placed again
// from where it would be without its shift.
import { describe, expect, it } from 'vitest';
import { HELP_MARGIN, helpRoom, helpShift } from '../../src/ui/help-place';

describe('helpShift', () => {
  it('leaves a help that ends inside the page where it is', () => {
    expect(helpShift(300, 700, 0, 1280)).toBe(0);
    expect(helpShift(41, 304, 0, 320)).toBe(0);
  });

  it(`moves a help left so that it ends ${HELP_MARGIN} px inside the right edge`, () => {
    expect(helpShift(660, 1300, 0, 1280)).toBe(36);
    expect(helpShift(660, 1300.2, 0, 1280)).toBe(37);
  });

  it(`never moves a help so far that it starts nearer the left edge than ${HELP_MARGIN} px`, () => {
    expect(helpShift(20, 400, 0, 320)).toBe(4);
    expect(helpShift(10, 400, 0, 320)).toBe(0);
  });

  it('places a help again from where it would be without its shift when the window narrows and widens', () => {
    // lcase_masking's help opened at 1280 px (660..1023), the window narrowed to 320 px: 41..329 before.
    expect(helpShift(41, 329, 0, 320)).toBe(25);
    // Measured with its 25 px of shift (16..304): the same place.
    expect(helpShift(16, 304, 25, 320)).toBe(25);
    // Narrowed further to 280 px (the help is 248 px wide there, at 41..289 without its shift).
    expect(helpShift(16, 264, 25, 280)).toBe(25);
    // Widened again: it goes back under its control.
    expect(helpShift(16, 304, 25, 1280)).toBe(0);
  });
});

describe('helpRoom', () => {
  it(`is the page less ${HELP_MARGIN} px at each side`, () => {
    expect(helpRoom(320)).toBe(288);
    expect(helpRoom(1280)).toBe(1248);
  });
});
