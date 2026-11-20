// The search page at the widths of phones, tablets and desktops (WCAG 2.1 SC 1.4.10 Reflow: 320 CSS
// px with no sideways scrolling). Line WQA's random walk found the page 2 px wider than a 320 px
// window in Firefox and WebKit: the mode tabs' fieldset was as wide as its content's longest line
// (a fieldset's min-inline-size is min-content), 16 px wider than the search panel, and the names
// "TBLASTN" and "TBLASTX" passed their tabs' borders in the three engines (2026-11-19).
import { expect } from '@playwright/test';
import { openSearch } from './support/search';
import { test } from './support/test';

const WIDTHS = [320, 330, 340, 350, 359, 360, 390, 768, 1280] as const;

for (const width of WIDTHS) {
  test(`at ${width} px the search page is no wider than the window, and each mode's name fits its tab`, async ({ page }) => {
    await page.setViewportSize({ width, height: 800 });
    await openSearch(page);
    await page.evaluate(() => document.fonts.ready);
    const layout = await page.evaluate(() => {
      const box = (element: Element) => element.getBoundingClientRect();
      const panel = box(document.querySelector('[data-testid="search-panel"]')!);
      const row = box(document.querySelector('[data-testid="program-tabs"]')!);
      const tabs = [...document.querySelectorAll('[data-testid="program-tabs"] label')].map((label) => {
        const text = [...label.childNodes].find((node) => node.nodeType === Node.TEXT_NODE && (node.textContent ?? '').trim() !== '')!;
        const range = document.createRange();
        range.selectNodeContents(text);
        const words = range.getBoundingClientRect();
        const tab = box(label);
        return { name: (text.textContent ?? '').trim(), left: words.left - tab.left, right: tab.right - words.right, size: parseFloat(getComputedStyle(label).fontSize) };
      });
      return {
        page: document.documentElement.scrollWidth,
        window: document.documentElement.clientWidth,
        rowInPanel: row.left >= panel.left - 0.5 && row.right <= panel.right + 0.5,
        tabs,
      };
    });
    expect(layout.page).toBeLessThanOrEqual(layout.window);
    expect(layout.rowInPanel).toBe(true);
    // Each name lies inside its tab's 1 px border.
    for (const tab of layout.tabs) {
      expect(tab.left, `${tab.name} at ${width} px`).toBeGreaterThanOrEqual(0.5);
      expect(tab.right, `${tab.name} at ${width} px`).toBeGreaterThanOrEqual(0.5);
    }
    // From 360 px the tabs keep the type of the narrow layout (0.9em of the panel's 14 px; Firefox keeps
    // sizes in 1/64 px) or the wide one.
    if (width >= 360) for (const tab of layout.tabs) expect(tab.size).toBeCloseTo(width > 800 ? 14 : 12.6, 1);
  });
}
