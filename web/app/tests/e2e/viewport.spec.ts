// The search page at the widths of phones, tablets and desktops (WCAG 2.1 SC 1.4.10 Reflow: 320 CSS
// px with no sideways scrolling). Line WQA's random walk found the page 2 px wider than a 320 px
// window in Firefox and WebKit: the mode tabs' fieldset was as wide as its content's longest line
// (a fieldset's min-inline-size is min-content), 16 px wider than the search panel, and the names
// "TBLASTN" and "TBLASTX" passed their tabs' borders in the three engines (2026-11-19).
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { REPOSITORY } from './support/harness-server';
import { openParameters, openSearch, pickRun, program, waitStatus } from './support/search';
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

const NARROW_WIDTHS = [320, 390, 768] as const;
const TIP_WIDTHS = [320, 390, 768, 850, 900] as const;
const COUNTS_FIXTURE = 'LOSAT/tests/winmask/counts/mg1655_u9.obinary';

interface HelpTipMeasurement {
  readonly id: string;
  readonly btnLeft: number;
  readonly btnTop: number;
  readonly btnBottom: number;
  readonly prevRight: number;
  readonly prevTop: number;
  readonly prevBottom: number;
  readonly overlapsY: boolean;
  readonly followsX: boolean;
}

/**
 * Checks every visible "?" (`.help-toggle`) on the current view: it shares a line with the content
 * right before it (for a sentence, the last rect of `Range.getClientRects()` over its text; for a
 * control, the control's bounding rect), so no "?" ever starts a line alone.
 */
async function measureAllHelpTips(page: Page): Promise<{ readonly scrollWidth: number; readonly clientWidth: number; readonly tips: readonly HelpTipMeasurement[] }> {
  return page.evaluate(() => {
    const isVisibleElement = (el: Element): boolean => {
      if (el.classList.contains('visually-hidden')) return false;
      const r = el.getBoundingClientRect();
      return r.width > 0.5 && r.height > 0.5;
    };

    const lastRangeRect = (start: Node, endBefore?: Node): DOMRect | undefined => {
      const range = document.createRange();
      if (endBefore !== undefined) {
        range.setStartBefore(start);
        range.setEndBefore(endBefore);
      } else {
        range.selectNodeContents(start);
      }
      return [...range.getClientRects()].filter((r) => r.width > 0.5 && r.height > 0.5).at(-1);
    };

    const findPrecedingElement = (node: Element): Element | null => {
      let cur: Element | null = node.previousElementSibling;
      while (cur !== null && !isVisibleElement(cur)) cur = cur.previousElementSibling;
      if (cur !== null) return cur;
      const parent = node.parentElement;
      if (parent !== null && parent.classList.contains('control-tip')) {
        return findPrecedingElement(parent);
      }
      return null;
    };

    const tips: HelpTipMeasurement[] = [];
    for (const btn of document.querySelectorAll('.help-toggle')) {
      const br = btn.getBoundingClientRect();
      if (br.width < 0.5 || br.height < 0.5) continue;
      const tip = btn.closest('.help-tip') ?? btn;
      const wrapper = tip.closest('.hint-line, .param-hint-line');

      let prevRect: DOMRect | undefined;
      if (wrapper !== null && wrapper.firstChild !== null && wrapper.firstChild !== tip) {
        prevRect = lastRangeRect(wrapper.firstChild, tip);
      } else {
        const prev = findPrecedingElement(tip);
        if (prev !== null) {
          if (prev.matches('.param-hint, .hint, .file-missing, .file-name')) {
            prevRect = lastRangeRect(prev);
          } else if (prev.matches('.file-pick')) {
            const children = [...prev.children].filter(isVisibleElement);
            const lastChild = children.at(-1);
            if (lastChild !== undefined) {
              prevRect = lastChild.matches('.file-missing, .file-name, .hint')
                ? lastRangeRect(lastChild)
                : lastChild.getBoundingClientRect();
            }
          } else {
            prevRect = prev.getBoundingClientRect();
          }
        }
      }

      if (prevRect === undefined) {
        tips.push({
          id: btn.getAttribute('data-testid') ?? btn.getAttribute('aria-label') ?? 'help-tip',
          btnLeft: br.left,
          btnTop: br.top,
          btnBottom: br.bottom,
          prevRight: 0,
          prevTop: 0,
          prevBottom: 0,
          overlapsY: false,
          followsX: false,
        });
        continue;
      }

      tips.push({
        id: btn.getAttribute('data-testid') ?? btn.getAttribute('aria-label') ?? 'help-tip',
        btnLeft: br.left,
        btnTop: br.top,
        btnBottom: br.bottom,
        prevRight: prevRect.right,
        prevTop: prevRect.top,
        prevBottom: prevRect.bottom,
        overlapsY: Math.min(br.bottom, prevRect.bottom) > Math.max(br.top, prevRect.top),
        followsX: br.left >= prevRect.right - 1,
      });
    }
    return {
      scrollWidth: document.documentElement.scrollWidth,
      clientWidth: document.documentElement.clientWidth,
      tips,
    };
  });
}

function expectAllHelpTipsShareLine(
  measured: { readonly scrollWidth: number; readonly clientWidth: number; readonly tips: readonly HelpTipMeasurement[] },
  context: string,
  requireAtLeastOne = true,
): void {
  expect(measured.scrollWidth, `${context}: page should not overflow horizontally`).toBeLessThanOrEqual(measured.clientWidth);
  if (requireAtLeastOne) {
    expect(measured.tips.length, `${context}: expected at least one visible HelpTip`).toBeGreaterThan(0);
  }
  for (const tip of measured.tips) {
    expect(
      tip.overlapsY,
      `${context}: ${tip.id} should overlap its preceding content vertically (btn ${tip.btnTop.toFixed(1)}..${tip.btnBottom.toFixed(1)} vs prev ${tip.prevTop.toFixed(1)}..${tip.prevBottom.toFixed(1)})`,
    ).toBe(true);
    expect(
      tip.followsX,
      `${context}: ${tip.id} should sit after its preceding content (btnLeft ${tip.btnLeft.toFixed(1)} vs prevRight ${tip.prevRight.toFixed(1)})`,
    ).toBe(true);
  }
}

/** Checks that the WindowMasker checkbox shares a line with the first line of its label. */
async function expectWmCheckboxSharesLine(page: Page, context: string): Promise<void> {
  const m = await page.evaluate(() => {
    const box = document.getElementById('param-window_masker_db')!.getBoundingClientRect();
    const label = document.getElementById('param-window_masker_db-label')!;
    const range = document.createRange();
    range.selectNodeContents(label);
    const first = [...range.getClientRects()].find((r) => r.width > 0.5 && r.height > 0.5)!;
    return {
      boxTop: box.top,
      boxBottom: box.bottom,
      boxRight: box.right,
      firstLeft: first.left,
      firstTop: first.top,
      firstBottom: first.bottom,
      overlapsY: Math.min(box.bottom, first.bottom) > Math.max(box.top, first.top),
      followsX: first.left >= box.right - 1,
    };
  });
  expect(m.overlapsY, `${context}: WindowMasker checkbox (${m.boxTop.toFixed(1)}..${m.boxBottom.toFixed(1)}) should share a line with its label (${m.firstTop.toFixed(1)}..${m.firstBottom.toFixed(1)})`).toBe(true);
  expect(m.followsX, `${context}: WindowMasker label (${m.firstLeft.toFixed(1)}) should start to the right of its checkbox (${m.boxRight.toFixed(1)})`).toBe(true);
}

async function openWmTab(page: Page): Promise<void> {
  await page.getByTestId('tab-windowmasker').click();
  await expect(page.getByTestId('wm-panel')).toHaveAttribute('data-describe', 'ready', { timeout: 60_000 });
}

async function openWmOptions(page: Page): Promise<void> {
  const wmToggle = page.getByTestId('wm-options-toggle');
  if ((await wmToggle.getAttribute('aria-expanded')) !== 'true') await wmToggle.click();
}

async function saveMadeCounts562(page: Page, runNumber: number): Promise<void> {
  await openWmTab(page);
  if (runNumber === 1) {
    await page.getByTestId('wm-fasta-input').fill('>g1\nACGTTGCAAAGCTTAGGCTAACGTTGCAAAGCTTAGGCTAACGTTGCAAAGCTTAGGCTA\n');
  }
  await expect(page.getByTestId('wm-argv-validation')).toHaveAttribute('data-state', 'ok', { timeout: 30_000 });
  await page.getByTestId('wm-run-losat').click();
  await waitStatus(page, runNumber, 'completed', 120_000);
  await pickRun(page, runNumber);
  await page.getByTestId('tool-tab-file-0').click();
  await page.getByTestId('tool-save-taxid-input').fill('562');
  await page.getByTestId('tool-save-taxid-button').click();
}

async function startWithSavedCounts562(page: Page): Promise<void> {
  await page.setViewportSize({ width: 1280, height: 900 });
  await openSearch(page);
  await saveMadeCounts562(page, 1);
  await expect(page.getByTestId('tool-save-status')).toBeVisible({ timeout: 30_000 });
}

async function forWidths(page: Page, widths: readonly number[], check: (width: number) => Promise<void>): Promise<void> {
  for (const width of widths) {
    await page.setViewportSize({ width, height: 844 });
    await page.evaluate(() => document.fonts.ready);
    await check(width);
  }
}

async function attachCountsFixture(page: Page): Promise<void> {
  await page.getByTestId('param-window_masker_db-by').selectOption('file');
  await page.getByTestId('param-window_masker_db-file').setInputFiles({
    name: 'mg1655_u9.obinary',
    mimeType: 'application/octet-stream',
    buffer: readFileSync(join(REPOSITORY, COUNTS_FIXTURE)),
  });
  await expect(page.getByTestId('param-window_masker_db-name')).toBeVisible();
}

async function openBlastnParameters(page: Page): Promise<void> {
  await page.getByTestId('tab-search').click();
  await program(page, 'blastn');
  await openParameters(page);
}

async function selectTaxidMade562(page: Page): Promise<void> {
  await page.getByTestId('param-window_masker_db').check();
  await page.getByTestId('param-window_masker_db-taxid').fill('562');
  await expect(page.getByTestId('counts-note')).toContainText(/Made in this (?:browser|tab only)/);
}

const checkAllTipsAtWidths = (page: Page, label: string, requireAtLeastOne = true): Promise<void> =>
  forWidths(page, TIP_WIDTHS, async (width) =>
    expectAllHelpTipsShareLine(await measureAllHelpTips(page), `${label} at ${width} px`, requireAtLeastOne),
  );

const checkWmCheckboxAtNarrowWidths = (page: Page, label: string): Promise<void> =>
  forWidths(page, NARROW_WIDTHS, (width) => expectWmCheckboxSharesLine(page, `${label} at ${width} px`));

test('at 320, 390 and 768 px the Replace dialog leaves at least 12 px on each side', async ({ page }) => {
  test.setTimeout(180_000);
  await startWithSavedCounts562(page);
  await saveMadeCounts562(page, 2);
  await expect(page.getByTestId('replace-counts-dialog')).toBeVisible();

  await forWidths(page, NARROW_WIDTHS, async (width) => {
    const gaps = await page.evaluate(() => {
      const r = document.querySelector('[data-testid="replace-counts-dialog"]')!.getBoundingClientRect();
      return {
        left: r.left,
        rightInner: window.innerWidth - r.right,
        rightClient: document.documentElement.clientWidth - r.right,
      };
    });
    expect(gaps.left, `Replace dialog left gap at ${width} px`).toBeGreaterThanOrEqual(11.5);
    expect(gaps.rightInner, `Replace dialog right gap (innerWidth) at ${width} px`).toBeGreaterThanOrEqual(11.5);
    expect(gaps.rightClient, `Replace dialog right gap (clientWidth) at ${width} px`).toBeGreaterThanOrEqual(11.5);
  });
  await page.getByTestId('replace-counts-cancel').click();
  await expect(page.getByTestId('replace-counts-dialog')).not.toBeVisible();
});

test('at 320, 390, 768, 850 and 900 px every visible "?" shares a line with its preceding content', async ({ page }) => {
  test.setTimeout(240_000);
  await startWithSavedCounts562(page);
  await checkAllTipsAtWidths(page, 'WindowMasker -mk_counts run page', false);

  await openWmTab(page);
  await page.getByTestId('wm-stage-mk_counts').check();
  await openWmOptions(page);
  await checkAllTipsAtWidths(page, 'WindowMasker Make counts');

  await page.getByTestId('wm-stage-mask').check();
  await openWmOptions(page);
  await checkAllTipsAtWidths(page, 'WindowMasker Mask');

  await page.getByTestId('wm-stage-convert').check();
  await openWmOptions(page);
  await checkAllTipsAtWidths(page, 'WindowMasker Convert counts');

  await openBlastnParameters(page);
  await page.getByTestId('query-input').fill('NM_000001.1');
  await expect(page.getByTestId('query-fetch-row')).toBeVisible();
  await page.getByTestId('param-extra').fill('-conffile test.ini');
  await expect(page.getByTestId('param-conffile')).toBeVisible();
  await checkAllTipsAtWidths(page, 'BLASTN WindowMasker unchecked');

  await page.getByTestId('param-window_masker_db').check();
  await expect(page.getByTestId('param-window_masker_db-missing')).toBeVisible();
  await checkAllTipsAtWidths(page, 'BLASTN Taxonomy ID empty');

  await page.getByTestId('param-window_masker_db-taxid').fill('562');
  await expect(page.getByTestId('counts-note')).toContainText(/Made in this (?:browser|tab only)/);
  await checkAllTipsAtWidths(page, 'BLASTN Taxonomy ID made');

  await page.getByTestId('param-window_masker_db-taxid').fill('9606');
  await expect(page.getByTestId('counts-note')).toContainText('Downloads wmasker.obinary for 9606');
  await checkAllTipsAtWidths(page, 'BLASTN Taxonomy ID unsaved');

  await attachCountsFixture(page);
  await checkAllTipsAtWidths(page, 'BLASTN Counts file');
});

test('at 320, 390 and 768 px the WindowMasker row checkbox shares a line with its name', async ({ page }) => {
  test.setTimeout(180_000);
  await startWithSavedCounts562(page);
  await openBlastnParameters(page);
  await checkWmCheckboxAtNarrowWidths(page, 'BLASTN WindowMasker unchecked');

  await selectTaxidMade562(page);
  await checkWmCheckboxAtNarrowWidths(page, 'BLASTN Taxonomy ID made');

  await attachCountsFixture(page);
  await checkWmCheckboxAtNarrowWidths(page, 'BLASTN Counts file');
});
