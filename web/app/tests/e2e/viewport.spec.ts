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
const COUNTS_FIXTURE = 'LOSAT/tests/winmask/counts/mg1655_u9.obinary';

interface SentenceTipMeasurement {
  readonly id: string;
  readonly btnLeft: number;
  readonly btnTop: number;
  readonly btnBottom: number;
  readonly lastRight: number;
  readonly lastTop: number;
  readonly lastBottom: number;
  readonly overlapsY: boolean;
  readonly followsX: boolean;
}

/**
 * Checks every visible "?" (HelpTip) that follows a sentence on the current view: the button sits on
 * the same line box as the last line of the sentence (or its trailing inline link, such as Delete),
 * never orphaned alone on the next line.
 */
async function measureSentenceHelpTips(page: Page): Promise<SentenceTipMeasurement[]> {
  return page.evaluate(() => {
    const tips = new Set<Element>();
    for (const selector of [
      '.hint-line .help-tip',
      '.param-hint-line .help-tip',
      '.param-control .param-hint ~ .help-tip',
      '.counts-line .help-tip',
      '.fetch-row .help-tip',
      '[data-kind="conffile"] .help-tip',
    ]) {
      for (const el of document.querySelectorAll(selector)) tips.add(el);
    }
    const results: SentenceTipMeasurement[] = [];
    for (const tip of tips) {
      const btn = tip.querySelector('button');
      if (btn === null) continue;
      const br = btn.getBoundingClientRect();
      if (br.width === 0 || br.height === 0) continue;

      const wrapper = tip.closest('.hint-line, .param-hint-line');
      let startNode: Node | null = null;
      if (wrapper !== null) {
        startNode = wrapper.firstChild;
      } else {
        const parent = tip.parentElement;
        if (parent !== null) {
          startNode =
            parent.querySelector('.param-hint') ??
            parent.querySelector('.hint') ??
            parent.querySelector('.file-missing') ??
            parent.querySelector('.file-name') ??
            tip.previousElementSibling;
        }
      }
      if (startNode === null || startNode === tip) continue;

      const range = document.createRange();
      range.setStartBefore(startNode);
      range.setEndBefore(tip);
      const last = [...range.getClientRects()].filter((r) => r.width > 0.5 && r.height > 0.5).at(-1);
      if (last === undefined) continue;

      results.push({
        id: btn.getAttribute('data-testid') ?? btn.getAttribute('aria-label') ?? 'help-tip',
        btnLeft: br.left,
        btnTop: br.top,
        btnBottom: br.bottom,
        lastRight: last.right,
        lastTop: last.top,
        lastBottom: last.bottom,
        overlapsY: Math.min(br.bottom, last.bottom) > Math.max(br.top, last.top),
        followsX: br.left >= last.right - 1,
      });
    }
    return results;
  });
}

function expectSentenceHelpTipsInline(tips: readonly SentenceTipMeasurement[], context: string): void {
  expect(tips.length, `${context}: expected at least one sentence HelpTip`).toBeGreaterThan(0);
  for (const tip of tips) {
    expect(tip.overlapsY, `${context}: ${tip.id} should overlap the sentence's last line vertically (btn ${tip.btnTop.toFixed(1)}..${tip.btnBottom.toFixed(1)} vs last ${tip.lastTop.toFixed(1)}..${tip.lastBottom.toFixed(1)})`).toBe(true);
    expect(tip.followsX, `${context}: ${tip.id} should sit after the sentence's last line (btnLeft ${tip.btnLeft.toFixed(1)} vs lastRight ${tip.lastRight.toFixed(1)})`).toBe(true);
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

async function forNarrowWidths(page: Page, check: (width: (typeof NARROW_WIDTHS)[number]) => Promise<void>): Promise<void> {
  for (const width of NARROW_WIDTHS) {
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

const checkTipsAtNarrowWidths = (page: Page, label: string): Promise<void> =>
  forNarrowWidths(page, async (width) => expectSentenceHelpTipsInline(await measureSentenceHelpTips(page), `${label} at ${width} px`));

const checkWmCheckboxAtNarrowWidths = (page: Page, label: string): Promise<void> =>
  forNarrowWidths(page, (width) => expectWmCheckboxSharesLine(page, `${label} at ${width} px`));

test('at 320, 390 and 768 px the Replace dialog leaves at least 12 px on each side', async ({ page }) => {
  test.setTimeout(180_000);
  await startWithSavedCounts562(page);
  await saveMadeCounts562(page, 2);
  await expect(page.getByTestId('replace-counts-dialog')).toBeVisible();

  await forNarrowWidths(page, async (width) => {
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

test('at 320, 390 and 768 px every "?" after a sentence stays on the sentence\'s last line', async ({ page }) => {
  test.setTimeout(180_000);
  await startWithSavedCounts562(page);

  await openWmTab(page);
  await page.getByTestId('wm-stage-mask').check();
  const wmToggle = page.getByTestId('wm-options-toggle');
  if ((await wmToggle.getAttribute('aria-expanded')) !== 'true') await wmToggle.click();
  await checkTipsAtNarrowWidths(page, 'WindowMasker Mask');

  await openBlastnParameters(page);
  await page.getByTestId('query-input').fill('NM_000001.1');
  await expect(page.getByTestId('query-fetch-row')).toBeVisible();
  await page.getByTestId('param-extra').fill('-conffile test.ini');
  await expect(page.getByTestId('param-conffile')).toBeVisible();
  await selectTaxidMade562(page);
  await checkTipsAtNarrowWidths(page, 'BLASTN Taxonomy ID made');

  await page.getByTestId('param-window_masker_db-taxid').fill('9606');
  await expect(page.getByTestId('counts-note')).toContainText('Downloads wmasker.obinary for 9606');
  await checkTipsAtNarrowWidths(page, 'BLASTN Taxonomy ID unsaved');

  await attachCountsFixture(page);
  await checkTipsAtNarrowWidths(page, 'BLASTN Counts file');
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
