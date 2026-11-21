// WA2-AX-12 (W7): the tabs that look like tabs are tabs in ARIA's pattern: Search / Results, the
// results' views (Hits, Dot Plot, Run details, Outputs) and the Outputs' formats. A tablist of
// tabs, the chosen one aria-selected and the only Tab stop of its list, each controlling a
// tabpanel that the chosen tab labels; Left and Right (around the ends), Home and End choose (the
// output formats: move the focus, and Enter or Space chooses).
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { focused } from './support/live';
import { BUILD_HAS_ENGINE } from './support/browser';
import { chooseRun, openSearch, paste, program, submit, waitStatus } from './support/search';

test.setTimeout(180_000);

/** Checks a tablist: its tabs (by test ID) in order, the chosen one, and the chosen one's panel. */
async function expectTabs(page: Page, name: string, testids: readonly string[], chosen: string): Promise<void> {
  const list = page.getByRole('tablist', { name });
  await expect(list.getByRole('tab')).toHaveCount(testids.length);
  for (const [i, testid] of testids.entries()) {
    const tab = list.getByRole('tab').nth(i);
    await expect(tab).toHaveAttribute('data-testid', testid);
    await expect(tab).toHaveAttribute('aria-selected', String(testid === chosen));
    await expect(tab).toHaveAttribute('tabindex', testid === chosen ? '0' : '-1');
  }
  const tab = page.getByTestId(chosen);
  const [id, panel] = [await tab.getAttribute('id'), await tab.getAttribute('aria-controls')];
  expect(id).toBeTruthy();
  const shown = page.locator(`[id="${panel}"]`);
  await expect(shown).toHaveAttribute('role', 'tabpanel');
  await expect(shown).toHaveAttribute('aria-labelledby', id!);
  await expect(shown).toBeVisible();
}

/** Presses a key on the focused tab; the tab chosen then has the focus. */
async function press(page: Page, key: string, chosen: string): Promise<void> {
  await page.keyboard.press(key);
  await expect(page.getByTestId(chosen)).toHaveAttribute('aria-selected', 'true');
  expect(await focused(page)).toBe(chosen);
}

const MAIN = ['tab-search', 'tab-results'];
const VIEWS = ['results-view-hits', 'pane-dotplot', 'results-view-details', 'results-view-outputs'];
const FORMATS = ['format-0', 'format-6', 'format-7'];
const PROTEIN = 'MKLVVLAAGGHHKLMKLVVLAAGGHHKLMSTNPKPQRKTKRNTNRRPQDVKFPGG';

test('Search and Results, the results views and the output formats are tabs that the arrow keys reach', async ({ page }) => {
  await openSearch(page);
  await expectTabs(page, 'Main', MAIN, 'tab-search');
  await page.getByTestId('tab-search').focus();
  await press(page, 'ArrowRight', 'tab-results');
  await expectTabs(page, 'Main', MAIN, 'tab-results');
  await expect(page.getByTestId('results-panel')).toBeVisible();
  await press(page, 'ArrowRight', 'tab-search');
  await press(page, 'End', 'tab-results');
  await press(page, 'Home', 'tab-search');
  await expect(page.getByTestId('search-panel')).toBeVisible();

  await program(page, 'blastp');
  await paste(page, 'query', `>q1\n${PROTEIN}\n`);
  await paste(page, 'subject', `>s1\n${PROTEIN}\n`);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await expectTabs(page, 'Main', MAIN, 'tab-results');
  await expectTabs(page, 'Results view', VIEWS, 'results-view-hits');
  await page.getByTestId('results-view-hits').focus();
  await press(page, 'ArrowRight', 'pane-dotplot');
  await expect(page.getByTestId('results-dotplot')).toBeVisible();
  await expectTabs(page, 'Results view', VIEWS, 'pane-dotplot');
  await press(page, 'End', 'results-view-outputs');
  await expectTabs(page, 'Results view', VIEWS, 'results-view-outputs');
  await press(page, 'ArrowRight', 'results-view-hits');
  await expect(page.getByTestId('results-hits')).toBeVisible();
  await press(page, 'ArrowLeft', 'results-view-outputs');
  await press(page, 'ArrowLeft', 'results-view-details');
  await expectTabs(page, 'Results view', VIEWS, 'results-view-details');
  await press(page, 'Home', 'results-view-hits');
  // One Tab stop: Tab leaves the list from the chosen tab, and Shift+Tab comes back to it.
  await page.keyboard.press('Tab');
  expect(VIEWS).not.toContain(await focused(page));
  await page.keyboard.press('Shift+Tab');
  expect(await focused(page)).toBe('results-view-hits');

  await page.getByTestId('results-view-outputs').click();
  await expectTabs(page, 'Output format', FORMATS, 'format-6');
  await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', '1:6');
  // The formats are chosen by Enter or Space (each output is read from storage): the keys move the focus only.
  const output = page.getByTestId('result-output');
  await page.getByTestId('format-6').focus();
  await page.keyboard.press('ArrowLeft');
  expect(await focused(page)).toBe('format-0');
  await expectTabs(page, 'Output format', FORMATS, 'format-6');
  await expect(output).toHaveAttribute('data-shown', '1:6');
  await page.keyboard.press('Enter');
  await expectTabs(page, 'Output format', FORMATS, 'format-0');
  await expect(output).toHaveAttribute('data-shown', '1:0');
  await page.keyboard.press('ArrowLeft');
  expect(await focused(page)).toBe('format-7');
  await page.keyboard.press('Home');
  expect(await focused(page)).toBe('format-0');
  await page.keyboard.press('End');
  expect(await focused(page)).toBe('format-7');
  await page.keyboard.press('Space');
  await expect(output).toHaveAttribute('data-shown', '1:7');
  await expectTabs(page, 'Output format', FORMATS, 'format-7');
});

// Line WQA, W-2 (2026-11-20): in the random walk at 320 px, WebKit's press on the Search tab came
// as the run on the Results tab ended, and the tab stayed unchosen. The page became shorter than
// where it was scrolled to ("Reading the results…" for a moment in place of the run's status), the
// browser moved it up by as much (106 px there), and the release came on what had been under the
// tab: a click goes where the press and the release both were. A cancelled or failed run, whose
// message is shorter than its status, and another run chosen did the same, in every browser. The
// page now keeps its height from a press until its click (ui/press-hold.ts).

/** Opens a page of 320 px for a run held in progress (the FakeEngine's phases wait on the page's timers). */
async function heldRun(page: Page): Promise<void> {
  await page.setViewportSize({ width: 320, height: 480 });
  await page.clock.install();
  await openSearch(page);
  await program(page, 'blastp');
  await paste(page, 'query', `>q1\n${PROTEIN}\n`);
  await paste(page, 'subject', `>s1\n${PROTEIN}\n>s2\n${PROTEIN}\n`);
  await page.clock.pauseAt((await page.evaluate(() => Date.now())) + 1000);
}

/**
 * Presses the Search tab while `change` makes the page shorter, then releases: the page has not
 * moved under the press, and the press was a click. The page is scrolled `scrollY`, or to its end;
 * it is laid out after each change, as a browser may do between any two of them, so a state that
 * lasts only until the next change moves the page all the same.
 */
async function pressSearchDuring(page: Page, change: () => Promise<void>, scrollY?: number): Promise<void> {
  if (scrollY === undefined) {
    // A window in which the page, scrolled to its end, still shows the tabs: scrolled half as far
    // as the tabs are from the top.
    const now = await page.evaluate(() => ({
      height: document.documentElement.scrollHeight,
      tab: document.querySelector('[data-testid="tab-search"]')!.getBoundingClientRect().top + window.scrollY,
    }));
    await page.setViewportSize({ width: 320, height: now.height - Math.round(now.tab / 2) });
  }
  const scrolled = await page.evaluate((y) => {
    window.scrollTo(0, y ?? document.documentElement.scrollHeight);
    const seen = new Set<number>();
    (window as unknown as { seenScroll: Set<number> }).seenScroll = seen;
    new MutationObserver(() => {
      void document.documentElement.scrollHeight;
      seen.add(window.scrollY);
    }).observe(document.body, { subtree: true, childList: true, attributes: true, characterData: true });
    return window.scrollY;
  }, scrollY);
  expect(scrolled).toBeGreaterThan(0);
  const tab = page.getByTestId('tab-search');
  const box = (await tab.boundingBox())!;
  await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
  await page.mouse.down();
  await change();
  expect(await page.evaluate(() => [...(window as unknown as { seenScroll: Set<number> }).seenScroll])).toEqual([scrolled]);
  await page.mouse.up();
  await expect(tab).toHaveAttribute('aria-selected', 'true');
}

test('a press on the Search tab while a run’s results are read is a click: the page does not move', async ({ page }) => {
  test.skip(BUILD_HAS_ENGINE, "the FakeEngine's runs wait on the page's timers, which the test holds");
  await heldRun(page);
  await submit(page);
  await expect(page.getByTestId('run-status')).toBeVisible();
  await pressSearchDuring(page, async () => {
    await page.clock.resume();
    await waitStatus(page, 1, 'completed');
    await expect(page.getByTestId('results-view-hits')).toBeVisible();
  });
});

test('a press on the Search tab while the run is cancelled is a click: its shorter message does not move the page', async ({ page }) => {
  test.skip(BUILD_HAS_ENGINE, "the FakeEngine's runs wait on the page's timers, which the test holds");
  await heldRun(page);
  await submit(page);
  await expect(page.getByTestId('run-status')).toBeVisible();
  await pressSearchDuring(page, async () => {
    await page.getByTestId('run-status-cancel').focus();
    await page.keyboard.press('Enter');
    await page.clock.resume();
    await waitStatus(page, 1, 'cancelled');
    await expect(page.getByTestId('results-status')).toHaveAttribute('data-run-status', 'cancelled');
  });
});

test('a press on the Search tab while another run is shown is a click: the shorter page does not move', async ({ page }) => {
  test.skip(BUILD_HAS_ENGINE, "the FakeEngine's runs wait on the page's timers, which the test holds");
  await heldRun(page);
  await submit(page);
  await page.getByTestId('run-status-cancel').click();
  await page.getByTestId('tab-search').click();
  await submit(page);
  await page.clock.resume();
  await waitStatus(page, 1, 'cancelled');
  await waitStatus(page, 2, 'completed');
  await expect(page.getByTestId('results-view-hits')).toBeVisible();
  // The results of run 2 are taller than the window; run 1 shows one line in their place.
  await page.setViewportSize({ width: 320, height: 1000 });
  const tab = await page.getByTestId('tab-search').evaluate((element) => element.getBoundingClientRect().top + window.scrollY);
  const run1 = await page.getByTestId('run-1-option').getAttribute('value');
  await pressSearchDuring(
    page,
    async () => {
      await page.getByTestId('results-run').selectOption(run1!);
      await expect(page.getByTestId('results-status')).toHaveAttribute('data-run-status', 'cancelled');
    },
    Math.round(tab / 2),
  );
});

// Review r2 of W-2, M1 (2026-11-21): a press whose release the page never saw (a right press, whose release a context
// menu takes) held until the next press: the Results tab, chosen with the keyboard after it, kept the taller Search
// form's height below it. Only a press of the main button holds now, and a lost release ends it at the next move.
test('a right press whose release never comes leaves the page its height', async ({ page }) => {
  await openSearch(page);
  const main = page.locator('main');
  const height = async () => (await main.boundingBox())!.height;
  await page.getByTestId('tab-results').click();
  const results = await height();
  await page.getByTestId('tab-search').click();
  expect(await height(), 'the Search form is taller than the Results tab with no run').toBeGreaterThan(results + 100);
  const form = (await page.getByTestId('search-panel').boundingBox())!;
  await page.mouse.move(form.x + form.width / 2, form.y + form.height / 2);
  await page.mouse.down({ button: 'right' });
  // The release does not come; the keyboard goes on.
  await page.getByTestId('tab-search').focus();
  await press(page, 'ArrowRight', 'tab-results');
  await expect.poll(height).toBeCloseTo(results, 0);
});
