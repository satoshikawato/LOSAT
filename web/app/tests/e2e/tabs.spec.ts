// WA2-AX-12 (W7): the tabs that look like tabs are tabs in ARIA's pattern: Search / Results, the
// results' views (Hits, Dot Plot, Run details, Outputs) and the Outputs' formats. A tablist of
// tabs, the chosen one aria-selected and the only Tab stop of its list, each controlling a
// tabpanel that the chosen tab labels; Left and Right (around the ends), Home and End choose (the
// output formats: move the focus, and Enter or Space chooses).
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { focused } from './support/live';
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
