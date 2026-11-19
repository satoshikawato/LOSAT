// Back and Forward between the Search and Results tabs (line WQA's sweep, D6, 2026-11-19). NCBI
// BLAST's flow, which "Run LOSAT" follows, returns from the results to the form on Back. Before,
// the tabs made no history entry: Back on the results left LOSAT ("Leave site?"; leaving lost
// every run and input) and Forward opened an empty LOSAT. The Results tab now has one entry
// above the Search tab's, at the same URL (ui/tab-history.ts).
import { expect, type Dialog, type Page } from '@playwright/test';
import { test } from './support/test';
import { navigate, reload } from './support/navigate';
import { described, openSearch, paste, program, submit, waitStatus } from './support/search';

test.setTimeout(120_000);

const PROTEIN = 'MKLVVLAAGGHHKLMKLVVLAAGGHHKLMSTNPKPQRKTKRNTNRRPQDVKFPGG';
const QUERY = `>q1\n${PROTEIN}\n`;
const TITLE = 'Back and Forward';

/** The dialogs that the page asked for (the leave guard's `beforeunload`), each accepted. */
function dialogs(page: Page): string[] {
  const seen: string[] = [];
  page.on('dialog', (dialog: Dialog) => {
    seen.push(dialog.type());
    void dialog.accept();
  });
  return seen;
}

async function expectTab(page: Page, tab: 'search' | 'results'): Promise<void> {
  await expect(page.getByTestId(`tab-${tab}`)).toHaveAttribute('aria-selected', 'true');
  await expect(page.getByTestId(`${tab}-panel`)).toBeVisible();
}

/** A page of the site before LOSAT, so that Back from LOSAT has somewhere to go. */
async function openAfterAnotherPage(page: Page): Promise<void> {
  await navigate(page, '/favicon.svg');
  await openSearch(page);
}

const scrollY = (page: Page) => page.evaluate(() => Math.round(window.scrollY));

/** Runs one search; returns how far the form was scrolled when "Run LOSAT" was pressed. */
async function runOnce(page: Page): Promise<number> {
  await program(page, 'blastp');
  await page.getByTestId('job-title').fill(TITLE);
  await paste(page, 'query', QUERY);
  await paste(page, 'subject', QUERY);
  await page.getByTestId('run-losat').scrollIntoViewIfNeeded();
  const formScroll = await scrollY(page);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await expectTab(page, 'results');
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
  return formScroll;
}

test('Back on the results returns to the form as it was, Forward shows the results again', async ({ page }) => {
  const asked = dialogs(page);
  await openAfterAnotherPage(page);
  const url = page.url();
  const formScroll = await runOnce(page);
  expect(formScroll).toBeGreaterThan(0);
  expect(page.url()).toBe(url);
  // Far down the results: Forward shows them from their heading again (their views read their content
  // after they mount, so a position far down could not be met).
  await page.evaluate(() => window.scrollTo(0, document.documentElement.scrollHeight));

  // The form as it was: its fields, where it was scrolled, and the focus in it or on its tab (never lost).
  await page.goBack();
  await expectTab(page, 'search');
  await expect(page.getByTestId('query-input')).toHaveValue(QUERY);
  await expect(page.getByTestId('job-title')).toHaveValue(TITLE);
  await expect(page.getByTestId('program-blastp')).toBeChecked();
  await expect.poll(() => scrollY(page)).toBe(formScroll);
  expect(
    await page.evaluate(() => {
      const active = document.activeElement;
      return active?.closest('[data-testid=search-panel]') !== null || active?.getAttribute('data-testid') === 'tab-search';
    }),
  ).toBe(true);
  expect(page.url()).toBe(url);
  expect(asked).toEqual([]);

  await page.goForward();
  await expectTab(page, 'results');
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
  await expect(page.getByTestId('results-heading')).toBeFocused();
  await expect(page.getByTestId('results-heading')).toBeInViewport();
  expect(asked).toEqual([]);

  // From the Search tab, Back leaves LOSAT, and the page asks first (the leave guard), as before.
  await page.goBack();
  await expectTab(page, 'search');
  await page.goBack();
  await expect.poll(() => asked).toEqual(['beforeunload']);
  await expect(page).toHaveURL(/\/favicon\.svg$/);
});

test('the tab buttons and the arrow keys keep one Results entry, and a reload opens the Search tab', async ({ page }) => {
  const asked = dialogs(page);
  await openAfterAnotherPage(page);
  await runOnce(page);
  const length = await page.evaluate(() => history.length);

  // The buttons and the arrow keys switch the tabs as before, without piling up entries, and the
  // Search tab's button, which goes back to the form's entry, does not scroll the page.
  await page.getByTestId('tab-search').scrollIntoViewIfNeeded();
  const before = await scrollY(page);
  await page.getByTestId('tab-search').click();
  await expectTab(page, 'search');
  await expect.poll(() => page.evaluate(() => history.state)).toBeNull();
  expect(await scrollY(page)).toBe(before);
  await page.getByTestId('tab-results').click();
  await expectTab(page, 'results');
  await page.getByTestId('tab-results').focus();
  await page.keyboard.press('ArrowLeft');
  await expectTab(page, 'search');
  await expect(page.getByTestId('tab-search')).toBeFocused();
  await page.keyboard.press('ArrowRight');
  await expectTab(page, 'results');
  await expect(page.getByTestId('tab-results')).toBeFocused();
  await page.keyboard.press('ArrowRight');
  await expectTab(page, 'search');
  await expect.poll(() => page.evaluate(() => history.length)).toBe(length);
  // The Search tab chosen by its button: Forward shows the results, Back the form.
  await page.goForward();
  await expectTab(page, 'results');
  await page.goBack();
  await expectTab(page, 'search');
  await expect(page.getByTestId('query-input')).toHaveValue(QUERY);
  expect(asked).toEqual([]);

  // A reload on the Results tab opens the page as a first visit does: the Search tab, nothing kept (REQ-20).
  await page.goForward();
  await expectTab(page, 'results');
  const url = page.url();
  await reload(page);
  await described(page);
  expect(page.url()).toBe(url);
  expect(asked).toEqual(['beforeunload']);
  await expectTab(page, 'search');
  await expect(page.getByTestId('query-input')).toHaveValue('');
  await page.getByTestId('tab-results').click();
  await expectTab(page, 'results');
  // "Run a search" of the empty Results tab goes to the form, and Forward comes back.
  await page.getByTestId('results-run-a-search').click();
  await expectTab(page, 'search');
  await page.goForward();
  await expectTab(page, 'results');
  await expect(page.getByTestId('results-run-a-search')).toBeVisible();
  expect(asked).toEqual(['beforeunload']);
});

// Review of 6fff89265: a reload kept the Results entry as a second Search entry, so the first Back
// after it did nothing, and Forward after Back and a reload showed an empty Results tab.
test('after a reload, Back leaves at once and Forward never opens the results that the reload took', async ({ page }) => {
  const asked = dialogs(page);
  await openAfterAnotherPage(page);
  await runOnce(page);
  await reload(page);
  await described(page);
  expect(asked).toEqual(['beforeunload']);
  await expectTab(page, 'search');
  await expect.poll(() => page.evaluate(() => history.state)).toBeNull();
  await page.goBack();
  await expect(page).toHaveURL(/\/favicon\.svg$/);
  expect(asked).toEqual(['beforeunload']);

  await openSearch(page);
  await runOnce(page);
  await page.goBack();
  await expectTab(page, 'search');
  await reload(page);
  await described(page);
  await page.goForward();
  await expect.poll(() => page.evaluate(() => history.state)).toBeNull();
  await expectTab(page, 'search');
  await expect(page.getByTestId('results-run-a-search')).toHaveCount(0);
  await page.goBack();
  await expect(page).toHaveURL(/\/favicon\.svg$/);
});
