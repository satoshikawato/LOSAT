// A press of a control is a click on it even if the page moves under it (line WQA's D27, 2026-11-27). A click goes where
// the press and the release both were. Content that came, grew or went above a control between its press and its
// release moved the control from under the pointer: the release came beside it, and the press did nothing and said
// nothing. The query card's warning that its check adds as it ends moved Run LOSAT down 34 px; so did any other line
// above a control, in every browser (WebKit has no scroll anchoring, and Chromium's and Firefox's do not act at the
// page's top or for a change between the line they keep and the control). From a press until its click, the page now
// keeps the pressed control where it was under the pointer (ui/press-hold.ts). The tests force the timing: a warning
// comes, or goes, above the control between the press and the release.
import { expect, type Locator, type Page } from '@playwright/test';
import { test } from './support/test';
import { chooseRun, openSearch, paste, program, submit, waitStatus } from './support/search';

test.setTimeout(120_000);

const PROTEIN = 'MKLVVLAAGGHHKLMKLVVLAAGGHHKLMSTNPKPQRKTKRNTNRRPQDVKFPGG';
const WIDTHS = [320, 1280] as const;

/** The search form with a query and a subject, the window `width` px wide and 720 px high. */
async function form(page: Page, width: number): Promise<void> {
  await page.setViewportSize({ width, height: 720 });
  await openSearch(page);
  await program(page, 'blastp');
  await paste(page, 'query', `>q1\n${PROTEIN}\n`);
  await paste(page, 'subject', `>s1\n${PROTEIN}\n`);
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
}

/**
 * Puts a warning, as an input's card shows one, just before the element of `selector` (`before`) or as its last line
 * (`end`), and returns its height.
 */
async function warn(page: Page, selector: string, where: 'before' | 'end'): Promise<number> {
  return page.evaluate(
    ([selector, where]) => {
      const warning = document.createElement('p');
      warning.className = 'warning';
      warning.dataset.testid = 'press-warning';
      warning.textContent = 'This input ends with a line that the command line cannot read.';
      const at = document.querySelector(selector)!;
      if (where === 'before') at.before(warning);
      else at.append(warning);
      return warning.getBoundingClientRect().height;
    },
    [selector, where] as const,
  );
}

/** Takes the warning of `warn` away. */
async function unwarn(page: Page): Promise<void> {
  await page.evaluate(() => document.querySelector('[data-testid="press-warning"]')?.remove());
}

/** Whether the point is on `control`, as the page is laid out now. */
async function under(control: Locator, at: { x: number; y: number }): Promise<boolean> {
  return control.evaluate((element, at) => element.contains(document.elementFromPoint(at.x, at.y)), at);
}

/**
 * Presses the middle of `control`, has `change` move it, and releases there. Before the release, the control is still
 * under the pointer. Returns the window's scroll at the press and at the release.
 */
async function pressWhile(page: Page, control: Locator, change: () => Promise<void>): Promise<{ press: number; release: number }> {
  const box = (await control.boundingBox())!;
  const at = { x: box.x + box.width / 2, y: box.y + box.height / 2 };
  await page.mouse.move(at.x, at.y);
  await page.mouse.down();
  const press = await page.evaluate(() => window.scrollY);
  await change();
  expect.soft(await under(control, at), 'the control is under the pointer at the release').toBe(true);
  const release = await page.evaluate(() => window.scrollY);
  await page.mouse.up();
  return { press, release };
}

/** Run LOSAT started a run: the page went to its results. */
async function started(page: Page): Promise<void> {
  await expect(page.getByTestId('tab-results')).toHaveAttribute('aria-selected', 'true');
  await expect(page.getByTestId('run-1-option')).toHaveCount(1);
}

for (const width of WIDTHS) {
  test.describe(`${width} px`, () => {
    test('a press on Run LOSAT is a click while a warning comes just above it', async ({ page }) => {
      await form(page, width);
      const button = page.getByTestId('run-losat');
      await button.evaluate((element) => element.scrollIntoView({ block: 'center' }));
      const height = await pressWhile(page, button, async () => {
        expect(await warn(page, '.blast-row:has([data-testid="run-losat"])', 'before')).toBeGreaterThan(20);
      });
      expect(height.press, 'the page is scrolled to the button').toBeGreaterThan(0);
      await started(page);
    });

    test('a press on Run LOSAT at the page’s top is a click while the query card’s warning comes', async ({ page }) => {
      await form(page, width);
      // A window high enough to show the button with the page at its top.
      const bottom = await page.getByTestId('run-losat').evaluate((element) => element.getBoundingClientRect().bottom + window.scrollY);
      await page.setViewportSize({ width, height: Math.ceil(bottom) + 40 });
      await page.evaluate(() => window.scrollTo(0, 0));
      const scroll = await pressWhile(page, page.getByTestId('run-losat'), async () => {
        await warn(page, '[data-testid="query-source-0"]', 'end');
      });
      expect(scroll.press).toBe(0);
      await started(page);
    });

    test('a press on Run LOSAT is a click while a warning comes above it in a window that holds the whole page', async ({ page }) => {
      await form(page, width);
      const whole = await page.evaluate(() => document.documentElement.scrollHeight);
      await page.setViewportSize({ width, height: whole + 100 });
      expect(await page.evaluate(() => document.documentElement.scrollHeight - window.innerHeight), 'the page does not scroll').toBeLessThanOrEqual(0);
      await pressWhile(page, page.getByTestId('run-losat'), async () => {
        await warn(page, '[data-testid="query-source-0"]', 'end');
      });
      await started(page);
    });

    test('a press on Run LOSAT at the page’s top is a click while a warning above it goes', async ({ page }) => {
      await form(page, width);
      await warn(page, '[data-testid="query-source-0"]', 'end');
      const bottom = await page.getByTestId('run-losat').evaluate((element) => element.getBoundingClientRect().bottom + window.scrollY);
      await page.setViewportSize({ width, height: Math.ceil(bottom) + 40 });
      await page.evaluate(() => window.scrollTo(0, 0));
      const scroll = await pressWhile(page, page.getByTestId('run-losat'), () => unwarn(page));
      expect(scroll.press).toBe(0);
      await started(page);
    });

    test('a press on an output format is a click while a warning comes above the formats', async ({ page }) => {
      await form(page, width);
      await submit(page);
      await waitStatus(page, 1, 'completed');
      await chooseRun(page, 1);
      await page.getByTestId('results-view-outputs').click();
      await expect(page.getByTestId('format-6')).toHaveAttribute('aria-selected', 'true');
      const format = page.getByTestId('format-0');
      await format.evaluate((element) => element.scrollIntoView({ block: 'center' }));
      await pressWhile(page, format, async () => {
        await warn(page, '[role="tablist"][aria-label="Output format"]', 'before');
      });
      await expect(format).toHaveAttribute('aria-selected', 'true');
      await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', '1:0');
    });

    test('a press on the Results tab is a click while a warning comes above the tabs', async ({ page }) => {
      await form(page, width);
      await page.evaluate(() => window.scrollTo(0, 0));
      const tab = page.getByTestId('tab-results');
      await pressWhile(page, tab, async () => {
        await warn(page, 'nav', 'before');
      });
      await expect(tab).toHaveAttribute('aria-selected', 'true');
    });

    // Review r1 of D27, M1: the page keeps the pressed element's top in place, not its middle, so a block pressed on its
    // blank space that grows at its end (the query card's warning comes in it) moves nothing, as without the keeping.
    test('a press on the query block’s blank space does not move the page while the query card’s warning comes in it', async ({ page }) => {
      await form(page, width);
      // A window high enough to show the query card with the page at its top.
      const bottom = await page.getByTestId('query-source-0').evaluate((card) => card.getBoundingClientRect().bottom + window.scrollY);
      await page.setViewportSize({ width, height: Math.ceil(bottom) + 80 });
      await page.evaluate(() => window.scrollTo(0, 0));
      const panel = page.getByTestId('query-panel');
      // A point in the block's space at its left, beside its card: the press is on the block itself.
      const y = await page.getByTestId('query-source-0').evaluate((card) => card.getBoundingClientRect().top + card.getBoundingClientRect().height / 2);
      const at = { x: (await panel.evaluate((element) => element.getBoundingClientRect().left)) + 4, y };
      expect(await panel.evaluate((element, at) => document.elementFromPoint(at.x, at.y) === element, at)).toBe(true);
      const top = await panel.evaluate((element) => element.getBoundingClientRect().top);
      await page.mouse.move(at.x, at.y);
      await page.mouse.down();
      await warn(page, '[data-testid="query-source-0"]', 'end');
      expect(await page.evaluate(() => window.scrollY), 'the page does not move').toBe(0);
      expect(await panel.evaluate((element) => element.getBoundingClientRect().top)).toBe(top);
      await page.mouse.up();
    });

    // Review r1 of D27, M2: the space put above at the page's top for a warning that went is taken back when it comes
    // back, and nothing moves at the click.
    test('a press at the page’s top is a click, and nothing moves at it, while a warning above goes and comes back', async ({ page }) => {
      await form(page, width);
      await warn(page, '[data-testid="query-source-0"]', 'end');
      const toggle = page.getByTestId('algorithm-parameters-toggle');
      const bottom = await toggle.evaluate((element) => element.getBoundingClientRect().bottom + window.scrollY);
      await page.setViewportSize({ width, height: Math.ceil(bottom) + 40 });
      await page.evaluate(() => window.scrollTo(0, 0));
      // Its middle, which its pressed look (a scale about its middle) does not move.
      const middle = () => toggle.evaluate((element) => element.getBoundingClientRect().top + element.getBoundingClientRect().height / 2);
      const before = await middle();
      const scroll = await pressWhile(page, toggle, async () => {
        await unwarn(page);
        await warn(page, '[data-testid="query-source-0"]', 'end');
      });
      expect(scroll).toEqual({ press: 0, release: 0 });
      await expect(toggle).toHaveAttribute('aria-expanded', 'true');
      expect(await page.evaluate(() => window.scrollY)).toBe(0);
      expect(await middle(), 'the control is where it was').toBeCloseTo(before, 1);
    });

    // A press with the keys goes to the focused control wherever it is: the page does not move for it, as it does not
    // when nothing is pressed.
    test('Space on Run LOSAT starts the run while a warning comes above it, and the page does not move', async ({ page }) => {
      await form(page, width);
      const button = page.getByTestId('run-losat');
      await button.evaluate((element) => element.scrollIntoView({ block: 'center' }));
      const before = await page.evaluate(() => window.scrollY);
      await warn(page, '.blast-row:has([data-testid="run-losat"])', 'before');
      expect(await page.evaluate(() => window.scrollY), 'no press: the page does not move').toBe(before);
      await unwarn(page);
      await button.evaluate((element: HTMLElement) => element.focus({ preventScroll: true }));
      await page.keyboard.down(' ');
      await warn(page, '.blast-row:has([data-testid="run-losat"])', 'before');
      expect(await page.evaluate(() => window.scrollY), 'a press with the keys: the page does not move').toBe(before);
      await page.keyboard.up(' ');
      await started(page);
    });
  });
}
