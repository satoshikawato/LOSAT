// The Dot Plot under a finger (W9: WA2-AX-10).
//
// The touch screen is a phone's: hasTouch and a phone's width. A finger that moves is
// Chromium's DevTools Protocol touch events, which also scroll the page as a phone does; Firefox
// and WebKit have no such events here, so there the finger is pointer events of type "touch" sent
// to the plot, which check the plot's own handling (the browser's scrolling is Chromium's test).
// A pinch of two fingers cannot be made here (CDP's pinch gesture did not zoom the page even over
// text, W5): the plot leaves it to the page with touch-action, which the tests read.
import { expect, type CDPSession, type Locator, type Page } from '@playwright/test';
import { test } from './support/test';
import { chooseRun, openSearch, paste, program, submit, waitStatus } from './support/search';

test.setTimeout(180_000);

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

/** Pseudo-random nucleotides (a linear congruential generator, as results.spec.ts's): the same letters in every run. */
function dna(seed: number, length: number): string {
  let state = seed;
  let out = '';
  for (let i = 0; i < length; i++) {
    state = (Math.imul(state, 1103515245) + 12345) >>> 0;
    out += 'ACGT'[(state >>> 16) % 4];
  }
  return out;
}

/** Runs the search set in the form and opens its results. */
async function runAndOpen(page: Page): Promise<void> {
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
}

async function showTab(page: Page, id: 'pane-dotplot' | 'results-view-hits'): Promise<void> {
  const tab = page.getByTestId(id);
  await tab.click();
  await expect(tab).toHaveAttribute('aria-selected', 'true');
}

/** Puts an element in the middle of the window. */
async function centre(locator: Locator): Promise<void> {
  await locator.evaluate((element) => element.scrollIntoView({ block: 'center' }));
  await locator.page().evaluate(() => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve))));
}

const view = async (canvas: Locator): Promise<number[]> => (await canvas.getAttribute('data-view'))!.split(',').map(Number);
const scrollY = (page: Page): Promise<number> => page.evaluate(() => window.scrollY);

/** The query and subject of the dot plot tests: one long diagonal HSP that a zoomed-in view can move along. */
const LONG = dna(51, 600);

/** A BLASTN plot of LONG against itself, zoomed in twice, in the middle of the window. */
async function zoomedPlot(page: Page): Promise<Locator> {
  await program(page, 'blastn');
  await paste(page, 'query', `>q\n${LONG}\n`);
  await paste(page, 'subject', `>s\n${LONG}\n`);
  await runAndOpen(page);
  await showTab(page, 'pane-dotplot');
  const canvas = page.getByTestId('dotplot-canvas');
  await expect(canvas).toHaveAttribute('data-segments', /^[1-9]/);
  const full = (await canvas.getAttribute('data-view'))!;
  await page.getByTestId('dotplot-zoom-in').click();
  await page.getByTestId('dotplot-zoom-in').click();
  await expect(canvas).not.toHaveAttribute('data-view', full);
  await centre(canvas);
  return canvas;
}

// --- WA2-AX-10: a finger moving up or down scrolls the page, not the plot ------------------------------

test.describe('a phone', () => {
  test.use({ hasTouch: true, viewport: { width: 412, height: 839 } });

  test('a finger moving up or down over the dot plot scrolls the page and leaves the plot; one moving sideways moves the plot (touch events)', async ({
    page,
    browserName,
  }) => {
    test.skip(browserName !== 'chromium', "touch events that scroll the page are Chromium's DevTools Protocol's");
    const canvas = await zoomedPlot(page);
    const cdp = await page.context().newCDPSession(page);

    // Up, in steps of 10 px and of 3 px (W5: both moved the plot by 7 to 10 % of its view), and down.
    for (const [dy, steps] of [
      [-150, 10],
      [-150, 50],
      [150, 50],
    ] as const) {
      await centre(canvas);
      const [before, top] = [await view(canvas), await scrollY(page)];
      const box = (await canvas.boundingBox())!;
      await touchDrag(page, cdp, box.x + box.width / 2, box.y + box.height / 2, 0, dy, steps);
      // The page scrolled; the plot did not move.
      await expect.poll(() => scrollY(page)).not.toBe(top);
      expect(Math.sign((await scrollY(page)) - top), `${dy} px in ${steps} steps`).toBe(-Math.sign(dy));
      expect(await view(canvas), `${dy} px in ${steps} steps`).toEqual(before);
      await expect(page.getByTestId('dotplot-popup')).toHaveCount(0);
    }

    // Sideways: the plot moves along the query, and the page stays.
    await centre(canvas);
    const [before, top] = [await view(canvas), await scrollY(page)];
    const box = (await canvas.boundingBox())!;
    await touchDrag(page, cdp, box.x + box.width / 2, box.y + box.height / 2, -100, 0, 10);
    await expect.poll(async () => (await view(canvas))[0]).toBeGreaterThan(before[0]!);
    const after = await view(canvas);
    expect(after.slice(2)).toEqual(before.slice(2));
    expect(await scrollY(page)).toBe(top);
    // The browser takes a finger moving up or down, and two fingers (a pinch zooms the page).
    await expect(canvas).toHaveCSS('touch-action', 'pan-y pinch-zoom');
  });

  test('a finger moving up or down over the dot plot leaves the plot; one moving sideways moves it (pointer events of a finger)', async ({
    page,
  }) => {
    const canvas = await zoomedPlot(page);
    const id = await activePointer(page, canvas);
    const start = async () => {
      const box = (await canvas.boundingBox())!;
      return { x: box.width / 2, y: box.height / 2 };
    };
    // Up and down, straight and slanted (a finger is never straight): the browser scrolls the page
    // (touch-action pan-y); the plot does not move, before the browser takes the finger or after.
    for (const [dx, dy] of [
      [0, -3],
      [0, 3],
      [1, -3],
      [-2, 3],
    ] as const) {
      const before = await view(canvas);
      const { x, y } = await start();
      await finger(canvas, id, Array.from({ length: 31 }, (_, i) => [x + i * dx, y + i * dy]));
      await page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve))));
      expect(await view(canvas), `${dx}, ${dy} px a step`).toEqual(before);
      await expect(page.getByTestId('dotplot-popup')).toHaveCount(0);
    }
    // Sideways, also slanted: the plot moves with the finger.
    for (const [dx, dy] of [
      [-3, 0],
      [-3, 1],
    ] as const) {
      const before = await view(canvas);
      const { x, y } = await start();
      await finger(canvas, id, Array.from({ length: 31 }, (_, i) => [x + i * dx, y + i * dy]));
      await expect.poll(async () => (await view(canvas))[0], `${dx}, ${dy} px a step`).toBeGreaterThan(before[0]!);
    }
  });

});

/** A finger of the DevTools Protocol: down at (x, y) of the window, moved by (dx, dy) in `steps` steps a frame apart, and up. */
async function touchDrag(page: Page, cdp: CDPSession, x: number, y: number, dx: number, dy: number, steps: number): Promise<void> {
  await cdp.send('Input.dispatchTouchEvent', { type: 'touchStart', touchPoints: [{ x, y }] });
  for (let i = 1; i <= steps; i++) {
    await cdp.send('Input.dispatchTouchEvent', { type: 'touchMove', touchPoints: [{ x: x + (dx * i) / steps, y: y + (dy * i) / steps }] });
    await page.waitForTimeout(16);
  }
  await cdp.send('Input.dispatchTouchEvent', { type: 'touchEnd', touchPoints: [] });
  await page.waitForTimeout(400);
}

/**
 * The ID of a pointer that the browser has (its mouse): a pointer event sent by a script for a
 * pointer the browser does not have makes setPointerCapture throw.
 */
async function activePointer(page: Page, canvas: Locator): Promise<number> {
  const seen = canvas.evaluate(
    (element) => new Promise<number>((resolve) => element.addEventListener('pointermove', (event) => resolve((event as PointerEvent).pointerId), { once: true })),
  );
  const box = (await canvas.boundingBox())!;
  await page.mouse.move(box.x + 2, box.y + 2);
  await page.mouse.move(box.x + 4, box.y + 3);
  return seen;
}

/** A finger as pointer events of type "touch" sent to the canvas: down at the first point (in CSS pixels of the canvas), moved through the others, up at the last. */
async function finger(canvas: Locator, pointerId: number, points: readonly (readonly [number, number])[]): Promise<void> {
  await canvas.evaluate(
    (element, [id, path]) => {
      const rect = element.getBoundingClientRect();
      const send = (type: string, [x, y]: readonly [number, number], buttons: number) =>
        element.dispatchEvent(
          new PointerEvent(type, {
            pointerId: id,
            pointerType: 'touch',
            isPrimary: true,
            clientX: rect.left + x,
            clientY: rect.top + y,
            buttons,
            bubbles: true,
            cancelable: true,
          }),
        );
      send('pointerdown', path[0]!, 1);
      for (const point of path.slice(1)) send('pointermove', point, 1);
      send('pointerup', path.at(-1)!, 0);
    },
    [pointerId, points] as const,
  );
}
