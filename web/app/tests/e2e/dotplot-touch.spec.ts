// The Dot Plot and the Graphic Summary under a finger (W9: WA2-AX-10, WA2-AX-11), the HSP table's
// Identity and Orientation columns (WA2-OP-01) and the numbers on the axes of a short span
// (WA2-OP-02).
//
// The touch screen is a phone's: hasTouch and a phone's width. A tap is Playwright's touchscreen
// tap, which gives pointer events of type "touch" in all three browsers. A finger that moves is
// Chromium's DevTools Protocol touch events, which also scroll the page as a phone does; Firefox
// and WebKit have no such events here, so there the finger is pointer events of type "touch" sent
// to the plot, which check the plot's own handling (the browser's scrolling is Chromium's test).
// A pinch of two fingers cannot be made here (CDP's pinch gesture did not zoom the page even over
// text, W5): the plot leaves it to the page with touch-action, which the tests read.
import { expect, type CDPSession, type Locator, type Page } from '@playwright/test';
import { test } from './support/test';
import { axeViolations } from './support/axe';
import { BUILD_HAS_ENGINE } from './support/browser';
import { chooseRun, openParameters, openSearch, paste, program, submit, task, waitStatus } from './support/search';

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

  // --- WA2-AX-11: touch targets ----------------------------------------------------------------------

  test('a tap 12 px from a dot plot line selects it; a mouse click there does not', async ({ page }) => {
    // Two HSPs far apart: A is twice in the subject (as results.spec.ts's dot plot test).
    const a = dna(11, 60);
    await program(page, 'blastn');
    await paste(page, 'query', `>q1\n${dna(21, 20)}\n>q2\n${a}\n`);
    await paste(page, 'subject', `>s1\n${a}${dna(12, 60)}${a}\n`);
    await runAndOpen(page);
    await page.getByTestId('query-row-1').click();
    await showTab(page, 'pane-dotplot');
    const canvas = page.getByTestId('dotplot-canvas');
    await expect(canvas).toHaveAttribute('data-segments', '2');
    const selected = (await canvas.getAttribute('data-selected'))!;
    const targets = JSON.parse((await canvas.getAttribute('data-targets'))!) as { hsp: string; ends: number[] }[];
    const target = targets.find((t) => t.hsp !== selected)!;
    // A point 12 px from the line, across it, in the plot (the frame inside the axes' 50 px bands,
    // DotPlot's MARGIN), 30 px or more from every other line (the FakeEngine's two lines cross).
    const [left, top] = [50, 50];
    const [width, height] = (await canvas.getAttribute('data-plot'))!.split('x').map(Number) as [number, number];
    const [ax, ay, bx, by] = target.ends as [number, number, number, number];
    const length = Math.hypot(bx - ax, by - ay);
    const at = [0.5, 0.4, 0.6, 0.3, 0.7, 0.2, 0.8]
      .flatMap((t) => [1, -1].map((side) => ({ x: ax + t * (bx - ax) + (side * 12 * -(by - ay)) / length, y: ay + t * (by - ay) + (side * 12 * (bx - ax)) / length })))
      .find(
        (point) =>
          point.x > left + 2 &&
          point.x < left + width - 2 &&
          point.y > top + 2 &&
          point.y < top + height - 2 &&
          targets.every((other) => other === target || distance(point, other.ends) >= 30),
      )!;
    expect(at, 'a point 12 px from the line and 30 px from the others').toBeDefined();
    expect(distance(at, target.ends)).toBeCloseTo(12, 0);
    await centre(canvas);
    let box = (await canvas.boundingBox())!;

    // The mouse: as before, 8 px.
    await page.mouse.click(box.x + at.x, box.y + at.y);
    await page.waitForTimeout(200);
    await expect(canvas).toHaveAttribute('data-selected', selected);
    await expect(page.getByTestId('dotplot-popup')).toHaveCount(0);

    // A finger: the line is chosen, and its popup opens.
    await page.touchscreen.tap(box.x + at.x, box.y + at.y);
    await expect(canvas).toHaveAttribute('data-selected', target.hsp);
    await expect(page.getByTestId('dotplot-popup')).toHaveAttribute('data-hsp', target.hsp);

    // A tap far from every line closes the popup and keeps the selection.
    await centre(canvas);
    box = (await canvas.boundingBox())!;
    const far = farPoint(box, targets);
    expect(far, 'a point of the canvas 40 px from every line').toBeDefined();
    await page.touchscreen.tap(box.x + far!.x, box.y + far!.y);
    await expect(page.getByTestId('dotplot-popup')).toHaveCount(0);
    await expect(canvas).toHaveAttribute('data-selected', target.hsp);
  });

  test("a tap between two rows of the Graphic Summary chooses the nearer row's HSP; the mouse keeps its row and 4 px", async ({ page }) => {
    // The query's thirds in three subjects, shorter and shorter, so that the rows come in this
    // order: in the engine build each row's bar starts where the one above it ends (the
    // FakeEngine's bars all start at the query's start).
    const [p, q, r] = [dna(61, 80), dna(62, 80), dna(63, 80)];
    await program(page, 'blastn');
    await paste(page, 'query', `>q\n${p}${q}${r}\n`);
    await paste(page, 'subject', `>s1\n${dna(64, 30)}${p}\n>s2\n${q.slice(0, 70)}${dna(65, 30)}\n>s3\n${dna(66, 30)}${r.slice(0, 60)}${dna(67, 30)}\n`);
    await runAndOpen(page);
    const canvas = page.getByTestId('graphic-canvas');
    await expect(canvas).toHaveAttribute('data-rows', '3');
    const bars = JSON.parse((await canvas.getAttribute('data-targets'))!) as Bar[];
    const rowYs = [...new Set(bars.map((bar) => bar.y))].sort((m, n) => m - n);
    expect(rowYs).toHaveLength(3);
    const first = bars.filter((bar) => bar.y === rowYs[0]).sort((m, n) => n.x1 - m.x1)[0]!;
    // Points in the gap between the first and the second row, on both sides of the line between
    // their bands: 6 and 9 px past the right end of the first row's bar, and under its middle.
    const points = [
      { x: first.x1 + 6, y: first.y + 4 },
      { x: first.x1 + 6, y: first.y + 7 },
      { x: first.x1 + 9, y: first.y + 2 },
      { x: first.x, y: first.y + 7 },
    ];
    /**
     * Selects another HSP than `hsp` (with a click on the middle of the last row's bar, which no
     * other bar overlaps), then presses at `point` by `how`; returns the selection before that press.
     */
    const press = async (point: { x: number; y: number }, how: 'tap' | 'click', hsp: string): Promise<string> => {
      const last = bars.filter((bar) => bar.y === rowYs[2]);
      expect(last).toHaveLength(1);
      const other = last[0]!.hsp === hsp ? first : last[0]!;
      for (const [at, by] of [
        [other, 'click'],
        [point, how],
      ] as const) {
        await page.evaluate(() => window.scrollTo(0, 0));
        await centre(canvas);
        const box = (await canvas.boundingBox())!;
        if (by === 'tap') await page.touchscreen.tap(box.x + at.x, box.y + at.y);
        else await page.mouse.click(box.x + at.x, box.y + at.y);
        if (at === other) await expect(canvas).toHaveAttribute('data-selected', other.hsp);
      }
      return other.hsp;
    };
    let across = 0;
    for (const point of points) {
      // A finger: the bar nearest to it (12 px at most), in its row or the one above or below.
      const near = nearestBar(bars, point, 12);
      expect(near, `a bar within 12 px of ${JSON.stringify(point)}`).toBeDefined();
      if (near!.y !== bandOf(rowYs, point.y)) across++;
      await press(point, 'tap', near!.hsp);
      await expect(canvas, `a tap at ${JSON.stringify(point)}`).toHaveAttribute('data-selected', near!.hsp);
    }
    // In the engine build, a point in one row's band is nearer to the bar of the row below.
    if (BUILD_HAS_ENGINE) expect(across).toBeGreaterThan(0);
    // The mouse picks in the row under it, 4 px beside a bar at most, as before: a click at the
    // points that are further from the bars of their rows chooses nothing.
    let missed = 0;
    for (const point of points) {
      if (bars.some((bar) => bar.y === bandOf(rowYs, point.y) && Math.max(bar.x0 - point.x, point.x - bar.x1, 0) <= 4)) continue;
      missed++;
      const before = await press(point, 'click', nearestBar(bars, point, 12)!.hsp);
      await page.waitForTimeout(300);
      await expect(canvas, `a click at ${JSON.stringify(point)}`).toHaveAttribute('data-selected', before);
    }
    expect(missed).toBeGreaterThanOrEqual(2);
  });

  test('the "First lines" and "Records" of an input are 24 px touch targets', async ({ page }) => {
    await program(page, 'blastn');
    await paste(page, 'query', `>q\n${dna(71, 120)}\n`);
    await paste(page, 'subject', `>s\n${dna(72, 120)}\n`);
    for (const summary of await page.locator('[data-testid$="-source-0"] summary').all()) {
      expect((await summary.boundingBox())!.height, (await summary.textContent()) ?? '').toBeGreaterThanOrEqual(24);
    }
    expect(await page.locator('[data-testid$="-source-0"] summary').count()).toBeGreaterThanOrEqual(2);
    const found = await axeViolations(page, ['[data-testid="query-source-0"]', '[data-testid="subject-source-0"]']);
    expect(found.filter((violation) => violation.id === 'target-size')).toEqual([]);
  });
});

/** A point of the canvas 40 px or more from every line, or undefined. */
function farPoint(box: { width: number; height: number }, targets: readonly { ends: number[] }[]): { x: number; y: number } | undefined {
  for (let y = 10; y < box.height - 10; y += 5) {
    for (let x = 10; x < box.width - 10; x += 5) if (targets.every((t) => distance({ x, y }, t.ends) >= 40)) return { x, y };
  }
  return undefined;
}

/** A bar of the Graphic Summary in CSS pixels of its canvas: its HSP, the middle, the ends and the row's y. */
interface Bar {
  readonly hsp: string;
  readonly x: number;
  readonly y: number;
  readonly x0: number;
  readonly x1: number;
}

/** The y of the row whose band (12 px, its bar in the middle) holds `y`. */
function bandOf(rowYs: readonly number[], y: number): number | undefined {
  return rowYs.find((row) => y >= row - 6 && y < row + 6);
}

/** The bar (6 px high) nearest to a point, `reach` px away at most. */
function nearestBar(bars: readonly Bar[], point: { x: number; y: number }, reach: number): Bar | undefined {
  let best: Bar | undefined;
  let bestDistance = reach;
  for (const bar of bars) {
    const dx = Math.max(bar.x0 - point.x, point.x - bar.x1, 0);
    const dy = Math.max(Math.abs(point.y - bar.y) - 3, 0);
    const d = Math.hypot(dx, dy);
    if (d <= bestDistance) [best, bestDistance] = [bar, d];
  }
  return best;
}

/** The distance from a point to a line from (ax, ay) to (bx, by). */
function distance(point: { x: number; y: number }, [ax, ay, bx, by]: readonly number[]): number {
  const [dx, dy] = [bx! - ax!, by! - ay!];
  const t = dx === 0 && dy === 0 ? 0 : Math.max(0, Math.min(1, ((point.x - ax!) * dx + (point.y - ay!) * dy) / (dx * dx + dy * dy)));
  return Math.hypot(point.x - (ax! + t * dx), point.y - (ay! + t * dy));
}

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

// --- WA2-OP-01: the HSP table's Identity and Orientation ---------------------------------------------

for (const width of [1280, 375]) {
  test(`at ${width} px, the HSP tables show Identity and Orientation whole: 100.000 and "Not in the record"`, async ({ page }) => {
    await page.setViewportSize({ width, height: 900 });
    // A one-letter BLASTN HSP (results.spec.ts's): its record has no orientation, and its identity is 100.000.
    const qIdx = BUILD_HAS_ENGINE ? 0 : 2;
    await program(page, 'blastn');
    await task(page, 'blastn');
    await openParameters(page);
    await page.getByTestId('param-word_size').selectOption('other');
    await page.getByTestId('param-word_size-other').fill('4');
    await paste(page, 'query', BUILD_HAS_ENGINE ? '>q\nTAGGACGG\n' : '>q1\nACGTACGTACGT\n>q2\nACGTACGTACGT\n>q\nTAGGACGG\n');
    await paste(page, 'subject', '>s\nYCAYAANTNCRGYACT\n');
    await runAndOpen(page);
    if (!BUILD_HAS_ENGINE) await page.getByTestId(`query-row-${qIdx}`).click();
    for (const tab of ['pane-dotplot', 'results-view-hits'] as const) {
      await showTab(page, tab);
      const table = page.locator(tab === 'pane-dotplot' ? '[data-testid="results-dotplot"]' : '[data-testid="results-hits"]').getByTestId('hsp-table');
      const single = table.locator('[data-testid^="hsp-row-"][data-orientation="unknown"]');
      await expect(single).toHaveCount(1);
      await expect(single.locator('[data-field="orientation"]')).toHaveText('Not in the record');
      if (BUILD_HAS_ENGINE) await expect(single.locator('[data-field="pident"]')).toHaveText('100.000');
      // The cells of the two columns, and a copy of the row with BLAST's widest identity (the
      // FakeEngine writes FAKE): no text is cut.
      const cut = await single.evaluate((row) => {
        const copy = row.cloneNode(true) as HTMLElement;
        copy.removeAttribute('data-testid');
        copy.querySelector('[data-field="pident"]')!.textContent = '100.000';
        row.after(copy);
        const cells = [row, copy].flatMap((r) => [...r.querySelectorAll<HTMLElement>('[data-field="pident"], [data-field="orientation"]')]);
        const out = cells.filter((cell) => cell.scrollWidth > cell.clientWidth).map((cell) => `${cell.dataset['field']} "${cell.textContent}" ${cell.scrollWidth} > ${cell.clientWidth}`);
        copy.remove();
        return out;
      });
      expect(cut, tab).toEqual([]);
    }
  });
}

// --- WA2-OP-02: numbers on the axes of a short span -------------------------------------------------

/** The tick labels that the dot plot drew on its query and subject axes. */
async function tickLabels(canvas: Locator): Promise<[string[], string[]]> {
  return JSON.parse((await canvas.getAttribute('data-ticks')) ?? '[[],[]]') as [string[], string[]];
}

test('the dot plot numbers both axes of a view shorter than 100 letters: a short pair, a zoomed view, an HSP of one letter', async ({ page }) => {
  // Two short sequences (8 against 16 letters, results.spec.ts's HSP of one letter), whole.
  const qIdx = BUILD_HAS_ENGINE ? 0 : 2;
  await program(page, 'blastn');
  await task(page, 'blastn');
  await openParameters(page);
  await page.getByTestId('param-word_size').selectOption('other');
  await page.getByTestId('param-word_size-other').fill('4');
  await paste(page, 'query', BUILD_HAS_ENGINE ? '>q\nTAGGACGG\n' : '>q1\nACGTACGTACGT\n>q2\nACGTACGTACGT\n>q\nTAGGACGG\n');
  await paste(page, 'subject', '>s\nYCAYAANTNCRGYACT\n');
  await runAndOpen(page);
  if (!BUILD_HAS_ENGINE) await page.getByTestId(`query-row-${qIdx}`).click();
  await showTab(page, 'pane-dotplot');
  const canvas = page.getByTestId('dotplot-canvas');
  await expect(canvas).toHaveAttribute('data-view', '0,8,0,16');
  // Steps of 2 letters on the query's 8, of 5 on the subject's 16 (blast2dotplot.py's 100 gave each axis its 0 only).
  await expect.poll(() => tickLabels(canvas)).toEqual([
    ['0', '2', '4', '6', '8'],
    ['0', '5', '10', '15'],
  ]);
  // "Zoom to HSP" on the HSP of one letter.
  const single = page.getByTestId('results-dotplot').getByTestId('hsp-list').locator('[data-testid^="hsp-row-"][data-orientation="unknown"]');
  await single.click();
  await page.getByTestId('dotplot-zoom-hsp').click();
  await page.waitForTimeout(200);
  for (const labels of await tickLabels(canvas)) expect(labels.length).toBeGreaterThanOrEqual(2);
});

test('the dot plot numbers both axes of a view zoomed to under 100 letters, and keeps blast2dotplot.py\'s ticks above 100', async ({ page }) => {
  await program(page, 'blastn');
  await paste(page, 'query', `>q\n${LONG}\n`);
  await paste(page, 'subject', `>s\n${LONG}\n`);
  await runAndOpen(page);
  await showTab(page, 'pane-dotplot');
  const canvas = page.getByTestId('dotplot-canvas');
  await expect(canvas).toHaveAttribute('data-view', '0,600,0,600');
  // The whole 600 letters: blast2dotplot.py's 100 a major tick.
  await expect.poll(async () => (await tickLabels(canvas))[0]).toEqual(['0', '100', '200', '300', '400', '500', '600']);
  // Zoomed in to 600 / 1.5^n letters: under 100 from n = 5 (79 letters) on, under 20 at n = 9 (16)
  // and 10 (10.4). Above 100 letters, the steps stay blast2dotplot.py's (100 a major tick up to 1000).
  for (let n = 1; n <= 10; n++) {
    const before = (await canvas.getAttribute('data-view'))!;
    await page.getByTestId('dotplot-zoom-in').click();
    await expect(canvas).not.toHaveAttribute('data-view', before);
    const [x0, x1, y0, y1] = await view(canvas);
    const labels = await tickLabels(canvas);
    const span = x1! - x0!;
    const says = (axis: number) => `${n}: ${axis === 0 ? `query ${x0}–${x1}` : `subject ${y0}–${y1}`}: ${labels[axis]!.join(' ')}`;
    for (const axis of [0, 1]) {
      if (span <= 100) expect(labels[axis]!.length, says(axis)).toBeGreaterThanOrEqual(2);
      else expect(labels[axis]!.map((label) => Number(label.replace(/,/g, '')) % 100), says(axis)).toEqual(labels[axis]!.map(() => 0));
    }
  }
});
