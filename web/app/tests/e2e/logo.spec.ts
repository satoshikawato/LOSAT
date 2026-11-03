// DW-43: the logo in the header, the tab's icon, and the loader (the logo's A drawn in one stroke)
// before the words of a state in progress. An accession fetch whose answer is held keeps one such
// state on the page.
import { expect, test } from '@playwright/test';
import { openSearch } from './support/search';

test('the header shows the logo, and the tab icon is served', async ({ page }) => {
  await openSearch(page);
  await expect(page.getByTestId('losat-logo')).toBeVisible();
  await expect(page.getByRole('heading', { level: 1 })).toHaveAccessibleName('LOSAT');
  const href = await page.locator('link[rel="icon"]').getAttribute('href');
  const icon = await page.request.get(new URL(href!, page.url()).href);
  expect(icon.ok()).toBe(true);
  expect(icon.headers()['content-type']).toContain('image/svg+xml');
});

test('a state in progress shows the loader; it draws, and stands still when the system asks for less motion', async ({ page }) => {
  let release!: () => void;
  const held = new Promise<void>((resolve) => (release = resolve));
  await page.route(/^https:\/\/eutils\.ncbi\.nlm\.nih\.gov\//, async (route) => {
    await held;
    await route.fulfill({
      status: 200,
      contentType: 'text/plain',
      headers: { 'Access-Control-Allow-Origin': '*' },
      body: '>NR_074540.1 16S ribosomal RNA\nTTATTGGAGAGTTTGATCCTGGCTCAGGATGAACGCTGGCGGCGTGCCTAATACATGCAAGTCGAGCGAA\n',
    });
  });
  await openSearch(page);
  await page.getByTestId('query-input').fill('NR_074540.1');
  await page.getByTestId('query-fetch').click();

  const status = page.getByTestId('query-fetch-status');
  await expect(status).toHaveText('Fetching 1 of 1 from NCBI…');
  const loader = status.getByTestId('losat-loader');
  await expect(loader).toBeVisible();
  // The pens' mask is the loader's own (each loader has one).
  const mask = await loader.locator('mask').getAttribute('id');
  await expect(loader.locator('g[mask]')).toHaveAttribute('mask', `url(#${mask})`);
  expect(await loader.evaluate((el) => el.getAnimations().map((a) => (a as CSSAnimation).animationName))).toEqual(['losat-loader-stroke']);
  // It moves: no part of the loop holds the first pen's dashes for 0.3 s.
  const pen = loader.locator('.losat-loader-pen-a');
  const dashes = () => pen.evaluate((el) => getComputedStyle(el).strokeDasharray);
  const first = await dashes();
  await expect.poll(dashes, { intervals: [300] }).not.toBe(first);

  await page.emulateMedia({ reducedMotion: 'reduce' });
  expect(await loader.evaluate((el) => el.getAnimations().length)).toBe(0);
  expect(await loader.locator('.losat-loader-pen').evaluateAll((els) => els.map((el) => getComputedStyle(el).strokeDasharray))).toEqual(['none', 'none']);

  release();
  await expect(status).toHaveCount(0);
});
