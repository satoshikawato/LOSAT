// DW-43: the logo in the header, the tab's icon (PNG icons DW-44, the web app manifest DW-45), and the loader (the logo's A drawn in one stroke)
// before the words of a state in progress. An accession fetch whose answer is held keeps one such
// state on the page.
import { expect } from '@playwright/test';
import { test } from './support/test';
import { openSearch } from './support/search';

test('the header shows the logo, and the tab and home-screen icons are served', async ({ page }) => {
  await openSearch(page);
  await expect(page.getByTestId('losat-logo')).toBeVisible();
  await expect(page.getByRole('heading', { level: 1 })).toHaveAccessibleName('LOSAT');
  const links = await page
    .locator('link[rel="icon"], link[rel="apple-touch-icon"]')
    .evaluateAll((els) => els.map((el) => [el.getAttribute('rel'), el.getAttribute('href')] as const));
  expect(links).toEqual([
    ['icon', '/favicon-32.png'],
    ['icon', '/favicon.svg'],
    ['apple-touch-icon', '/apple-touch-icon.png'],
  ]);
  for (const [, href] of links) {
    const icon = await page.request.get(new URL(href!, page.url()).href);
    expect(icon.ok(), href!).toBe(true);
    expect(icon.headers()['content-type'], href!).toContain(href!.endsWith('.svg') ? 'image/svg+xml' : 'image/png');
  }
});

test('the web app manifest and its icons are served', async ({ page, browserName }) => {
  await openSearch(page);
  await expect(page.locator('link[rel="manifest"]')).toHaveAttribute('href', '/manifest.webmanifest');
  const response = await page.request.get(new URL('/manifest.webmanifest', page.url()).href);
  expect(response.ok()).toBe(true);
  expect(response.headers()['content-type']).toContain('application/manifest+json');
  const manifest = (await response.json()) as { icons: { src: string }[] };
  for (const { src } of manifest.icons) {
    const icon = await page.request.get(new URL(src, page.url()).href);
    expect(icon.ok(), src).toBe(true);
    expect(icon.headers()['content-type'], src).toContain('image/png');
  }
  // Chromium's own reading (only Chromium answers): it finds the manifest from the page and reads
  // every member without an error (an invalid colour, for one, is an error) and every icon.
  if (browserName !== 'chromium') return;
  const read = await (await page.context().newCDPSession(page)).send('Page.getAppManifest');
  expect(read.url).toBe(response.url());
  expect(read.errors).toEqual([]);
  expect(read.manifest.icons?.map((icon) => icon.url)).toEqual(manifest.icons.map(({ src }) => new URL(src, page.url()).href));
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
  await page.getByTestId('query-fetch-source').selectOption('ncbi');
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

test('the race for the fastest site shows the loader too', async ({ page }) => {
  let release!: () => void;
  const held = new Promise<void>((resolve) => (release = resolve));
  await page.route(/^https:\/\/(?:eutils\.ncbi\.nlm\.nih\.gov|www\.ebi\.ac\.uk|rest\.uniprot\.org|ngdc\.cncb\.ac\.cn)\//, async (route) => {
    await held;
    await route.fulfill({ status: 404, headers: { 'Access-Control-Allow-Origin': '*' }, body: 'not found' });
  });
  await openSearch(page);
  await page.getByTestId('query-input').fill('NR_074540.1');
  await expect(page.getByTestId('query-fetch-source')).toHaveValue('fastest');
  await page.getByTestId('query-fetch').click();
  const status = page.getByTestId('query-fetch-status');
  await expect(status).toHaveText('Finding the fastest site…');
  await expect(status.getByTestId('losat-loader')).toBeVisible();
  release();
});
