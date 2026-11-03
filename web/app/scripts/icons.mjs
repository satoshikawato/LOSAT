// Draws the PNG icons from public/favicon.svg (DW-44) with Playwright's Chromium:
//   public/favicon-32.png        the tab icon for browsers that take no SVG icon
//   public/apple-touch-icon.png  180 px for an iOS home screen, which takes no SVG icon. iOS rounds
//                                the corners itself and shows transparent pixels as black, so the
//                                tile fills the square and the A keeps a wider margin.
// Run after a change of favicon.svg: node scripts/icons.mjs
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { chromium } from '@playwright/test';

const PUBLIC = resolve(import.meta.dirname, '../public');
const favicon = readFileSync(resolve(PUBLIC, 'favicon.svg'), 'utf8');
const [, viewBox] = /viewBox="([^"]+)"/.exec(favicon);
const [, a] = /(<g fill="#fff">.*<\/g>)/.exec(favicon);
const [x, y, side] = viewBox.split(' ').map(Number);

// The A is 10.27 high in the favicon's 13.53; here 62 % of the side.
const homeSide = 10.267 / 0.62;
const homeX = x + side / 2 - homeSide / 2;
const homeY = y + side / 2 - homeSide / 2;
const home =
  `<svg xmlns="http://www.w3.org/2000/svg" viewBox="${homeX.toFixed(3)} ${homeY.toFixed(3)} ${homeSide.toFixed(3)} ${homeSide.toFixed(3)}">` +
  `<rect x="${homeX.toFixed(3)}" y="${homeY.toFixed(3)}" width="${homeSide.toFixed(3)}" height="${homeSide.toFixed(3)}" fill="#0b3a5e"/>${a}</svg>`;

const browser = await chromium.launch();
const page = await browser.newPage();
for (const [file, svg, size] of [
  ['favicon-32.png', favicon, 32],
  ['apple-touch-icon.png', home, 180],
]) {
  await page.setViewportSize({ width: size, height: size });
  const src = `data:image/svg+xml;base64,${Buffer.from(svg).toString('base64')}`;
  await page.setContent(`<body style="margin:0"><img src="${src}" width="${size}" height="${size}" style="display:block">`);
  await page.locator('img').evaluate((img) => img.decode());
  await page.screenshot({ path: resolve(PUBLIC, file), omitBackground: true, clip: { x: 0, y: 0, width: size, height: size } });
}
await browser.close();
