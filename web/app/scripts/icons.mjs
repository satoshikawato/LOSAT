// Draws the PNG icons from public/favicon.svg (DW-44) with Playwright's Chromium:
//   public/favicon-32.png        the tab icon for browsers that take no SVG icon
//   public/apple-touch-icon.png  180 px for an iOS home screen, which takes no SVG icon. iOS rounds
//                                the corners itself and shows transparent pixels as black, so the
//                                tile fills the square and the A keeps a wider margin.
//   public/icon-192.png, icon-512.png, icon-maskable-512.png
//                                the icons of manifest.webmanifest (DW-45): the tab icon's tile, and for
//                                Android's home screen, which cuts its own shape out of the middle, a
//                                full square whose A stays inside the middle circle of 80 % (the
//                                maskable safe zone).
// Run after a change of favicon.svg: node scripts/icons.mjs
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { chromium } from '@playwright/test';

const PUBLIC = resolve(import.meta.dirname, '../public');
const favicon = readFileSync(resolve(PUBLIC, 'favicon.svg'), 'utf8');
const [, viewBox] = /viewBox="([^"]+)"/.exec(favicon);
const [, a] = /(<g fill="#fff">.*<\/g>)/.exec(favicon);
const [x, y, side] = viewBox.split(' ').map(Number);

/** A full square of the tile's navy around the A, which is 10.27 high (in the favicon's 13.53) and here `share` of the side. */
function square(share) {
  const full = 10.267 / share;
  const [sx, sy, s] = [x + side / 2 - full / 2, y + side / 2 - full / 2, full].map((v) => v.toFixed(3));
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="${sx} ${sy} ${s} ${s}"><rect x="${sx}" y="${sy}" width="${s}" height="${s}" fill="#0b3a5e"/>${a}</svg>`;
}

const browser = await chromium.launch();
const page = await browser.newPage();
for (const [file, svg, size] of [
  ['favicon-32.png', favicon, 32],
  ['apple-touch-icon.png', square(0.62), 180],
  ['icon-192.png', favicon, 192],
  ['icon-512.png', favicon, 512],
  // The A's lower corners, the farthest from the middle (7.12 in the favicon), at 0.37 of the side.
  ['icon-maskable-512.png', square(0.54), 512],
]) {
  await page.setViewportSize({ width: size, height: size });
  const src = `data:image/svg+xml;base64,${Buffer.from(svg).toString('base64')}`;
  await page.setContent(`<body style="margin:0"><img src="${src}" width="${size}" height="${size}" style="display:block">`);
  await page.locator('img').evaluate((img) => img.decode());
  await page.screenshot({ path: resolve(PUBLIC, file), omitBackground: true, clip: { x: 0, y: 0, width: size, height: size } });
}
await browser.close();
