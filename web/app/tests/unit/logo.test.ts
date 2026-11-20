import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { LOADER_LENGTH, LOADER_PENS, LOGO_A, LOGO_LETTERS, LOGO_VIEWBOX } from '../../src/ui/logo';

// styles.css repeats the loader pens' numbers (the markup carries no style attributes, which the CSP refuses).
const css = readFileSync(resolve(import.meta.dirname, '../../src/ui/styles.css'), 'utf8');

/** The pen's length along its path, from its `d` (straight segments only). */
function measured(d: string): number {
  const points = [...d.matchAll(/(-?[\d.]+),(-?[\d.]+)/g)].map((m) => [Number(m[1]), Number(m[2])] as const);
  let length = 0;
  for (let i = 1; i < points.length; i++) length += Math.hypot(points[i]![0] - points[i - 1]![0], points[i]![1] - points[i - 1]![1]);
  return length;
}

describe('the loader (DW-43)', () => {
  it('each pen in styles.css starts and runs as logo.ts says, and its length is its path', () => {
    for (const pen of LOADER_PENS) {
      const rule = new RegExp(`\\.losat-loader-pen-${pen.id} \\{\\s*--pen-start: ([\\d.]+);\\s*--pen-length: ([\\d.]+);`).exec(css);
      expect(rule, pen.id).not.toBeNull();
      expect(Number(rule![1])).toBe(pen.start);
      expect(Number(rule![2])).toBe(pen.length);
      expect(measured(pen.d)).toBeCloseTo(pen.length, 2);
    }
  });

  it('runs to the end of the last pen, and the second pen starts where the first lands on the base', () => {
    expect(LOADER_LENGTH).toBeCloseTo(Math.max(...LOADER_PENS.map((pen) => pen.start + pen.length)), 3);
    const [a, b] = LOADER_PENS;
    expect(b!.d.startsWith(`M ${a!.d.split(' L ')[2]}`)).toBe(true);
    expect(measured(a!.d.split(' L ').slice(0, 3).join(' L '))).toBeCloseTo(b!.start, 2);
    const ends = [...css.matchAll(/--losat-loader-(?:from|to): ([\d.]+);|initial-value: ([\d.]+);/g)]
      .map((m) => Number(m[1] ?? m[2]))
      .filter((value) => value !== 0);
    expect(ends.length).toBeGreaterThan(0);
    expect(new Set(ends)).toEqual(new Set([LOADER_LENGTH]));
  });
});

/** The `d` of every path in an SVG file, in order. */
function paths(file: string): string[] {
  return [...readFileSync(file, 'utf8').matchAll(/<path d="([^"]+)"/g)].map((m) => m[1]!);
}

describe('the logo outside the app (DW-44)', () => {
  it("the README's logos and the tab icon carry logo.ts's paths", () => {
    for (const name of ['losat-logo.svg', 'losat-logo-white.svg']) {
      const file = resolve(import.meta.dirname, '../../../../docs/assets', name);
      expect(paths(file), name).toEqual([...LOGO_LETTERS, ...LOGO_A]);
      expect(readFileSync(file, 'utf8'), name).toContain(`viewBox="${LOGO_VIEWBOX}"`);
    }
    expect(paths(resolve(import.meta.dirname, '../../public/favicon.svg'))).toEqual([...LOGO_A]);
  });

  it('the PNG icons have their sizes (scripts/icons.mjs draws them from favicon.svg)', () => {
    for (const [name, size] of [['favicon-32.png', 32], ['apple-touch-icon.png', 180]] as const) {
      expect(pngSize(name), name).toEqual([size, size]);
    }
  });
});

/** The width and the height of a PNG in public/: IHDR's follow the signature (8 bytes) and the chunk's length and type. */
function pngSize(name: string): [number, number] {
  const png = readFileSync(resolve(import.meta.dirname, '../../public', name));
  return [png.readUInt32BE(16), png.readUInt32BE(20)];
}

describe('the web app manifest (DW-50)', () => {
  const manifest = JSON.parse(readFileSync(resolve(import.meta.dirname, '../../public/manifest.webmanifest'), 'utf8')) as {
    theme_color: string;
    background_color: string;
    icons: { src: string; sizes: string; type: string; purpose?: string }[];
  };

  it('opens LOSAT at the root in a window of its own', () => {
    expect(manifest).toMatchObject({ id: '/', name: 'LOSAT', start_url: '/', display: 'standalone' });
  });

  it("takes the header band's navy and the page's ground from styles.css", () => {
    expect(manifest.theme_color).toBe(/--header: (#[0-9a-f]{6});/.exec(css)![1]);
    expect(manifest.background_color).toBe(/--page: (#[0-9a-f]{6});/.exec(css)![1]);
  });

  it('lists icons of 192 and 512 px and a maskable one, each PNG of the size it says', () => {
    for (const icon of manifest.icons) {
      const [width, height] = icon.sizes.split('x').map(Number);
      expect(icon.type, icon.src).toBe('image/png');
      expect(pngSize(icon.src.replace(/^\//, '')), icon.src).toEqual([width, height]);
    }
    const any = manifest.icons.filter((icon) => (icon.purpose ?? 'any') === 'any').map((icon) => icon.sizes);
    expect(any).toEqual(['192x192', '512x512']);
    expect(manifest.icons.filter((icon) => icon.purpose === 'maskable').map((icon) => icon.sizes)).toEqual(['512x512']);
  });
});
