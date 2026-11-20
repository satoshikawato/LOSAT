// The text colours of styles.css :root on the grounds that the page draws them on are 4.5:1 or more
// (WCAG 2.1 SC 1.4.3, axe's color-contrast). Line WQA's random walk found the error colour #b84f00
// at 4.3:1 on the page ground, under a refused Run's message and the Additional options' check, in
// the three engines (2026-11-19): the token had been chosen on the white group alone.
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';

const css = readFileSync(resolve(import.meta.dirname, '../../src/ui/styles.css'), 'utf8');
const root = css.slice(css.indexOf(':root {'), css.indexOf('\n}', css.indexOf(':root {')));
const TOKENS = new Map([...root.matchAll(/--([a-z0-9-]+):\s*(#[0-9a-f]{6});/gi)].map((m) => [m[1]!, m[2]!.toLowerCase()]));

/** The page's own grounds, under any text (the comment on :root). */
const GROUNDS = ['page', 'group', 'band', 'bg'] as const;

/**
 * Each text colour and the grounds it is drawn on: the page's grounds, and the tinted ones of its
 * rules. The error colour also marks a reverse orientation in a selected row (--accent-bg) and a
 * missing file in an error box (--error-bg).
 */
const TEXT_ON: ReadonlyArray<readonly [string, readonly string[]]> = [
  ['fg', [...GROUNDS, 'legend', 'panel']],
  ['fg-2', [...GROUNDS, 'panel']],
  ['muted', [...GROUNDS, 'panel']],
  ['accent', ['group', 'band', 'bg', 'accent-bg', 'panel']],
  ['ok', [...GROUNDS, 'ok-bg']],
  ['error', [...GROUNDS, 'error-bg', 'accent-bg', 'panel']],
  ['warn-fg', ['warn-bg']],
  ['tab-fg', ['tab', 'tab-hover']],
  ['header-fg', ['header', 'accent']],
  ['header-muted', ['header']],
];

function luminance(hex: string): number {
  const [r, g, b] = [1, 3, 5].map((at) => {
    const c = parseInt(hex.slice(at, at + 2), 16) / 255;
    return c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4;
  });
  return 0.2126 * r! + 0.7152 * g! + 0.0722 * b!;
}

function contrast(a: string, b: string): number {
  const [light, dark] = [luminance(a), luminance(b)].sort((x, y) => y - x);
  return (light! + 0.05) / (dark! + 0.05);
}

describe('colour contrast of the text colours', () => {
  it('measures as WCAG does', () => {
    expect(contrast('#000000', '#ffffff')).toBeCloseTo(21, 5);
    expect(contrast('#b84f00', '#e9eef3')).toBeLessThan(4.5);
  });

  for (const [text, grounds] of TEXT_ON) {
    it(`--${text} is 4.5:1 or more on ${grounds.map((ground) => `--${ground}`).join(', ')}`, () => {
      const colour = TOKENS.get(text);
      expect(colour, `--${text} in :root`).toBeDefined();
      for (const ground of grounds) {
        const behind = TOKENS.get(ground);
        expect(behind, `--${ground} in :root`).toBeDefined();
        expect(contrast(colour!, behind!), `--${text} ${colour} on --${ground} ${behind}`).toBeGreaterThanOrEqual(4.5);
      }
    });
  }

  it('the orientation of a reverse HSP is the error colour', () => {
    expect(css).toMatch(/\.orientation\[data-orientation='reverse'\] \{\s*color: var\(--error\);/);
  });
});
