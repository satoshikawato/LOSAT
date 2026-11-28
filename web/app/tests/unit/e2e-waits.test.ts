// The end-to-end tests wait for the answer to a tool tab's check (the `data-state` of `wm-`, `ms-`, `mk-` or
// `bf-argv-validation`) through `checked` alone (tests/e2e/support/search.ts), whose time limit allows for the start
// of the Engine worker that the page's first check makes (ENGINE_CHECK_TIMEOUT). A spec that read the state with its
// own expect waited Playwright's 5 s, and failed under load in WebKit and Chromium while its check waited for the
// engine's two files to be fetched and compiled (line WW, WW4j).
import { readdirSync, readFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { describe, expect, it } from 'vitest';

const E2E = resolve(import.meta.dirname, '../e2e');
const HELPER = join(E2E, 'support/search.ts');

/** A tool tab's check (a literal prefix or one given by a template) whose `data-state` the same statement reads. */
const OWN_WAIT = /(?:\b(?:wm|ms|mk|bf)|\$\{[^}]*\})-argv-validation[^;]*?\bdata-state\b/gs;

function files(dir: string): string[] {
  return readdirSync(dir, { withFileTypes: true }).flatMap((entry) =>
    entry.isDirectory() ? files(join(dir, entry.name)) : entry.name.endsWith('.ts') ? [join(dir, entry.name)] : [],
  );
}

/** The statements of `text` that read the state of a tool tab's check themselves, with their line numbers. */
function ownWaits(text: string): string[] {
  return [...text.matchAll(OWN_WAIT)].map((m) => `${text.slice(0, m.index).split('\n').length}: ${m[0].replace(/\s+/g, ' ')}`);
}

describe('the waits of the end-to-end tests for a tool tab’s check', () => {
  it('finds the forms that a spec could write', () => {
    const found = [
      "await expect(page.getByTestId('wm-argv-validation')).toHaveAttribute('data-state', 'ok');",
      "await expect(page.getByTestId('mk-argv-validation')).not.toHaveAttribute('data-state', 'checking');",
      "await expect(page.getByTestId('ms-argv-validation'))\n    .toHaveAttribute('data-state', state, { timeout: 30_000 });",
      "await expect(page.getByTestId(`${prefix}-argv-validation`)).toHaveAttribute('data-state', 'ok');",
      "await expect(page.locator('[data-testid=\"bf-argv-validation\"]')).toHaveAttribute('data-state', 'invalid');",
      "const state = await page.getByTestId('wm-argv-validation').getAttribute('data-state');",
    ];
    for (const text of found) expect(ownWaits(text), text).toHaveLength(1);
    const kept = [
      // The search form's check is the Data worker's, which the describe already started.
      "await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');",
      // The text of an answer, read after `checked`.
      "await expect(page.getByTestId('wm-argv-validation')).toHaveText(refused);\nawait expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');",
      "await checked(page, 'wm-panel', 'ok');",
    ];
    for (const text of kept) expect(ownWaits(text), text).toEqual([]);
  });

  it('go through `checked` in every spec', () => {
    const all = files(E2E);
    expect(all.length).toBeGreaterThan(40);
    const own = all
      .filter((file) => file !== HELPER)
      .flatMap((file) => ownWaits(readFileSync(file, 'utf8')).map((wait) => `${file.slice(E2E.length + 1)}:${wait}`));
    expect(own).toEqual([]);
  });
});
