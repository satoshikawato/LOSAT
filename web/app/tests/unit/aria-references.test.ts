// The application gives ID references (aria-labelledby, aria-controls, ...) as attributes only, never
// through ARIA element reflection (`element.ariaLabelledByElements` and the like): the axe helper of
// the end-to-end tests reads the attributes alone (tests/e2e/support/axe.ts, because reading the
// reflected property crashed Firefox 159 build 1572, 2026-11-19), so a reference set the other way
// would go unchecked.
import { readdirSync, readFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { describe, expect, it } from 'vitest';

const SRC = resolve(import.meta.dirname, '../../src');

function files(dir: string): string[] {
  return readdirSync(dir, { withFileTypes: true }).flatMap((entry) =>
    entry.isDirectory() ? files(join(dir, entry.name)) : /\.(ts|vue)$/.test(entry.name) ? [join(dir, entry.name)] : [],
  );
}

describe('ID references of the application', () => {
  it('are attributes: no reflected ARIA element property in src', () => {
    const used = files(SRC).flatMap((file) =>
      [...readFileSync(file, 'utf8').matchAll(/\baria[A-Z][A-Za-z]*Elements?\b/g)].map((m) => `${file.slice(SRC.length + 1)}: ${m[0]}`),
    );
    expect(files(SRC).length).toBeGreaterThan(20);
    expect(used).toEqual([]);
  });
});
