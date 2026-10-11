import { describe, expect, it } from 'vitest';
import { RELEASE_ID_PATTERN, releaseId } from '../../build/release';

describe('release ID', () => {
  it('takes LOSAT_WEB_RELEASE_ID when it is set', () => {
    expect(releaseId(undefined, { LOSAT_WEB_RELEASE_ID: 'test-a.1' })).toBe('test-a.1');
  });

  it('refuses a name that cannot be a path segment', () => {
    for (const name of ['A', '../x', 'a/b', '-a', 'a b', 'x'.repeat(65)]) {
      expect(() => releaseId(undefined, { LOSAT_WEB_RELEASE_ID: name })).toThrow(/LOSAT_WEB_RELEASE_ID/);
    }
  });

  it('is made from the commit and the engine otherwise', () => {
    const id = releaseId(undefined, {});
    expect(id).toMatch(RELEASE_ID_PATTERN);
    expect(id).toMatch(/-fake(-dirty)?$/);
  });

  it('is defined for the application', () => {
    expect(__LOSAT_RELEASE_ID__).toMatch(RELEASE_ID_PATTERN);
  });
});
