import { execFileSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { afterEach, describe, expect, it } from 'vitest';
import { BUILD_INPUTS, CERTIFICATION_RECORDS, SHARED_MEMORY_GUARD, VERIFICATION_CELLS } from '../../build/inputs';
import { appBuild, buildEnvironmentProblem, RELEASE_ID_PATTERN, releaseId, sourceState } from '../../build/release';

const COMMIT = '0123456789abcdef0123456789abcdef01234567';

describe('the build environment', () => {
  it('refuses a build with NODE_ENV set to anything but production: that would be a development build, with source paths in the bundle', () => {
    expect(buildEnvironmentProblem({})).toBeUndefined();
    expect(buildEnvironmentProblem({ NODE_ENV: '' })).toBeUndefined();
    expect(buildEnvironmentProblem({ NODE_ENV: 'production' })).toBeUndefined();
    for (const value of ['test', 'development', 'Production']) {
      expect(buildEnvironmentProblem({ NODE_ENV: value })).toMatch(new RegExp(`^NODE_ENV=${value}: this would be a development build.*Unset NODE_ENV`));
    }
  });
});

describe('release ID', () => {
  it('takes LOSAT_WEB_RELEASE_ID when it is set', () => {
    expect(releaseId(undefined, { env: { LOSAT_WEB_RELEASE_ID: 'test-a.1' } })).toBe('test-a.1');
    expect(releaseId(undefined, { env: { LOSAT_WEB_RELEASE_ID: 'test-a.1' }, source: undefined })).toBe('test-a.1');
  });

  it('refuses a name that cannot be a path segment', () => {
    for (const name of ['A', '../x', 'a/b', '-a', 'a b', 'x'.repeat(65)]) {
      expect(() => releaseId(undefined, { env: { LOSAT_WEB_RELEASE_ID: name } })).toThrow(/LOSAT_WEB_RELEASE_ID/);
    }
  });

  it('is the commit cut to 12 hex digits, the engine, and -dirty when an input differs from the commit', () => {
    expect(releaseId(undefined, { env: {}, source: { commit: COMMIT, dirty: false } })).toBe('0123456789ab-fake');
    expect(releaseId(undefined, { env: {}, source: { commit: COMMIT, dirty: true } })).toBe('0123456789ab-fake-dirty');
    // The app's build is the commit cut to 8 hex digits: the same in every clone, whatever length git abbreviates to there.
    expect(appBuild({ commit: COMMIT, dirty: false })).toBe('01234567');
    expect(appBuild(undefined)).toBe('unknown');
  });

  it('refuses a build without a commit unless it is named; the development server and Vitest run without one', () => {
    expect(() => releaseId(undefined, { env: {}, source: undefined })).toThrow(/no git commit.*LOSAT_WEB_RELEASE_ID/);
    expect(() => releaseId(undefined, { env: {}, source: undefined, build: true })).toThrow(/LOSAT_WEB_RELEASE_ID/);
    expect(releaseId(undefined, { env: {}, source: undefined, build: false })).toBe('unknown-fake');
  });

  it('is made from this repository’s full commit, not from git’s abbreviation', () => {
    const source = sourceState();
    expect(source?.commit).toMatch(/^[0-9a-f]{40}$/);
    const id = releaseId(undefined, { env: {} });
    expect(id).toMatch(RELEASE_ID_PATTERN);
    expect(id).toBe(`${source!.commit.slice(0, 12)}-fake${source!.dirty ? '-dirty' : ''}`);
  });

  it('is defined for the application', () => {
    expect(__LOSAT_RELEASE_ID__).toMatch(RELEASE_ID_PATTERN);
    expect(__LOSAT_APP_BUILD__).toMatch(/^(?:[0-9a-f]{8}|unknown)$/);
  });
});

describe('the inputs that the dirty check covers', () => {
  let dirs: string[] = [];
  afterEach(() => {
    for (const dir of dirs) rmSync(dir, { recursive: true, force: true });
    dirs = [];
  });

  /** A scratch repository with one commit of `files`. */
  function repository(files: Record<string, string>): string {
    const dir = mkdtempSync(join(tmpdir(), 'losat-release-id-'));
    dirs.push(dir);
    for (const [path, text] of Object.entries(files)) {
      mkdirSync(join(dir, path, '..'), { recursive: true });
      writeFileSync(join(dir, path), text);
    }
    const git = (...args: string[]) =>
      execFileSync('git', ['-C', dir, '-c', 'user.name=t', '-c', 'user.email=t@example.invalid', '-c', 'commit.gpgsign=false', ...args], {
        stdio: ['ignore', 'pipe', 'pipe'],
      });
    git('init', '-q');
    git('add', '.');
    git('commit', '-q', '-m', 'inputs');
    return dir;
  }

  it('lists web/app and every file of the repository that the build plugins read', () => {
    expect(BUILD_INPUTS).toEqual(['web/app', VERIFICATION_CELLS, ...Object.values(CERTIFICATION_RECORDS), SHARED_MEMORY_GUARD]);
    expect(BUILD_INPUTS).toContain('docs/web/verification_cells.tsv');
    expect(BUILD_INPUTS).toContain('LOSAT/tests/wasi_shared_memory.js');
  });

  it('marks a change to a listed input, an untracked file among them, and nothing else', () => {
    const inputs = ['web/app', 'docs/web/cells.tsv'];
    const dir = repository({ 'web/app/a.ts': 'a\n', 'docs/web/cells.tsv': 'x\n', 'docs/other.md': 'o\n' });
    const clean = sourceState(dir, inputs);
    expect(clean).toEqual({ commit: expect.stringMatching(/^[0-9a-f]{40}$/), dirty: false });
    writeFileSync(join(dir, 'docs/other.md'), 'changed\n');
    expect(sourceState(dir, inputs)?.dirty).toBe(false);
    writeFileSync(join(dir, 'docs/web/cells.tsv'), 'y\n');
    expect(sourceState(dir, inputs)).toEqual({ commit: clean!.commit, dirty: true });
    writeFileSync(join(dir, 'docs/web/cells.tsv'), 'x\n');
    expect(sourceState(dir, inputs)?.dirty).toBe(false);
    writeFileSync(join(dir, 'web/app/new.ts'), 'n\n');
    expect(sourceState(dir, inputs)?.dirty).toBe(true);
  });

  it('gives no commit outside a repository', () => {
    const dir = mkdtempSync(join(tmpdir(), 'losat-release-id-'));
    dirs.push(dir);
    expect(sourceState(dir, ['.'])).toBeUndefined();
  });
});
