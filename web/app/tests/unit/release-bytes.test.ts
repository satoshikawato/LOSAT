// The release ID fixes the build's bytes (S16 decision 1, review M3): a release directory is
// immutable and cached for a year, so two builds of one commit must write the same site byte for
// byte, under the same ID, wherever the commit is checked out (a CI build and a local build of one
// commit; code review 2, L5); a build that has no commit to name its release stops.
import { spawnSync } from 'node:child_process';
import { cpSync, existsSync, mkdtempSync, readFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { afterAll, describe, expect, it } from 'vitest';
import { APP, BUILD_INPUTS, REPOSITORY } from '../../build/inputs';
import { appBuild, releaseId, sourceState } from '../../build/release';
import { findReactors } from '../../build/reactors';
import { listFiles } from '../../build/site-files';

const APP_DIR = fileURLToPath(new URL('../..', import.meta.url));
const work = mkdtempSync(join(tmpdir(), 'losat-release-bytes-'));
/** The second checkout's folder while it is a worktree of the repository. */
let worktree: string | undefined;

const git = (args: readonly string[]) => spawnSync('git', ['-C', REPOSITORY, ...args], { encoding: 'utf8' });

afterAll(() => {
  if (worktree !== undefined && git(['worktree', 'remove', '--force', worktree]).status !== 0) {
    rmSync(worktree, { recursive: true, force: true });
    git(['worktree', 'prune']);
  }
  rmSync(work, { recursive: true, force: true });
});

/**
 * What Vitest adds to the environment (and VITE_USER_NODE_ENV, which a shell may export).
 * NODE_ENV=test above all would make the build a development build (Vue's plugin then writes each
 * component's absolute path into `__file`), which is neither the site nor the same in two
 * checkouts; vite.config.ts stops such a build.
 */
const VITEST_ENV = ['NODE_ENV', 'VITE_USER_NODE_ENV', 'TEST', 'VITEST', 'VITEST_MODE', 'VITEST_POOL_ID', 'VITEST_WORKER_ID', 'FORCE_TTY', 'BASE_URL', 'MODE', 'DEV', 'PROD', 'SSR'];

/**
 * `vite build` of the application in `app` (this checkout's by default) into `<work>/<name>` (the
 * same build as `npm run build` in a shell that sets neither NODE_ENV nor LOSAT_WEB_RELEASE_ID),
 * with `env` added.
 */
function build(name: string, env: NodeJS.ProcessEnv = {}, app = APP_DIR, args: readonly string[] = []): { dir: string; status: number | null; stderr: string } {
  const dir = join(work, name);
  const base: NodeJS.ProcessEnv = { ...process.env };
  for (const key of ['LOSAT_WEB_RELEASE_ID', ...VITEST_ENV]) delete base[key];
  const result = spawnSync(process.execPath, [join(app, 'node_modules/vite/bin/vite.js'), 'build', '--outDir', dir, '--logLevel', 'error', ...args], {
    cwd: app,
    env: { ...base, ...env },
    encoding: 'utf8',
  });
  return { dir, status: result.status, stderr: result.stderr };
}

let firstBuild: ReturnType<typeof build> | undefined;
/** The first build of this checkout, made once for the tests that compare with it. */
const first = () => (firstBuild ??= build('first'));

/** Every file of a site directory with its SHA-256, by path. */
const bytesOf = (dir: string) => Object.fromEntries(listFiles(dir).map((file) => [file.path, file.sha256]));

/**
 * A second checkout of the commit at another path: a worktree of HEAD under the system's temporary
 * folder with only the build's inputs (build/inputs.ts) checked out, so that its build also shows
 * that the build reads nothing else of the repository, and a copy of this checkout's node_modules
 * (what `npm ci` installs from the same package-lock.json). The reason why there is none when the
 * inputs have uncommitted changes (a checkout of the commit would not hold them), or when git cannot
 * add a worktree.
 */
function secondCheckout(): { app: string } | { skip: string } {
  const source = sourceState();
  if (source === undefined) return { skip: 'no git commit to check out' };
  if (source.dirty) return { skip: "the build's inputs (build/inputs.ts) have changes that are not committed: a checkout of the commit would not hold them" };
  const root = join(work, 'second-checkout');
  const added = git(['worktree', 'add', '--no-checkout', '--detach', root, source.commit]);
  if (added.status !== 0) return { skip: `git cannot add a worktree: ${added.stderr.trim()}` };
  worktree = root;
  const checkout = spawnSync('git', ['-C', root, 'checkout', source.commit, '--', ...BUILD_INPUTS], { encoding: 'utf8' });
  expect(checkout.status, checkout.stderr).toBe(0);
  const app = join(root, APP);
  cpSync(join(APP_DIR, 'node_modules'), join(app, 'node_modules'), { recursive: true, verbatimSymlinks: true });
  return { app };
}

describe('the bytes of a build', () => {
  it('are the same for two builds of one tree, under the ID made from the commit', () => {
    const second = build('second');
    expect(first().status, first().stderr).toBe(0);
    expect(second.status, second.stderr).toBe(0);
    const source = sourceState();
    const release = releaseId(findReactors(), { env: {} });
    const files = bytesOf(first().dir);
    expect(Object.keys(files).some((path) => path.startsWith(`r/${release}/assets/`))).toBe(true);
    expect(bytesOf(second.dir)).toEqual(files);
    const manifest = JSON.parse(readFileSync(join(first().dir, 'r', release, 'release-manifest.json'), 'utf8')) as { app: { build: string } };
    expect(manifest.app.build).toBe(appBuild(source));
    expect(manifest.app.build).toMatch(/^[0-9a-f]{8}$/);
  }, 180_000);

  it('are the same from a second checkout of the commit at another path that holds only the build inputs', (context) => {
    const checkout = secondCheckout();
    if ('skip' in checkout) return context.skip(checkout.skip);
    expect(checkout.app).not.toBe(APP_DIR);
    const there = build('second-checkout-build', {}, checkout.app);
    expect(there.status, there.stderr).toBe(0);
    expect(first().status, first().stderr).toBe(0);
    const release = releaseId(findReactors(), { env: {} });
    expect(release).not.toMatch(/-dirty$/);
    const files = bytesOf(first().dir);
    expect(Object.keys(files).some((path) => path.startsWith(`r/${release}/assets/`))).toBe(true);
    expect(bytesOf(there.dir)).toEqual(files);
  }, 180_000);

  // vite.config.ts reads what Vite resolved (code review 3 L1): VITE_USER_NODE_ENV (or NODE_ENV in a
  // .env file, which Vite turns into it) is read after the config function has run, and still stops
  // the build before it writes anything.
  it('stops a development build or another mode however Vite decided it, before writing anything', () => {
    for (const [name, env, args, problem] of [
      ['user-node-env', { VITE_USER_NODE_ENV: 'development' }, [], /NODE_ENV=development, VITE_USER_NODE_ENV=development .*: this would be a development build/],
      ['node-env-test', { NODE_ENV: 'test' }, [], /NODE_ENV=test: this would be a development build/],
      ['mode', {}, ['--mode', 'staging'], /vite build --mode staging: a release is built in production mode only/],
    ] as const) {
      const refused = build(name, { ...env, LOSAT_WEB_RELEASE_ID: 'refused-1' }, APP_DIR, args);
      expect(refused.status, name).not.toBe(0);
      expect(refused.stderr, name).toMatch(problem);
      expect(existsSync(refused.dir), name).toBe(false);
    }
  }, 180_000);

  it('stops without a commit unless LOSAT_WEB_RELEASE_ID names the release', () => {
    // GIT_DIR names no repository: git finds no commit, as in a tree without git.
    const noGit = { GIT_DIR: join(work, 'no-repository') };
    const refused = build('refused', noGit);
    expect(refused.status).not.toBe(0);
    expect(refused.stderr).toMatch(/no git commit to name its release: set LOSAT_WEB_RELEASE_ID/);
    const named = build('named', { ...noGit, LOSAT_WEB_RELEASE_ID: 'named-1' });
    expect(named.status, named.stderr).toBe(0);
    const manifest = JSON.parse(readFileSync(join(named.dir, 'r', 'named-1', 'release-manifest.json'), 'utf8')) as { app: { build: string } };
    expect(manifest.app.build).toBe('unknown');
  }, 180_000);
});
