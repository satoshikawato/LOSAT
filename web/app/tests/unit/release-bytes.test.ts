// The release ID fixes the build's bytes (S16 decision 1, review M3): a release directory is
// immutable and cached for a year, so two builds of one tree must write the same site byte for
// byte, under the same ID; a build that has no commit to name its release stops.
import { spawnSync } from 'node:child_process';
import { mkdtempSync, readFileSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { afterAll, describe, expect, it } from 'vitest';
import { appBuild, releaseId, sourceState } from '../../build/release';
import { findReactors } from '../../build/reactors';
import { listFiles } from '../../build/site-files';

const APP_DIR = fileURLToPath(new URL('../..', import.meta.url));
const VITE = join(APP_DIR, 'node_modules/vite/bin/vite.js');
const work = mkdtempSync(join(tmpdir(), 'losat-release-bytes-'));

afterAll(() => rmSync(work, { recursive: true, force: true }));

/** `vite build` into `<work>/<name>` (the same build as `npm run build`), with `env` added and LOSAT_WEB_RELEASE_ID removed. */
function build(name: string, env: NodeJS.ProcessEnv = {}): { dir: string; status: number | null; stderr: string } {
  const dir = join(work, name);
  const base: NodeJS.ProcessEnv = { ...process.env };
  delete base['LOSAT_WEB_RELEASE_ID'];
  const result = spawnSync(process.execPath, [VITE, 'build', '--outDir', dir, '--logLevel', 'error'], {
    cwd: APP_DIR,
    env: { ...base, ...env },
    encoding: 'utf8',
  });
  return { dir, status: result.status, stderr: result.stderr };
}

/** Every file of a site directory with its SHA-256, by path. */
const bytesOf = (dir: string) => Object.fromEntries(listFiles(dir).map((file) => [file.path, file.sha256]));

describe('the bytes of a build', () => {
  it('are the same for two builds of one tree, under the ID made from the commit', () => {
    const first = build('first');
    const second = build('second');
    expect(first.status, first.stderr).toBe(0);
    expect(second.status, second.stderr).toBe(0);
    const source = sourceState();
    const release = releaseId(findReactors(), { env: {} });
    const files = bytesOf(first.dir);
    expect(Object.keys(files).some((path) => path.startsWith(`r/${release}/assets/`))).toBe(true);
    expect(bytesOf(second.dir)).toEqual(files);
    const manifest = JSON.parse(readFileSync(join(first.dir, 'r', release, 'release-manifest.json'), 'utf8')) as { app: { build: string } };
    expect(manifest.app.build).toBe(appBuild(source));
    expect(manifest.app.build).toMatch(/^[0-9a-f]{8}$/);
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
