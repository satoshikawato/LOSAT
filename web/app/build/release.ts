// The release ID of a build (plan §5.9): the name of the directory `/r/<release-id>/` that
// holds the build's assets, the ID that a browser pins, and the ID that measurement may
// send. A release directory is immutable and cached for a year, so the ID must fix the
// build's bytes. `LOSAT_WEB_RELEASE_ID` names it; otherwise it is made from what decides
// them: the commit (its first 12 hex digits, whatever length git would abbreviate it to in
// this clone), the engine modules (or `fake` for a FakeEngine build), and `-dirty` when any
// input of the build (build/inputs.ts) differs from the commit. A build without a commit must
// be named; only the development server and Vitest run without one.
import { execFileSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { RELEASE_ID_PATTERN } from '../src/domain/site-index.ts';
import { BUILD_INPUTS, REPOSITORY } from './inputs.ts';
import type { Reactors } from './reactors.ts';

export { RELEASE_ID_PATTERN };

export const RELEASE_ENV = 'LOSAT_WEB_RELEASE_ID';

/** The hex digits of the commit in a release ID. */
export const RELEASE_COMMIT_LENGTH = 12;
/** The hex digits of the commit in the manifest's `app.build` and in session files (`__LOSAT_APP_BUILD__`). */
export const APP_BUILD_LENGTH = 8;

/** The commit of the working tree and whether a build input differs from it. */
export interface SourceState {
  /** The full hex object name of HEAD. */
  readonly commit: string;
  readonly dirty: boolean;
}

function git(repository: string, args: readonly string[]): string | undefined {
  try {
    return execFileSync('git', ['-C', repository, ...args], { stdio: ['ignore', 'pipe', 'ignore'] }).toString().trim();
  } catch {
    return undefined;
  }
}

/**
 * The commit of `repository` and whether any of `inputs` (paths relative to it) has a change
 * that is not committed, untracked files included; undefined when there is no commit (no git).
 * A status that cannot be read counts as dirty.
 */
export function sourceState(repository: string = REPOSITORY, inputs: readonly string[] = BUILD_INPUTS): SourceState | undefined {
  const commit = git(repository, ['rev-parse', '--verify', '--quiet', 'HEAD^{commit}']);
  if (commit === undefined || !/^(?:[0-9a-f]{40}|[0-9a-f]{64})$/.test(commit)) return undefined;
  const status = git(repository, ['status', '--porcelain', '--', ...inputs]);
  return { commit, dirty: status !== '' };
}

/** The engine part of a release ID: `fake`, or a digest of the two served modules. */
export function engineTag(reactors: Reactors | undefined): string {
  if (reactors === undefined) return 'fake';
  return createHash('sha256').update(`${reactors.serial.sha256}\n${reactors.threads.sha256}\n`).digest('hex').slice(0, 8);
}

export interface ReleaseIdOptions {
  readonly env?: NodeJS.ProcessEnv;
  /** The commit and the inputs' state (default: this repository's); undefined without git. */
  readonly source?: SourceState | undefined;
  /** False for the development server and Vitest, whose pages are not a site: they may run without a commit. */
  readonly build?: boolean;
}

/**
 * The release ID of this build. Throws when `LOSAT_WEB_RELEASE_ID` does not match
 * RELEASE_ID_PATTERN, and when a build has neither that name nor a commit.
 */
export function releaseId(reactors: Reactors | undefined, options: ReleaseIdOptions = {}): string {
  const env = options.env ?? process.env;
  const named = env[RELEASE_ENV];
  if (named !== undefined && named !== '') {
    if (!RELEASE_ID_PATTERN.test(named)) throw new Error(`${RELEASE_ENV}=${named}: not ${String(RELEASE_ID_PATTERN)}`);
    return named;
  }
  const source = 'source' in options ? options.source : sourceState();
  if (source === undefined) {
    if (options.build === false) return `unknown-${engineTag(reactors)}`;
    throw new Error(`this build has no git commit to name its release: set ${RELEASE_ENV} (docs/web/release_manifest.md)`);
  }
  return `${source.commit.slice(0, RELEASE_COMMIT_LENGTH)}-${engineTag(reactors)}${source.dirty ? '-dirty' : ''}`;
}

/** The manifest's `app.build` and `__LOSAT_APP_BUILD__`: the commit's first 8 hex digits, or `unknown` without git. */
export function appBuild(source: SourceState | undefined): string {
  return source === undefined ? 'unknown' : source.commit.slice(0, APP_BUILD_LENGTH);
}

/** What Vite resolved for a command (the config of a `configResolved` hook), as far as the build rule reads it. */
export interface ResolvedBuild {
  readonly command: 'build' | 'serve';
  readonly mode: string;
  readonly isProduction: boolean;
}

/**
 * Why this `vite build` must stop, or undefined: a release is a production build in production
 * mode. Vite decides that only after the config function has run: from NODE_ENV, else from
 * VITE_USER_NODE_ENV, which it also sets from a NODE_ENV in a `.env*` file of web/app (Vite 8.3.1
 * `resolveConfig`, code review 3 L1); so the rule reads the final word, `config.isProduction`, in
 * a `configResolved` hook (vite.config.ts). A development build (Vitest exports NODE_ENV=test; a
 * shell or a `.env` file may say development) makes the Vue plugin write each component's
 * absolute source path into the bundle (`__file`): a path of the build machine in a served file
 * (TD-11), and bytes that differ between clones, so the release ID would no longer fix the bytes.
 * Another `--mode` reads other `.env.<mode>` files and changes `import.meta.env.MODE` under the
 * same ID. The development server and Vitest (`serve`) run as they are.
 */
export function productionBuildProblem(resolved: ResolvedBuild, env: NodeJS.ProcessEnv = process.env): string | undefined {
  if (resolved.command !== 'build') return undefined;
  if (resolved.mode !== 'production') {
    return `vite build --mode ${resolved.mode}: a release is built in production mode only. Build without --mode (docs/web/deploy_cloudflare.md).`;
  }
  if (resolved.isProduction) return undefined;
  const userNodeEnv = env['VITE_USER_NODE_ENV'];
  const causes = [
    `NODE_ENV=${env['NODE_ENV'] ?? ''}`,
    ...(userNodeEnv === undefined || userNodeEnv === '' ? [] : [`VITE_USER_NODE_ENV=${userNodeEnv} (exported, or NODE_ENV in a .env file)`]),
  ];
  return (
    `${causes.join(', ')}: this would be a development build, with the build machine's paths in the served files. ` +
    'Unset NODE_ENV and VITE_USER_NODE_ENV (or set NODE_ENV to production), take NODE_ENV out of any .env file of web/app, ' +
    'and build again (docs/web/deploy_cloudflare.md).'
  );
}
