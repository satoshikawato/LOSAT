// The release ID of a build (plan §5.9): the name of the directory `/r/<release-id>/` that
// holds the build's assets, the ID that a browser pins, and the ID that measurement may
// send. `LOSAT_WEB_RELEASE_ID` names it; otherwise it is made from what decides the
// build's bytes: the commit, the engine modules (or `fake` for a FakeEngine build), and
// `-dirty` when web/app has uncommitted changes. The same inputs give the same ID.
import { execSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { fileURLToPath } from 'node:url';
import { RELEASE_ID_PATTERN } from '../src/domain/site-index.ts';
import type { Reactors } from './reactors.ts';

export { RELEASE_ID_PATTERN };

export const RELEASE_ENV = 'LOSAT_WEB_RELEASE_ID';

const APP_DIR = fileURLToPath(new URL('..', import.meta.url));

function git(args: string): string | undefined {
  try {
    return execSync(`git ${args}`, { cwd: APP_DIR, stdio: ['ignore', 'pipe', 'ignore'] }).toString().trim();
  } catch {
    return undefined;
  }
}

/** The engine part of a release ID: `fake`, or a digest of the two served modules. */
export function engineTag(reactors: Reactors | undefined): string {
  if (reactors === undefined) return 'fake';
  return createHash('sha256').update(`${reactors.serial.sha256}\n${reactors.threads.sha256}\n`).digest('hex').slice(0, 8);
}

/** The release ID of this build; throws when `LOSAT_WEB_RELEASE_ID` does not match RELEASE_ID_PATTERN. */
export function releaseId(reactors: Reactors | undefined, env: NodeJS.ProcessEnv = process.env): string {
  const named = env[RELEASE_ENV];
  if (named !== undefined && named !== '') {
    if (!RELEASE_ID_PATTERN.test(named)) throw new Error(`${RELEASE_ENV}=${named}: not ${String(RELEASE_ID_PATTERN)}`);
    return named;
  }
  const commit = git('rev-parse --short=12 HEAD') ?? 'unknown';
  const dirty = (git('status --porcelain -- .') ?? '') === '' ? '' : '-dirty';
  return `${commit}-${engineTag(reactors)}${dirty}`;
}
