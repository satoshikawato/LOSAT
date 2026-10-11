// The boot page's decisions (src/boot/plan.ts): which release starts, when it is pinned, and
// what the page says when it cannot start. A pinned release is never switched silently.
import { describe, expect, it } from 'vitest';
import {
  chooseRelease,
  entryUrls,
  failureMessage,
  judgeManifest,
  PIN_KEY,
  readPin,
  recheckManifest,
  writePin,
} from '../../src/boot/plan';
import { parseReleaseManifest } from '../../src/domain/release-manifest';
import { parseSiteIndex } from '../../src/domain/site-index';

const SHA = 'b'.repeat(64);
const index = parseSiteIndex({
  format: 1,
  latest: 'new',
  releases: [
    { id: 'old', version: '0.0.0', build: 'aaaaaaa' },
    { id: 'new', version: '0.0.1', build: 'bbbbbbb' },
  ],
});
const manifestOf = (release: string) => ({
  format: 1,
  release,
  app: { version: '0.0.1', build: 'bbbbbbb' },
  engine: null,
  entry: { script: 'assets/index-1.js', styles: ['assets/index-2.css'] },
  assets: [
    { path: 'assets/index-1.js', size: 1, sha256: SHA },
    { path: 'assets/index-2.css', size: 1, sha256: SHA },
  ],
});

/** A Storage double; `fail` makes the access itself throw, as a blocked localStorage does. */
function storage(values: Record<string, string>, fail: 'none' | 'access' | 'write' = 'none') {
  return () => {
    if (fail === 'access') throw new DOMException('blocked', 'SecurityError');
    return {
      getItem: (key: string) => values[key] ?? null,
      setItem: (key: string, value: string) => {
        if (fail === 'write') throw new DOMException('full', 'QuotaExceededError');
        values[key] = value;
      },
    };
  };
}

describe('boot: the pin', () => {
  it('reads a pinned release ID, and nothing else', () => {
    expect(readPin(storage({ [PIN_KEY]: 'old' }))).toBe('old');
    expect(readPin(storage({}))).toBeUndefined();
    expect(readPin(storage({ [PIN_KEY]: '../x' }))).toBeUndefined();
    expect(readPin(storage({ [PIN_KEY]: 'old' }, 'access'))).toBeUndefined();
  });

  it('writes the pin, and says when the browser refuses', () => {
    const values: Record<string, string> = {};
    expect(writePin(storage(values), 'new')).toBe(true);
    expect(values[PIN_KEY]).toBe('new');
    expect(writePin(storage({}, 'write'), 'new')).toBe(false);
    expect(writePin(storage({}, 'access'), 'new')).toBe(false);
  });
});

describe('boot: the release to start', () => {
  it('starts the pinned release even when the site has a newer one', () => {
    expect(chooseRelease('old', index, true)).toEqual({ kind: 'start', release: 'old', pin: false });
    expect(chooseRelease('old', undefined, false)).toEqual({ kind: 'start', release: 'old', pin: false });
  });

  it("pins the site's latest release on a first visit", () => {
    expect(chooseRelease(undefined, index, true)).toEqual({ kind: 'start', release: 'new', pin: true });
  });

  it('needs the network on a first visit', () => {
    expect(chooseRelease(undefined, undefined, false)).toMatchObject({
      kind: 'first-visit-offline',
      text: 'LOSAT Web needs a network connection the first time it is opened.',
    });
    expect(chooseRelease(undefined, undefined, true)).toMatchObject({ kind: 'no-index' });
  });
});

describe('boot: the manifest', () => {
  it("accepts the release's own manifest", () => {
    const outcome = judgeManifest('new', { status: 200, body: manifestOf('new'), json: true });
    expect(outcome.kind).toBe('ok');
  });

  it('calls a release not served on 404 or 410, on an answer that is not JSON, and on a wrong manifest', () => {
    expect(judgeManifest('old', { status: 404, body: undefined, json: false })).toMatchObject({ kind: 'not-served' });
    expect(judgeManifest('old', { status: 410, body: undefined, json: false })).toMatchObject({ kind: 'not-served' });
    // A host that answers every path with the boot page (a single-page fallback).
    expect(judgeManifest('old', { status: 200, body: undefined, json: false })).toMatchObject({ kind: 'not-served' });
    expect(judgeManifest('old', { status: 200, body: manifestOf('new'), json: true })).toMatchObject({ kind: 'not-served' });
  });

  it('tells a network failure and a server error apart', () => {
    expect(judgeManifest('old', undefined)).toEqual({ kind: 'unreachable' });
    expect(judgeManifest('old', { status: 503, body: undefined, json: false })).toMatchObject({ kind: 'failed' });
  });
});

describe('boot: an answer that means "not served" is read once more past the HTTP cache', () => {
  it('for 404, 410 and an answer that is not JSON; not for a manifest, a server error or no answer', () => {
    expect(recheckManifest({ status: 404, body: undefined, json: false })).toBe(true);
    expect(recheckManifest({ status: 410, body: undefined, json: false })).toBe(true);
    expect(recheckManifest({ status: 200, body: undefined, json: false })).toBe(true);
    expect(recheckManifest({ status: 200, body: {}, json: true })).toBe(false);
    expect(recheckManifest({ status: 503, body: undefined, json: false })).toBe(false);
    expect(recheckManifest(undefined)).toBe(false);
  });
});

describe('boot: the messages', () => {
  it('offers the latest release for a pinned release that the site no longer serves', () => {
    expect(failureMessage('old', { kind: 'not-served', reason: '404' }, index, true)).toEqual({
      kind: 'not-served',
      text: 'This browser uses LOSAT Web release old, which this site no longer serves.',
      useLatest: 'new',
    });
  });

  it('offers nothing when the site index was not read or names the same release', () => {
    expect(failureMessage('old', { kind: 'not-served', reason: '404' }, undefined, true).useLatest).toBeUndefined();
    expect(failureMessage('new', { kind: 'not-served', reason: '404' }, index, true).useLatest).toBeUndefined();
  });

  it('says that a release is not stored when the page is offline', () => {
    expect(failureMessage('old', { kind: 'unreachable' }, undefined, false)).toMatchObject({ kind: 'offline-not-stored' });
    expect(failureMessage('old', { kind: 'unreachable' }, undefined, true)).toMatchObject({ kind: 'failed' });
  });

  it('loads the entry from the release directory', () => {
    expect(entryUrls(parseReleaseManifest(manifestOf('new')))).toEqual({
      script: '/r/new/assets/index-1.js',
      styles: ['/r/new/assets/index-2.css'],
    });
  });
});
