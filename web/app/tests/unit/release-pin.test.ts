// The pinned release in localStorage (src/infra/browser/offline.ts), and the names that every
// release must keep: a release that renamed them would lose the pin, or delete what an older
// tab still uses (plan §5.9; design §7.2 and S10 decision 8 for tmp/).
import { describe, expect, it } from 'vitest';
import { SESSION_LOCK_PREFIX, TMP_DIRECTORY } from '../../src/infra/data/session';
import { CHANNEL_NAME, PIN_KEY, releasePin } from '../../src/infra/browser/offline';
import { RELEASE_CACHE_PREFIX, RELEASE_LOCK_PREFIX, SHELL_CACHE, SITE_INDEX_URL, SW_SCOPE, SW_URL, MANIFEST_NAME } from '../../src/sw/protocol';

class MemoryStorage {
  readonly values = new Map<string, string>();
  getItem(key: string): string | null {
    return this.values.get(key) ?? null;
  }
  setItem(key: string, value: string): void {
    this.values.set(key, value);
  }
}

const throwing = (what: string) => () => {
  throw new DOMException(what, 'SecurityError');
};

describe('the pinned release', () => {
  it('reads and writes losat-web:release', () => {
    const storage = new MemoryStorage();
    const pin = releasePin(() => storage as unknown as Storage);
    expect(pin.read()).toBeUndefined();
    expect(pin.write('0a1b2c3d4e5f-fake')).toEqual({ ok: true });
    expect(storage.values.get('losat-web:release')).toBe('0a1b2c3d4e5f-fake');
    expect(pin.read()).toBe('0a1b2c3d4e5f-fake');
  });

  it('reads no pin from a value that is not a release ID', () => {
    const storage = new MemoryStorage();
    const pin = releasePin(() => storage as unknown as Storage);
    for (const value of ['', 'Upper', '../r/x', 'a'.repeat(65), '<script>']) {
      storage.values.set(PIN_KEY, value);
      expect(pin.read()).toBeUndefined();
    }
  });

  it('survives storage that throws, is missing, or does not keep the value', () => {
    const denied = releasePin(throwing('The operation is insecure.'));
    expect(denied.read()).toBeUndefined();
    expect(denied.write('b')).toEqual({ ok: false, reason: 'The operation is insecure.' });

    const missing = releasePin(() => undefined);
    expect(missing.read()).toBeUndefined();
    expect(missing.write('b')).toEqual({ ok: false, reason: 'this browser gives the page no storage' });

    const full = new MemoryStorage();
    full.setItem = () => {
      throw new DOMException('The quota has been exceeded.', 'QuotaExceededError');
    };
    expect(releasePin(() => full as unknown as Storage).write('b')).toEqual({ ok: false, reason: 'The quota has been exceeded.' });

    const forgetful = new MemoryStorage();
    forgetful.setItem = () => undefined;
    expect(releasePin(() => forgetful as unknown as Storage).write('b')).toEqual({ ok: false, reason: 'the browser did not keep it' });

    const readOnly = new MemoryStorage();
    readOnly.getItem = throwing('denied');
    expect(releasePin(() => readOnly as unknown as Storage).read()).toBeUndefined();
  });
});

describe('names shared by every release', () => {
  it('keeps the pin, the channel, the release locks and caches, and the worker', () => {
    expect(PIN_KEY).toBe('losat-web:release');
    expect(CHANNEL_NAME).toBe('losat-web:release');
    expect(RELEASE_LOCK_PREFIX).toBe('losat-web:release:');
    expect(RELEASE_CACHE_PREFIX).toBe('losat-web-release-');
    expect(SHELL_CACHE).toBe('losat-web-shell');
    expect(MANIFEST_NAME).toBe('release-manifest.json');
    expect(SITE_INDEX_URL).toBe('/releases.json');
    expect([SW_URL, SW_SCOPE]).toEqual(['/sw.js', '/']);
  });

  it('keeps the tmp/ directory and its lock names, so that a newer release never deletes an older tab data', () => {
    expect(TMP_DIRECTORY).toBe('tmp');
    expect(SESSION_LOCK_PREFIX).toBe('losat-web:tmp:');
  });
});
