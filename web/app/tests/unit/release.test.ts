// The release manager (src/application/release.ts): offline use of this tab's release, the
// update check, and the explicit update.
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import type { AppState, RunView } from '../../src/application/coordinator';
import { ReleaseManager, switchImpact, type ReleaseDeps } from '../../src/application/release';
import { Store } from '../../src/application/store';
import type { RunStatus } from '../../src/domain/run';
import type { SiteIndex } from '../../src/domain/site-index';
import type { OfflinePort, PinWrite, PrepareProgress, PrepareResult, ReleaseChannel, ReleasePage, ReleasePin } from '../../src/ports/offline';
import { releaseCachesToRemove } from '../../src/sw/cleanup';
import { RELEASE_CACHE_PREFIX, RELEASE_LOCK_PREFIX } from '../../src/sw/protocol';

const run = (number: number, status: RunStatus): RunView =>
  ({ snapshot: { runId: `run-${number}`, number }, status, record: {} }) as unknown as RunView;

const index = (latest: string): SiteIndex => ({
  format: 1,
  latest,
  releases: [
    { id: 'a', version: '0.0.0', build: 'aaa' },
    { id: latest, version: '0.0.1', build: 'bbb' },
  ],
});

class FakeOffline implements OfflinePort {
  supported = true;
  connectError: Error | undefined;
  held: string[] = [];
  /** The release locks held now, by name, as navigator.locks.query() lists them. */
  locks: string[] = [];
  /** The releases stored whole. */
  storedReleases = new Set<string>();
  prepared: string[] = [];
  cleanups: string[][] = [];
  activations = 0;
  /** What each prepare reports and returns; by default every file arrives. */
  results: PrepareResult[] = [];
  private pending: Array<() => void> = [];
  /** Hold the next preparations until `release()` is called. */
  hold = false;
  connect = vi.fn(async () => {
    if (this.connectError !== undefined) throw this.connectError;
  });
  async prepare(release: string, onProgress: (progress: PrepareProgress) => void): Promise<PrepareResult> {
    this.prepared.push(release);
    onProgress({ files: 0, filesTotal: 2, bytes: 0, bytesTotal: 300 });
    if (this.hold) await new Promise<void>((resolve) => this.pending.push(resolve));
    onProgress({ files: 1, filesTotal: 2, bytes: 100, bytesTotal: 300 });
    const result = this.results.shift() ?? { ok: true };
    if (result.ok) this.storedReleases.add(release);
    return result;
  }
  async stored(release: string): Promise<boolean> {
    return this.storedReleases.has(release);
  }
  /** The cleanup that another tab asks for, by the worker's rule (src/sw/cleanup.ts); returns the releases removed. */
  cleanupElsewhere(keep: readonly string[]): string[] {
    const names = [...this.storedReleases].map((release) => RELEASE_CACHE_PREFIX + release);
    const removed = releaseCachesToRemove(names, { held: this.locks.map((name) => ({ name })) }, keep, []).map((name) => name.slice(RELEASE_CACHE_PREFIX.length));
    for (const release of removed) this.storedReleases.delete(release);
    return removed;
  }
  releaseHeld(): void {
    for (const resolve of this.pending.splice(0)) resolve();
  }
  async cleanup(keep: readonly string[]): Promise<void> {
    this.cleanups.push([...keep]);
  }
  async activateWaiting(): Promise<boolean> {
    this.activations++;
    return true;
  }
  holdRelease(release: string): () => void {
    this.held.push(release);
    const name = RELEASE_LOCK_PREFIX + release;
    this.locks.push(name);
    let holding = true;
    return () => {
      if (!holding) return;
      holding = false;
      this.locks.splice(this.locks.indexOf(name), 1);
    };
  }
}

class FakePin implements ReleasePin {
  value: string | undefined = 'a';
  writeError: string | undefined;
  writes: string[] = [];
  read(): string | undefined {
    return this.value;
  }
  write(release: string): PinWrite {
    this.writes.push(release);
    if (this.writeError !== undefined) return { ok: false, reason: this.writeError };
    this.value = release;
    return { ok: true };
  }
}

class FakeChannel implements ReleaseChannel {
  announced: string[] = [];
  private listeners = new Set<(release: string) => void>();
  announce(release: string): void {
    this.announced.push(release);
  }
  onSwitch(listener: (release: string) => void): () => void {
    this.listeners.add(listener);
    return () => this.listeners.delete(listener);
  }
  /** Another tab switched. */
  receive(release: string): void {
    for (const listener of this.listeners) listener(release);
  }
}

class FakeHost implements ReleasePage {
  isOnline = true;
  reloads = 0;
  private idle: () => void = () => undefined;
  private onlineListeners = new Set<() => void>();
  private idlePromise = new Promise<void>((resolve) => (this.idle = resolve));
  online(): boolean {
    return this.isOnline;
  }
  onOnline(listener: () => void): () => void {
    this.onlineListeners.add(listener);
    return () => this.onlineListeners.delete(listener);
  }
  whenIdle(): Promise<void> {
    return this.idlePromise;
  }
  becomeIdle(): void {
    this.idle();
  }
  goOnline(): void {
    this.isOnline = true;
    for (const listener of this.onlineListeners) listener();
  }
  reload(): void {
    this.reloads++;
  }
}

class FakePage {
  visible = true;
  private listeners = new Set<(visible: boolean) => void>();
  isVisible(): boolean {
    return this.visible;
  }
  onVisibilityChange(listener: (visible: boolean) => void): () => void {
    this.listeners.add(listener);
    return () => this.listeners.delete(listener);
  }
  setVisible(visible: boolean): void {
    this.visible = visible;
    for (const listener of this.listeners) listener(visible);
  }
}

const flush = async (): Promise<void> => {
  for (let i = 0; i < 20; i++) await Promise.resolve();
};

function setup(options: { latest?: string; served?: boolean; runs?: RunView[] } = {}) {
  const offline = new FakeOffline();
  const pin = new FakePin();
  const channel = new FakeChannel();
  const host = new FakeHost();
  const page = new FakePage();
  const site = { latest: options.latest ?? 'a', reads: 0, index: vi.fn(async () => undefined as SiteIndex | undefined) };
  site.index.mockImplementation(async () => {
    site.reads++;
    return index(site.latest);
  });
  const runs = new Store<AppState>({ runs: options.runs ?? [] });
  const deps: ReleaseDeps = {
    release: 'a',
    served: options.served ?? true,
    offline,
    pin,
    channel,
    site,
    page,
    host,
    runs,
    checkEveryMs: 1000,
    onStored: vi.fn(),
  };
  const manager = new ReleaseManager(deps);
  const state = () => manager.state.get();
  /** Starts and lets the page become idle. */
  const start = async () => {
    manager.start();
    host.becomeIdle();
    await flush();
  };
  return { offline, pin, channel, host, page, site, runs, manager, state, start, deps };
}

beforeEach(() => {
  vi.useFakeTimers();
});
afterEach(() => {
  vi.useRealTimers();
});

describe('offline use', () => {
  it('prepares this tab release once the page is idle, then says it is ready and cleans up', async () => {
    const { offline, manager, state, host } = setup();
    expect(state()).toEqual({ release: 'a', offline: { kind: 'waiting' }, update: { kind: 'none' } });
    manager.start();
    expect(offline.held).toEqual(['a']);
    await flush();
    expect(offline.prepared).toEqual([]);
    host.becomeIdle();
    await flush();
    expect(offline.prepared).toEqual(['a']);
    expect(state().offline).toEqual({ kind: 'ready' });
    expect(offline.cleanups).toEqual([['a']]);
  });

  it('tells when a release was stored, so that the storage estimate counts it: at "Offline ready", after the cleanup, after an update’s download', async () => {
    const { deps, offline, manager, state, start } = setup({ latest: 'b' });
    offline.hold = true;
    await start();
    expect(deps.onStored).not.toHaveBeenCalled();
    offline.hold = false;
    offline.releaseHeld();
    await flush();
    expect(state().offline).toEqual({ kind: 'ready' });
    expect(deps.onStored).toHaveBeenCalledTimes(2);
    manager.review();
    offline.results = [{ ok: false, reason: 'a file did not match its checksum (engine/serial.wasm)' }];
    await manager.download();
    expect(deps.onStored).toHaveBeenCalledTimes(2);
    await manager.download();
    expect(state().update).toEqual({ kind: 'downloaded', latest: 'b' });
    expect(deps.onStored).toHaveBeenCalledTimes(3);
  });

  it('shows the progress while the files are checked', async () => {
    const { offline, state, start } = setup();
    offline.hold = true;
    await start();
    expect(state().offline).toEqual({ kind: 'preparing', progress: { files: 0, filesTotal: 2, bytes: 0, bytesTotal: 300 } });
    offline.releaseHeld();
    await flush();
    expect(state().offline).toEqual({ kind: 'ready' });
  });

  it('waits while a run is being prepared, so as not to slow the first search', async () => {
    const { offline, runs, start } = setup({ runs: [run(1, 'preparing')] });
    await start();
    expect(offline.prepared).toEqual([]);
    runs.set({ runs: [run(1, 'running')] });
    await flush();
    expect(offline.prepared).toEqual(['a']);
  });

  it('keeps the pinned release too when another tab pinned another one', async () => {
    const { offline, pin, start } = setup();
    pin.value = 'b';
    await start();
    expect(offline.cleanups).toEqual([['a', 'b']]);
  });

  it('says "Offline ready" only while the pin read back names this tab release', async () => {
    // The browser did not keep the boot page's pin (storage that throws or is full).
    const unpinned = setup();
    unpinned.pin.value = undefined;
    await unpinned.start();
    expect(unpinned.state().offline).toEqual({ kind: 'stored' });

    // Another tab switched before this tab's release was stored, and switches back later.
    const { channel, pin, state, start } = setup();
    pin.value = 'b';
    await start();
    expect(state().offline).toEqual({ kind: 'stored', pin: 'b' });
    pin.value = 'a';
    channel.receive('a');
    expect(state().offline).toEqual({ kind: 'ready' });
    channel.receive('c');
    expect(state().offline).toEqual({ kind: 'stored', pin: 'c' });
  });

  it('says why a preparation failed, and tries again on request', async () => {
    const { offline, manager, state, start } = setup();
    offline.results = [{ ok: false, reason: 'engine/serial.wasm: the server answered 404' }];
    await start();
    expect(state().offline).toEqual({ kind: 'failed', reason: 'engine/serial.wasm: the server answered 404' });
    expect(offline.cleanups).toEqual([]);
    manager.retryOffline();
    await flush();
    expect(state().offline).toEqual({ kind: 'ready' });
    expect(offline.connect).toHaveBeenCalledTimes(1);
  });

  it('is unsupported without service workers', async () => {
    const without = setup();
    without.offline.supported = false;
    const manager = new ReleaseManager({ ...without.deps, offline: without.offline });
    expect(manager.state.get().offline).toEqual({ kind: 'unsupported' });
  });

  it('fails, and tries again on request, when the worker cannot be registered or does not start where service workers exist', async () => {
    const refused = setup();
    refused.offline.connectError = new Error('the offline worker did not start');
    await refused.start();
    expect(refused.state().offline).toEqual({ kind: 'failed', reason: 'the offline worker could not be started (the offline worker did not start)' });
    expect(refused.offline.prepared).toEqual([]);
    refused.offline.connectError = undefined;
    refused.manager.retryOffline();
    await flush();
    expect(refused.offline.connect).toHaveBeenCalledTimes(2);
    expect(refused.state().offline).toEqual({ kind: 'ready' });
  });

  it('does nothing on a page that is not served as a release', async () => {
    const { offline, site, state, start } = setup({ served: false, latest: 'b' });
    await start();
    expect(state().offline).toEqual({ kind: 'unsupported', reason: 'this page is not served as a release of the site' });
    expect(offline.held).toEqual([]);
    expect(offline.prepared).toEqual([]);
    expect(site.index).not.toHaveBeenCalled();
    expect(state().update).toEqual({ kind: 'none' });
  });
});

describe('update check', () => {
  it('offers a release other than this tab', async () => {
    const { state, start } = setup({ latest: 'b' });
    await start();
    expect(state().update).toEqual({ kind: 'available', latest: 'b' });
  });

  it('offers nothing while the site serves this tab release', async () => {
    const { state, start } = setup();
    await start();
    expect(state().update).toEqual({ kind: 'none' });
  });

  it('checks hourly while visible and online, and catches up when the page is shown or online again', async () => {
    const { host, page, site, start, state } = setup();
    await start();
    expect(site.reads).toBe(1);
    await vi.advanceTimersByTimeAsync(1000);
    expect(site.reads).toBe(2);
    page.visible = false;
    site.latest = 'b';
    await vi.advanceTimersByTimeAsync(1000);
    expect(site.reads).toBe(2);
    page.setVisible(true);
    await flush();
    expect(site.reads).toBe(3);
    expect(state().update).toEqual({ kind: 'available', latest: 'b' });
    host.isOnline = false;
    await vi.advanceTimersByTimeAsync(3000);
    expect(site.reads).toBe(3);
    host.goOnline();
    await flush();
    expect(site.reads).toBe(4);
  });

  it('never reads the site index while offline', async () => {
    const { host, site, start } = setup({ latest: 'b' });
    host.isOnline = false;
    await start();
    await vi.advanceTimersByTimeAsync(5000);
    expect(site.index).not.toHaveBeenCalled();
  });

  it('"Not now" holds for this release only', async () => {
    const { manager, site, state, start } = setup({ latest: 'b' });
    await start();
    manager.notNow();
    expect(state().update).toEqual({ kind: 'dismissed', latest: 'b' });
    await manager.check();
    expect(state().update).toEqual({ kind: 'dismissed', latest: 'b' });
    site.latest = 'c';
    await manager.check();
    expect(state().update).toEqual({ kind: 'available', latest: 'c' });
    site.latest = 'a';
    await manager.check();
    expect(state().update).toEqual({ kind: 'none' });
  });

  it('does not offer the release that another tab switched to, and says so', async () => {
    const { channel, state, start } = setup({ latest: 'b' });
    await start();
    expect(state().update).toEqual({ kind: 'available', latest: 'b' });
    channel.receive('b');
    expect(state()).toMatchObject({ elsewhere: 'b', update: { kind: 'none' } });
  });

  it('reads a switch of another tab from the pin when its message was missed', async () => {
    const { pin, state, start } = setup({ latest: 'b' });
    pin.value = 'b';
    await start();
    expect(state()).toMatchObject({ elsewhere: 'b', update: { kind: 'none' } });
  });
});

describe('explicit update', () => {
  async function reviewed(options: Parameters<typeof setup>[0] = {}) {
    const s = setup({ latest: 'b', ...options });
    await s.start();
    s.manager.review();
    return s;
  }

  it('downloads and checks, then switches only when the user chooses', async () => {
    const { offline, pin, channel, host, manager, state } = await reviewed();
    expect(state().update).toEqual({ kind: 'reviewing', latest: 'b' });
    offline.hold = true;
    const downloading = manager.download();
    await flush();
    expect(state().update).toEqual({ kind: 'downloading', latest: 'b', progress: { files: 0, filesTotal: 2, bytes: 0, bytesTotal: 300 } });
    offline.releaseHeld();
    await downloading;
    expect(state().update).toEqual({ kind: 'downloaded', latest: 'b' });
    expect(offline.prepared).toEqual(['a', 'b']);
    // Nothing has changed yet.
    expect(pin.writes).toEqual([]);
    expect(host.reloads).toBe(0);
    await manager.switchNow();
    expect(pin.value).toBe('b');
    expect(channel.announced).toEqual(['b']);
    expect(offline.activations).toBe(1);
    expect(host.reloads).toBe(1);
    expect(state().update).toEqual({ kind: 'switched', latest: 'b' });
    // The browser kept the page (the user chose to stay): it can reload later.
    manager.reloadNow();
    expect(host.reloads).toBe(2);
  });

  it('stays on this release when the download fails, and can try again', async () => {
    const { offline, pin, host, manager, state } = await reviewed();
    offline.results = [{ ok: false, reason: 'engine/serial.wasm: its SHA-256 is not the one that the list gives' }];
    await manager.download();
    expect(state().update).toEqual({ kind: 'failed', latest: 'b', reason: 'engine/serial.wasm: its SHA-256 is not the one that the list gives' });
    await manager.switchNow();
    expect(pin.writes).toEqual([]);
    expect(host.reloads).toBe(0);
    await manager.download();
    expect(state().update).toEqual({ kind: 'downloaded', latest: 'b' });
  });

  it('does not reload when the pin cannot be written', async () => {
    const { pin, channel, host, offline, manager, state } = await reviewed();
    await manager.download();
    pin.writeError = 'QuotaExceededError';
    await manager.switchNow();
    expect(state().update).toEqual({ kind: 'not-kept', latest: 'b', reason: 'QuotaExceededError' });
    expect(pin.value).toBe('a');
    expect(channel.announced).toEqual([]);
    expect(offline.activations).toBe(0);
    expect(host.reloads).toBe(0);
    manager.cancel();
    expect(state().update).toEqual({ kind: 'available', latest: 'b' });
  });

  it('cannot switch before the download is checked', async () => {
    const { pin, host, manager } = await reviewed();
    await manager.switchNow();
    expect(pin.writes).toEqual([]);
    expect(host.reloads).toBe(0);
  });

  it('switches without a download in a browser without service workers', async () => {
    const s = setup({ latest: 'b' });
    s.offline.supported = false;
    const manager = new ReleaseManager(s.deps);
    manager.start();
    s.host.becomeIdle();
    await flush();
    manager.review();
    expect(manager.switchesDirectly()).toBe(true);
    await manager.download();
    expect(s.offline.prepared).toEqual([]);
    await manager.switchNow();
    expect(s.pin.value).toBe('b');
    expect(s.offline.activations).toBe(0);
    expect(s.host.reloads).toBe(1);
  });

  it('holds the downloaded release until the switch, so that another tab’s cleanup keeps it', async () => {
    const { offline, pin, host, manager, state } = await reviewed();
    offline.hold = true;
    const downloading = manager.download();
    await flush();
    // Held from the start of the download, not only while the worker prepares it.
    expect(offline.locks).toContain(`${RELEASE_LOCK_PREFIX}b`);
    offline.releaseHeld();
    await downloading;
    expect(state().update).toEqual({ kind: 'downloaded', latest: 'b' });
    // Another tab of release a starts, is Offline ready, and asks the worker to clean up (keep: a).
    expect(offline.cleanupElsewhere(['a'])).toEqual([]);
    await manager.switchNow();
    expect(pin.value).toBe('b');
    expect(host.reloads).toBe(1);
  });

  it('lets the downloaded release go when the user cancels or says "Not now"', async () => {
    const { offline, manager } = await reviewed();
    await manager.download();
    manager.cancel();
    expect(offline.locks).toEqual([`${RELEASE_LOCK_PREFIX}a`]);
    expect(offline.cleanupElsewhere(['a'])).toEqual(['b']);

    manager.review();
    await manager.download();
    expect(offline.locks).toEqual([`${RELEASE_LOCK_PREFIX}a`, `${RELEASE_LOCK_PREFIX}b`]);
    manager.notNow();
    expect(offline.locks).toEqual([`${RELEASE_LOCK_PREFIX}a`]);
  });

  it('does not switch to a release that is no longer stored, and can download it again', async () => {
    const { offline, pin, channel, host, manager, state } = await reviewed();
    await manager.download();
    // The browser removed it (storage pressure, site data cleared): no lock prevents that.
    offline.storedReleases.delete('b');
    await manager.switchNow();
    expect(state().update).toEqual({ kind: 'removed', latest: 'b' });
    expect(pin.writes).toEqual([]);
    expect(channel.announced).toEqual([]);
    expect(host.reloads).toBe(0);
    await manager.download();
    expect(state().update).toEqual({ kind: 'downloaded', latest: 'b' });
    await manager.switchNow();
    expect(pin.value).toBe('b');
    expect(host.reloads).toBe(1);
  });

  it('treats a browser whose registration shows service workers blocked as one without them', async () => {
    // A run is being prepared, so this tab's release waits; the user asks for the download first.
    const { offline, pin, host, manager, runs, state, start } = setup({ latest: 'b', runs: [run(1, 'preparing')] });
    offline.connect.mockImplementation(async () => {
      offline.supported = false;
      throw new Error('service workers are blocked in this browser');
    });
    await start();
    manager.review();
    await manager.download();
    expect(state().update).toEqual({ kind: 'reviewing', latest: 'b' });
    expect(manager.switchesDirectly()).toBe(true);
    expect(state().offline).toEqual({ kind: 'unsupported', reason: 'service workers are blocked in this browser' });
    runs.set({ runs: [run(1, 'running')] });
    await flush();
    expect(state().offline).toEqual({ kind: 'unsupported', reason: 'service workers are blocked in this browser' });
    await manager.switchNow();
    expect(pin.value).toBe('b');
    expect(host.reloads).toBe(1);
  });

  it('offers no switch until the download is checked when the worker cannot be registered where service workers exist', async () => {
    const { offline, pin, host, manager, state, start } = setup({ latest: 'b' });
    offline.connectError = new Error('the operation is insecure');
    await start();
    manager.review();
    expect(state().update).toEqual({ kind: 'reviewing', latest: 'b' });
    expect(manager.switchesDirectly()).toBe(false);
    await manager.switchNow();
    expect(pin.writes).toEqual([]);
    await manager.download();
    expect(state().update).toEqual({ kind: 'failed', latest: 'b', reason: 'the offline worker could not be started (the operation is insecure)' });
    await manager.switchNow();
    expect(pin.writes).toEqual([]);
    expect(host.reloads).toBe(0);
  });

  it('"Cancel" goes back to the offer, and a check does not disturb a review', async () => {
    const { manager, site, state } = await reviewed();
    site.latest = 'c';
    await manager.check();
    expect(state().update).toEqual({ kind: 'reviewing', latest: 'b' });
    manager.cancel();
    expect(state().update).toEqual({ kind: 'available', latest: 'b' });
  });

  it('a switch in another tab after this tab switched: the pin names the other tab release', async () => {
    const { channel, manager, state } = await reviewed();
    await manager.download();
    await manager.switchNow();
    channel.receive('c');
    expect(state()).toMatchObject({ elsewhere: 'c', update: { kind: 'none' } });
    channel.receive('a');
    expect(state().elsewhere).toBeUndefined();
  });
});

describe('what a switch does to the work of this tab', () => {
  it('counts the runs that the reload stops and what only a session file keeps', () => {
    const runs = [run(1, 'completed'), run(2, 'running'), run(3, 'queued'), run(4, 'queued'), run(5, 'failed'), run(6, 'preparing')];
    expect(switchImpact(runs, [{ note: '' }, { note: 'check this' }, { note: 'and this' }])).toEqual({
      running: 2,
      queued: 2,
      completed: 1,
      candidates: 3,
      notes: 2,
    });
    expect(switchImpact([], [])).toEqual({ running: 0, queued: 0, completed: 0, candidates: 0, notes: 0 });
  });
});
