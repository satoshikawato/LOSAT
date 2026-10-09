// The browser's page lifecycle (ports/page.ts): Screen Wake Lock, `beforeunload` and
// `visibilitychange`.
import type { PagePort, WakeLockHandle } from '../../ports/page';

interface WakeLockSentinelLike extends EventTarget {
  release(): Promise<void>;
}

interface WakeLockLike {
  request(type: 'screen'): Promise<WakeLockSentinelLike>;
}

function wakeLock(): WakeLockLike | undefined {
  return (navigator as Navigator & { wakeLock?: WakeLockLike }).wakeLock;
}

function guard(event: BeforeUnloadEvent): void {
  // Browsers show their own text; the page cannot set it.
  event.preventDefault();
  event.returnValue = '';
}

export const browserPage: PagePort = {
  wakeLockSupported: wakeLock() !== undefined,

  async requestWakeLock(): Promise<WakeLockHandle> {
    const api = wakeLock();
    if (api === undefined) throw new Error('this browser cannot keep the screen on');
    const sentinel = await api.request('screen');
    return {
      release: () => sentinel.release(),
      onRelease: (listener) => sentinel.addEventListener('release', () => listener(), { once: true }),
    };
  },

  setLeaveGuard(on: boolean): void {
    if (on) window.addEventListener('beforeunload', guard);
    else window.removeEventListener('beforeunload', guard);
  },

  isVisible(): boolean {
    return document.visibilityState === 'visible';
  },

  onVisibilityChange(listener: (visible: boolean) => void): () => void {
    const handler = () => listener(document.visibilityState === 'visible');
    document.addEventListener('visibilitychange', handler);
    return () => document.removeEventListener('visibilitychange', handler);
  },
};
