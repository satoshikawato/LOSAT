// Page lifecycle port (design §9.4, REQ-22): the screen wake lock, the warning before the
// page is left, and the page's visibility. The browser implementation is
// src/infra/browser/page.ts; tests use a fake.

export interface WakeLockHandle {
  release(): Promise<void>;
  /** Called once when the lock ends, by `release` or by the browser (for example when the page is hidden). */
  onRelease(listener: () => void): void;
}

export interface PagePort {
  /** Whether the browser has the Screen Wake Lock API. */
  readonly wakeLockSupported: boolean;
  /** Rejects when the browser refuses the lock (for example on low battery). */
  requestWakeLock(): Promise<WakeLockHandle>;
  /** While on, the browser asks before the page is closed or reloaded. */
  setLeaveGuard(on: boolean): void;
  isVisible(): boolean;
  /** Returns a function that removes the listener. */
  onVisibilityChange(listener: (visible: boolean) => void): () => void;
}
