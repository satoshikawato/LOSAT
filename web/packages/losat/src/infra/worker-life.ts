// The life of a worker as a Web Lock (line WQA's sweep, D2, 2026-11-19). A worker that ends without an error event
// (Worker.terminate() and self.close() look so to the page, as does a browser that ends a worker by itself) left the
// page waiting for its answers for good: a search or its results spun forever. Such a worker holds a Web Lock of its
// own for as long as it lives, and the page waits for that lock: the browser releases it when the worker ends, also
// while the worker is busy, and never while it lives, however long a synchronous call runs (Chromium, Firefox and
// WebKit: losat-baselines/web-qa-20261119/fix-workers/probe/lock-release.mjs). Nothing passes between them while the
// worker lives, so a search pays nothing for it, and a worker that cannot take a message for minutes (the Engine
// worker during a search, whose call into the engine is synchronous) is never taken for ended, as an answer with a
// time limit would take it.

/** The part of the Web Locks API (`navigator.locks`) that the life of a worker uses. */
export interface LifeLocks {
  request(name: string, options: { readonly signal?: AbortSignal }, callback: () => unknown): Promise<unknown>;
}

/** The Web Locks of this page or worker, or undefined where the browser has none. */
export function webLocks(): LifeLocks | undefined {
  const locks = (globalThis.navigator as { locks?: LifeLocks } | undefined)?.locks;
  return typeof locks?.request === 'function' ? locks : undefined;
}

/**
 * In the worker: takes the lock `name` and holds it until the worker ends. Resolves with whether it holds it: not
 * without Web Locks, nor when the browser refused the request. Chromium and Firefox refuse it with SecurityError where
 * the site's data is blocked (line WQA's sweep, D8, 2026-11-20). The worker works then, but nothing tells the page that
 * it ended without an error event: no browser signals the end of a worker on a MessagePort either (no close event, no
 * error when posting to it; the review of D8, losat-baselines/web-qa-20261119/wqa2/review-d8-probe/close.log). Such an
 * end leaves the page as before D2: a search of the Engine worker waits until it is cancelled, and the calls of
 * the Data worker (its lock is that of its session, data-worker.ts in the application) wait for good.
 */
export function holdLife(locks: LifeLocks | undefined, name: string): Promise<boolean> {
  if (locks === undefined) return Promise.resolve(false);
  return new Promise<boolean>((resolve) => {
    locks
      .request(name, {}, () => {
        resolve(true);
        return new Promise<void>(() => {});
      })
      .catch(() => resolve(false));
  });
}

/**
 * On the page: calls `ended` once the worker that holds the lock `name` (holdLife, or another lock that it holds for as
 * long as it lives) has ended, unless `signal` aborts first, as it does when the page ends the worker itself. Call it
 * only once the worker holds the lock, else the page takes the lock first.
 */
export function watchLife(locks: LifeLocks, name: string, signal: AbortSignal | undefined, ended: () => void): void {
  locks
    .request(name, signal === undefined ? {} : { signal }, () => {
      if (signal?.aborted !== true) ended();
    })
    // Aborted: the page ended the worker itself.
    .catch(() => undefined);
}
