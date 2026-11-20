// Main-thread side of the Data worker: a DataGateway, and the engine's describe/start/validate,
// whose calls run in the worker. A worker that could not start, or that ended, rejects every call
// with an UnavailableError, and the page is told (line WQA's sweep, D1 and D2, 2026-11-19): it
// ends without an error event when the browser ends it, which the release of the Web Lock of its
// session shows (worker-life.ts in the package). Before, its calls waited for good: a run stayed
// "Organizing results", a pasted text "Reading the records…" and Run LOSAT "Starting…". Where the
// browser refuses that lock (a site whose data is blocked, in Chromium and Firefox), nothing shows
// such an end, and its calls wait for good as before (worker-life.ts, holdLife).
import { UnavailableError, watchLife, webLocks } from 'losat';
import { DATA_WORKER_METHODS, type DataWorkerApi } from './methods';
import { rpcClient } from './rpc';

/** Why every call fails when the worker's script could not be loaded or did not start (line WQA's sweep, D1). */
export const DATA_NOT_STARTED = 'LOSAT could not start the part that keeps the results';
/** Why every call fails once the worker ended (line WQA's sweep, D2); the name is ResumeNotice's. */
export const DATA_ENDED = 'The part of LOSAT that keeps the results stopped';

/**
 * Starts a Data worker for this tab. It lives as long as the tab (plan §3.1): `onEnd` is told
 * when it could not start or ended, with the error that its calls get from then on.
 */
export function startDataWorker(onEnd?: (error: Error) => void): DataWorkerApi {
  const worker = new Worker(new URL('./data-worker.ts', import.meta.url), { type: 'module', name: 'losat-data' });
  let end: (error: Error) => void = () => undefined;
  const ended = new Promise<Error>((resolve) => (end = resolve));
  const data = rpcClient<DataWorkerApi>(worker, DATA_WORKER_METHODS, {
    notStarted: () => new UnavailableError(DATA_NOT_STARTED),
    ended,
    ...(onEnd === undefined ? {} : { onEnd }),
  });
  const locks = webLocks();
  if (locks !== undefined) {
    data.lifeLock().then(
      (name) => {
        if (name !== undefined) watchLife(locks, name, undefined, () => end(new UnavailableError(DATA_ENDED)));
      },
      // A worker that could not start has said so.
      () => undefined,
    );
  }
  return data;
}
