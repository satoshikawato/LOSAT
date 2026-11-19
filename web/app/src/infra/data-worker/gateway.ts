// Main-thread side of the Data worker: a DataGateway, and the engine's describe/start/validate,
// whose calls run in the worker. A worker that could not start rejects every call with an
// UnavailableError, and the page is told (line WQA's sweep, D1, 2026-11-19): before, each input
// card said "This input cannot be read: the worker stopped".
import { UnavailableError } from 'losat';
import { DATA_WORKER_METHODS, type DataWorkerApi } from './methods';
import { rpcClient } from './rpc';

/** Why every call fails when the worker's script could not be loaded or did not start (line WQA's sweep, D1). */
export const DATA_NOT_STARTED = 'LOSAT could not start the part that keeps the results';

/**
 * Starts a Data worker for this tab. It lives as long as the tab (plan §3.1): `onEnd` is told
 * when it could not start, with the error that its calls get from then on.
 */
export function startDataWorker(onEnd?: (error: Error) => void): DataWorkerApi {
  const worker = new Worker(new URL('./data-worker.ts', import.meta.url), { type: 'module', name: 'losat-data' });
  return rpcClient<DataWorkerApi>(worker, DATA_WORKER_METHODS, {
    notStarted: () => new UnavailableError(DATA_NOT_STARTED),
    ...(onEnd === undefined ? {} : { onEnd }),
  });
}
