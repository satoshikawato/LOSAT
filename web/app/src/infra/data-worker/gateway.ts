// Main-thread side of the Data worker: a DataGateway, and the engine's describe/validate,
// whose calls run in the worker.
import { DATA_WORKER_METHODS, type DataWorkerApi } from './methods';
import { rpcClient } from './rpc';

/** Starts a Data worker for this tab. It lives as long as the tab (plan §3.1). */
export function startDataWorker(): DataWorkerApi {
  const worker = new Worker(new URL('./data-worker.ts', import.meta.url), { type: 'module', name: 'losat-data' });
  return rpcClient<DataWorkerApi>(worker, DATA_WORKER_METHODS);
}
