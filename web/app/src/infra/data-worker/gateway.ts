// Main-thread side of the Data worker: a DataGateway whose calls run in the worker.
import type { DataGateway } from '../../ports/data';
import { DATA_GATEWAY_METHODS } from './methods';
import { rpcClient } from './rpc';

/** Starts a Data worker for this tab. It lives as long as the tab (plan §3.1). */
export function startDataWorker(): DataGateway {
  const worker = new Worker(new URL('./data-worker.ts', import.meta.url), { type: 'module', name: 'losat-data' });
  return rpcClient<DataGateway>(worker, DATA_GATEWAY_METHODS);
}
