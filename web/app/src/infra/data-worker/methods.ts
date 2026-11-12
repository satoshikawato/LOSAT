// The methods that the Data worker serves: the DataGateway, the saved Example sets
// (ports/example.ts), the inflation of a downloaded Example set (ports/inflate.ts), and
// `describe` / `validate` of the engine, which its serial reactor answers
// (src/infra/reactor/control.ts). The compile-time checks below keep the lists complete.
import type { DataGateway } from '../../ports/data';
import type { ExampleStore } from '../../ports/example';
import type { WholeFileInflater } from '../../ports/inflate';
import type { EngineControl } from '../reactor/control';

export const DATA_GATEWAY_METHODS = [
  'addSource',
  'indexSource',
  'reviseDataset',
  'buildRunInput',
  'checkInput',
  'previewSource',
  'dropSource',
  'openRun',
  'commitRun',
  'discardRun',
  'readOutput',
  'readOutputRange',
  'readHits',
  'readHitTable',
  'readDiagnostics',
  'deleteRun',
  'storageInfo',
  'checkStoredResults',
] as const satisfies readonly (keyof DataGateway)[];

export const EXAMPLE_STORE_METHODS = ['readExample', 'saveExample', 'deleteExample'] as const satisfies readonly (keyof ExampleStore)[];

export const INFLATE_METHODS = ['inflateFile'] as const satisfies readonly (keyof WholeFileInflater)[];

export const ENGINE_CONTROL_METHODS = ['describe', 'validate'] as const satisfies readonly (keyof EngineControl)[];

/** Everything the Data worker serves. */
export type DataWorkerApi = DataGateway & ExampleStore & WholeFileInflater & EngineControl;
export const DATA_WORKER_METHODS = [
  ...DATA_GATEWAY_METHODS,
  ...EXAMPLE_STORE_METHODS,
  ...INFLATE_METHODS,
  ...ENGINE_CONTROL_METHODS,
] as const;

type Missing = Exclude<keyof DataWorkerApi, (typeof DATA_WORKER_METHODS)[number]>;
const complete: [Missing] extends [never] ? true : Missing = true;
void complete;
