// The methods that the Data worker serves: the DataGateway, the saved Example sets
// (ports/example.ts), the saved WindowMasker counts (ports/counts.ts), the inflation of a downloaded Example set (losat/ports: inflate.ts), and
// `describe` / `validate` of the engine, which its serial reactor answers
// (web/packages/losat/src/infra/reactor/control.ts). The compile-time checks below keep the lists complete.
import type { EngineControl } from 'losat';
import type { WholeFileInflater } from 'losat/ports';
import type { CountsStore } from '../../ports/counts';
import type { DataGateway } from '../../ports/data';
import type { ExampleStore } from '../../ports/example';

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

export const COUNTS_STORE_METHODS = [
  'readCounts',
  'beginCounts',
  'appendCounts',
  'finishCounts',
  'closeCounts',
  'deleteCounts',
] as const satisfies readonly (keyof CountsStore)[];

export const INFLATE_METHODS = ['inflateFile'] as const satisfies readonly (keyof WholeFileInflater)[];

export const ENGINE_CONTROL_METHODS = ['describe', 'validate'] as const satisfies readonly (keyof EngineControl)[];

/** Everything the Data worker serves. */
export type DataWorkerApi = DataGateway & ExampleStore & CountsStore & WholeFileInflater & EngineControl;
export const DATA_WORKER_METHODS = [
  ...DATA_GATEWAY_METHODS,
  ...EXAMPLE_STORE_METHODS,
  ...COUNTS_STORE_METHODS,
  ...INFLATE_METHODS,
  ...ENGINE_CONTROL_METHODS,
] as const;

type Missing = Exclude<keyof DataWorkerApi, (typeof DATA_WORKER_METHODS)[number]>;
const complete: [Missing] extends [never] ? true : Missing = true;
void complete;
