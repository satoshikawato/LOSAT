// The methods that the Data worker serves: the DataGateway, the saved Example sets
// (ports/example.ts), the saved WindowMasker counts (ports/counts.ts), the inflation of a downloaded Example set (losat/ports: inflate.ts), and
// `describe` / `start` / `validate` of the engine, which its serial reactor answers
// (web/packages/losat/src/infra/reactor/control.ts), and the name of the lock by which the page learns that it ended.
// The compile-time checks below keep the lists complete.
import type { EngineControl, ToolControl } from 'losat';
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
  'keepFailedRun',
  'readOutput',
  'readOutputRange',
  'readHits',
  'readHitTable',
  'readDiagnostics',
  'deleteRun',
  'storageInfo',
  'checkStoredResults',
  'openToolRun',
  'readToolOutputRange',
  'getToolOutputFile',
  'inflateToSessionFile',
  'dropSessionFile',
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

export const ENGINE_CONTROL_METHODS = [
  'describe',
  'describeTool',
  'start',
  'validate',
] as const satisfies readonly (keyof (EngineControl & ToolControl))[];

/** How the page learns that the Data worker ended (gateway.ts; line WQA's sweep, D2). */
export interface DataWorkerLife {
  /** The name of the Web Lock that the worker holds for as long as it lives (its session's, session.ts), or undefined when it holds none. */
  lifeLock(): Promise<string | undefined>;
}

export const LIFE_METHODS = ['lifeLock'] as const satisfies readonly (keyof DataWorkerLife)[];

/** Everything the Data worker serves. */
export type DataWorkerApi = DataGateway &
  ExampleStore &
  CountsStore &
  WholeFileInflater &
  EngineControl &
  ToolControl &
  DataWorkerLife;
export const DATA_WORKER_METHODS = [
  ...DATA_GATEWAY_METHODS,
  ...EXAMPLE_STORE_METHODS,
  ...COUNTS_STORE_METHODS,
  ...INFLATE_METHODS,
  ...ENGINE_CONTROL_METHODS,
  ...LIFE_METHODS,
] as const;

type Missing = Exclude<keyof DataWorkerApi, (typeof DATA_WORKER_METHODS)[number]>;
const complete: [Missing] extends [never] ? true : Missing = true;
void complete;
