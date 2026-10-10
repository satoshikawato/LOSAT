// The methods that the Data worker serves: the DataGateway, and `describe` / `validate`
// of the engine, which its serial reactor answers (src/infra/reactor/control.ts). The
// compile-time checks below keep the lists complete.
import type { DataGateway } from '../../ports/data';
import type { EngineControl } from '../reactor/control';

export const DATA_GATEWAY_METHODS = [
  'addSource',
  'indexSource',
  'reviseDataset',
  'buildRunInput',
  'checkInput',
  'previewSource',
  'readResidues',
  'openRun',
  'commitRun',
  'discardRun',
  'readOutput',
  'readOutputRange',
  'readHits',
  'readHspRecords',
  'readHitTable',
  'readDiagnostics',
  'deleteRun',
  'storageInfo',
] as const satisfies readonly (keyof DataGateway)[];

export const ENGINE_CONTROL_METHODS = ['describe', 'validate'] as const satisfies readonly (keyof EngineControl)[];

/** Everything the Data worker serves. */
export type DataWorkerApi = DataGateway & EngineControl;
export const DATA_WORKER_METHODS = [...DATA_GATEWAY_METHODS, ...ENGINE_CONTROL_METHODS] as const;

type Missing = Exclude<keyof DataWorkerApi, (typeof DATA_WORKER_METHODS)[number]>;
const complete: [Missing] extends [never] ? true : Missing = true;
void complete;
