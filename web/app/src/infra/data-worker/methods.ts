// The DataGateway methods that the Data worker serves. The compile-time check below keeps
// the list complete.
import type { DataGateway } from '../../ports/data';

export const DATA_GATEWAY_METHODS = [
  'addSource',
  'indexSource',
  'reviseDataset',
  'buildRunInput',
  'openRun',
  'commitRun',
  'discardRun',
  'readOutput',
  'readHits',
  'readDiagnostics',
  'deleteRun',
  'storageInfo',
] as const satisfies readonly (keyof DataGateway)[];

type Missing = Exclude<keyof DataGateway, (typeof DATA_GATEWAY_METHODS)[number]>;
const complete: [Missing] extends [never] ? true : Missing = true;
void complete;
