// A worker thread for tests/tool-credit.test.ts: it posts "bytes" as the Engine worker does (tool-credit.ts) and reports
// each step to the main thread, so that the test sees where it blocks. Run by Node itself (type stripping), so it imports
// the module under test by its file name and uses no syntax that needs a transform.
import { parentPort, workerData } from 'node:worker_threads';
import { ToolCreditWorker } from '../../src/infra/engine-worker/tool-credit.ts';

const { buffer, limit, sliceMs, sizes } = workerData as { buffer: SharedArrayBuffer; limit?: number; sliceMs?: number; sizes: number[] };
const credit = limit === undefined ? new ToolCreditWorker(buffer) : new ToolCreditWorker(buffer, limit, sliceMs);
sizes.forEach((size, index) => {
  credit.add(size);
  parentPort!.postMessage(`posted ${index + 1}`);
  credit.wait();
  parentPort!.postMessage(`after ${index + 1}`);
});
parentPort!.postMessage('done');
